"""Run the whole analysis end to end and write the results out."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import mzid_parser, peptides, scoring
from .structures import StructureProvider


@dataclass
class Settings:
    """Everything the run needs, already resolved to real paths and values."""

    mzid_path: Path
    output_dir: Path
    treatment_samples: list[str]
    control_samples: list[str]
    scaffold_size: float = scoring.DEFAULT_SCAFFOLD_SIZE
    rigidity: float = scoring.DEFAULT_RIGIDITY
    structure_dir: Path | None = None
    cache_dir: Path | None = None
    allow_download: bool = True
    allow_esmfold: bool = False
    min_psms: int = 1
    min_sites: int = scoring.MIN_SITES
    plddt_floor: float = scoring.DEFAULT_PLDDT_FLOOR
    permutations: int = scoring.DEFAULT_PERMUTATIONS
    max_proteins: int | None = None
    plot_top: int = 5
    make_plots: bool = True
    random_seed: int = 0

    @property
    def has_control(self) -> bool:
        return bool(self.control_samples)


def run(settings: Settings, log=print) -> pd.DataFrame:
    """Execute the pipeline; returns the ranked results table."""
    settings.output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = settings.cache_dir or (settings.output_dir / "structures")

    log(f"Reading {settings.mzid_path.name} ...")
    psms = mzid_parser.parse_mzid(settings.mzid_path)
    psms = mzid_parser.drop_contaminants_and_decoys(psms)
    protein_info = mzid_parser.parse_proteins(settings.mzid_path)
    log(f"  {len(psms)} peptide-spectrum matches across "
        f"{psms['accession'].nunique()} proteins")

    table = peptides.build_peptide_table(
        psms,
        settings.treatment_samples,
        settings.control_samples,
        min_psms=settings.min_psms,
    )
    log(peptides.summarize(table, has_control=settings.has_control))

    candidates = peptides.proteins_with_enough_signal(table, settings.min_sites)
    if settings.max_proteins:
        candidates = candidates[: settings.max_proteins]
    log(f"\nProteins with at least {settings.min_sites} capture sites: {len(candidates)}")
    if not candidates:
        log("Nothing to score. Try lowering --min-sites or --min-psms.")
        return pd.DataFrame()

    provider = StructureProvider(
        cache_dir,
        local_dir=settings.structure_dir,
        allow_download=settings.allow_download,
        allow_esmfold=settings.allow_esmfold,
    )

    genes = protein_info.set_index("accession")["gene"].to_dict()
    rng = np.random.default_rng(settings.random_seed)
    scores, skipped = [], []

    log(f"\nScoring with scaffold size {settings.scaffold_size:g} A, "
        f"rigidity {settings.rigidity:g}")
    for index, accession in enumerate(candidates, start=1):
        gene = genes.get(accession)
        label = f"{gene} ({accession})" if gene else accession
        positions = (
            table[(table["accession"] == accession) & table["is_signal"]]["start"]
            .astype(int)
            .tolist()
        )
        coords = provider.get(accession, focus_positions=positions)
        if coords is None or coords.empty:
            skipped.append((accession, provider.notes.get(accession, "no structure")))
            log(f"  [{index}/{len(candidates)}] {label}: skipped, no structure")
            continue

        score = scoring.score_protein(
            accession,
            table,
            coords,
            scaffold_size=settings.scaffold_size,
            rigidity=settings.rigidity,
            plddt_floor=settings.plddt_floor,
            n_permutations=settings.permutations,
            min_sites=settings.min_sites,
            gene=gene,
            rng=rng,
        )
        if score is None:
            skipped.append((accession, "too few sites in the modelled region"))
            log(f"  [{index}/{len(candidates)}] {label}: skipped, sites fall outside the model")
            continue

        scores.append(score)
        log(
            f"  [{index}/{len(candidates)}] {label}: "
            f"{score.n_signal_sites} sites, {score.fraction_reachable:.0%} in reach, "
            f"z = {score.z_score:.2f}, p = {score.p_value:.3f}"
        )

    results = scoring.score_to_frame(scores)
    if results.empty:
        log("\nNo protein could be scored.")
        return results

    results_path = settings.output_dir / "binding_site_scores.csv"
    results.to_csv(results_path, index=False)
    log(f"\nWrote {results_path}")

    peptide_path = settings.output_dir / "peptide_sites.csv"
    table[table["is_signal"] | table["is_background"]].to_csv(peptide_path, index=False)
    log(f"Wrote {peptide_path}")

    if skipped:
        skipped_path = settings.output_dir / "skipped_proteins.csv"
        pd.DataFrame(skipped, columns=["accession", "reason"]).to_csv(
            skipped_path, index=False
        )
        log(f"Wrote {skipped_path} ({len(skipped)} proteins)")

    if settings.make_plots and settings.plot_top > 0:
        _write_plots(settings, results, scores, table, provider, log)

    return results


def _write_plots(settings, results, scores, table, provider, log) -> None:
    try:
        import plotly  # noqa: F401
    except ImportError:
        log("\nSkipping plots: plotly is not installed (pip install plotly).")
        return

    from .visualize import plot_protein

    by_accession = {s.accession: s for s in scores}
    plot_dir = settings.output_dir / "plots"
    top = results.head(settings.plot_top)["accession"].tolist()
    log(f"\nWriting {len(top)} interactive plots ...")

    for accession in top:
        score = by_accession.get(accession)
        coords = provider.get(accession)
        if score is None or coords is None:
            continue
        path = plot_protein(
            score, table, coords, plot_dir, plddt_floor=settings.plddt_floor
        )
        if path:
            log(f"  {path}")
