"""Collapse per-spectrum matches into one row per peptide, and label them.

A negative control (scrambled scaffold, beads-only, DMSO, ...) is optional. When
one is supplied, a peptide counts as signal only if it never appears in the
control. When none is supplied, every confidently detected peptide is signal and
the spatial test alone decides whether a protein looks like a real binder.
"""

from __future__ import annotations

import pandas as pd

SIGNAL_COLUMNS = [
    "accession",
    "peptide_sequence",
    "start",
    "stop",
    "treatment_psms",
    "control_psms",
    "enrichment",
    "is_signal",
    "is_background",
]


def build_peptide_table(
    psms: pd.DataFrame,
    treatment_samples: list[str],
    control_samples: list[str] | None = None,
    *,
    min_psms: int = 1,
) -> pd.DataFrame:
    """One row per (protein, peptide, position) with treatment/control counts.

    `min_psms` is the number of spectra a peptide needs across the treatment
    samples before it is trusted as signal. It matters most when there is no
    control, since it is then the only abundance filter in play.
    """
    control_samples = control_samples or []
    overlap = set(treatment_samples) & set(control_samples)
    if overlap:
        raise ValueError(
            f"Samples cannot be both treatment and control: {sorted(overlap)}"
        )
    if not treatment_samples:
        raise ValueError("At least one treatment sample is required.")

    usable = psms.dropna(subset=["accession", "start", "stop"]).copy()
    usable["start"] = usable["start"].astype(int)
    usable["stop"] = usable["stop"].astype(int)

    known = set(treatment_samples) | set(control_samples)
    usable = usable[usable["sample"].isin(known)]
    if usable.empty:
        raise ValueError(
            "No peptide-spectrum matches remain after restricting to the "
            "selected samples."
        )

    counts = (
        usable.groupby(
            ["accession", "peptide_sequence", "start", "stop", "sample"]
        )
        .size()
        .reset_index(name="psm_count")
    )

    wide = counts.pivot_table(
        index=["accession", "peptide_sequence", "start", "stop"],
        columns="sample",
        values="psm_count",
        fill_value=0,
    ).reset_index()
    wide.columns.name = None

    for sample in known:
        if sample not in wide.columns:
            wide[sample] = 0

    wide["treatment_psms"] = wide[treatment_samples].sum(axis=1)
    wide["control_psms"] = (
        wide[control_samples].sum(axis=1) if control_samples else 0
    )
    wide["enrichment"] = wide["treatment_psms"] / (wide["control_psms"] + 1)

    detected = wide["treatment_psms"] >= min_psms
    if control_samples:
        wide["is_signal"] = detected & (wide["control_psms"] == 0)
        wide["is_background"] = wide["control_psms"] > 0
    else:
        wide["is_signal"] = detected
        wide["is_background"] = False

    return wide.reset_index(drop=True)


def summarize(table: pd.DataFrame, *, has_control: bool) -> str:
    """A short human-readable description of what the labelling produced."""
    lines = [
        f"Peptide-position rows: {len(table)}",
        f"Signal peptides: {int(table['is_signal'].sum())}",
    ]
    if has_control:
        lines.append(f"Control (background) peptides: {int(table['is_background'].sum())}")
        both = int(((table["treatment_psms"] > 0) & (table["control_psms"] > 0)).sum())
        lines.append(f"Seen in both treatment and control: {both}")
    else:
        lines.append("No control supplied - all detected peptides treated as signal.")
    proteins = table.loc[table["is_signal"], "accession"].nunique()
    lines.append(f"Proteins with at least one signal peptide: {proteins}")
    return "\n".join(lines)


def proteins_with_enough_signal(
    table: pd.DataFrame, min_sites: int = 4
) -> list[str]:
    """Accessions carrying enough distinct signal positions to fit a sphere."""
    signal = table[table["is_signal"]]
    per_protein = signal.groupby("accession")["start"].nunique()
    return sorted(per_protein[per_protein >= min_sites].index.tolist())
