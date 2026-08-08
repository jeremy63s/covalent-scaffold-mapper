"""Command-line entry point.

Anything not given as a flag is asked for, so the tool works both as a one-liner
in a script and as a guided session for someone running it the first time.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import mzid_parser, scoring
from .io_utils import PathError, confirm, resolve_input_path, resolve_output_dir
from .pipeline import Settings, run


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="scaffold-binding",
        description=(
            "Find candidate covalent binding pockets by testing whether a "
            "scaffold's captured peptides cluster on the 3D structure."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    source = parser.add_argument_group("input and output")
    source.add_argument("--mzid", help="Scaffold / mzIdentML export (.mzid or .mzid.gz)")
    source.add_argument("--output", help="directory to write results into")
    source.add_argument(
        "--structures",
        help="directory of PDB files named by accession, e.g. P08670.pdb",
    )
    source.add_argument("--cache", help="where downloaded structures are kept")

    samples = parser.add_argument_group("samples")
    samples.add_argument(
        "--treatment", nargs="+", help="sample names from the real scaffold"
    )
    samples.add_argument(
        "--control",
        nargs="+",
        default=None,
        help="optional negative control samples (scrambled scaffold, beads-only)",
    )
    samples.add_argument(
        "--no-control",
        action="store_true",
        help="run without a negative control and do not ask for one",
    )
    samples.add_argument("--list-samples", action="store_true",
                         help="print the sample names in the export and exit")

    model = parser.add_argument_group("scaffold model")
    model.add_argument(
        "--scaffold-size",
        type=float,
        help="reach of the scaffold in angstroms; the radius of the fitted shell",
    )
    model.add_argument(
        "--rigidity",
        type=float,
        help=(
            "0 = fully flexible (nothing inside the shell is penalised), "
            "1 = fully rigid (only the shell surface is unpenalised)"
        ),
    )

    tuning = parser.add_argument_group("filters and tuning")
    tuning.add_argument("--min-psms", type=int, default=1,
                        help="spectra a peptide needs before it counts as signal")
    tuning.add_argument("--min-sites", type=int, default=scoring.MIN_SITES,
                        help="capture sites a protein needs to be scored")
    tuning.add_argument("--plddt-floor", type=float, default=scoring.DEFAULT_PLDDT_FLOOR,
                        help="lowest model confidence a residue may have")
    tuning.add_argument("--permutations", type=int, default=scoring.DEFAULT_PERMUTATIONS,
                        help="random draws used for the significance test")
    tuning.add_argument("--max-proteins", type=int, default=None,
                        help="stop after this many proteins (useful for a trial run)")
    tuning.add_argument("--seed", type=int, default=0, help="random seed")

    behaviour = parser.add_argument_group("behaviour")
    behaviour.add_argument("--no-download", action="store_true",
                           help="use only local structures, never fetch AlphaFold")
    behaviour.add_argument("--esmfold", action="store_true",
                           help="fold proteins AlphaFold lacks with ESMFold")
    behaviour.add_argument("--no-plots", action="store_true",
                           help="skip the interactive 3D output")
    behaviour.add_argument("--plot-top", type=int, default=5,
                           help="how many top-ranked proteins to plot")

    return parser


def _ask_float(prompt: str, default: float, low: float, high: float) -> float:
    """Read a number in range, re-asking on bad input, defaulting when piped."""
    if not sys.stdin.isatty():
        return default
    while True:
        raw = input(f"{prompt} [{default:g}]: ").strip()
        if not raw:
            return default
        try:
            value = float(raw)
        except ValueError:
            print("  Please enter a number.")
            continue
        if not low <= value <= high:
            print(f"  Please enter a value between {low:g} and {high:g}.")
            continue
        return value


def _choose_samples(
    available: list[str], treatment: list[str] | None, control: list[str] | None,
    no_control: bool,
) -> tuple[list[str], list[str]]:
    """Work out which runs are treatment and which are control."""
    unknown = set(treatment or []) | set(control or [])
    unknown -= set(available)
    if unknown:
        raise SystemExit(
            f"These sample names are not in the export: {sorted(unknown)}\n"
            f"Available: {available}"
        )

    if treatment:
        return treatment, ([] if no_control else (control or []))

    print("\nSamples found in this export:")
    for index, name in enumerate(available, start=1):
        print(f"  {index}. {name}")

    if not sys.stdin.isatty():
        raise SystemExit(
            "Pass --treatment (and optionally --control) when running "
            "non-interactively."
        )

    def pick(prompt: str, allow_empty: bool) -> list[str]:
        while True:
            raw = input(prompt).strip()
            if not raw:
                if allow_empty:
                    return []
                print("  At least one sample is required.")
                continue
            chosen = []
            for token in raw.replace(",", " ").split():
                if token.isdigit() and 1 <= int(token) <= len(available):
                    chosen.append(available[int(token) - 1])
                elif token in available:
                    chosen.append(token)
                else:
                    print(f"  Not a listed sample: {token}")
                    break
            else:
                return chosen

    treatment = pick("Treatment samples (numbers or names, space separated): ", False)
    if no_control:
        return treatment, []

    remaining = [s for s in available if s not in treatment]
    if not remaining:
        return treatment, []
    if not confirm("\nDo you have a negative control sample?", default=False):
        return treatment, []
    control = pick("Control samples (blank for none): ", True)
    return treatment, [c for c in control if c not in treatment]


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    try:
        mzid_path = resolve_input_path(
            args.mzid, "Path to the Scaffold / mzIdentML export (.mzid or .mzid.gz)"
        )
    except PathError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    available = mzid_parser.list_samples(mzid_parser.parse_mzid(mzid_path))
    if args.list_samples:
        print("\n".join(available))
        return 0

    treatment, control = _choose_samples(
        available, args.treatment, args.control, args.no_control
    )

    try:
        output_dir = resolve_output_dir(
            args.output,
            "Directory to write results into",
            default=Path.cwd() / "scaffold_binding_results",
        )
        structure_dir = (
            resolve_input_path(args.structures, "Directory of PDB structures",
                               expect_dir=True)
            if args.structures
            else None
        )
    except PathError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    scaffold_size = args.scaffold_size
    if scaffold_size is None:
        print(
            "\nScaffold size is how far the scaffold reaches from its anchor "
            "point, in angstroms.\nIt sets the radius of the shell fitted "
            "around each candidate pocket."
        )
        scaffold_size = _ask_float(
            "Scaffold size (A)", scoring.DEFAULT_SCAFFOLD_SIZE, 1.0, 200.0
        )

    rigidity = args.rigidity
    if rigidity is None:
        print(
            "\nRigidity is how much the scaffold can fold back on itself.\n"
            "  0.0  fully flexible (PEG-like): any site within reach is fine\n"
            "  0.7  sites from 0.7x to 1x the scaffold size are unpenalised\n"
            "  1.0  fully rigid: only sites at exactly the scaffold size fit"
        )
        rigidity = _ask_float("Rigidity (0-1)", scoring.DEFAULT_RIGIDITY, 0.0, 1.0)

    settings = Settings(
        mzid_path=mzid_path,
        output_dir=output_dir,
        treatment_samples=treatment,
        control_samples=control,
        scaffold_size=scaffold_size,
        rigidity=rigidity,
        structure_dir=structure_dir,
        cache_dir=Path(args.cache).expanduser() if args.cache else None,
        allow_download=not args.no_download,
        allow_esmfold=args.esmfold,
        min_psms=args.min_psms,
        min_sites=args.min_sites,
        plddt_floor=args.plddt_floor,
        permutations=args.permutations,
        max_proteins=args.max_proteins,
        plot_top=args.plot_top,
        make_plots=not args.no_plots,
        random_seed=args.seed,
    )

    print()
    results = run(settings)
    if not results.empty:
        columns = ["gene", "accession", "n_signal_sites", "fraction_reachable",
                   "z_score", "p_value"]
        if settings.has_control:
            columns.append("signal_vs_control")
        print("\nTop candidates:")
        print(results[columns].head(15).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
