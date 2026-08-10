"""Collect a scaffold's chemistry from the user, one backbone bond at a time.

Three numbers per bond, all of them table lookups: how long the bond is, what
angle it makes, and how it twists. Only the twist needs judgement, and only in
one case -- an ordinary single bond is RIS and wants rotamer energies, while a
double bond, an amide, or an aromatic bond is rigid and wants nothing.

If the rotamer energies cannot be found, the persistence-length route asks for
two numbers instead of a form and is tabulated for anything common.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from . import reach
from .reach import ReachModel

# Bond lengths in angstroms, for the prompt's suggestions.
COMMON_BONDS = "C-C 1.53, C-O 1.43, C-N 1.47, amide C-N 1.33, C=C 1.34, aromatic 1.39"
COMMON_ANGLES = "sp3 (single-bonded C or O) 109-112, sp2 (double/aromatic) 120"


def _interactive() -> bool:
    return sys.stdin is not None and sys.stdin.isatty()


def _ask(prompt: str, default: str | None = None) -> str:
    suffix = f" [{default}]" if default is not None else ""
    while True:
        raw = input(f"{prompt}{suffix}: ").strip()
        if raw:
            return raw
        if default is not None:
            return default
        print("  This one is required.")


def _ask_float(prompt: str, default: float | None = None,
               low: float = -1e9, high: float = 1e9) -> float:
    while True:
        raw = _ask(prompt, None if default is None else f"{default:g}")
        try:
            value = float(raw)
        except ValueError:
            print("  Please enter a number.")
            continue
        if not low <= value <= high:
            print(f"  Please enter a value between {low:g} and {high:g}.")
            continue
        return value


def _ask_int(prompt: str, default: int | None = None, low: int = 1) -> int:
    while True:
        raw = _ask(prompt, None if default is None else str(default))
        try:
            value = int(raw)
        except ValueError:
            print("  Please enter a whole number.")
            continue
        if value < low:
            print(f"  Please enter a value of at least {low}.")
            continue
        return value


def _ask_choice(prompt: str, options: dict[str, str], default: str) -> str:
    keys = "/".join(options)
    while True:
        raw = _ask(f"{prompt} ({keys})", default).lower()
        for key in options:
            if raw == key or raw == key[0]:
                return key
        print(f"  Please answer one of: {keys}")


def ask_polymer() -> ReachModel:
    """Walk the user through defining their scaffold, and build a reach model."""
    if not _interactive():
        raise RuntimeError(
            "No scaffold definition was given and there is no terminal to ask on. "
            "Pass --polymer, or --persistence-length together with --scaffold-size."
        )

    print("\n" + "=" * 70)
    print("Scaffold definition")
    print("=" * 70)
    print(
        "The scaffold is a polymer, so where its free end can land is set by its\n"
        "chemistry. Two ways to describe it:\n\n"
        "  bonds  Fill in one line per backbone bond of the repeat unit. Bond\n"
        "         lengths and angles come straight off a chemistry table; the\n"
        "         only judgement is whether each bond twists freely or is locked.\n"
        "         Most faithful, especially for short chains.\n\n"
        "  lp     Give the polymer's persistence length and total length instead.\n"
        "         Both are tabulated for essentially every common polymer, so it\n"
        "         is one lookup. Slightly less faithful for very short chains."
    )
    route = _ask_choice("\nWhich route", {"bonds": "", "lp": ""}, "bonds")

    if route == "lp":
        print(
            "\nPersistence length is how far along the chain it stays pointing the\n"
            "same way -- the standard measure of stiffness. Contour length is the\n"
            "chain's full stretched-out length."
        )
        lp = _ask_float("Persistence length (A)", low=0.1, high=1e4)
        contour = _ask_float("Contour length of the whole scaffold (A)",
                             low=0.1, high=1e5)
        name = _ask("Name for this scaffold", "scaffold")
        return reach.from_persistence_length(contour, lp, label=name)

    print(
        "\nOne entry per backbone bond in a single repeat unit.\n"
        f"  Typical lengths (A): {COMMON_BONDS}\n"
        f"  Typical angles (deg): {COMMON_ANGLES}\n"
        "  Twist: 'ris' for an ordinary single bond, 'rigid' for a double bond,\n"
        "         amide, or aromatic bond, 'free' for no barrier at all.\n"
        "  RIS rotamer energies are the cost of a kink versus staying straight,\n"
        "  in kcal/mol, from a polymer reference (Flory, the Polymer Handbook).\n"
        "  A typical alkane gauche penalty is 0.5-0.9; PEG's O-C-C-O gauche\n"
        "  effect goes the other way and favours the kink."
    )

    name = _ask("\nName of the repeat unit", "monomer")
    n_bonds = _ask_int("How many backbone bonds in one repeat unit", low=1)

    bonds = []
    for i in range(1, n_bonds + 1):
        print(f"\n-- bond {i} of {n_bonds} --")
        length = _ask_float("  Bond length (A)", low=0.1, high=10.0)
        angle = _ask_float("  Bond angle (deg)", 109.5, low=1.0, high=179.9)
        kind = _ask_choice("  Twist", {"ris": "", "rigid": "", "free": ""}, "ris")

        if kind == "free":
            torsion: dict = {"type": "free"}
        elif kind == "rigid":
            fixed = _ask_float("  Locked dihedral (deg)", 180.0, low=-180.0, high=360.0)
            torsion = {"type": "rigid", "fixed_deg": fixed}
        else:
            print("  Rotamer energies relative to trans, in kcal/mol.")
            trans = _ask_float("    trans (deg 180) energy", 0.0)
            gauche_plus = _ask_float("    gauche+ (deg 60) energy", 0.5)
            gauche_minus = _ask_float("    gauche- (deg -60) energy", 0.5)
            torsion = {
                "type": "ris",
                "states_deg": [180, 60, -60],
                "energies_kcal": [trans, gauche_plus, gauche_minus],
            }
        bonds.append({"length_A": length, "angle_deg": angle, "torsion": torsion})

    temperature = _ask_float("\nTemperature (K)", 298.0, low=1.0, high=1000.0)
    n_units = _ask_int(
        "How many repeat units make up the scaffold (PEG8 would be 8)", low=1
    )

    spec = {
        "name": name,
        "temperature_K": temperature,
        "bonds": bonds,
        "n_units": n_units,
    }

    save = _ask("\nSave this definition to a file for reuse? (path, or blank to skip)", "")
    if save:
        path = Path(save).expanduser()
        path.write_text(json.dumps(spec, indent=2))
        print(f"  Saved {path} -- reuse it with --polymer {path}")

    model = reach.from_spec(spec, n_units)
    print(f"\n  {model.describe()}")
    return model


def build_from_args(args) -> ReachModel:
    """Resolve the scaffold from command-line flags, asking only if none given."""
    if args.polymer:
        return reach.from_json(args.polymer, args.monomer_units)

    if args.persistence_length is not None:
        if args.scaffold_size is None:
            raise RuntimeError(
                "--persistence-length needs --scaffold-size too: the contour "
                "length of the whole scaffold, in angstroms"
            )
        return reach.from_persistence_length(
            args.scaffold_size, args.persistence_length
        )

    return ask_polymer()
