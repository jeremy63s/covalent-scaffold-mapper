"""Score how well a protein's capture sites fit a scaffold-shaped shell.

The model is a sphere of fixed radius centred on a candidate binding pocket. The
radius is the scaffold's reach (`scaffold_size`, in angstroms), so it is set by
the chemistry rather than fitted to the data -- only the pocket centre is fitted.

How far a site may sit *inside* that shell is governed by `rigidity`:

    rigidity 0.0   a floppy scaffold (PEG8 and similar) can fold back on itself,
                   so a site anywhere inside the shell is equally consistent with
                   the pocket and nothing inside is penalised.
    rigidity 0.7   sites between 0.7R and R are free; closer in than 0.7R the
                   penalty grows with how far past that inner ring the site sits.
    rigidity 1.0   the scaffold cannot bend at all, so anything off the shell
                   surface is penalised.

Sites *outside* the shell are penalised by their distance beyond it regardless of
rigidity -- no amount of flexibility lets a scaffold reach further than its length.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import minimize

DEFAULT_SCAFFOLD_SIZE = 25.0
DEFAULT_RIGIDITY = 0.5
DEFAULT_PLDDT_FLOOR = 50.0
DEFAULT_PERMUTATIONS = 200
MIN_SITES = 4

# A site this close to the shell counts as "on" it when reporting fractions.
ON_SHELL_TOLERANCE = 0.0

# A scaffold can only reach residues on the outside of a protein, so the null has
# to draw from those too. Buried residues sit much closer together than surface
# ones and would fit any small sphere easily, which would make chance look better
# than it is. Burial is estimated by counting neighbouring alpha-carbons within
# this radius; residues above the cutoff percentile are treated as interior.
BURIAL_RADIUS = 12.0
BURIAL_PERCENTILE = 60.0


def exposed_residues(
    coords: pd.DataFrame,
    *,
    plddt_floor: float = DEFAULT_PLDDT_FLOOR,
    burial_radius: float = BURIAL_RADIUS,
    burial_percentile: float = BURIAL_PERCENTILE,
) -> np.ndarray:
    """Coordinates of confidently modelled residues near the protein surface."""
    confident = coords[coords["plddt"] >= plddt_floor]
    points = confident[["x", "y", "z"]].to_numpy()
    if len(points) < 3:
        return coords[["x", "y", "z"]].to_numpy()

    distances = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    neighbours = (distances < burial_radius).sum(axis=1) - 1
    cutoff = np.percentile(neighbours, burial_percentile)
    exposed = points[neighbours <= cutoff]
    return exposed if len(exposed) >= 3 else points


def penalty(
    distances: np.ndarray, scaffold_size: float, rigidity: float
) -> np.ndarray:
    """Per-site penalty in angstroms for the given distances to the pocket centre.

    Zero means the site is somewhere the scaffold could plausibly reach.
    """
    distances = np.asarray(distances, dtype=float)
    inner_edge = rigidity * scaffold_size
    too_far = np.maximum(0.0, distances - scaffold_size)
    too_close = np.maximum(0.0, inner_edge - distances)
    return too_far + too_close


def _mean_penalty(center: np.ndarray, points: np.ndarray, size: float, rigidity: float) -> float:
    distances = np.linalg.norm(points - center, axis=1)
    return float(np.mean(penalty(distances, size, rigidity)))


def _seed_centers(points: np.ndarray, scaffold_size: float) -> list[np.ndarray]:
    """Starting guesses for the pocket centre.

    Capture sites sit on the outside of a protein, so the pocket that anchored
    them is usually offset from their centroid by roughly one scaffold length --
    along the spread of the sites, or straight out along the axis they vary least
    in. Both directions are tried, as is the plain centroid.
    """
    centroid = points.mean(axis=0)
    seeds = [centroid]

    if len(points) >= 3:
        centred = points - centroid
        # Principal axes of the site cloud; the last one is the shell normal.
        _, _, axes = np.linalg.svd(centred, full_matrices=False)
        for axis in axes:
            seeds.append(centroid + scaffold_size * axis)
            seeds.append(centroid - scaffold_size * axis)

    return seeds


def fit_center(
    points: np.ndarray,
    scaffold_size: float,
    rigidity: float,
    *,
    max_seeds: int | None = None,
) -> tuple[np.ndarray, float]:
    """Find the pocket centre that minimises the mean penalty over `points`."""
    seeds = _seed_centers(points, scaffold_size)
    if max_seeds is not None:
        seeds = seeds[:max_seeds]

    best_center, best_score = None, np.inf
    for seed in seeds:
        result = minimize(
            _mean_penalty,
            seed,
            args=(points, scaffold_size, rigidity),
            method="Nelder-Mead",
            options={"maxiter": 2000, "xatol": 0.05, "fatol": 0.001},
        )
        if result.fun < best_score:
            best_center, best_score = result.x, float(result.fun)

    return best_center, best_score


@dataclass
class ProteinScore:
    """Everything the pipeline learned about one protein."""

    accession: str
    gene: str | None = None
    n_signal_sites: int = 0
    n_control_sites: int = 0
    center: tuple[float, float, float] | None = None
    scaffold_size: float = DEFAULT_SCAFFOLD_SIZE
    rigidity: float = DEFAULT_RIGIDITY
    mean_penalty: float = float("nan")
    fraction_reachable: float = float("nan")
    null_mean: float = float("nan")
    null_std: float = float("nan")
    z_score: float = float("nan")
    p_value: float = float("nan")
    control_mean_penalty: float | None = None
    control_fraction_reachable: float | None = None
    signal_vs_control: float | None = None
    note: str = ""

    def as_row(self) -> dict:
        row = {
            "accession": self.accession,
            "gene": self.gene,
            "n_signal_sites": self.n_signal_sites,
            "n_control_sites": self.n_control_sites,
            "center_x": self.center[0] if self.center is not None else None,
            "center_y": self.center[1] if self.center is not None else None,
            "center_z": self.center[2] if self.center is not None else None,
            "scaffold_size": self.scaffold_size,
            "rigidity": self.rigidity,
            "mean_penalty": self.mean_penalty,
            "fraction_reachable": self.fraction_reachable,
            "null_mean_penalty": self.null_mean,
            "null_std_penalty": self.null_std,
            "z_score": self.z_score,
            "p_value": self.p_value,
            "control_mean_penalty": self.control_mean_penalty,
            "control_fraction_reachable": self.control_fraction_reachable,
            "signal_vs_control": self.signal_vs_control,
            "note": self.note,
        }
        return row


def site_coordinates(
    peptides: pd.DataFrame,
    coords: pd.DataFrame,
    *,
    plddt_floor: float = DEFAULT_PLDDT_FLOOR,
    search_offsets: tuple[int, ...] = (0, 1, 2, -1),
) -> tuple[np.ndarray, list[dict]]:
    """Map each peptide to the 3D position of its N-terminal residue.

    The N-terminus is where the scaffold's reactive group attached, so that is
    the point the geometry is about. Residues modelled with low confidence are
    skipped, and a couple of neighbouring positions are tried when the exact
    residue is missing from the structure.
    """
    lookup = coords.set_index("residue_num")
    points, metadata = [], []

    for _, row in peptides.iterrows():
        start = int(row["start"])
        for offset in search_offsets:
            residue = start + offset
            if residue not in lookup.index:
                continue
            entry = lookup.loc[residue]
            if isinstance(entry, pd.DataFrame):
                entry = entry.iloc[0]
            if float(entry["plddt"]) < plddt_floor:
                break
            points.append([float(entry["x"]), float(entry["y"]), float(entry["z"])])
            metadata.append(
                {
                    "peptide_sequence": row.get("peptide_sequence"),
                    "start": start,
                    "stop": row.get("stop"),
                    "treatment_psms": row.get("treatment_psms"),
                    "control_psms": row.get("control_psms"),
                    "plddt": float(entry["plddt"]),
                }
            )
            break

    return np.array(points) if points else np.empty((0, 3)), metadata


def score_protein(
    accession: str,
    peptide_table: pd.DataFrame,
    coords: pd.DataFrame,
    *,
    scaffold_size: float = DEFAULT_SCAFFOLD_SIZE,
    rigidity: float = DEFAULT_RIGIDITY,
    plddt_floor: float = DEFAULT_PLDDT_FLOOR,
    n_permutations: int = DEFAULT_PERMUTATIONS,
    min_sites: int = MIN_SITES,
    gene: str | None = None,
    rng: np.random.Generator | None = None,
) -> ProteinScore | None:
    """Fit a pocket to one protein's signal sites and test it against chance.

    The null asks: if the same number of capture sites had landed anywhere on
    this protein's modelled surface, how well would the best pocket have fitted?
    Each random draw is fitted from scratch with the same objective, so the
    comparison is like-for-like.
    """
    rng = rng or np.random.default_rng(0)
    protein_peptides = peptide_table[peptide_table["accession"] == accession]
    signal = protein_peptides[protein_peptides["is_signal"]]
    control = protein_peptides[protein_peptides["is_background"]]

    signal_points, _ = site_coordinates(signal, coords, plddt_floor=plddt_floor)
    if len(signal_points) < min_sites:
        return None

    center, observed = fit_center(signal_points, scaffold_size, rigidity)
    distances = np.linalg.norm(signal_points - center, axis=1)
    penalties = penalty(distances, scaffold_size, rigidity)

    score = ProteinScore(
        accession=accession,
        gene=gene,
        n_signal_sites=len(signal_points),
        center=tuple(float(v) for v in center),
        scaffold_size=scaffold_size,
        rigidity=rigidity,
        mean_penalty=observed,
        fraction_reachable=float(np.mean(penalties <= ON_SHELL_TOLERANCE)),
    )

    # Null distribution: same number of sites, drawn from surface residues
    # anywhere on this protein, each refitted independently.
    surface = exposed_residues(coords, plddt_floor=plddt_floor)
    if len(surface) < len(signal_points):
        surface = coords[["x", "y", "z"]].to_numpy()

    if len(surface) >= len(signal_points) and n_permutations > 0:
        null_scores = []
        for _ in range(n_permutations):
            index = rng.choice(len(surface), size=len(signal_points), replace=False)
            _, null_value = fit_center(
                surface[index], scaffold_size, rigidity, max_seeds=3
            )
            null_scores.append(null_value)
        null_scores = np.asarray(null_scores)
        score.null_mean = float(null_scores.mean())
        score.null_std = float(null_scores.std())
        if score.null_std > 0:
            # Lower penalty is better, so a positive z means tighter than chance.
            score.z_score = float((score.null_mean - observed) / score.null_std)
        else:
            score.z_score = 0.0
        score.p_value = float((np.sum(null_scores <= observed) + 1) / (len(null_scores) + 1))
    else:
        score.note = "structure too small for a permutation null"

    if len(control) > 0:
        control_points, _ = site_coordinates(control, coords, plddt_floor=plddt_floor)
        score.n_control_sites = len(control_points)
        if len(control_points) > 0:
            control_distances = np.linalg.norm(control_points - center, axis=1)
            control_penalties = penalty(control_distances, scaffold_size, rigidity)
            score.control_mean_penalty = float(np.mean(control_penalties))
            score.control_fraction_reachable = float(
                np.mean(control_penalties <= ON_SHELL_TOLERANCE)
            )
            # How much more of the signal lands in reach than the control does.
            score.signal_vs_control = float(
                score.fraction_reachable / (score.control_fraction_reachable + 0.01)
            )

    return score


def score_to_frame(scores: list[ProteinScore]) -> pd.DataFrame:
    """Collect scores into a table ranked by how convincing the pocket is."""
    if not scores:
        return pd.DataFrame()
    frame = pd.DataFrame([s.as_row() for s in scores])
    return frame.sort_values(
        ["z_score", "fraction_reachable"], ascending=[False, False]
    ).reset_index(drop=True)
