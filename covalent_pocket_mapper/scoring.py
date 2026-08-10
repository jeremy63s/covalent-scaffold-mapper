"""Locate the pocket a scaffold was anchored to, and test it against chance.

Each capture site is one end of the scaffold pinned to a known atom. The reach
model says where the other end can be, so every site casts a cloud of
probability into the space around it. Where a protein really does have one
binding pocket, all those clouds overlap on it.

The clouds are combined by multiplying -- in log space, since the numbers get
small fast. Multiplying is what makes the reasoning strict rather than
suggestive: the reach density is exactly zero past the scaffold's contour
length, so a single site that cannot reach a location rules that location out
however well the rest agree. No averaging can talk it back in.

The answer is read off where that combined cloud is highest *on the protein
surface*, since that is where a covalent partner has to be.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .reach import ReachModel

DEFAULT_PLDDT_FLOOR = 50.0
DEFAULT_PERMUTATIONS = 200
MIN_SITES = 4

# Grid spacing for the search over candidate pocket locations, in angstroms.
DEFAULT_GRID_SPACING = 1.5
# A pocket sits in contact with the protein but not inside it: near enough to a
# backbone atom to touch, far enough not to overlap one.
CONTACT_MIN = 3.0
CONTACT_MAX = 8.0
# Buried points have backbone all around them; exposed ones have it on one side.
# This is the minimum directional asymmetry a candidate must show. Kept low on
# purpose: binding pockets are concave, so demanding a strongly convex surface
# throws away exactly the locations worth finding. Raising it past ~0.2 starts
# losing real pockets; lowering it below ~0.1 only adds interior voids.
EXPOSURE_MIN = 0.15
EXPOSURE_RADIUS = 12.0
MAX_CANDIDATES = 80_000

# Burial proxy for choosing which residues the null may draw its sites from.
BURIAL_RADIUS = 12.0
BURIAL_PERCENTILE = 60.0
# Share of the posterior a credible region has to contain.
CREDIBLE_MASS = 0.5


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


def surface_candidates(
    coords: pd.DataFrame,
    *,
    grid_spacing: float = DEFAULT_GRID_SPACING,
    plddt_floor: float = DEFAULT_PLDDT_FLOOR,
    contact_min: float = CONTACT_MIN,
    contact_max: float = CONTACT_MAX,
    exposure_min: float = EXPOSURE_MIN,
    max_candidates: int = MAX_CANDIDATES,
) -> np.ndarray:
    """Grid points lying in the solvent shell just outside the protein.

    Two filters. The first keeps points in contact range of the backbone -- not
    overlapping an atom, not out in bulk solvent. The second throws out interior
    cavities, which pass the first test perfectly well: a buried point has
    backbone in every direction, so the mean direction to its neighbours nearly
    cancels, while a point on the outside sees backbone only on one side.
    """
    atoms = coords[coords["plddt"] >= plddt_floor][["x", "y", "z"]].to_numpy()
    if len(atoms) < 4:
        atoms = coords[["x", "y", "z"]].to_numpy()
    if len(atoms) < 4:
        return np.empty((0, 3))

    low = atoms.min(axis=0) - contact_max
    high = atoms.max(axis=0) + contact_max
    axes = [
        np.arange(low[i], high[i] + grid_spacing, grid_spacing) for i in range(3)
    ]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)

    tree = cKDTree(atoms)
    nearest, _ = tree.query(grid, k=1)
    shell = grid[(nearest >= contact_min) & (nearest <= contact_max)]
    if len(shell) == 0:
        return shell

    keep = np.zeros(len(shell), dtype=bool)
    neighbourhoods = tree.query_ball_point(shell, EXPOSURE_RADIUS)
    for i, idx in enumerate(neighbourhoods):
        if not idx:
            continue
        vectors = atoms[idx] - shell[i]
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        keep[i] = np.linalg.norm((vectors / norms).mean(axis=0)) >= exposure_min

    candidates = shell[keep]
    if len(candidates) > max_candidates:
        step = int(np.ceil(len(candidates) / max_candidates))
        candidates = candidates[::step]
    return candidates


def combined_log_density(
    candidates: np.ndarray, sites: np.ndarray, model: ReachModel
) -> np.ndarray:
    """Log of the product of every site's reach cloud, at each candidate point.

    -inf marks a location at least one site cannot reach, which is the zero that
    the multiplication is there to propagate.
    """
    if len(candidates) == 0 or len(sites) == 0:
        return np.full(len(candidates), -np.inf)

    total = np.zeros(len(candidates))
    for site in sites:
        distances = np.linalg.norm(candidates - site, axis=1)
        total += model.log_density_at(distances)
    return total


def _credible_region(candidates: np.ndarray, log_density: np.ndarray) -> float:
    """Width of the smallest set of points holding CREDIBLE_MASS of the posterior."""
    finite = np.isfinite(log_density)
    if finite.sum() < 2:
        return 0.0
    values = log_density[finite]
    points = candidates[finite]

    weights = np.exp(values - values.max())
    weights /= weights.sum()
    order = np.argsort(weights)[::-1]
    cumulative = np.cumsum(weights[order])
    cutoff = int(np.searchsorted(cumulative, CREDIBLE_MASS)) + 1
    chosen = points[order[:cutoff]]
    if len(chosen) < 2:
        return 0.0
    return float(np.linalg.norm(chosen.max(axis=0) - chosen.min(axis=0)))


@dataclass
class ProteinScore:
    """Everything the pipeline learned about one protein."""

    accession: str
    gene: str | None = None
    n_signal_sites: int = 0
    n_control_sites: int = 0
    binding_site: tuple[float, float, float] | None = None
    scaffold: str = ""
    max_reach: float = float("nan")
    log_likelihood: float = float("-inf")
    site_probability: float = float("nan")
    site_posterior: float = float("nan")
    credible_region: float = float("nan")
    mean_site_distance: float = float("nan")
    null_mean: float = float("nan")
    null_std: float = float("nan")
    null_feasible_fraction: float = float("nan")
    z_score: float = float("nan")
    p_value: float = float("nan")
    control_reachable: int | None = None
    control_fraction_reachable: float | None = None
    note: str = ""

    def as_row(self) -> dict:
        site = self.binding_site
        return {
            "accession": self.accession,
            "gene": self.gene,
            "n_signal_sites": self.n_signal_sites,
            "n_control_sites": self.n_control_sites,
            "binding_site_x": site[0] if site else None,
            "binding_site_y": site[1] if site else None,
            "binding_site_z": site[2] if site else None,
            "site_probability": self.site_probability,
            "site_posterior": self.site_posterior,
            "credible_region_A": self.credible_region,
            "log_likelihood": self.log_likelihood,
            "mean_site_distance_A": self.mean_site_distance,
            "z_score": self.z_score,
            "p_value": self.p_value,
            "null_feasible_fraction": self.null_feasible_fraction,
            "control_sites_reachable": self.control_reachable,
            "control_fraction_reachable": self.control_fraction_reachable,
            "scaffold": self.scaffold,
            "max_reach_A": self.max_reach,
            "note": self.note,
        }


def score_protein(
    accession: str,
    peptide_table: pd.DataFrame,
    coords: pd.DataFrame,
    model: ReachModel,
    *,
    plddt_floor: float = DEFAULT_PLDDT_FLOOR,
    grid_spacing: float = DEFAULT_GRID_SPACING,
    n_permutations: int = DEFAULT_PERMUTATIONS,
    min_sites: int = MIN_SITES,
    gene: str | None = None,
    rng: np.random.Generator | None = None,
    candidates: np.ndarray | None = None,
) -> ProteinScore | None:
    """Find the most likely pocket for one protein and test it against chance.

    The null asks: had the same number of capture sites landed anywhere on this
    protein's surface, how well would the best pocket have explained them? Every
    draw is scored the same way, including the possibility -- common, and
    informative -- that no single pocket can reach a random scatter at all.
    """
    rng = rng or np.random.default_rng(0)
    protein_peptides = peptide_table[peptide_table["accession"] == accession]
    signal = protein_peptides[protein_peptides["is_signal"]]
    control = protein_peptides[protein_peptides["is_background"]]

    signal_points, _ = site_coordinates(signal, coords, plddt_floor=plddt_floor)
    if len(signal_points) < min_sites:
        return None

    if candidates is None:
        candidates = surface_candidates(
            coords, grid_spacing=grid_spacing, plddt_floor=plddt_floor
        )
    if len(candidates) == 0:
        return None

    field = combined_log_density(candidates, signal_points, model)
    best_index = int(np.argmax(field))
    observed = float(field[best_index])

    score = ProteinScore(
        accession=accession,
        gene=gene,
        n_signal_sites=len(signal_points),
        scaffold=model.label,
        max_reach=model.max_reach,
        log_likelihood=observed,
    )

    if not np.isfinite(observed):
        # Every candidate location is out of reach of at least one site, so no
        # single anchor point explains this protein's capture pattern.
        score.note = "no surface point is within reach of all capture sites"
        score.p_value = 1.0
        return score

    best = candidates[best_index]
    score.binding_site = tuple(float(v) for v in best)
    score.credible_region = _credible_region(candidates, field)
    score.mean_site_distance = float(
        np.mean(np.linalg.norm(signal_points - best, axis=1))
    )
    # Per-site geometric mean of the reach density, so the number stays
    # comparable between proteins carrying different numbers of sites.
    score.site_probability = float(np.exp(observed / len(signal_points)))

    finite = np.isfinite(field)
    weights = np.exp(field[finite] - observed)
    score.site_posterior = float(1.0 / weights.sum())

    surface = exposed_residues(coords, plddt_floor=plddt_floor)
    if len(surface) >= len(signal_points) and n_permutations > 0:
        null_scores = np.empty(n_permutations)
        for i in range(n_permutations):
            index = rng.choice(len(surface), size=len(signal_points), replace=False)
            null_field = combined_log_density(candidates, surface[index], model)
            null_scores[i] = null_field.max()

        feasible = np.isfinite(null_scores)
        score.null_feasible_fraction = float(feasible.mean())
        if feasible.sum() >= 2:
            score.null_mean = float(null_scores[feasible].mean())
            score.null_std = float(null_scores[feasible].std())
            if score.null_std > 0:
                score.z_score = float((observed - score.null_mean) / score.null_std)
            else:
                score.z_score = 0.0
        # -inf never beats a finite observation, so unreachable draws count as
        # losses rather than being dropped.
        score.p_value = float(
            (np.sum(null_scores >= observed) + 1) / (len(null_scores) + 1)
        )
    else:
        score.note = "structure too small for a permutation null"

    if len(control) > 0:
        control_points, _ = site_coordinates(
            control, coords, plddt_floor=plddt_floor
        )
        score.n_control_sites = len(control_points)
        if len(control_points) > 0:
            control_distances = np.linalg.norm(control_points - best, axis=1)
            reachable = int(np.sum(control_distances <= model.max_reach))
            score.control_reachable = reachable
            score.control_fraction_reachable = float(reachable / len(control_points))

    return score


def score_to_frame(scores: list[ProteinScore]) -> pd.DataFrame:
    """Collect scores into a table, most convincing pocket first."""
    if not scores:
        return pd.DataFrame()
    frame = pd.DataFrame([s.as_row() for s in scores])
    return frame.sort_values(
        ["p_value", "z_score"], ascending=[True, False]
    ).reset_index(drop=True)
