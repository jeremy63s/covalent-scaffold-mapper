"""Turn a scaffold's chemistry into a map of where its far end can be.

`polymer_pr` gives P(r), the chance the two ends of the chain are a distance r
apart. That is a *radial* distribution: it already contains the 4*pi*r^2 factor
for the surface area of the shell at radius r. To ask the question this tool
actually needs -- "if one end is bonded here, how likely is the other end to be
at that particular point in space?" -- the shell factor has to come back out:

    density(r) = P(r) / (4 * pi * r^2)          units: 1/A^3

That density, evaluated at every point around a capture site, is the reach
cloud for one anchor. Beyond the chain's contour length it is exactly zero, and
zero is the whole point: a site the scaffold cannot reach rules a location out
no matter how well the other sites agree with it.

Everything here works in angstroms. `polymer_pr` works in nanometres, and the
conversion happens once, on the way in.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .polymer_pr import Monomer, _grow_endpoints, monomer_from_spec, wlc_pr

NM_TO_A = 10.0
DEFAULT_CHAINS = 120_000
DEFAULT_BINS = 120


@dataclass
class ReachModel:
    """Where the free end of a scaffold can be, relative to its anchored end."""

    radii: np.ndarray        # bin centres, angstroms
    log_density: np.ndarray  # log of the 3D density at that radius, log(1/A^3)
    max_reach: float         # contour length, angstroms; beyond this, zero
    label: str = "scaffold"
    contour_length: float = float("nan")
    most_likely_r: float = float("nan")
    mean_r: float = float("nan")

    def log_density_at(self, distances: np.ndarray) -> np.ndarray:
        """Log density at each distance; -inf where the scaffold cannot reach."""
        distances = np.asarray(distances, dtype=float)
        out = np.interp(
            distances,
            self.radii,
            self.log_density,
            left=self.log_density[0],
            right=-np.inf,
        )
        return np.where(distances > self.max_reach, -np.inf, out)

    def describe(self) -> str:
        return (
            f"{self.label}: contour {self.contour_length:.1f} A, "
            f"most likely end-to-end {self.most_likely_r:.1f} A, "
            f"mean {self.mean_r:.1f} A"
        )


def _density_from_radial(
    r_a: np.ndarray, p_r_a: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Strip the 4*pi*r^2 shell factor off a radial distribution.

    Bins the sampler never reached carry no information rather than a true zero,
    so they are dropped and the log density is interpolated across them.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        density = p_r_a / (4.0 * np.pi * r_a ** 2)
    usable = np.isfinite(density) & (density > 0) & (r_a > 0)
    if usable.sum() < 2:
        raise ValueError(
            "the sampled end-to-end distribution is too sparse to build a reach "
            "model; raise the chain count"
        )
    return r_a[usable], np.log(density[usable])


def _density_from_samples(
    r_a: np.ndarray, bins: int
) -> tuple[np.ndarray, np.ndarray]:
    """Estimate the 3D density directly from sampled end-to-end distances.

    Equal-count bins rather than equal-width ones. A uniform histogram starves
    its small-r bins -- there is very little volume down there -- and dividing
    those few counts by r^2 turns a handful of samples into a wild density. With
    equal occupancy every bin carries the same relative error, and each one is
    divided by the exact volume of the shell it spans.
    """
    r_sorted = np.sort(r_a)
    total = len(r_sorted)
    edges = np.quantile(r_sorted, np.linspace(0.0, 1.0, bins + 1))
    edges = np.unique(edges)
    if len(edges) < 3:
        raise ValueError("end-to-end distances are too degenerate to bin")

    counts, _ = np.histogram(r_sorted, bins=edges)
    shell_volume = (4.0 * np.pi / 3.0) * (edges[1:] ** 3 - edges[:-1] ** 3)

    with np.errstate(divide="ignore", invalid="ignore"):
        density = (counts / total) / shell_volume
    centers = 0.5 * (edges[:-1] + edges[1:])

    usable = np.isfinite(density) & (density > 0)
    if usable.sum() < 2:
        raise ValueError(
            "the sampled end-to-end distribution is too sparse to build a reach "
            "model; raise the chain count"
        )
    return centers[usable], _smooth(np.log(density[usable]))


def _smooth(values: np.ndarray, window: int = 5) -> np.ndarray:
    """Damp bin-to-bin sampling noise without moving the curve.

    The underlying density is smooth in r, so residual jitter between adjacent
    bins is Monte Carlo error. Edges are handled by shrinking the window rather
    than padding, which would drag the ends toward zero.
    """
    if window < 3 or len(values) < window:
        return values
    half = window // 2
    out = np.empty_like(values)
    for i in range(len(values)):
        lo, hi = max(0, i - half), min(len(values), i + half + 1)
        out[i] = values[lo:hi].mean()
    return out


def _extend_to_contour(
    radii: np.ndarray, log_density: np.ndarray, contour: float, points: int = 16
) -> tuple[np.ndarray, np.ndarray]:
    """Carry the density out to the chain's true limit.

    Sampling never produces a fully stretched chain: that configuration is real
    but vanishingly rare, so the longest distance the simulation happens to
    reach is not the longest the scaffold can span. Cutting off there would rule
    out locations that are merely improbable, and the difference matters -- the
    gap between the two is several angstroms even with a hundred thousand
    chains. The tail is extrapolated in log space out to the contour length,
    which is the only hard limit there is.
    """
    if contour <= radii[-1] or len(radii) < 4:
        return radii, log_density

    tail = max(2, len(radii) // 10)
    slope = float(np.polyfit(radii[-tail:], log_density[-tail:], 1)[0])
    slope = min(slope, -1e-3)  # the tail can only fall away

    extra_radii = np.linspace(radii[-1], contour, points + 1)[1:]
    extra = log_density[-1] + slope * (extra_radii - radii[-1])
    return np.concatenate([radii, extra_radii]), np.concatenate([log_density, extra])


def from_monomer(
    monomer: Monomer,
    n_units: int,
    *,
    chains: int = DEFAULT_CHAINS,
    bins: int = DEFAULT_BINS,
    seed: int = 0,
    label: str | None = None,
) -> ReachModel:
    """Build a reach model by simulating the chain, bond by bond."""
    if n_units < 1:
        raise ValueError("a polymer needs at least one repeat unit")

    rng = np.random.default_rng(seed)
    ends_nm = _grow_endpoints(monomer, n_units, chains, rng)
    r_a = np.linalg.norm(ends_nm, axis=1) * NM_TO_A

    radii, log_density = _density_from_samples(r_a, bins)
    contour = monomer.contour_per_monomer() * n_units * NM_TO_A
    radii, log_density = _extend_to_contour(radii, log_density, contour)

    return ReachModel(
        radii=radii,
        log_density=log_density,
        # The contour length is the one true limit: the chain physically cannot
        # span further, and everything shorter is only a matter of probability.
        max_reach=float(contour),
        label=label or f"{monomer.name} x{n_units}",
        contour_length=float(contour),
        # Peak of the radial distribution -- the distance the two ends are most
        # often found apart, not the peak of the 3D density (which for a coil
        # always sits at the origin).
        most_likely_r=float(radii[int(np.argmax(np.exp(log_density) * radii ** 2))]),
        mean_r=float(r_a.mean()),
    )


def from_spec(spec: dict, n_units: int, **kwargs) -> ReachModel:
    """Build a reach model from the hand-entered per-bond form."""
    monomer = monomer_from_spec(spec)
    return from_monomer(monomer, n_units, label=spec.get("name"), **kwargs)


def from_json(path: str | Path, n_units: int | None = None, **kwargs) -> ReachModel:
    """Load the per-bond form from a JSON file.

    The file may carry `n_units` itself, so a saved polymer is one argument.
    """
    import json

    spec = json.loads(Path(path).read_text())
    units = n_units if n_units is not None else spec.get("n_units")
    if units is None:
        raise ValueError(
            f"{path} does not set 'n_units' and none was given on the command line"
        )
    return from_spec(spec, int(units), **kwargs)


def from_persistence_length(
    contour_length: float,
    persistence_length: float,
    *,
    bins: int = 2000,
    label: str | None = None,
) -> ReachModel:
    """Build a reach model from just contour length and persistence length.

    The shortcut for when rotamer energies are not available: both numbers are
    tabulated for essentially every common polymer. Slightly less faithful than
    simulating the bonds for very short chains, but it costs one lookup.
    """
    if contour_length <= 0 or persistence_length <= 0:
        raise ValueError("contour length and persistence length must both be positive")

    r_a, p_r_a = wlc_pr(contour_length, persistence_length, bins=bins)
    radii, log_density = _density_from_radial(r_a, p_r_a)

    peak = float(r_a[int(np.argmax(p_r_a))])
    mean = float(np.trapezoid(p_r_a * r_a, r_a)) if hasattr(np, "trapezoid") else float(
        np.trapz(p_r_a * r_a, r_a)
    )
    return ReachModel(
        radii=radii,
        log_density=log_density,
        max_reach=float(contour_length),
        label=label or f"worm-like chain (L={contour_length:.0f} A, lp={persistence_length:.0f} A)",
        contour_length=float(contour_length),
        most_likely_r=peak,
        mean_r=mean,
    )
