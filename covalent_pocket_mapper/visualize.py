"""Interactive 3D views of the inferred pocket, written as standalone HTML files.

Each file opens in a browser with no server and no internet, so results can be
handed to a collaborator as-is.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .reach import ReachModel
from .scoring import (
    CREDIBLE_MASS,
    ProteinScore,
    combined_log_density,
    site_coordinates,
    surface_candidates,
)

SITE_COLOR = "#2eb872"
CONTROL_COLOR = "#d64545"
POCKET_COLOR = "#1f4e79"


def plot_protein(
    score: ProteinScore,
    peptide_table: pd.DataFrame,
    coords: pd.DataFrame,
    model: ReachModel,
    output_dir: Path,
    *,
    plddt_floor: float = 50.0,
    grid_spacing: float = 1.5,
) -> Path | None:
    """Write one protein's pocket view; returns the file path."""
    import plotly.graph_objects as go

    if score.binding_site is None:
        return None

    best = np.array(score.binding_site)
    label = score.gene or score.accession

    protein_peptides = peptide_table[peptide_table["accession"] == score.accession]
    signal = protein_peptides[protein_peptides["is_signal"]]
    control = protein_peptides[protein_peptides["is_background"]]

    signal_points, signal_meta = site_coordinates(
        signal, coords, plddt_floor=plddt_floor
    )
    control_points, control_meta = site_coordinates(
        control, coords, plddt_floor=plddt_floor
    )

    traces = []

    backbone = coords.sort_values("residue_num")
    traces.append(
        go.Scatter3d(
            x=backbone["x"], y=backbone["y"], z=backbone["z"],
            mode="lines",
            line=dict(color="#8a8a8a", width=3),
            opacity=0.3,
            name="Protein backbone",
            hoverinfo="skip",
        )
    )

    # The credible region: surface points holding most of the posterior. This is
    # the honest answer to "where is the pocket" -- a cloud, not a dot.
    candidates = surface_candidates(
        coords, grid_spacing=grid_spacing, plddt_floor=plddt_floor
    )
    field = combined_log_density(candidates, signal_points, model)
    finite = np.isfinite(field)
    if finite.any():
        values = field[finite]
        points = candidates[finite]
        weights = np.exp(values - values.max())
        weights /= weights.sum()
        order = np.argsort(weights)[::-1]
        cutoff = int(np.searchsorted(np.cumsum(weights[order]), CREDIBLE_MASS)) + 1
        region = points[order[:cutoff]]
        traces.append(
            go.Scatter3d(
                x=region[:, 0], y=region[:, 1], z=region[:, 2],
                mode="markers",
                marker=dict(size=2.5, color=values[order[:cutoff]],
                            colorscale="Plasma", opacity=0.55,
                            colorbar=dict(title="log likelihood", x=1.02)),
                name=f"Most likely {CREDIBLE_MASS:.0%} of pocket locations",
                hoverinfo="skip",
            )
        )

    if len(signal_points):
        distances = np.linalg.norm(signal_points - best, axis=1)
        sizes = [
            6 + 3 * min(int(m.get("treatment_psms") or 1), 5) for m in signal_meta
        ]
        text = [
            f"{m['peptide_sequence']}<br>residues {m['start']}-{m['stop']}"
            f"<br>spectra: {int(m.get('treatment_psms') or 0)}"
            f"<br>distance to pocket: {d:.1f} A"
            for m, d in zip(signal_meta, distances)
        ]
        traces.append(
            go.Scatter3d(
                x=signal_points[:, 0], y=signal_points[:, 1], z=signal_points[:, 2],
                mode="markers",
                marker=dict(size=sizes, color=SITE_COLOR,
                            line=dict(color="#1a1a1a", width=1)),
                name="Capture sites (treatment)",
                customdata=text,
                hovertemplate="%{customdata}<extra></extra>",
            )
        )
        # One line per site: the scaffold, drawn as the tether it represents.
        for point in signal_points:
            traces.append(
                go.Scatter3d(
                    x=[point[0], best[0]], y=[point[1], best[1]],
                    z=[point[2], best[2]],
                    mode="lines",
                    line=dict(color=SITE_COLOR, width=2, dash="dot"),
                    opacity=0.5, showlegend=False, hoverinfo="skip",
                )
            )

    if len(control_points):
        control_distances = np.linalg.norm(control_points - best, axis=1)
        control_text = [
            f"{m['peptide_sequence']}<br>residues {m['start']}-{m['stop']}"
            f"<br>control spectra: {int(m.get('control_psms') or 0)}"
            f"<br>distance to pocket: {d:.1f} A"
            f"<br>{'within' if d <= model.max_reach else 'beyond'} scaffold reach"
            for m, d in zip(control_meta, control_distances)
        ]
        traces.append(
            go.Scatter3d(
                x=control_points[:, 0], y=control_points[:, 1], z=control_points[:, 2],
                mode="markers",
                marker=dict(size=6, color=CONTROL_COLOR,
                            line=dict(color="#5c1f1f", width=1)),
                name="Control sites",
                customdata=control_text,
                hovertemplate="%{customdata}<extra></extra>",
            )
        )

    traces.append(
        go.Scatter3d(
            x=[best[0]], y=[best[1]], z=[best[2]],
            mode="markers",
            marker=dict(size=10, color=POCKET_COLOR, symbol="cross"),
            name="Most likely pocket",
            hovertemplate=(
                f"({best[0]:.1f}, {best[1]:.1f}, {best[2]:.1f})"
                f"<br>posterior {score.site_posterior:.2%}"
                f"<extra></extra>"
            ),
        )
    )

    subtitle = (
        f"{model.label} | {score.n_signal_sites} sites | "
        f"credible region {score.credible_region:.0f} A | "
        f"p = {score.p_value:.3f}"
    )
    if score.control_fraction_reachable is not None:
        subtitle += f" | control in reach {score.control_fraction_reachable:.0%}"

    figure = go.Figure(data=traces)
    figure.update_layout(
        title=f"{label} ({score.accession})<br><sub>{subtitle}</sub>",
        scene=dict(
            xaxis_title="X (A)", yaxis_title="Y (A)", zaxis_title="Z (A)",
            aspectmode="data",
        ),
        width=1000, height=750,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{label}_{score.accession}_pocket.html"
    figure.write_html(str(out_path), include_plotlyjs="inline")
    return out_path
