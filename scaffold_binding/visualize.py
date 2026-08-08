"""Interactive 3D views of a fitted pocket, written as standalone HTML files.

Each file opens in a browser with no server and no internet, so results can be
handed to a collaborator as-is.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .scoring import ProteinScore, penalty, site_coordinates

# Green means the scaffold could reach this site from the fitted pocket; amber
# means it could not, and the number is how many angstroms off it was.
REACHABLE_COLOR = "#2eb872"
UNREACHABLE_COLOR = "#f0a202"
CONTROL_COLOR = "#d64545"
SHELL_COLOR = "#4c9bd5"


def _shell_surface(center: np.ndarray, radius: float, resolution: int = 40):
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radius * np.cos(v)[np.newaxis, :]
    return x, y, z


def plot_protein(
    score: ProteinScore,
    peptide_table: pd.DataFrame,
    coords: pd.DataFrame,
    output_dir: Path,
    *,
    plddt_floor: float = 50.0,
) -> Path | None:
    """Write one protein's pocket view; returns the file path."""
    import plotly.graph_objects as go

    if score.center is None:
        return None

    center = np.array(score.center)
    radius = score.scaffold_size
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
            x=backbone["x"],
            y=backbone["y"],
            z=backbone["z"],
            mode="lines",
            line=dict(color="#8a8a8a", width=3),
            opacity=0.35,
            name="Protein backbone",
            hoverinfo="skip",
        )
    )

    if len(signal_points):
        distances = np.linalg.norm(signal_points - center, axis=1)
        penalties = penalty(distances, radius, score.rigidity)
        colors = [
            REACHABLE_COLOR if p <= 0 else UNREACHABLE_COLOR for p in penalties
        ]
        sizes = [
            6 + 3 * min(int(m.get("treatment_psms") or 1), 5) for m in signal_meta
        ]
        text = [
            f"{m['peptide_sequence']}<br>residues {m['start']}-{m['stop']}"
            f"<br>spectra: {int(m.get('treatment_psms') or 0)}"
            f"<br>distance to pocket: {d:.1f} A"
            f"<br>penalty: {p:.1f} A"
            for m, d, p in zip(signal_meta, distances, penalties)
        ]
        traces.append(
            go.Scatter3d(
                x=signal_points[:, 0],
                y=signal_points[:, 1],
                z=signal_points[:, 2],
                mode="markers",
                marker=dict(
                    size=sizes, color=colors, line=dict(color="#1a1a1a", width=1)
                ),
                name="Capture sites (treatment)",
                customdata=text,
                hovertemplate="%{customdata}<extra></extra>",
            )
        )

    if len(control_points):
        control_distances = np.linalg.norm(control_points - center, axis=1)
        control_text = [
            f"{m['peptide_sequence']}<br>residues {m['start']}-{m['stop']}"
            f"<br>control spectra: {int(m.get('control_psms') or 0)}"
            f"<br>distance to pocket: {d:.1f} A"
            for m, d in zip(control_meta, control_distances)
        ]
        traces.append(
            go.Scatter3d(
                x=control_points[:, 0],
                y=control_points[:, 1],
                z=control_points[:, 2],
                mode="markers",
                marker=dict(
                    size=7, color=CONTROL_COLOR, line=dict(color="#5c1f1f", width=1)
                ),
                name="Control sites",
                customdata=control_text,
                hovertemplate="%{customdata}<extra></extra>",
            )
        )

    x, y, z = _shell_surface(center, radius)
    traces.append(
        go.Surface(
            x=x,
            y=y,
            z=z,
            opacity=0.12,
            colorscale=[[0, SHELL_COLOR], [1, SHELL_COLOR]],
            showscale=False,
            name=f"Scaffold reach ({radius:.0f} A)",
            hoverinfo="skip",
        )
    )

    # With a flexible scaffold the reachable zone is a filled ball rather than a
    # shell; draw the inner limit so the free band is visible.
    if 0 < score.rigidity < 1:
        ix, iy, iz = _shell_surface(center, radius * score.rigidity)
        traces.append(
            go.Surface(
                x=ix,
                y=iy,
                z=iz,
                opacity=0.08,
                colorscale=[[0, "#9b59b6"], [1, "#9b59b6"]],
                showscale=False,
                name=f"Inner limit ({radius * score.rigidity:.0f} A)",
                hoverinfo="skip",
            )
        )

    traces.append(
        go.Scatter3d(
            x=[center[0]],
            y=[center[1]],
            z=[center[2]],
            mode="markers",
            marker=dict(size=10, color="#1f4e79", symbol="cross"),
            name="Inferred pocket centre",
            hovertemplate=(
                f"({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})<extra></extra>"
            ),
        )
    )

    subtitle = (
        f"scaffold {radius:.0f} A, rigidity {score.rigidity:g} | "
        f"{score.n_signal_sites} sites, {score.fraction_reachable:.0%} in reach | "
        f"z = {score.z_score:.2f}, p = {score.p_value:.3f}"
    )
    if score.signal_vs_control is not None:
        subtitle += f" | signal/control = {score.signal_vs_control:.1f}"

    figure = go.Figure(data=traces)
    figure.update_layout(
        title=f"{label} ({score.accession})<br><sub>{subtitle}</sub>",
        scene=dict(
            xaxis_title="X (A)", yaxis_title="Y (A)", zaxis_title="Z (A)",
            aspectmode="data",
        ),
        width=1000,
        height=750,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{label}_{score.accession}_pocket.html"
    figure.write_html(str(out_path), include_plotlyjs="inline")
    return out_path
