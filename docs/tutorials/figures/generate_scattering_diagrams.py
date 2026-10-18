"""Rebuild the vector diagrams for scattering-from-phases.qmd.

Run with the documentation Python environment from any directory. Geometry is
illustrative: only the six marked points scatter, not the surrounding outline.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch
import numpy as np

POINTS = np.array([[-2., -.65], [-1.2, 1.], [-.4, -1.25],
                   [.4, 1.25], [1.2, -1.], [2., .65]])
CENTER = np.array([0., 0.])
EDGE = np.array([-3., 0.])
INK = "#172f46"
ORANGE = "#b95017"
BLUE = "#176b96"
GOLD = "#f3bb45"
OUTPUT = Path(__file__).resolve().parent
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "svg.fonttype": "none", "svg.hashsalt": "pyseb-scattering"})


def arrow(ax, start, end, color=ORANGE, width=1.8, alpha=1.):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
        mutation_scale=13, linewidth=width, color=color, alpha=alpha,
        shrinkA=5, shrinkB=5, zorder=2))


def geometry(ax, title):
    ax.add_patch(Ellipse((0, 0), 6, 3.6, facecolor="#edf4f8",
                         edgecolor=INK, linewidth=1.7))
    ax.scatter(*POINTS.T, s=40, color=INK, zorder=4)
    ax.set(xlim=(-3.65, 3.65), ylim=(-2.55, 2.55), aspect="equal")
    ax.set_title(title, loc="left", color=INK, fontsize=13, weight="bold", pad=6)
    ax.axis("off")


def reference(ax, position, label, offset):
    ax.scatter(*position, s=95, color=GOLD, edgecolor=INK, linewidth=.8, zorder=5)
    ax.annotate(label, position, xytext=offset, textcoords="offset points",
                color=INK, weight="bold", zorder=6)


def main(preview_directory=None):
    if preview_directory is not None:
        preview_directory.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12.5, 4.2), facecolor="white")
    grid = fig.add_gridspec(1, 3, width_ratios=(1, .82, 1))
    midpoint_ax = fig.add_subplot(grid[0, 0])
    edge_ax = fig.add_subplot(grid[0, 2])
    change_ax = fig.add_subplot(grid[0, 1])
    for ax, title in (
            (midpoint_ax, "Midpoint reference"),
            (edge_ax, "Edge reference")):
        geometry(ax, title)
    for ax, position, label, offset in (
            (midpoint_ax, CENTER, r"$R_m$", (8, 5)),
            (edge_ax, EDGE, r"$R_e$", (-6, 12))):
        for point in POINTS:
            arrow(ax, position, point)
        reference(ax, position, label, offset)
    midpoint_ax.text(0, -2.25, r"$a_m(\mathbf{q})$: midpoint → scatterers",
                     ha="center", color=INK)
    edge_ax.text(0, -2.25, r"$a_e(\mathbf{q})$: edge → scatterers",
                 ha="center", color=INK)
    change_ax.axis("off")
    change_ax.set(xlim=(0, 1), ylim=(0, 1))
    change_ax.add_patch(FancyArrowPatch((.02, .5), (.98, .5),
        arrowstyle="-|>", mutation_scale=22, linewidth=2.5, color=BLUE))
    change_ax.text(.5, .62,
        r"$\times\; e^{i\mathbf{q}\cdot(\mathbf{R}_m-\mathbf{R}_e)}$",
        ha="center", va="center", color=BLUE, fontsize=14)
    change_ax.text(.5, .38, "Move the reference\n(before averaging)",
                   ha="center", va="center", color=INK, fontsize=10)
    fig.subplots_adjust(left=.04, right=.98, top=.90, bottom=.12, wspace=.05)
    fig.text(.5, .025,
        "The same six scatterers. Moving the reference multiplies each unaveraged phase contribution.",
        ha="center", color=INK, fontsize=10)
    fig.savefig(OUTPUT / "scattering-distance-vectors.svg", facecolor="white",
                metadata={"Date": None})
    if preview_directory is not None:
        fig.savefig(preview_directory / "scattering-distance-vectors.png",
                    dpi=150, facecolor="white")
    plt.close(fig)

    # Pair distances belong with the definitions, separate from the reference shift.
    fig, (pair_ax, phase_ax) = plt.subplots(1, 2, figsize=(10, 4.2), facecolor="white")
    geometry(pair_ax, "Form factor: all scatterer pairs")
    geometry(phase_ax, "Phase factor: two references")
    for i, j in ((0, 3), (1, 4), (2, 5)):
        arrow(pair_ax, POINTS[i], POINTS[j], color=BLUE)
    pair_ax.text(0, -2.25,
        r"$F$: pairs $\mathbf{r}_i-\mathbf{r}_j$ (a few shown)",
        ha="center", color=INK)
    reference(phase_ax, CENTER, r"$R_m$", (8, 6))
    reference(phase_ax, EDGE, r"$R_e$", (-6, 12))
    arrow(phase_ax, EDGE, CENTER, color=BLUE, width=2.5)
    phase_ax.text(-1.5, -.48, r"$\mathbf{R}_m-\mathbf{R}_e$",
                  ha="center", color=BLUE)
    phase_ax.text(0, -2.25, r"$\Psi_{m,e}$: reference → reference",
                  ha="center", color=INK)
    fig.subplots_adjust(left=.04, right=.98, top=.90, bottom=.10, wspace=.15)
    fig.savefig(OUTPUT / "scattering-pair-vectors.svg", facecolor="white",
                metadata={"Date": None})
    if preview_directory is not None:
        fig.savefig(preview_directory / "scattering-pair-vectors.png",
                    dpi=150, facecolor="white")
    plt.close(fig)

    # A reference shift rotates every complex contribution by the same phase.
    q = np.array([.9, .5])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), facecolor="white")
    for ax, ref, label in zip(axes, (CENTER, EDGE), ("Midpoint reference", "Edge reference")):
        terms = np.exp(1j*((POINTS-ref)@q)) / len(POINTS)
        total = terms.sum()
        for z in terms:
            arrow(ax, (0, 0), (z.real, z.imag), width=1.6, alpha=.7)
        arrow(ax, (0, 0), (total.real, total.imag), color=BLUE, width=3)
        ax.add_patch(plt.Circle((0, 0), abs(total), fill=False, color=BLUE,
                                linestyle=":", linewidth=1))
        ax.axhline(0, color="#cbd5df", lw=.8)
        ax.axvline(0, color="#cbd5df", lw=.8)
        ax.set(xlim=(-.38, .38), ylim=(-.38, .38), aspect="equal",
               xlabel="real part", ylabel="imaginary part", title=label)
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle("Before orientational averaging: moving the reference rotates the phase", color=INK, fontsize=12)
    fig.text(.5, .025, "Orange: individual contributions / 6. Blue: their sum. The blue length stays the same.",
             ha="center", color=INK, fontsize=10)
    fig.tight_layout(rect=(0, .06, 1, .94))
    fig.savefig(OUTPUT / "scattering-reference-phase.svg", facecolor="white",
                metadata={"Date": None})
    if preview_directory is not None:
        fig.savefig(preview_directory / "scattering-reference-phase.png",
                    dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview-dir", type=Path,
                        help="Optionally also write PNG previews in this directory")
    main(parser.parse_args().preview_dir)
