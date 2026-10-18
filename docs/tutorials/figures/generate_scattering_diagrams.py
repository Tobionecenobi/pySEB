"""Rebuild the vector diagrams for scattering-from-phases.qmd.

Run with the documentation Python environment from any directory. The fan sketch is illustrative; the rod sketch depicts the continuous material
used in the worked derivation.
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

    # Four vector constructions on the continuous rod used in the tutorial.
    fig, axes = plt.subplots(2, 2, figsize=(11, 5.8), facecolor="white")
    titles = ("1. Form factor: scatterer pairs", "2. Amplitude: easiest reference",
              "3. Amplitude: subtract the shift", "4. Phase factor: reference pairs")
    for ax, title in zip(axes.flat, titles):
        ax.plot([-2, 2], [0, 0], color=INK, lw=7, solid_capstyle="round")
        ax.set(xlim=(-2.8, 2.8), ylim=(-1.4, 1.2))
        ax.set_title(title, loc="left", color=INK, weight="bold", fontsize=12)
        ax.axis("off")
        ax.text(-2, -.27, r"$-L/2$", ha="center", color=INK)
        ax.text(2, -.27, r"$L/2$", ha="center", color=INK)
    ax = axes[0, 0]
    ax.scatter([-.9, 1.2], [0, 0], color=ORANGE, s=45, zorder=4)
    ax.text(-.9, .2, r"$t$", ha="center", color=INK)
    ax.text(1.2, .2, r"$s$", ha="center", color=INK)
    arrow(ax, (-.9, .65), (1.2, .65), color=BLUE, width=2.2)
    ax.text(.15, .88, r"$(s-t)\hat{\mathbf{u}}$", ha="center", color=BLUE)
    ax.text(0, -1.05, "Integrate over both scattering elements", ha="center", color=INK)
    ax = axes[0, 1]
    reference(ax, (0, 0), r"$R_m=0$", (-22, -32))
    for endpoint, height in ((-1.4, .48), (1.4, .72)):
        ax.plot([0, 0], [0, height], color="#cbd5df", lw=1)
        ax.plot([endpoint, endpoint], [0, height], color="#cbd5df", lw=1)
        arrow(ax, (0, height), (endpoint, height))
    ax.text(1.4, .94, r"$s\hat{\mathbf{u}}$", ha="center", color=ORANGE)
    ax.text(0, -1.05, "Integrate from the midpoint to the material", ha="center", color=INK)
    ax = axes[1, 0]
    reference(ax, (-2, 0), r"$R_{new}$", (-15, -34))
    ax.scatter(0, 0, color=GOLD, edgecolor=INK, s=65, zorder=4)
    ax.plot([-2, -2], [0, .7], color="#cbd5df", lw=1)
    ax.plot([1.2, 1.2], [0, .7], color="#cbd5df", lw=1)
    arrow(ax, (-2, .7), (1.2, .7))
    ax.text(-.4, .94, r"$s\hat{\mathbf{u}}-\mathbf{r}$", ha="center", color=ORANGE)
    arrow(ax, (0, -.77), (-2, -.77), color=BLUE)
    ax.text(-1, -1.12, r"$\mathbf{r}=\mathbf{R}_{new}-\mathbf{R}_m$", ha="center", color=BLUE)
    ax.text(.22, -.35, r"$R_m$", color=INK)
    ax = axes[1, 1]
    reference(ax, (-2, 0), "end1", (-12, -32))
    reference(ax, (0, 0), "middle", (-18, -32))
    reference(ax, (2, 0), "end2", (-12, -32))
    arrow(ax, (-2, .65), (2, .65), color=BLUE, width=2.2)
    ax.text(0, .87, r"$|\mathbf{R}_2-\mathbf{R}_1|=L$", ha="center", color=BLUE)
    ax.text(0, -1.05, "Use reference separations; no material integral", ha="center", color=INK)
    fig.subplots_adjust(left=.05, right=.98, bottom=.04, top=.93, hspace=.32, wspace=.18)
    fig.savefig(OUTPUT / "rod-derivation-vectors.svg", facecolor="white",
                metadata={"Date": None})
    svg_path = OUTPUT / "rod-derivation-vectors.svg"
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n")
    if preview_directory is not None:
        fig.savefig(preview_directory / "rod-derivation-vectors.png", dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview-dir", type=Path,
                        help="Optionally also write PNG previews in this directory")
    main(parser.parse_args().preview_dir)
