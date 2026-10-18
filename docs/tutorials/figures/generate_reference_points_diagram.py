"""Draw the reference-point primer without adding material at the references."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTPUT = Path(__file__).resolve().parent
INK = "#172f46"
ORANGE = "#b95017"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "svg.fonttype": "none", "svg.hashsalt": "pyseb-references"})


def main(preview_directory=None):
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.4), facecolor="white")
    for ax, title in zip(axes, ("Scattering material", "The same material, with reference labels")):
        ax.scatter([-1, 1], [0, 0], color=INK, s=100, zorder=3)
        for x in (-1, 0, 1, 2):
            ax.text(x, -.48, str(x), ha="center", color=INK)
        ax.set(xlim=(-1.8, 2.8), ylim=(-1.25, 1.05))
        ax.set_title(title, color=INK, loc="left", fontsize=12, weight="bold")
        ax.axis("off")
    axes[0].text(0, .38, "Each filled point has scattering weight 1", ha="center", color=INK, fontsize=10)
    axes[0].text(.5, -1.04, "2 scatterers", ha="center", color=INK, weight="bold")
    for x, label in ((-1, "left"), (0, "middle"), (2, "outside")):
        axes[1].scatter([x], [0], s=280, marker="D", facecolors="none",
                        edgecolors=ORANGE, linewidths=2, zorder=4)
        axes[1].text(x, .35, label, ha="center", color=ORANGE, fontsize=11, weight="bold")
    axes[1].text(.5, -1.04, "Still 2 scatterers; 3 reference locations", ha="center", color=INK, weight="bold")
    fig.text(.5, .04, "Filled dot: scattering material     Hollow diamond: reference location (no extra scattering weight)",
             ha="center", color=INK, fontsize=10)
    fig.subplots_adjust(left=.05, right=.97, top=.84, bottom=.17, wspace=.2)
    path = OUTPUT / "reference-points-vs-scatterers.svg"
    fig.savefig(path, facecolor="white", metadata={"Date": None})
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    if preview_directory:
        preview_directory.mkdir(parents=True, exist_ok=True)
        fig.savefig(preview_directory / "reference-points-vs-scatterers.png", dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview-dir", type=Path)
    main(parser.parse_args().preview_dir)
