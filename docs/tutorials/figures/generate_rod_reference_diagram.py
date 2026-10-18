"""Draw the symmetric third references used in add-subunit-reference.qmd."""
import argparse
from pathlib import Path

from generate_scattering_diagrams import plt, INK, BLUE, GOLD, ORANGE


def main(preview_directory=None):
    fig, ax = plt.subplots(figsize=(10, 2.7), facecolor="white")
    ax.plot([0, 1], [0, 0], color=INK, lw=5, solid_capstyle="round")
    for x, label, coordinate, color, label_y in (
            (0, "end1", "0", INK, -.30),
            (1/3, "third1 (new)", r"$L/3$", GOLD, -.30),
            (.5, "middle", r"$L/2$", ORANGE, -.57),
            (2/3, "third2 (new)", r"$2L/3$", GOLD, -.30),
            (1, "end2", r"$L$", INK, -.30)):
        ax.scatter(x, 0, s=90, color=color, edgecolor=INK, zorder=4)
        ax.plot([x, x], [-.08, .08], color=INK, lw=1)
        ax.text(x, label_y, label, ha="center", color=INK, fontsize=11)
        ax.text(x, label_y-.20, coordinate, ha="center", color=INK, fontsize=12)
    for start, end, label in ((0, 1/3, r"$L/3$"), (1/3, 2/3, r"$L/3$"), (2/3, 1, r"$L/3$")):
        ax.annotate("", (end, .45), (start, .45),
                    arrowprops={"arrowstyle": "<->", "color": BLUE, "lw": 1.7})
        ax.text((start+end)/2, .53, label, ha="center", color=BLUE, fontsize=13)
    ax.set(xlim=(-.08, 1.08), ylim=(-1., .95))
    ax.set_title("The rod is unchanged; add two symmetric attachment locations", color=INK)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(Path(__file__).with_name("rod-third-reference.svg"),
                metadata={"Date": None}, facecolor="white")
    if preview_directory is not None:
        preview_directory.mkdir(parents=True, exist_ok=True)
        fig.savefig(preview_directory / "rod-third-reference.png", dpi=150,
                    facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview-dir", type=Path)
    main(parser.parse_args().preview_dir)
