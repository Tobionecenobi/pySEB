"""Ordered calculation source for the linked-titin tutorial sequence.

Only tagged calculation cells contribute to the provenance hash. Moving a cell
between these pages must preserve its contents and the sequence order.
"""
from pathlib import Path
import re

TUTORIAL_NAMES = (
    "fit-linked-titin.qmd",
    "titin-coordinate-domains.qmd",
    "titin-hydration.qmd",
    "titin-relative-weights.qmd",
)


def generation_cells(directory):
    cells = []
    for name in TUTORIAL_NAMES:
        source = (Path(directory) / name).read_text()
        cells.extend(block for block in re.findall(
            r"```\{python\}\n(.*?)\n```", source, re.S)
            if block.splitlines()[0] == "#| tags: [titin-generation]")
    return cells
