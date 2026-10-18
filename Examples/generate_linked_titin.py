#!/usr/bin/env python3
"""Regenerate the linked-titin tutorial's fits and OpenUSD scenes.

The displayed, tagged Python cells in the tutorial sequence are the
calculation source.
No Quarto/Jupyter process is needed. Publication happens only after every
stage succeeds; an interrupted fit leaves the previous bundle intact.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = ROOT / "docs/tutorials/fit-linked-titin.qmd"
sys.path.insert(0, str(TUTORIAL.parent))
from _titin_sources import generation_cells as tutorial_generation_cells
DESTINATION = ROOT / "docs/tutorials/data/titin-results"
INPUT_NAMES = ("SASDKC2-titin.csv", "1TIU-heavy-ensemble.pdb")


def generation_cells():
    return tutorial_generation_cells(TUTORIAL.parent)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    sys.path.insert(0, str(ROOT))
    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib.pyplot as plt
    import plotly.graph_objects as go
    from pxr import Usd

    # The QMD makes every plot and interactive view from the saved data.
    # Generation exports real USDA scenes, without opening browser windows.
    go.Figure.show = lambda self, *args, **kwargs: None
    # The controlled HTML viewer subclasses Figure and supplies its own show.
    from pyseb._viewer import figure_class
    figure_class().show = go.Figure.show
    plt.show = lambda *args, **kwargs: plt.close("all")
    cells = generation_cells()
    if not cells:
        raise RuntimeError("No titin-generation cells found in the tutorial")
    started = time.monotonic()
    source_directory = ROOT / "docs/tutorials/data"
    with tempfile.TemporaryDirectory(prefix="pyseb-titin-") as workspace:
        working = Path(workspace)
        (working / "data").mkdir()
        for name in INPUT_NAMES:
            shutil.copyfile(source_directory / name, working / "data" / name)
        previous_directory = Path.cwd()
        os.chdir(working)
        try:
            namespace = {"__name__": "__titin_generation__"}
            for index, cell in enumerate(cells, 1):
                print(f"Titin generation cell {index}/{len(cells)} "
                      f"({time.monotonic()-started:.1f}s)", flush=True)
                exec(compile(cell, f"{TUTORIAL.name}:generation-{index}", "exec"), namespace)
            bundle = working / "data/titin-results"
            scenes = sorted((bundle / "scenes").glob("*.usda"))
            if len(scenes) != 6:
                raise RuntimeError(f"Expected six OpenUSD scenes, got {len(scenes)}")
            for path in scenes:
                stage = Usd.Stage.Open(str(path))
                if stage is None or not list(stage.Traverse()):
                    raise RuntimeError(f"Invalid OpenUSD scene: {path}")
            versions = {}
            for package in ("numpy", "scipy", "matplotlib", "plotly", "usd-core", "pyseb"):
                try:
                    versions[package] = importlib.metadata.version(package)
                except importlib.metadata.PackageNotFoundError:
                    versions[package] = "local build"
            manifest = {
                "schema_version": 1,
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "generation_sha256": hashlib.sha256("\n".join(cells).encode()).hexdigest(),
                "generator_sha256": digest(Path(__file__).resolve()),
                "git_commit": subprocess.check_output(
                    ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
                "python": sys.version,
                "versions": versions,
                "inputs": {name: digest(working / "data" / name) for name in INPUT_NAMES},
                "files": {str(path.relative_to(bundle)): digest(path)
                          for path in sorted(bundle.rglob("*")) if path.is_file()},
                "generation_seconds": time.monotonic()-started,
                "note": "Generation-cell hash identifies the working-tree calculation; git commit alone does not.",
            }
            (bundle / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
        finally:
            os.chdir(previous_directory)
        DESTINATION.mkdir(parents=True, exist_ok=True)
        for path in sorted(bundle.rglob("*")):
            if path.is_file() and path.name != "manifest.json":
                target = DESTINATION / path.relative_to(bundle)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, target)
        # Copy the manifest last: the render loader detects incomplete updates.
        shutil.copyfile(bundle / "manifest.json", DESTINATION / "manifest.json")
    print(f"Saved {DESTINATION} ({time.monotonic()-started:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
