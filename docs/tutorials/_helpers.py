"""Helpers used while rendering the executable tutorial pages."""

import os
from pathlib import Path
import subprocess


def tutorial_binary(name):
    configured = os.environ.get("PYSEB_TUTORIAL_BIN_DIR")
    if configured:
        directory = Path(configured)
    else:
        repository = Path.cwd().parents[1]
        directory = repository / "build-docs" / "docs-tutorials"
    return directory / f"tutorial-{name}"


def run_tutorial(name):
    """Run a compiled documentation tutorial and display its stdout."""
    result = subprocess.run(
        [tutorial_binary(name)],
        check=True,
        capture_output=True,
        text=True,
    )
    print(result.stdout, end="")
    return result.stdout
