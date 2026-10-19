"""Rebuild preceding AOT stages from their displayed calculation cells.

Only explicitly tagged cells execute. Plotting, visualization and prose never
run through this helper. Paths are relative to the tutorial checkout so each
continuation page works from a fresh notebook kernel.
"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import re


def previous_stage(stage):
    pages = {
        'baseline': ('fit-real-sans-data.qmd',),
        'distributions': ('fit-real-sans-data.qmd', 'sans-size-distributions.qmd'),
    }
    if stage not in pages:
        raise ValueError(f'Unknown AOT stage: {stage}')
    directory = Path(__file__).resolve().parent
    if Path.cwd().resolve() != directory:
        raise ValueError(f'Run the AOT examples from {directory}')
    namespace = {'__name__': '__aot_setup__'}
    with redirect_stdout(StringIO()):
        for name in pages[stage]:
            cells = re.findall(r'```\{python\}\n(.*?)\n```',
                               (directory / name).read_text(), re.S)
            tagged = [cell for cell in cells if cell.startswith('#| tags: [aot-setup]\n')]
            if not tagged:
                raise ValueError(f'No calculation setup found in {name}')
            for index, cell in enumerate(tagged, 1):
                exec(compile(cell, f'{name}:setup-{index}', 'exec'), namespace)
    print(f'Rebuilt AOT {stage} stage from the preceding tutorial calculations.')
    return namespace
