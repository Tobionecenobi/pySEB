"""Helpers used while rendering the executable tutorial pages."""

import os
from pathlib import Path
import subprocess


def repository_root():
    """Return the repository containing the documentation being rendered."""
    current = Path.cwd().resolve()
    return next(
        candidate
        for candidate in (current, *current.parents)
        if (candidate / "pyproject.toml").is_file()
    )


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


def run_usd_example(script_name, output_name):
    """Run a maintained USD example in a temporary documentation workspace."""
    import sys
    import tempfile

    repository = repository_root()
    workspace = tempfile.TemporaryDirectory()
    environment = os.environ.copy()
    existing_python_path = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = os.pathsep.join(
        value for value in (str(repository), existing_python_path) if value
    )
    subprocess.run(
        [
            sys.executable,
            repository / "pyseb" / "examples" / "python" / script_name,
        ],
        cwd=workspace.name,
        env=environment,
        check=True,
    )
    output = Path(workspace.name) / output_name
    if not output.is_file():
        raise FileNotFoundError(f"example did not create {output_name}")
    return workspace, output


from pyseb.visualization import _usd_traces, usd_to_plotly


def bundled_subunit_catalogue():
    """Load bundled YAML models and choose validated preview parameters.

    A bundled model appears automatically when a new ``models/*.pyseb.yaml``
    file is added. Its first validation case containing every structural
    parameter supplies the deterministic scale used by the homepage preview.
    """
    import pyseb
    import yaml

    catalogue = []
    for path in sorted((repository_root() / "models").glob("*.pyseb.yaml")):
        definition = pyseb.load_subunit_definition(str(path))
        if not definition.visualization.present:
            raise ValueError(
                f"bundled model {definition.id} needs visualization metadata "
                "for the generated homepage catalogue"
            )

        document = yaml.safe_load(path.read_text(encoding="utf-8"))
        required = set(definition.parameters)
        preview_parameters = None
        for case in document.get("validation", {}).get("cases", []):
            values = case.get("parameters", {})
            if required.issubset(values):
                preview_parameters = {
                    name: float(values[name]) for name in sorted(required)
                }
                break
        if preview_parameters is None:
            names = ", ".join(sorted(required)) or "no parameters"
            raise ValueError(
                f"bundled model {definition.id} needs a validation case covering "
                f"all preview parameters ({names})"
            )

        catalogue.append({
            "api_name": definition.api_name,
            "description": definition.metadata.description,
            "id": definition.id,
            "parameters": preview_parameters,
            "path": path,
            "title": definition.metadata.title,
        })

    return catalogue


def bundled_subunit_markdown(catalogue):
    """Return the generated homepage bullet for the bundled catalogue."""
    names = ", ".join(f"`{model['api_name']}`" for model in catalogue)
    return f"- **{len(catalogue)} bundled YAML subunits:** {names}."


def bundled_subunit_figure(catalogue, *, columns=3):
    """Render the catalogue with both layouts and the standard viewer controls."""
    import math
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import pyseb
    from pyseb._viewer import figure_class

    rows = math.ceil(len(catalogue) / columns)
    specs = [[{"type": "scene"} for _ in range(columns)] for _ in range(rows)]
    titles = [model["api_name"] for model in catalogue]
    titles.extend([""] * (rows * columns - len(titles)))
    galleries = {layout: make_subplots(
        rows=rows, cols=columns, specs=specs, subplot_titles=titles,
        horizontal_spacing=.015, vertical_spacing=.06)
        for layout in ("random", "readable")}
    palette = ((.12, .38, .72), (.10, .62, .66), (.32, .68, .38),
               (.91, .63, .13), (.82, .25, .18), (.56, .36, .73))
    for index, model in enumerate(catalogue):
        row, column = divmod(index, columns)
        color = palette[index % len(palette)]
        world = pyseb.World(f"catalogue-{model['api_name']}")
        graph = world.Add(model["id"], "sample", "preview")
        world.Add(graph, "previewStructure")
        parameters = {f"{name}_preview": value
                      for name, value in model["parameters"].items()}
        options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
        options.seed = 1000 + index
        options.curve_samples = 64
        options.surface_samples = 22
        options.color_overrides = {"sample": color}
        options.opacity_overrides = {"sample": .82}
        preview = pyseb.structure_figure(
            world, "previewStructure", parameters, options,
            layouts=["random", "readable"])
        for variant in preview._pyseb_variants:
            traces = [trace for trace in go.Figure(variant['data']).data
                      if trace.meta['pyseb']['role'] == 'geometry']
            if not traces:
                # Point is an invisible connector; show its location in the catalogue.
                rgb = "rgb({},{},{})".format(*(round(255*c) for c in color))
                traces = [go.Scatter3d(x=[0], y=[0], z=[0], mode="markers",
                    marker={"color": rgb, "size": 8, "symbol": "diamond"},
                    hovertemplate=f"{model['api_name']}: fixed connector reference<extra></extra>",
                    showlegend=False)]
            for trace in traces:
                trace.name = model['api_name']
                trace.meta = {'pyseb': {'instance': model['api_name'],
                                        'role': 'geometry', 'neighbors': []}}
                galleries[variant['layout']].add_trace(trace, row=row+1, col=column+1)
    figure = figure_class()(galleries['random'])
    scene_style = {"aspectmode": "data", "xaxis": {"visible": False},
                   "yaxis": {"visible": False}, "zaxis": {"visible": False},
                   "camera": {"eye": {"x": 1.45, "y": 1.45, "z": 1.0}},
                   "bgcolor": "rgba(0,0,0,0)"}
    for index in range(1, rows*columns+1):
        figure.layout['scene' if index == 1 else f'scene{index}'].update(scene_style)
    figure.update_layout(height=270*rows, margin={"l": 0, "r": 0, "t": 70, "b": 0},
                         paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                         showlegend=False)
    figure._pyseb_variants = [
        {'layout': layout, 'seed': 'per-subunit catalogue seeds',
         'label': 'Model' if layout == 'random' else 'Unfolded',
         'data': gallery.to_plotly_json()['data']}
        for layout, gallery in galleries.items()]
    return figure
