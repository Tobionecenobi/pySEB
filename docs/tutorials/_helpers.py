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


def run_usd_example(script_name, output_name):
    """Run a maintained USD example in a temporary documentation workspace."""
    import sys
    import tempfile

    current = Path.cwd().resolve()
    repository = next(
        candidate
        for candidate in (current, *current.parents)
        if (candidate / "pyproject.toml").is_file()
    )
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


def usd_to_plotly(path, *, height=650):
    """Convert the pySEB geometry in a USDA file to an interactive Plotly figure.

    This documentation-only adapter deliberately lives outside the pyseb
    package. OpenUSD resolves the scene hierarchy and transforms; Plotly
    supplies the browser renderer used by the generated tutorial site.
    """
    import plotly.graph_objects as go
    from pxr import Usd, UsdGeom

    stage = Usd.Stage.Open(str(path))
    if stage is None:
        raise ValueError(f"OpenUSD could not open {path}")
    transform_cache = UsdGeom.XformCache(Usd.TimeCode.Default())
    traces = []

    def transformed_points(prim, points):
        matrix = transform_cache.GetLocalToWorldTransform(prim)
        return [tuple(float(value) for value in matrix.Transform(point)) for point in points]

    def display_style(prim):
        gprim = UsdGeom.Gprim(prim)
        colors = gprim.GetDisplayColorPrimvar().ComputeFlattened()
        opacities = gprim.GetDisplayOpacityPrimvar().ComputeFlattened()
        color = colors[0] if colors else (0.55, 0.65, 0.8)
        opacity = float(opacities[0]) if opacities else 1.0
        rgb = "rgb({},{},{})".format(
            *(round(255 * float(channel)) for channel in color)
        )
        return rgb, opacity

    for prim in stage.Traverse():
        name = str(prim.GetPath())

        if prim.IsA(UsdGeom.Mesh):
            mesh = UsdGeom.Mesh(prim)
            points = transformed_points(prim, mesh.GetPointsAttr().Get() or [])
            counts = mesh.GetFaceVertexCountsAttr().Get() or []
            indices = mesh.GetFaceVertexIndicesAttr().Get() or []
            triangles = []
            cursor = 0
            for count in counts:
                polygon = indices[cursor:cursor + count]
                cursor += count
                triangles.extend(
                    (int(polygon[0]), int(polygon[offset]), int(polygon[offset + 1]))
                    for offset in range(1, count - 1)
                )
            if not points or not triangles:
                continue
            x, y, z = zip(*points)
            i, j, k = zip(*triangles)
            color, opacity = display_style(prim)
            traces.append(go.Mesh3d(
                x=x, y=y, z=z, i=i, j=j, k=k,
                name=name,
                color=color,
                opacity=opacity,
                flatshading=False,
                hovertemplate=f"{name}<extra></extra>",
                showlegend=False,
            ))

        elif prim.IsA(UsdGeom.BasisCurves):
            curves = UsdGeom.BasisCurves(prim)
            points = transformed_points(prim, curves.GetPointsAttr().Get() or [])
            counts = curves.GetCurveVertexCountsAttr().Get() or []
            coordinates = [[], [], []]
            cursor = 0
            for count in counts:
                for point in points[cursor:cursor + count]:
                    for axis in range(3):
                        coordinates[axis].append(point[axis])
                cursor += count
                for axis in range(3):
                    coordinates[axis].append(None)
            if not points:
                continue
            color, opacity = display_style(prim)
            traces.append(go.Scatter3d(
                x=coordinates[0], y=coordinates[1], z=coordinates[2],
                mode="lines",
                name=name,
                line={"color": color, "width": 5},
                opacity=opacity,
                hovertemplate=f"{name}<extra></extra>",
                showlegend=False,
            ))

        elif prim.IsA(UsdGeom.PointInstancer):
            instancer = UsdGeom.PointInstancer(prim)
            points = transformed_points(
                prim, instancer.GetPositionsAttr().Get() or [])
            if not points:
                continue
            x, y, z = zip(*points)
            colors = prim.GetAttribute("primvars:displayColor").Get() or []
            color = colors[0] if colors else (0.65, 0.65, 0.65)
            rgb = "rgb({},{},{})".format(
                *(round(255 * float(channel)) for channel in color)
            )
            opacities = prim.GetAttribute("primvars:displayOpacity").Get() or []
            opacity = float(opacities[0]) if opacities else 1.0
            traces.append(go.Scatter3d(
                x=x, y=y, z=z,
                mode="markers",
                name=name,
                marker={"color": rgb, "size": 5, "opacity": opacity},
                hovertemplate=f"{name}<extra></extra>",
                showlegend=False,
            ))

    figure = go.Figure(traces)
    figure.update_layout(
        height=height,
        margin={"l": 0, "r": 0, "t": 35, "b": 0},
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        scene={
            "aspectmode": "data",
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "zaxis": {"visible": False},
            "camera": {"eye": {"x": 1.5, "y": 1.5, "z": 1.1}},
            "bgcolor": "rgba(0,0,0,0)",
        },
        title={
            "text": "Interactive OpenUSD realization — drag to rotate, scroll to zoom",
            "x": 0.5,
            "xanchor": "center",
            "font": {"size": 15},
        },
        showlegend=False,
    )
    return figure
