"""Optional interactive 3D views of representative SEB structures.

Install with ``pip install 'pyseb[visualization]'``. Geometry and colors come
from the OpenUSD exporter; cloud points use display-sized markers.
"""

from pathlib import Path
from tempfile import TemporaryDirectory


def _dependencies():
    try:
        import plotly.graph_objects as go
        from pxr import Usd, UsdGeom
    except ImportError as exc:
        raise ImportError(
            "3D visualization requires plotly and usd-core. "
            "Install them with: pip install 'pyseb[visualization]'"
        ) from exc
    return go, Usd, UsdGeom


def _usd_traces(path):
    """Convert the pySEB geometry in a USDA file to Plotly traces.

    OpenUSD resolves the scene hierarchy and transforms.
    """
    go, Usd, UsdGeom = _dependencies()

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
        imageable = UsdGeom.Imageable(prim)
        if imageable and imageable.ComputeVisibility() == UsdGeom.Tokens.invisible:
            continue
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
            if opacity == 0:
                continue
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
            if opacity == 0:
                continue
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
            if opacity == 0:
                continue
            traces.append(go.Scatter3d(
                x=x, y=y, z=z,
                mode="markers",
                name=name,
                marker={"color": rgb, "size": 5, "opacity": opacity},
                hovertemplate=f"{name}<extra></extra>",
                showlegend=False,
            ))

    return traces


def usd_to_plotly(path, *, height=650):
    """Convert the pySEB geometry in a USDA file to an interactive Plotly figure."""
    go, _, _ = _dependencies()

    figure = go.Figure(_usd_traces(path))
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


def structure_figure(world, structure, parameters=None, options=None, *, height=650):
    """Return an interactive Plotly figure of a named structure.

    ``parameters`` is the parameter mapping accepted by ``world.export_usd``.
    Pass ``USDExportOptions`` as ``options`` to set units, colors, sampling,
    seed and layout. Defaults use the exporter's seeded, representative layout.
    Temporary scene files are removed before returning; the figure is standalone.
    """
    _dependencies()
    from . import LengthUnit, USDExportOptions

    if options is None:
        options = USDExportOptions(LengthUnit.Angstrom)
    with TemporaryDirectory(prefix="pyseb-view-") as directory:
        path = Path(directory) / "structure.usda"
        world.export_usd(structure, str(path), {} if parameters is None else parameters, options)
        return usd_to_plotly(path, height=height)


def visualize(world, structure, parameters=None, options=None, *, height=650, **show_kwargs):
    """Display a structure in 3D and return its Plotly figure.

    Extra keyword arguments go to ``Figure.show`` (for example ``renderer``
    or ``config``). In notebooks, use ``structure_figure`` as the last expression
    to display once using the notebook's normal figure rendering.
    """
    figure = structure_figure(world, structure, parameters, options, height=height)
    show_kwargs.setdefault("config", {"displaylogo": False, "scrollZoom": True})
    figure.show(**show_kwargs)
    return figure


show_structure = visualize
