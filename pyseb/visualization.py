"""Optional interactive 3D views of representative SEB structures.

Install with ``pip install 'pyseb[visualization]'``. Geometry and colors come
from the OpenUSD exporter; cloud points use display-sized markers.
"""

from pathlib import Path
import math
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


def _volume_values(volumes):
    result = {}
    for key, value in ({} if volumes is None else volumes).items():
        value = float(value)
        if not math.isfinite(value) or value <= 0:
            raise ValueError("material volumes must be finite and positive")
        result[str(key)] = value
    return result


def _usd_traces(path, *, references=False, parameters=None, volumes=None):
    """Convert the pySEB geometry in a USDA file to Plotly traces.

    OpenUSD resolves the scene hierarchy and transforms.
    """
    go, Usd, UsdGeom = _dependencies()

    stage = Usd.Stage.Open(str(path))
    if stage is None:
        raise ValueError(f"OpenUSD could not open {path}")
    transform_cache = UsdGeom.XformCache(Usd.TimeCode.Default())
    traces = []
    reference_groups = {}
    scattering = {}
    parameters = {} if parameters is None else parameters
    volumes = _volume_values(volumes)

    def number(value):
        try:
            value = float(value)
            return value if math.isfinite(value) else None
        except (TypeError, ValueError, OverflowError):
            return None

    def scattering_info(prim):
        instance = instance_name(prim)
        if instance in scattering:
            return scattering[instance]
        owner = prim
        while owner and not owner.IsPseudoRoot():
            if owner.GetAttribute("pyseb:model_id").Get():
                break
            owner = owner.GetParent()
        if owner and owner.IsPseudoRoot():
            owner = None
        tag = str(owner.GetAttribute("pyseb:tag").Get() or "") if owner else ""
        source = str(owner.GetAttribute("pyseb:source_name").Get() or "") if owner else ""
        beta = number(parameters.get(f"beta_{tag}", parameters.get("beta"))) if tag else None
        volume = None
        origin = None
        for key in (instance, instance.replace("/", ":"), tag, source):
            if key and key in volumes:
                volume, origin = volumes[key], "supplied material volume"
                break
        model = owner.GetAttribute("pyseb:model_id").Get() if owner else None
        if owner and volume is None:
            dims = {a.GetName()[12:]: number(a.Get()) for a in owner.GetAttributes()
                    if a.GetName().startswith("pyseb:param_")}
            try:
                if model == "pyseb/SolidSphere":
                    volume = 4 * math.pi * dims["R"]**3 / 3
                elif model == "pyseb/SolidSphericalShell":
                    volume = 4 * math.pi * (dims["Ro"]**3 - dims["Ri"]**3) / 3
                elif model == "pyseb/SolidCylinder":
                    volume = math.pi * dims["R"]**2 * dims["L"]
                elif model == "pyseb/Spheroid":
                    volume = 4 * math.pi * dims["a"]**2 * dims["c"] / 3
            except (KeyError, TypeError, OverflowError):
                volume = None
            if volume is not None:
                origin = "geometric material volume"
        if model == "pyseb/DebyeSphereCloud":
            values = [number(v) for child in Usd.PrimRange(owner)
                      for v in (child.GetAttribute("pyseb:beta").Get() or [])]
            if values and all(v is not None for v in values):
                beta = number(math.fsum(values))
        volume = number(volume)
        if volume is not None and volume <= 0:
            volume, origin = None, None
        result = {"beta": beta, "volume": volume, "volume_source": origin}
        scattering[instance] = result
        return result

    def instance_name(prim):
        current = prim
        while current and not current.IsPseudoRoot():
            name = current.GetAttribute("pyseb:instance_path").Get()
            if name:
                return str(name)
            current = current.GetParent()
        return str(prim.GetParent().GetPath()).strip("/") or str(prim.GetPath())

    def metadata(prim, role="geometry"):
        return {"pyseb": {"instance": instance_name(prim), "role": role,
                           **scattering_info(prim)}}

    def hover(prim, name):
        info = scattering_info(prim)
        lines = [name]
        if info["beta"] is not None:
            lines.append(f"β = {info['beta']:.5g}")
        if info["volume"] is not None:
            lines.append(f"V = {info['volume']:.5g} (model length units³)")
            if info["beta"] is not None:
                lines.append(f"β / V = {info['beta'] / info['volume']:.5g}")
        return "<br>".join(lines) + "<extra></extra>"

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
        reference = prim.GetAttribute("pyseb:reference").Get()
        if references and reference is not None:
            matrix = transform_cache.GetLocalToWorldTransform(prim)
            point = tuple(float(value) for value in matrix.ExtractTranslation())
            group = reference_groups.setdefault(instance_name(prim), [])
            group.append((point, str(reference)))
            continue
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
                meta=metadata(prim),
                color=color,
                opacity=opacity,
                flatshading=False,
                hovertemplate=hover(prim, name),
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
                meta=metadata(prim),
                line={"color": color, "width": 5},
                opacity=opacity,
                hovertemplate=hover(prim, name),
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
                meta=metadata(prim),
                marker={"color": rgb, "size": 5},
                opacity=opacity,
                hovertemplate=hover(prim, name),
                showlegend=False,
            ))

    for instance, points in reference_groups.items():
        x, y, z = zip(*(point for point, _ in points))
        traces.append(go.Scatter3d(
            x=x, y=y, z=z, mode="markers", visible=False,
            name=f"{instance}: references", text=[name for _, name in points],
            meta={"pyseb": {"instance": instance, "role": "reference"}},
            marker={"color": "#f5cf61", "size": 5, "symbol": "diamond"},
            hovertemplate=f"{instance}: %{{text}}<extra></extra>", showlegend=False,
        ))
    return traces


def usd_to_plotly(path, *, height=650, controls=True, parameters=None, volumes=None):
    """Read a USD scene into a Plotly figure with optional in-frame controls."""
    go, _, _ = _dependencies()

    if controls:
        from ._viewer import figure_class
        figure = figure_class()(_usd_traces(path, references=True, parameters=parameters, volumes=volumes))
    else:
        figure = go.Figure(_usd_traces(path, parameters=parameters, volumes=volumes))
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


def structure_figure(world, structure, parameters=None, options=None, *, height=650,
                     controls=True, seeds=None, layouts=None, volumes=None):
    """Return an interactive Plotly figure of a named structure.

    ``parameters`` is the parameter mapping accepted by ``world.export_usd``.
    Pass ``USDExportOptions`` as ``options`` to set units, colors, sampling,
    seed and layout. Defaults use the exporter's seeded, representative layout.
    Temporary scene files are removed before returning; the figure is standalone.
    HTML views have a hover-revealed settings panel. ``controls=False`` returns
    a plain Plotly figure. Optional ``seeds`` and ``layouts`` prepare realizations
    for switching in static HTML, without a running Python session. Layouts are
    ``"random"`` or ``"readable"`` (or the corresponding ``LayoutMode`` values).
    ``volumes`` optionally maps occurrence paths, tags, or source names to
    physical material volumes in the cube of the model's length unit. Solids
    use analytic volumes automatically; polymer coil envelopes are not volumes.
    """
    _dependencies()
    volumes = _volume_values(volumes)
    from . import LengthUnit, USDExportOptions, LayoutMode

    if options is None:
        options = USDExportOptions(LengthUnit.Angstrom)
    seed_values = list(dict.fromkeys([options.seed, *(() if seeds is None else seeds)]))
    layout_values = [options.layout_mode]
    for layout in (() if layouts is None else layouts):
        if isinstance(layout, str):
            try:
                layout = {"random": LayoutMode.Random, "readable": LayoutMode.Readable}[layout.lower()]
            except KeyError as exc:
                raise ValueError("layouts must contain 'random' or 'readable'") from exc
        if layout not in (LayoutMode.Random, LayoutMode.Readable):
            raise ValueError("layouts must contain 'random' or 'readable'")
        if layout not in layout_values:
            layout_values.append(layout)
    if not controls and (seeds is not None or layouts is not None):
        raise ValueError("seeds and layouts require controls=True")

    # pybind export options do not support deepcopy. Copy their public fields
    # so preparing another seed/layout never mutates the caller's options.
    export_options = USDExportOptions(options.meters_per_unit)
    for name in dir(options):
        if not name.startswith("_"):
            setattr(export_options, name, getattr(options, name))
    variants = []
    figure = None
    with TemporaryDirectory(prefix="pyseb-view-") as directory:
        path = Path(directory) / "structure.usda"
        for layout in layout_values:
            for seed in seed_values:
                export_options.seed = seed
                export_options.layout_mode = layout
                world.export_usd(structure, str(path), {} if parameters is None else parameters, export_options)
                current = usd_to_plotly(path, height=height, controls=controls,
                                        parameters=parameters, volumes=volumes)
                if figure is None:
                    figure = current
                if len(seed_values) * len(layout_values) > 1:
                    mode = "Readable" if layout == LayoutMode.Readable else "Random"
                    variants.append({"label": f"{mode} · seed {seed}",
                                     "layout": mode.lower(), "seed": str(seed),
                                     "data": current.to_plotly_json()["data"]})
    if controls:
        figure._pyseb_variants = variants
    return figure


def visualize(world, structure, parameters=None, options=None, *, height=650,
              controls=True, seeds=None, layouts=None, volumes=None, **show_kwargs):
    """Display a structure in 3D and return its Plotly figure.

    Hover inside the HTML view for settings in its upper-left corner.
    Optional ``seeds`` and ``layouts`` prepare switchable realizations;
    ``controls=False`` uses standard Plotly rendering.

    Extra keyword arguments go to ``Figure.show`` (for example ``renderer``
    or ``config``). In notebooks, use ``structure_figure`` as the last expression
    to display once using the notebook's normal figure rendering.
    """
    figure = structure_figure(world, structure, parameters, options, height=height,
                              controls=controls, seeds=seeds, layouts=layouts, volumes=volumes)
    show_kwargs.setdefault("config", {"displaylogo": False, "scrollZoom": True})
    figure.show(**show_kwargs)
    return figure


show_structure = visualize
