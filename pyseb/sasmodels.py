"""Optional, isotropic 1-D bridge to sasmodels (install ``pyseb[sasmodels]``).

The world factory returns either a World or an object with a ``world`` attribute.
``bind(context, values)`` must apply *all* mutable settings and return a complete
pySEB parameter dictionary. Values are plain floats, not sasmodels Parameters.
"""
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Callable, Mapping, Optional, Sequence
import keyword
import math
import re
import threading

import numpy as np


@dataclass(frozen=True)
class Parameter:
    name: str
    default: float
    bounds: tuple = (-math.inf, math.inf)
    units: str = ""
    description: str = ""
    polydisperse: bool = False


@dataclass(frozen=True)
class InteractionSpec:
    """Explicit opt-in to a reference-based decoupling approximation."""
    reference: str
    effective_radius: Callable
    approximation: str = "decoupling"


@dataclass(frozen=True)
class SymbolicSpec:
    """Map World parameter symbols into public parameter expressions.

    amplitude, when needed, is an explicitly unnormalized World expression.
    It is independently checked against native evaluation at default parameters.
    volume and effective_radius use public parameter symbols.
    """
    parameters: Mapping
    volume: object = None
    amplitude: object = None
    effective_radius: object = None


@dataclass(frozen=True)
class ModelSpec:
    name: str
    parameters: Sequence[Parameter]
    world_factory: Callable
    structure: str
    bind: Callable
    intensity_factor: float
    volume: Optional[Callable] = None
    interaction: Optional[InteractionSpec] = None
    symbolic: Optional[SymbolicSpec] = None
    title: str = ""
    description: str = ""
    # These names are passed to volume/effective_radius callbacks. None means all.
    volume_parameters: Optional[Sequence[str]] = None

    @property
    def capabilities(self):
        return dict(intensity=True, resolution=True,
                    dispersity=any(p.polydisperse for p in self.parameters),
                    absolute_scale=self.volume is not None,
                    decoupling=self.interaction is not None,
                    compiled=self.symbolic is not None)


def _sasmodels():
    try:
        from sasmodels.modelinfo import make_model_info
    except ImportError as exc:
        raise ImportError('Install the optional integration: pip install "pyseb[sasmodels]"') from exc
    return make_model_info


def _validate(spec):
    reserved = {"q", "scale", "background"}
    if not re.fullmatch(r"[a-z][a-z0-9_]*", spec.name):
        raise ValueError("Model name must be a lowercase sasmodels identifier")
    names = [p.name for p in spec.parameters]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate parameter names")
    for p in spec.parameters:
        if (not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", p.name)
                or keyword.iskeyword(p.name) or p.name in reserved):
            raise ValueError(f"Invalid or reserved parameter name: {p.name}")
        if not p.bounds[0] <= p.default <= p.bounds[1] or not np.isfinite(p.default):
            raise ValueError(f"Invalid default/bounds for {p.name}")
    if not np.isfinite(spec.intensity_factor) or spec.intensity_factor <= 0:
        raise ValueError("intensity_factor must be an explicit positive finite conversion")
    if spec.volume_parameters is not None and not set(spec.volume_parameters) <= set(names):
        raise ValueError("Unknown volume parameter")
    if spec.interaction is not None:
        if spec.volume is None or spec.interaction.approximation != "decoupling":
            raise ValueError("Decoupling requires a particle volume and an explicit approximation")
        if not spec.interaction.reference or not callable(spec.interaction.effective_radius):
            raise ValueError("Interactions require a reference and effective-radius callback")


def _positive(value, label):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be finite and positive")
    return value


def _restrict_distributions(info, spec):
    allowed = {p.name for p in spec.parameters if p.polydisperse}
    for p in info.parameters.call_parameters:
        if p.name in {v.name for v in spec.parameters}:
            p.polydisperse = p.name in allowed
    info.parameters.pd_1d.intersection_update(allowed)
    info.parameters.pd_2d.intersection_update(allowed)
    return info


def make_model_info(spec, backend="runtime"):
    """Return sasmodels ModelInfo; use its build_model/DirectModel interfaces.

    Intensity is ``scale * factor * <I_raw> / <volume> + background``.
    Without a declared volume, volume=1 and scale has no implied physical units.
    Numerical factories/callbacks require the runtime backend.
    """
    make_info = _sasmodels()
    _validate(spec)
    if backend not in ("runtime", "compiled"):
        raise ValueError("backend must be 'runtime' or 'compiled'")
    names = [p.name for p in spec.parameters]
    volume_names = list(spec.volume_parameters) if spec.volume_parameters is not None else names
    # All volume dependencies must reach form_volume, even when not dispersible.
    module = ModuleType(spec.name)
    module.__file__ = str(Path(__file__).with_name(spec.name + ".py"))
    module.name = spec.name
    module.title = spec.title or spec.name
    module.description = spec.description or "Connected pySEB particle"
    module.__doc__ = module.description
    module.parameters = [
        [p.name, p.units, p.default, p.bounds,
         "volume" if p.name in volume_names or p.polydisperse else "",
         p.description or p.name] for p in spec.parameters
    ]
    # sasmodels calls form_volume with all parameters tagged 'volume'.
    volume_names = [p.name for p in spec.parameters
                    if p.name in volume_names or p.polydisperse]
    module.opencl = False
    module.single = False
    module.have_Fq = spec.interaction is not None
    if spec.interaction:
        module.radius_effective_modes = ["declared pySEB approximation"]

    # Keep backend selection explicit: a numerical cloud cannot silently
    # become a symbolic model, and a failed compilation must remain visible.
    if backend == "compiled":
        _compiled(module, spec, names, volume_names)
    else:
        # Share the model definition, but never share its mutable World between
        # threads. Complete bindings reset that thread's state before each call.
        local = threading.local()

        def context():
            if not hasattr(local, "context"):
                local.context = spec.world_factory()
            return local.context

        def evaluate(q, args, amplitude=False):
            q = np.asarray(q, dtype=float)
            if q.ndim != 1 or not np.all(np.isfinite(q)) or np.any(q < 0):
                raise ValueError("Expected finite, nonnegative 1-D q in inverse angstroms")
            values = dict(zip(names, map(float, args)))
            ctx = context()
            parameters = dict(spec.bind(ctx, values))
            world = getattr(ctx, "world", ctx)
            if amplitude:
                # Interactions need both F and I. Keep numeric buffers and share
                # one native traversal context per q sample for these queries.
                intensity, f = world._evaluate_form_factor_numpy(
                    spec.structure, parameters, q, spec.interaction.reference)
            else:
                # Intensity-only models need no reference-amplitude query.
                intensity = np.asarray(world.EvaluateFormFactorUnnormalized(
                    spec.structure, parameters, q.tolist()), dtype=float)
            intensity *= spec.intensity_factor
            if intensity.shape != q.shape or not np.all(np.isfinite(intensity)):
                raise ValueError("pySEB returned invalid intensity")
            # Intensity comes from the complete graph independently of F. For
            # general particles, the mean reference amplitude squared is not I.
            if not amplitude:
                return intensity
            if f.shape != q.shape or not np.all(np.isfinite(f)):
                raise ValueError("pySEB returned invalid amplitude")
            return math.sqrt(spec.intensity_factor) * f, intensity

        def iq(q, *args):
            return evaluate(q, args)

        def fq(q, *args):
            return evaluate(q, args, amplitude=True)

        iq.vectorized = fq.vectorized = True
        if spec.interaction:
            module.Fq = fq
        else:
            module.Iq = iq

        def volume(*args):
            return _positive(spec.volume(dict(zip(volume_names, map(float, args)))), "volume") if spec.volume else 1.0

        module.form_volume = volume
        if spec.interaction:
            def radius(mode, *args):
                return _positive(spec.interaction.effective_radius(
                    dict(zip(volume_names, map(float, args)))), "effective radius")
            module.radius_effective = radius
    info = make_info(module)
    # sasmodels infers magnetic parameters from names such as sld_core.
    # This bridge is explicitly nonmagnetic, including its generated kernels.
    from sasmodels.modelinfo import make_parameter_table
    info.parameters = info.base = make_parameter_table(module.parameters, allow_magnetic=False)
    info = _restrict_distributions(info, spec)
    info.pyseb_capabilities = spec.capabilities
    info.pyseb_backend = backend
    return info


def _compiled(module, spec, names, volume_names):
    import sympy as sp
    from sympy.core.function import AppliedUndef
    if spec.symbolic is None:
        raise ValueError("Compiled backend requires explicit SymbolicSpec")
    ctx = spec.world_factory()
    world = getattr(ctx, "world", ctx)
    symbols = {name: sp.Symbol(name) for name in names}
    q = sp.Symbol("q")
    substitutions = {sp.Symbol(str(k)): sp.sympify(v) for k, v in spec.symbolic.parameters.items()}
    try:
        raw = sp.sympify(str(world.FormFactor_Unnormalized(spec.structure)))
    except Exception as exc:
        raise ValueError("This structure cannot provide a compilable analytical expression") from exc
    # Derive the intensity from pySEB itself. The optional signed amplitude is
    # a separate declaration used only when interaction coupling is enabled.
    expressions = [sp.simplify(raw.subs(substitutions, simultaneous=True)) * spec.intensity_factor]
    if spec.interaction:
        if spec.symbolic.amplitude is None or spec.symbolic.effective_radius is None:
            raise ValueError("Compiled interactions require symbolic amplitude and effective radius")
        amplitude = sp.sympify(spec.symbolic.amplitude).subs(substitutions, simultaneous=True)
        amplitude = sp.simplify(amplitude) * math.sqrt(spec.intensity_factor)
        # Use a proven identity to avoid a huge squared expression for spheres.
        if sp.simplify(expressions[0] - amplitude**2) == 0:
            expressions[0] = amplitude**2
        expressions.insert(0, amplitude)
    for expr in expressions:
        if (expr.free_symbols - set(symbols.values()) - {q}
                or expr.has(sp.Integral, sp.Derivative, sp.Sum, AppliedUndef)):
            raise ValueError("Unsupported symbols, callbacks or integrals in compiled expression")
    volume = sp.sympify(spec.symbolic.volume) if spec.volume else sp.Integer(1)
    if spec.volume and spec.symbolic.volume is None:
        raise ValueError("Compiled absolute scaling requires a symbolic volume")
    if volume.free_symbols - {symbols[n] for n in volume_names}:
        raise ValueError("Symbolic volume contains undeclared dependencies")
    declarations = ", ".join("double " + n for n in names)
    declarations = (", " + declarations) if declarations else ""
    # Generate the zero-q limit explicitly, rather than returning 0/0.
    try:
        limits = [sp.simplify(sp.limit(expr, q, 0, dir="+")) for expr in expressions]
        if any(v.has(sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Limit) for v in limits):
            raise ValueError("No finite zero-q limit")
        terms, reduced = sp.cse(expressions)
        body = "\n".join(f"double {s} = {sp.ccode(e)};" for s, e in terms)
        if spec.interaction:
            assignments = f"*F1={sp.ccode(reduced[0])}; *F2={sp.ccode(reduced[1])};"
            at_zero = f"*F1={sp.ccode(limits[0])}; *F2={sp.ccode(limits[1])}; return;"
            module.c_code = f"static void Fq(double q, double *F1, double *F2{declarations}) {{ if(q==0.) {{{at_zero}}} {body} {assignments} }}"
        else:
            module.Iq = f"if(q==0.) return {sp.ccode(limits[0])}; {body} return {sp.ccode(reduced[0])};"
        module.form_volume = "return " + sp.ccode(volume) + ";"
        if spec.interaction:
            radius = sp.sympify(spec.symbolic.effective_radius)
            if radius.free_symbols - {symbols[n] for n in volume_names}:
                raise ValueError("Symbolic effective radius contains undeclared dependencies")
            args = ", ".join("double " + n for n in volume_names)
            module.c_code += f"\nstatic double radius_effective(int mode{', ' if args else ''}{args}) {{return {sp.ccode(radius)};}}"
    except (NotImplementedError, TypeError, ValueError) as exc:
        raise ValueError(f"Unsupported compiled expression: {exc}") from exc
    # Verify declarations against native values before exposing a compiled model.
    values = {p.name: p.default for p in spec.parameters}
    parameters = spec.bind(ctx, values)
    grid = np.array([0., .01, .1])
    expected = np.asarray(world.EvaluateFormFactorUnnormalized(spec.structure, parameters, grid.tolist())) * spec.intensity_factor
    actual = sp.lambdify((q, *symbols.values()), expressions[-1], "numpy")
    # The nonzero points avoid singular forms in the Python expression printer.
    np.testing.assert_allclose(actual(grid[1:], *values.values()), expected[1:], rtol=1e-7, atol=1e-9)
    if spec.volume:
        np.testing.assert_allclose(float(volume.subs(values)), _positive(spec.volume(values), "volume"), rtol=1e-10)
    if spec.interaction:
        expected_f = np.asarray(world.EvaluateFormFactorAmplitudeUnnormalized(spec.interaction.reference, parameters, grid[1:].tolist())) * math.sqrt(spec.intensity_factor)
        actual_f = sp.lambdify((q, *symbols.values()), expressions[0], "numpy")
        np.testing.assert_allclose(actual_f(grid[1:], *values.values()), expected_f, rtol=1e-7, atol=1e-9)


def with_hardsphere(spec, backend="runtime"):
    """sasmodels P@hardsphere with explicit beta/decoupling default.

    sasmodels volfraction controls both interaction strength and intensity.
    For independent phi and eta, callers use scale=phi/eta and volfraction=eta;
    use the particle alone for eta=0. radius_effective_mode=0 permits an
    externally supplied effective radius (e.g. from mean particle volume).
    """
    if spec.interaction is None or spec.volume is None:
        raise ValueError("Hard-sphere coupling requires declared volume and InteractionSpec")
    _sasmodels()
    from sasmodels.core import load_model_info
    from sasmodels.product import make_product_info
    info = make_product_info(make_model_info(spec, backend), load_model_info("hardsphere"))
    from sasmodels.modelinfo import ParameterTable
    info.parameters = info.base = ParameterTable(info.parameters.kernel_parameters, allow_magnetic=False)
    _restrict_distributions(info, spec)
    info.parameters["structure_factor_mode"].default = 1
    info.parameters["radius_effective_mode"].default = 1
    info.parameters.defaults.update(structure_factor_mode=1, radius_effective_mode=1)
    return info


def positive_gaussian():
    """Register and return the smooth positive-Gaussian distribution name."""
    _sasmodels()
    from scipy.special import roots_legendre
    from sasmodels.weights import Dispersion, DISTRIBUTIONS
    name = "pyseb_positive_gaussian"
    if name not in DISTRIBUTIONS:
        class PositiveGaussian(Dispersion):
            type = name

            def _weights(self, center, sigma, lb, ub):
                nodes, weights = roots_legendre(self.npts)
                lower = max(lb, center - self.nsigmas * sigma)
                upper = min(ub, center + self.nsigmas * sigma)
                values = lower + (nodes + 1) * (upper - lower) / 2
                return values, weights * np.exp(-0.5 * ((values-center)/sigma)**2)
        DISTRIBUTIONS[name] = PositiveGaussian
    return name
