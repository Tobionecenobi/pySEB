"""Minimal smoke test for an installed pyseb wheel."""

import math
from importlib import resources
from pathlib import Path
import tempfile

import pyseb


def main():
    backends = pyseb.available_backends()
    if "portable" not in backends:
        raise AssertionError(f"portable backend missing from {backends!r}")

    pyseb.set_backend("portable")

    models = resources.files("pyseb").joinpath("models")
    point_model = models.joinpath("Point.pyseb.yaml")
    gaussian_model = models.joinpath("GaussianPolymer.pyseb.yaml")
    if not point_model.is_file() or not gaussian_model.is_file():
        raise AssertionError("bundled .pyseb.yaml model files are missing")

    definition = pyseb.load_subunit_definition(str(gaussian_model))
    if definition.id != "pyseb/GaussianPolymer":
        raise AssertionError(f"unexpected bundled model ID: {definition.id!r}")
    report = pyseb.validate_subunit_file(str(gaussian_model))
    if not report.ok:
        raise AssertionError(f"bundled model validation failed: {report.failures!r}")

    world = pyseb.World()
    graph_id = world.add_subunit_file(str(gaussian_model), "poly1")
    world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")
    world.Add(graph_id, "diblockcopolymer")

    form_factor = world.FormFactor("diblockcopolymer")
    params = {
        "beta_poly1": 1.0,
        "beta_poly2": 1.0,
        "Rg_poly1": 5.0,
        "Rg_poly2": 7.0,
    }
    value = pyseb.evaluate_expression(world, form_factor, params, 0.1)
    if not math.isfinite(value):
        raise AssertionError(f"non-finite form factor value: {value!r}")

    x = pyseb.symbol("x")
    expression = pyseb.sin(x * 2.0 + 3.0)
    sympy_expression = pyseb.to_sympy(expression)
    sympy_value = float(sympy_expression.subs("x", 4.0).evalf())

    if not math.isclose(sympy_value, math.sin(11.0), rel_tol=1e-12, abs_tol=1e-12):
        raise AssertionError(f"unexpected SymPy conversion value: {sympy_value!r}")

    with tempfile.TemporaryDirectory() as directory:
        invalid = Path(directory, "invalid.pyseb.yaml")
        invalid.write_text("format: wrong\n", encoding="utf-8")
        try:
            pyseb.load_subunit_definition(str(invalid))
        except RuntimeError:
            pass
        else:
            raise AssertionError("invalid subunit definition was accepted")

    print("pyseb installed wheel smoke test passed")


if __name__ == "__main__":
    main()
