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
    model_files = sorted(
        (path for path in models.iterdir() if path.name.endswith(".pyseb.yaml")),
        key=lambda path: path.name,
    )
    if not model_files:
        raise AssertionError("bundled .pyseb.yaml model files are missing")

    catalogue_world = pyseb.World("catalogue")
    catalogue = {
        model.id: model
        for model in catalogue_world.list_subunit_models()
        if model.source.startswith("<bundled:")
    }
    if len(catalogue) != len(model_files):
        raise AssertionError(
            f"packaged/embedded model count differs: {len(model_files)} != {len(catalogue)}"
        )
    definitions = []
    for index, model_file in enumerate(model_files):
        definition = pyseb.load_subunit_definition(str(model_file))
        definitions.append(definition)
        report = pyseb.validate_subunit_file(str(model_file))
        if not report.ok:
            raise AssertionError(
                f"bundled model validation failed for {model_file.name}: {report.failures!r}"
            )
        if definition.id not in catalogue:
            raise AssertionError(f"bundled model is not registered: {definition.id!r}")
        if catalogue[definition.id].api_name != definition.api_name:
            raise AssertionError(f"wrong API name for {definition.id!r}")
        catalogue_world.Add(definition.id, f"canonical{index}")
        catalogue_world.Add(definition.api_name, f"alias{index}")

    definitions_by_id = {definition.id: definition for definition in definitions}
    required_models = {"pyseb/Point", "pyseb/GaussianPolymer"}
    missing_models = required_models.difference(definitions_by_id)
    if missing_models:
        raise AssertionError(f"required bundled models are missing: {sorted(missing_models)!r}")
    gaussian_model = models.joinpath("GaussianPolymer.pyseb.yaml")

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
