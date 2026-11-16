# Tests

This folder contains the C++ Google Test suite and the Python unittest suite
in `python/`. Shared reference data lives in `data/`.

## Structure
- Each test file corresponds to a specific module or functionality in the project.
- Test cases are written to validate the behavior of individual components.

## Running Tests
1. Ensure Google Test and the optional GiNaC dependencies are available.
2. Build and run the tests from the repository root:
   ```
   make tests
   ```

Or run the CTest build manually:
   ```
   cmake -S . -B build-ninja-ginac -G Ninja -DBUILD_PYTHON=OFF -DSEB_ENABLE_GINAC_BACKEND=ON
   cmake --build build-ninja-ginac
   cmake -S Tests -B Tests/build-ninja -G Ninja -DSEB_BUILD_DIR=$PWD/build-ninja-ginac
   cmake --build Tests/build-ninja
   ctest --test-dir Tests/build-ninja --output-on-failure
   ```

## Python tests

After building or installing the Python extension, run from the repository root:

```sh
PYTHONPATH="$PWD/build-python/pyseb:$PWD" python -m unittest discover -s Tests/python -p 'test_*.py'
```

The CTest commands above run the C++ suite; the Python suite runs separately.
Optional integrations may require additional dependencies. Tests locate models,
fixtures, documentation, and `Examples/python/` relative to the repository root.

## Geometric sampling

`sampling/` contains the C++ geometric sampling drivers and shared helpers.
The `Sample_*_test.cpp` wrappers exercise them through CTest. Existing sampled
reference curves are preserved in `data/sampling/`. The related derivation
notebooks and PDFs live under `docs/derivations/mathematica/`.

## Reference-data validation

`validation_test` compares the legacy subunit expressions and their Guinier
limits with the numerical data in `Examples/Validation/`. A failed comparison
or exception fails the test. Reference files must contain two finite numbers
per nonblank line (q and the expected value); missing, empty, and malformed
files are errors. Regression tests exercise these failure paths.

New bundled YAML models are covered automatically by
`file_defined_subunit_test`, which runs each model's `validation.cases`.
Expected reference values must be supplied in the YAML; a model without cases
currently produces a warning rather than a failure.

## Adding New Tests
1. Create a new test file in this folder.
2. Write test cases using the Google Test framework.
3. Update the `CMakeLists.txt` file to include the new test file in the build process.

## Notes
- Tests should be modular and focus on individual components.
- Avoid modifying the `Examples/` folder directly; instead, refactor code into reusable components if needed.
