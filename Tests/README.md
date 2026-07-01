# Tests

This folder contains unit tests for the project, written using Google Test. The purpose of these tests is to ensure the correctness and reliability of the codebase.

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

## Adding New Tests
1. Create a new test file in this folder.
2. Write test cases using the Google Test framework.
3. Update the `CMakeLists.txt` file to include the new test file in the build process.

## Notes
- Tests should be modular and focus on individual components.
- Avoid modifying the `Examples/` folder directly; instead, refactor code into reusable components if needed.
