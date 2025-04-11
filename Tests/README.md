# Tests

This folder contains unit tests for the project, written using Google Test. The purpose of these tests is to ensure the correctness and reliability of the codebase.

## Structure
- Each test file corresponds to a specific module or functionality in the project.
- Test cases are written to validate the behavior of individual components.

## Running Tests
1. Ensure Google Test is set up in the project.
2. Build the project using CMake or the provided build system.
3. Run the tests using the following command:
   ```
   ctest
   ```

## Adding New Tests
1. Create a new test file in this folder.
2. Write test cases using the Google Test framework.
3. Update the `CMakeLists.txt` file to include the new test file in the build process.

## Notes
- Tests should be modular and focus on individual components.
- Avoid modifying the `Examples/` folder directly; instead, refactor code into reusable components if needed.