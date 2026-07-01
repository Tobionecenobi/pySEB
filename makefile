# Convenience wrapper around the canonical CMake + Ninja builds.
#
# Direct CMake commands remain the source of truth; these targets are shortcuts
# for common developer workflows.

GENERATOR ?= Ninja
CMAKE_BUILD_TYPE ?= Release
CORE_BUILD_DIR ?= build-core
GINAC_BUILD_DIR ?= build-ninja-ginac
PYTHON_BUILD_DIR ?= build-python
TEST_BUILD_DIR ?= Tests/build-ninja

.PHONY: all core ginac python tests clean clean-core clean-ginac clean-python clean-tests

all: core

core:
	cmake -S . -B $(CORE_BUILD_DIR) -G "$(GENERATOR)" -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) -DBUILD_PYTHON=OFF -DSEB_ENABLE_GINAC_BACKEND=OFF
	cmake --build $(CORE_BUILD_DIR) --target seb seb-symbolic

ginac:
	cmake -S . -B $(GINAC_BUILD_DIR) -G "$(GENERATOR)" -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) -DBUILD_PYTHON=OFF -DSEB_ENABLE_GINAC_BACKEND=ON
	cmake --build $(GINAC_BUILD_DIR)

python:
	cmake -S . -B $(PYTHON_BUILD_DIR) -G "$(GENERATOR)" -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) -DBUILD_PYTHON=ON -DSEB_ENABLE_GINAC_BACKEND=OFF
	cmake --build $(PYTHON_BUILD_DIR)

tests: ginac
	cmake -S Tests -B $(TEST_BUILD_DIR) -G "$(GENERATOR)" -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) -DSEB_BUILD_DIR=$(abspath $(GINAC_BUILD_DIR))
	cmake --build $(TEST_BUILD_DIR)
	ctest --test-dir $(TEST_BUILD_DIR) --output-on-failure

clean: clean-core clean-ginac clean-python clean-tests

clean-core:
	cmake -E rm -rf $(CORE_BUILD_DIR)

clean-ginac:
	cmake -E rm -rf $(GINAC_BUILD_DIR)

clean-python:
	cmake -E rm -rf $(PYTHON_BUILD_DIR)

clean-tests:
	cmake -E rm -rf $(TEST_BUILD_DIR)
