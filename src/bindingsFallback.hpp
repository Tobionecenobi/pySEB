#ifndef BINDINGS_FALLBACK_HPP
#define BINDINGS_FALLBACK_HPP

#include <pybind11/pybind11.h>

#include "World.hpp"

namespace py = pybind11;

// Function to register fallback bindings (string representation)
void register_fallback_bindings(py::class_<World>& world);

#endif // BINDINGS_FALLBACK_HPP
