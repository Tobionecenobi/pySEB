#ifndef BINDINGS_SYMPY_HPP
#define BINDINGS_SYMPY_HPP

#include <pybind11/pybind11.h>

#include "World.hpp"

namespace py = pybind11;

// Function to register SymPy-specific bindings
void register_sympy_bindings(py::class_<World>& world);

#endif // BINDINGS_SYMPY_HPP
