#ifndef BINDINGS_GINAC_HPP
#define BINDINGS_GINAC_HPP

#include <pybind11/pybind11.h>

#include "World.hpp"

namespace py = pybind11;

// Function to register GiNaC-specific bindings
void register_ginac_bindings(py::class_<World>& world);

#endif // BINDINGS_GINAC_HPP
