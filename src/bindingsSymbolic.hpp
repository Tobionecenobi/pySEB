#ifndef BINDINGS_SYMBOLIC_HPP
#define BINDINGS_SYMBOLIC_HPP

#include <pybind11/pybind11.h>

#include "Expression.hpp"
#include "World.hpp"

namespace py = pybind11;

py::object pyseb_expression_to_sympy(const Expression& expr);
py::object pyseb_zero_sympy_expression();
void register_symbolic_world_bindings(py::class_<World>& world);

#endif // BINDINGS_SYMBOLIC_HPP
