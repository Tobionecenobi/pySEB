#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>

#include "SEB.hpp"
#include "World.hpp"
#include "Subunit.hpp"
#include "Subunits/Subunits.hpp"
#include "Expression.hpp"
#include "SymbolicInterface.hpp"
#include "bindingsSymbolic.hpp"

namespace py = pybind11;

// Define constants if not already defined
#ifndef GENERIC
#define GENERIC 0
#endif

#ifndef XVAR
#define XVAR 1
#endif

#ifndef QVAR
#define QVAR 2
#endif

// This function registers all the common types that are used in multiple binding files
void register_common_types(py::module& m) {
    // refPoint/subName/structName are typedefs of std::string, and GraphID is typedef int.
    // Python should pass str/int directly for these aliases.

    py::class_<Expression>(m, "Expression")
        .def(py::init<>())
        .def(py::init<double>())
        .def("__str__", &Expression::to_string)
        .def("__repr__", [](const Expression& self) {
            return "Expression(" + self.to_string() + ")";
        })
        .def(py::self + py::self)
        .def(py::self + double())
        .def(double() + py::self)
        .def(py::self - py::self)
        .def(py::self - double())
        .def(double() - py::self)
        .def(py::self * py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / py::self)
        .def(py::self / double())
        .def(double() / py::self)
        .def(-py::self)
        .def("__pow__", [](const Expression& self, const Expression& exponent) {
            return sebsym::pow(self, exponent);
        }, py::is_operator())
        .def("exp", &Expression::exp)
        .def("log", &Expression::log)
        .def("sin", &Expression::sin)
        .def("cos", &Expression::cos)
        .def("tan", &Expression::tan)
        .def("asin", &Expression::asin)
        .def("acos", &Expression::acos)
        .def("atan", &Expression::atan)
        .def("sinh", &Expression::sinh)
        .def("cosh", &Expression::cosh)
        .def("tanh", &Expression::tanh)
        .def("sqrt", &Expression::sqrt)
        .def("abs", &Expression::abs)
        .def("bessel_j0", &Expression::bessel_j0)
        .def("bessel_j1", &Expression::bessel_j1)
        .def("dawson", &Expression::dawson)
        .def("erf", &Expression::erf)
        .def("erfc", &Expression::erfc)
        .def("subs", py::overload_cast<const std::string&, double>(&Expression::subs, py::const_),
             py::arg("symbol"), py::arg("value"))
        .def("subs", py::overload_cast<const std::string&, const Expression&>(&Expression::subs, py::const_),
             py::arg("symbol"), py::arg("value"))
        .def("eval", &Expression::eval)
        .def("evalf", &Expression::evalf)
        .def("is_numeric", &Expression::is_numeric)
        .def("to_double", &Expression::to_double)
        .def("to_string", &Expression::to_string)
        .def("to_latex", &Expression::to_latex)
        .def("to_python", &Expression::to_python)
        .def("to_cform", &Expression::to_cform);

    py::implicitly_convertible<double, Expression>();

    m.def("symbol", &symbol, py::arg("name"));
    m.def("constant", &constant, py::arg("value"));
    m.def("pi", &pi);
    m.def("e", &e);
    m.def("i", &i);
    m.def("pow", static_cast<Expression (*)(const Expression&, const Expression&)>(&sebsym::pow),
          py::arg("base"), py::arg("exponent"));
    m.def("sin", [](const Expression& expr) { return expr.sin(); }, py::arg("expr"));
    m.def("cos", [](const Expression& expr) { return expr.cos(); }, py::arg("expr"));
    m.def("tan", [](const Expression& expr) { return expr.tan(); }, py::arg("expr"));
    m.def("exp", [](const Expression& expr) { return expr.exp(); }, py::arg("expr"));
    m.def("log", [](const Expression& expr) { return expr.log(); }, py::arg("expr"));
    m.def("sqrt", [](const Expression& expr) { return expr.sqrt(); }, py::arg("expr"));
    m.def("abs", [](const Expression& expr) { return expr.abs(); }, py::arg("expr"));
    m.def("to_sympy", &pyseb_expression_to_sympy, py::arg("expr"));

    // Expose SubUnit class
    py::class_<SubUnit>(m, "SubUnit")
        .def("getType", &SubUnit::getType)
        .def("getName", &SubUnit::getName)
        .def("getTag", &SubUnit::getTag)
        .def("getReferencePoints", &SubUnit::getReferencePoints);

    // Expose constants
    m.attr("GENERIC") = GENERIC;
    m.attr("XVAR") = XVAR;
    m.attr("QVAR") = QVAR;
    m.attr("WORLDMAXDEPTH") = WORLDMAXDEPTH;
}
