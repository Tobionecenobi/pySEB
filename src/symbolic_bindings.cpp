#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "SymbolicInterface.hpp"
#include "GiNaCSymbolic.hpp"

namespace py = pybind11;

// Forward declaration of the Python symbolic expression class
class PySymbolicExpression : public SymbolicExpression {
public:
    using SymbolicExpression::SymbolicExpression;

    // Basic operations
    SymExprPtr add(const SymExprPtr& other) const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, add, other);
    }

    SymExprPtr subtract(const SymExprPtr& other) const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, subtract, other);
    }

    SymExprPtr multiply(const SymExprPtr& other) const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, multiply, other);
    }

    SymExprPtr divide(const SymExprPtr& other) const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, divide, other);
    }

    SymExprPtr negate() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, negate, );
    }

    SymExprPtr pow(const SymExprPtr& exponent) const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, pow, exponent);
    }

    // Functions
    SymExprPtr exp() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, exp, );
    }

    SymExprPtr log() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, log, );
    }

    SymExprPtr sin() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, sin, );
    }

    SymExprPtr cos() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, cos, );
    }

    SymExprPtr tan() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, tan, );
    }

    SymExprPtr asin() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, asin, );
    }

    SymExprPtr acos() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, acos, );
    }

    SymExprPtr atan() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, atan, );
    }

    SymExprPtr sinh() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, sinh, );
    }

    SymExprPtr cosh() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, cosh, );
    }

    SymExprPtr tanh() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, tanh, );
    }

    SymExprPtr sqrt() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, sqrt, );
    }

    SymExprPtr abs() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, abs, );
    }

    // Special functions needed for scattering
    SymExprPtr bessel_j0() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, bessel_j0, );
    }

    SymExprPtr bessel_j1() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, bessel_j1, );
    }

    SymExprPtr dawson() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, dawson, );
    }

    SymExprPtr erf() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, erf, );
    }

    SymExprPtr erfc() const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, erfc, );
    }

    // Substitution and evaluation
    SymExprPtr subs(const std::string& symbol, const SymExprPtr& value) const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, subs, symbol, value);
    }

    SymExprPtr subs(const ParameterMap& params) const override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, subs, params);
    }

    double eval() const override {
        PYBIND11_OVERRIDE_PURE(double, SymbolicExpression, eval, );
    }

    bool is_numeric() const override {
        PYBIND11_OVERRIDE_PURE(bool, SymbolicExpression, is_numeric, );
    }

    // String representation
    std::string to_string() const override {
        PYBIND11_OVERRIDE_PURE(std::string, SymbolicExpression, to_string, );
    }

    std::string to_latex() const override {
        PYBIND11_OVERRIDE_PURE(std::string, SymbolicExpression, to_latex, );
    }

    std::string to_python() const override {
        PYBIND11_OVERRIDE_PURE(std::string, SymbolicExpression, to_python, );
    }

    std::string to_cform() const override {
        PYBIND11_OVERRIDE_PURE(std::string, SymbolicExpression, to_cform, );
    }
};

// Forward declaration of the Python symbolic factory class
class PySymbolicFactory : public SymbolicFactory {
public:
    using SymbolicFactory::SymbolicFactory;

    SymExprPtr createSymbol(const std::string& name) override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicFactory, createSymbol, name);
    }

    SymExprPtr createConstant(double value) override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicFactory, createConstant, value);
    }

    SymExprPtr createPi() override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicFactory, createPi, );
    }

    SymExprPtr createE() override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicFactory, createE, );
    }

    SymExprPtr createI() override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicFactory, createI, );
    }
};

void init_symbolic(py::module& m) {
    // Bind the SymbolicExpression class
    py::class_<SymbolicExpression, PySymbolicExpression, std::shared_ptr<SymbolicExpression>>(m, "SymbolicExpression")
        .def(py::init<>())
        .def("add", &SymbolicExpression::add)
        .def("subtract", &SymbolicExpression::subtract)
        .def("multiply", &SymbolicExpression::multiply)
        .def("divide", &SymbolicExpression::divide)
        .def("negate", &SymbolicExpression::negate)
        .def("pow", &SymbolicExpression::pow)
        .def("exp", &SymbolicExpression::exp)
        .def("log", &SymbolicExpression::log)
        .def("sin", &SymbolicExpression::sin)
        .def("cos", &SymbolicExpression::cos)
        .def("tan", &SymbolicExpression::tan)
        .def("asin", &SymbolicExpression::asin)
        .def("acos", &SymbolicExpression::acos)
        .def("atan", &SymbolicExpression::atan)
        .def("sinh", &SymbolicExpression::sinh)
        .def("cosh", &SymbolicExpression::cosh)
        .def("tanh", &SymbolicExpression::tanh)
        .def("sqrt", &SymbolicExpression::sqrt)
        .def("abs", &SymbolicExpression::abs)
        .def("bessel_j0", &SymbolicExpression::bessel_j0)
        .def("bessel_j1", &SymbolicExpression::bessel_j1)
        .def("dawson", &SymbolicExpression::dawson)
        .def("erf", &SymbolicExpression::erf)
        .def("erfc", &SymbolicExpression::erfc)
        .def("subs", py::overload_cast<const std::string&, const SymExprPtr&>(&SymbolicExpression::subs, py::const_))
        .def("subs", py::overload_cast<const ParameterMap&>(&SymbolicExpression::subs, py::const_))
        .def("eval", &SymbolicExpression::eval)
        .def("is_numeric", &SymbolicExpression::is_numeric)
        .def("to_string", &SymbolicExpression::to_string)
        .def("to_latex", &SymbolicExpression::to_latex)
        .def("to_python", &SymbolicExpression::to_python)
        .def("to_cform", &SymbolicExpression::to_cform)
        .def("__str__", &SymbolicExpression::to_string)
        .def("__repr__", &SymbolicExpression::to_string);

    // Bind the SymbolicFactory class
    py::class_<SymbolicFactory, PySymbolicFactory>(m, "SymbolicFactory")
        .def(py::init<>())
        .def("createSymbol", &SymbolicFactory::createSymbol)
        .def("createConstant", &SymbolicFactory::createConstant)
        .def("createPi", &SymbolicFactory::createPi)
        .def("createE", &SymbolicFactory::createE)
        .def("createI", &SymbolicFactory::createI)
        .def_static("instance", &SymbolicFactory::instance)
        .def_static("setInstance", &SymbolicFactory::setInstance);

    // Bind the GiNaCExpression class
    py::class_<GiNaCExpression, SymbolicExpression, std::shared_ptr<GiNaCExpression>>(m, "GiNaCExpression")
        .def(py::init<const GiNaC::ex&>());

    // Bind the GiNaCFactory class
    py::class_<GiNaCFactory, SymbolicFactory>(m, "GiNaCFactory")
        .def(py::init<>());

    // Bind helper functions
    m.def("symbol", &symbol);
    m.def("constant", &constant);
    m.def("pi", &pi);
    m.def("e", &e);
    m.def("i", &i);
}
