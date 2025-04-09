#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "SymbolicInterface.hpp"
#ifdef USE_GINAC
#include "GiNaCSymbolic.hpp"
#endif

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
}
