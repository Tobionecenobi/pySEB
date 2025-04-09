// Functions
SymExprPtr PySymbolicExpression::exp() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, exp, );
}

SymExprPtr PySymbolicExpression::log() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, log, );
}

SymExprPtr PySymbolicExpression::sin() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, sin, );
}

SymExprPtr PySymbolicExpression::cos() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, cos, );
}

SymExprPtr PySymbolicExpression::tan() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, tan, );
}

SymExprPtr PySymbolicExpression::asin() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, asin, );
}

SymExprPtr PySymbolicExpression::acos() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, acos, );
}

SymExprPtr PySymbolicExpression::atan() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, atan, );
}

SymExprPtr PySymbolicExpression::sinh() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, sinh, );
}

SymExprPtr PySymbolicExpression::cosh() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, cosh, );
}

SymExprPtr PySymbolicExpression::tanh() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, tanh, );
}

SymExprPtr PySymbolicExpression::sqrt() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, sqrt, );
}

SymExprPtr PySymbolicExpression::abs() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, abs, );
}
