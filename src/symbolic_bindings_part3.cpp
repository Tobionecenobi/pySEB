// Special functions needed for scattering
SymExprPtr PySymbolicExpression::bessel_j0() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, bessel_j0, );
}

SymExprPtr PySymbolicExpression::bessel_j1() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, bessel_j1, );
}

SymExprPtr PySymbolicExpression::dawson() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, dawson, );
}

SymExprPtr PySymbolicExpression::erf() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, erf, );
}

SymExprPtr PySymbolicExpression::erfc() const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, erfc, );
}

// Substitution and evaluation
SymExprPtr PySymbolicExpression::subs(const std::string& symbol, const SymExprPtr& value) const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, subs, symbol, value);
}

SymExprPtr PySymbolicExpression::subs(const ParameterMap& params) const override {
    PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicExpression, subs, params);
}

double PySymbolicExpression::eval() const override {
    PYBIND11_OVERRIDE_PURE(double, SymbolicExpression, eval, );
}

bool PySymbolicExpression::is_numeric() const override {
    PYBIND11_OVERRIDE_PURE(bool, SymbolicExpression, is_numeric, );
}
