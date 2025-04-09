#include "SpecialFunctionsWrapper.hpp"
#include <cmath>

// Implementation of the special functions

Expression csc(const Expression& x) {
    return Expression(SymbolicExpression::constant(1)) / x.sin();
}

Expression sec(const Expression& x) {
    return Expression(SymbolicExpression::constant(1)) / x.cos();
}

Expression BesselJ0(const Expression& x) {
    return x.bessel_j0();
}

Expression BesselJ1(const Expression& x) {
    return x.bessel_j1();
}

Expression BesselJ2(const Expression& x) {
    return (Expression(SymbolicExpression::constant(2.0)) * x.bessel_j1() - x.bessel_j0()) / x;
}

Expression DawsonF(const Expression& x) {
    return x.dawson();
}

Expression Six(const Expression& x) {
    return x.sin() / x;
}

Expression Erf(const Expression& x) {
    return x.erf();
}

Expression Erfc(const Expression& x) {
    return x.erfc();
}

Expression Hypergeometric0F1Regularized(const Expression& a, const Expression& x) {
    // This is a placeholder implementation
    return Expression(SymbolicExpression::constant(1.0));
}

Expression StruveH0(const Expression& x) {
    // This is a placeholder implementation
    return Expression(SymbolicExpression::constant(1.0));
}

Expression StruveH1(const Expression& x) {
    // This is a placeholder implementation
    return Expression(SymbolicExpression::constant(1.0));
}
