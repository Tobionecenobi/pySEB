#include "SpecialFunctionsWrapper.hpp"
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_erf.h>

// Implementation of the special functions

Expression csc(const Expression& x) {
    return Expression(SymbolicExpression::constant(1)) / x.sin();
}

Expression sec(const Expression& x) {
    return Expression(SymbolicExpression::constant(1)) / x.cos();
}

Expression BesselJ0_wrapper(const Expression& x) {
    if (x.is_numeric()) {
        return Expression(SymbolicExpression::constant(gsl_sf_bessel_J0(x.eval())));
    }
    return Expression(x.bessel_j0());
}

// Ensure BesselJ0 is available regardless of symbolic backend
Expression BesselJ0(const Expression& x) {
    return BesselJ0_wrapper(x);
}

Expression BesselJ1_wrapper(const Expression& x) {
    if (x.is_numeric()) {
        return Expression(SymbolicExpression::constant(gsl_sf_bessel_J1(x.eval())));
    }
    return Expression(x.bessel_j1());
}

// Ensure BesselJ1 is available regardless of symbolic backend
Expression BesselJ1(const Expression& x) {
    return BesselJ1_wrapper(x);
}

Expression BesselJ2(const Expression& x) {
    return (Expression(SymbolicExpression::constant(2.0)) * x.bessel_j1() - x.bessel_j0()) / x;
}

Expression DawsonF_wrapper(const Expression& x) {
    if (x.is_numeric()) {
        return Expression(SymbolicExpression::constant(gsl_sf_dawson(x.eval())));
    }
    return Expression(x.dawson());
}

// Ensure DawsonF is available regardless of symbolic backend
Expression DawsonF(const Expression& x) {
    return DawsonF_wrapper(x);
}

Expression Six(const Expression& x) {
    return x.sin() / x;
}

Expression Erf_wrapper(const Expression& x) {
    if (x.is_numeric()) {
        return Expression(SymbolicExpression::constant(gsl_sf_erf(x.eval())));
    }
    return Expression(x.erf());
}

// Ensure Erf is available regardless of symbolic backend
Expression Erf(const Expression& x) {
    return Erf_wrapper(x);
}

Expression Erfc_wrapper(const Expression& x) {
    if (x.is_numeric()) {
        return Expression(SymbolicExpression::constant(gsl_sf_erfc(x.eval())));
    }
    return Expression(x.erfc());
}

// Ensure Erfc is available regardless of symbolic backend
Expression Erfc(const Expression& x) {
    return Erfc_wrapper(x);
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
