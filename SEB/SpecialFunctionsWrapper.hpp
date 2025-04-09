#ifndef INCLUDE_SPECIALFUNCTIONS_WRAPPER_HPP
#define INCLUDE_SPECIALFUNCTIONS_WRAPPER_HPP

#include "Expression.hpp"
#include "GiNaCSymbolic.hpp"

// Define our own versions of the GiNaC macros
#define REGISTER_FUNCTION(name, options) \
    namespace { \
        Expression name##_wrapper(const Expression& x) { \
            return name##_eval(x); \
        } \
    }

#define REGISTER_FUNCTION_2P(name, options) \
    namespace { \
        Expression name##_wrapper(const Expression& x, const Expression& y) { \
            return name##_eval(x, y); \
        } \
    }

// Helper functions for special functions
Expression csc(const Expression& x);
Expression sec(const Expression& x);

inline Expression power(const Expression& x, const Expression& a) {
    return pow(x, a);
}

Expression BesselJ0(const Expression& x);
Expression BesselJ1(const Expression& x);
Expression BesselJ2(const Expression& x);

Expression DawsonF(const Expression& x);
Expression Six(const Expression& x);
Expression Erf(const Expression& x);
Expression Erfc(const Expression& x);

Expression Hypergeometric0F1Regularized(const Expression& a, const Expression& x);
Expression StruveH0(const Expression& x);
Expression StruveH1(const Expression& x);

#endif // INCLUDE_SPECIALFUNCTIONS_WRAPPER_HPP
