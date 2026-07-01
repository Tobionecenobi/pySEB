#ifndef INCLUDE_EXPRESSION_FUNCTIONS_HPP
#define INCLUDE_EXPRESSION_FUNCTIONS_HPP

#include "Expression.hpp"

// Global functions that work with Expression objects

// Trigonometric functions
inline Expression sin(const Expression& x) {
    return x.sin();
}

// Make sure the sin function is visible in the global namespace
using std::sin;

inline Expression cos(const Expression& x) {
    return x.cos();
}

inline Expression tan(const Expression& x) {
    return x.tan();
}

inline Expression asin(const Expression& x) {
    return x.asin();
}

inline Expression acos(const Expression& x) {
    return x.acos();
}

inline Expression atan(const Expression& x) {
    return x.atan();
}

// Hyperbolic functions
inline Expression sinh(const Expression& x) {
    return x.sinh();
}

inline Expression cosh(const Expression& x) {
    return x.cosh();
}

inline Expression tanh(const Expression& x) {
    return x.tanh();
}

// These functions are not implemented in the Expression class yet
inline Expression asinh(const Expression& x) {
    // Placeholder implementation
    return x;
}

inline Expression acosh(const Expression& x) {
    // Placeholder implementation
    return x;
}

inline Expression atanh(const Expression& x) {
    // Placeholder implementation
    return x;
}

// Exponential and logarithmic functions
inline Expression exp(const Expression& x) {
    return x.exp();
}

// Make sure the exp function is visible in the global namespace
using std::exp;

inline Expression log(const Expression& x) {
    return x.log();
}

inline Expression log10(const Expression& x) {
    // Placeholder implementation
    return x.log() / Expression(SymbolicExpression::constant(std::log(10.0)));
}

// Power functions
// These are already defined in Expression.hpp, so we don't need to redefine them here

inline Expression sqrt(const Expression& x) {
    return x.sqrt();
}

// Bessel functions
inline Expression bessel_j0(const Expression& x) {
    return x.bessel_j0();
}

inline Expression bessel_j1(const Expression& x) {
    return x.bessel_j1();
}

// Error functions
inline Expression erf(const Expression& x) {
    return x.erf();
}

inline Expression erfc(const Expression& x) {
    return x.erfc();
}

// Constants as functions to avoid static initialization order issues
inline Expression Pi() {
    return Expression(SymbolicExpression::pi());
}

inline Expression SMALL() {
    return Expression(SymbolicExpression::constant(1e-6));
}

// Integration function (renamed from integral to avoid C++20 conflict)
inline Expression integrate(const Expression& var, const Expression& a, const Expression& b, const Expression& integrand) {
    // This is a placeholder implementation
    // In a real implementation, you would perform numerical integration
    return integrand * (b - a);
}

// Overloads for integrate function to handle integers
inline Expression integrate(const Expression& var, int a, const Expression& b, const Expression& integrand) {
    return integrate(var, Expression(SymbolicExpression::constant(a)), b, integrand);
}

inline Expression integrate(const Expression& var, const Expression& a, int b, const Expression& integrand) {
    return integrate(var, a, Expression(SymbolicExpression::constant(b)), integrand);
}

inline Expression integrate(const Expression& var, int a, int b, const Expression& integrand) {
    return integrate(var, Expression(SymbolicExpression::constant(a)), Expression(SymbolicExpression::constant(b)), integrand);
}

// Hypergeometric function
inline Expression Hypergeometric0F1Regularized(int a, const Expression& x) {
    // This is a placeholder implementation
    return Expression(SymbolicExpression::constant(1.0));
}

// GiNaC compatibility functions
inline bool is_a_numeric(const Expression& x) {
    return x.is_numeric();
}

inline double ex_to_numeric(const Expression& x) {
    return x.to_numeric();
}

#endif // INCLUDE_EXPRESSION_FUNCTIONS_HPP
