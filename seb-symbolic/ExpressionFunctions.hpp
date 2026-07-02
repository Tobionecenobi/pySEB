#ifndef INCLUDE_EXPRESSION_FUNCTIONS_HPP
#define INCLUDE_EXPRESSION_FUNCTIONS_HPP

#include "Expression.hpp"
#include "SymbolicPortable.hpp"

#ifdef USE_GINAC_IMPL
#include "GiNaCSymbolic.hpp"
#include <ginac/integral.h>
#endif

#include <cmath>
#include <stdexcept>
#include <vector>

Expression Hypergeometric0F1Regularized(const Expression& a, const Expression& x);

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
    (void)x;
    throw std::runtime_error("asinh is not part of the SEB expression construction interface");
}

inline Expression acosh(const Expression& x) {
    (void)x;
    throw std::runtime_error("acosh is not part of the SEB expression construction interface");
}

inline Expression atanh(const Expression& x) {
    (void)x;
    throw std::runtime_error("atanh is not part of the SEB expression construction interface");
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
#ifdef USE_GINAC_IMPL
    if (SymbolicFactory::activeBackendName() == "ginac") {
        auto g_var = std::dynamic_pointer_cast<GiNaCExpression>(var.get());
        auto g_a = std::dynamic_pointer_cast<GiNaCExpression>(a.get());
        auto g_b = std::dynamic_pointer_cast<GiNaCExpression>(b.get());
        auto g_integrand = std::dynamic_pointer_cast<GiNaCExpression>(integrand.get());
        if (!g_var || !g_a || !g_b || !g_integrand) {
            throw std::runtime_error("GiNaC integral requires GiNaC expressions");
        }
        return Expression(std::make_shared<GiNaCExpression>(
            GiNaC::integral(g_var->get_ginac_expr(),
                            g_a->get_ginac_expr(),
                            g_b->get_ginac_expr(),
                            g_integrand->get_ginac_expr())));
    }
#endif

    if (SymbolicFactory::activeBackendName() != "portable") {
        throw std::runtime_error("symbolic backend '" + SymbolicFactory::activeBackendName() +
                                 "' does not support unevaluated integrals");
    }

    return Expression(std::make_shared<sebsym::PortableExpression>(
        sebsym::PortableExpression::Kind::Function,
        "Integral",
        std::vector<SymExprPtr>{var.get(), a.get(), b.get(), integrand.get()}));
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
    return Hypergeometric0F1Regularized(Expression(SymbolicExpression::constant(a)), x);
}

// GiNaC compatibility functions
inline bool is_a_numeric(const Expression& x) {
    return x.is_numeric();
}

inline double ex_to_numeric(const Expression& x) {
    return x.to_numeric();
}

#endif // INCLUDE_EXPRESSION_FUNCTIONS_HPP
