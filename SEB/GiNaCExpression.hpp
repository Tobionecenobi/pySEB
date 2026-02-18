#ifndef INCLUDE_GINAC_EXPRESSION_HPP
#define INCLUDE_GINAC_EXPRESSION_HPP

#include "SymbolicInterface.hpp"

// Forward declaration to avoid circular dependencies
class GiNaCExpression;

#ifdef USE_GINAC_IMPL
// The real implementation is in GiNaCSymbolic.hpp when USE_GINAC_IMPL is defined
#include <ginac/ginac.h>

// Forward declaration only - implementation is in GiNaCSymbolic.hpp
class GiNaCFactory;

#else
// Dummy implementation when USE_GINAC_IMPL is not defined
// This is just to make the compiler happy when USE_GINAC_IMPL is OFF

class GiNaCExpression : public SymbolicExpression {
public:
    GiNaCExpression() {}

    // Basic operations - all throw exceptions
    SymExprPtr add(const SymExprPtr& other) const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr subtract(const SymExprPtr& other) const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr multiply(const SymExprPtr& other) const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr divide(const SymExprPtr& other) const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr negate() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr pow(const SymExprPtr& exponent) const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr pow(double exponent) const override {
        throw std::runtime_error("GiNaC not available");
    }

    // Basic functions
    SymExprPtr exp() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr log() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr sin() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr cos() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr tan() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr asin() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr acos() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr atan() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr sinh() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr cosh() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr tanh() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr sqrt() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr abs() const override {
        throw std::runtime_error("GiNaC not available");
    }

    // Special functions for scattering
    SymExprPtr bessel_j0() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr bessel_j1() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr dawson() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr erf() const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr erfc() const override {
        throw std::runtime_error("GiNaC not available");
    }

    // Substitution and evaluation
    SymExprPtr subs(const std::string& symbol, const SymExprPtr& value) const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr subs(const std::string& symbol, double value) const override {
        throw std::runtime_error("GiNaC not available");
    }

    double eval() const override {
        throw std::runtime_error("GiNaC not available");
    }

    bool is_numeric() const override {
        throw std::runtime_error("GiNaC not available");
    }

    // Series expansion
    SymExprPtr series(const SymExprPtr& var, const SymExprPtr& point, int order) const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr series_to_poly(const SymExprPtr& series) const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr coeff(const SymExprPtr& var, int power) const override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr evalf() const override {
        throw std::runtime_error("GiNaC not available");
    }

    bool is_zero() const override {
        throw std::runtime_error("GiNaC not available");
    }

    double to_double() const override {
        throw std::runtime_error("GiNaC not available");
    }

    // String representations
    std::string to_string() const override {
        return "GiNaC not available";
    }

    std::string to_latex() const override {
        return "GiNaC not available";
    }

    std::string to_python() const override {
        return "GiNaC not available";
    }

    std::string to_cform() const override {
        return "GiNaC not available";
    }
};

class GiNaCFactory : public SymbolicFactory {
public:
    GiNaCFactory() {}

    SymExprPtr createSymbol(const std::string& name) override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr createConstant(double value) override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr createPi() override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr createE() override {
        throw std::runtime_error("GiNaC not available");
    }

    SymExprPtr createI() override {
        throw std::runtime_error("GiNaC not available");
    }
};
#endif // USE_GINAC_IMPL

#endif // INCLUDE_GINAC_EXPRESSION_HPP
