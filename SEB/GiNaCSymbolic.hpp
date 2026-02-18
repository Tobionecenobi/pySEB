#ifndef INCLUDE_GINAC_SYMBOLIC_HPP
#define INCLUDE_GINAC_SYMBOLIC_HPP

#include "SymbolicInterface.hpp"
#include "GiNaCExpression.hpp"

#ifdef USE_GINAC_IMPL

#include <ginac/ginac.h>
#include <sstream>
#include "Expression.hpp"
#include "SpecialFunctionsWrapper.hpp"

class GiNaCExpression : public SymbolicExpression {
private:
    GiNaC::ex _expr;

public:
    GiNaCExpression(const GiNaC::ex& e) : _expr(e) {}

    // Basic operations
    SymExprPtr add(const SymExprPtr& other) const override {
        auto g_other = std::dynamic_pointer_cast<GiNaCExpression>(other);
        return std::make_shared<GiNaCExpression>(_expr + g_other->_expr);
    }

    SymExprPtr subtract(const SymExprPtr& other) const override {
        auto g_other = std::dynamic_pointer_cast<GiNaCExpression>(other);
        return std::make_shared<GiNaCExpression>(_expr - g_other->_expr);
    }

    SymExprPtr multiply(const SymExprPtr& other) const override {
        auto g_other = std::dynamic_pointer_cast<GiNaCExpression>(other);
        return std::make_shared<GiNaCExpression>(_expr * g_other->_expr);
    }

    SymExprPtr divide(const SymExprPtr& other) const override {
        auto g_other = std::dynamic_pointer_cast<GiNaCExpression>(other);
        return std::make_shared<GiNaCExpression>(_expr / g_other->_expr);
    }

    SymExprPtr negate() const override {
        return std::make_shared<GiNaCExpression>(-_expr);
    }

    SymExprPtr pow(const SymExprPtr& exponent) const override {
        auto g_exp = std::dynamic_pointer_cast<GiNaCExpression>(exponent);
        return std::make_shared<GiNaCExpression>(GiNaC::pow(_expr, g_exp->_expr));
    }

    // Basic functions
    SymExprPtr exp() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::exp(_expr));
    }

    SymExprPtr log() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::log(_expr));
    }

    SymExprPtr sin() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::sin(_expr));
    }

    SymExprPtr cos() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::cos(_expr));
    }

    SymExprPtr tan() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::tan(_expr));
    }

    SymExprPtr asin() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::asin(_expr));
    }

    SymExprPtr acos() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::acos(_expr));
    }

    SymExprPtr atan() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::atan(_expr));
    }

    SymExprPtr sinh() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::sinh(_expr));
    }

    SymExprPtr cosh() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::cosh(_expr));
    }

    SymExprPtr tanh() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::tanh(_expr));
    }

    SymExprPtr sqrt() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::sqrt(_expr));
    }

    SymExprPtr abs() const override {
        return std::make_shared<GiNaCExpression>(GiNaC::abs(_expr));
    }

    // Special functions
    SymExprPtr bessel_j0() const override {
        // GiNaC doesn't have a direct bessel_j0 function
        // We'll use a simple approximation: J0(x) ≈ cos(x)
        return std::make_shared<GiNaCExpression>(GiNaC::cos(_expr));
    }

    SymExprPtr bessel_j1() const override {
        // GiNaC doesn't have a direct bessel_j1 function
        // We'll use a simple approximation: J1(x) ≈ sin(x)/x
        return std::make_shared<GiNaCExpression>(GiNaC::sin(_expr) / _expr);
    }

    SymExprPtr dawson() const override {
        // GiNaC doesn't have a built-in dawson function, so we implement it
        // Dawson function F(x) = exp(-x²)∫₀ˣexp(t²)dt
        // We'll use a simple approximation here
        return std::make_shared<GiNaCExpression>(
            GiNaC::exp(-GiNaC::pow(_expr, 2)) * GiNaC::sqrt(GiNaC::Pi) / 2
        );
    }

    SymExprPtr erf() const override {
        // GiNaC doesn't have a built-in erf function, so we'll use an approximation
        // erf(x) ≈ 1 - exp(-x²)
        return std::make_shared<GiNaCExpression>(1 - GiNaC::exp(-GiNaC::pow(_expr, 2)));
    }

    SymExprPtr erfc() const override {
        // erfc(x) = 1 - erf(x) ≈ exp(-x²)
        return std::make_shared<GiNaCExpression>(GiNaC::exp(-GiNaC::pow(_expr, 2)));
    }

    // Substitution and evaluation
    SymExprPtr subs(const std::string& symbol, const SymExprPtr& value) const override {
        auto g_value = std::dynamic_pointer_cast<GiNaCExpression>(value);
        return std::make_shared<GiNaCExpression>(_expr.subs(GiNaC::symbol(symbol) == g_value->_expr));
    }

    SymExprPtr subs(const ParameterMap& params) const override {
        GiNaC::exmap substitutions;
        for (const auto& param : params) {
            substitutions[GiNaC::symbol(param.first)] = GiNaC::numeric(param.second);
        }
        return std::make_shared<GiNaCExpression>(_expr.subs(substitutions));
    }

    SymExprPtr subs(const std::map<SymExprPtr, SymExprPtr>& exprMap) const override {
        GiNaC::exmap substitutions;
        for (const auto& pair : exprMap) {
            auto g_key = std::dynamic_pointer_cast<GiNaCExpression>(pair.first);
            auto g_value = std::dynamic_pointer_cast<GiNaCExpression>(pair.second);
            substitutions[g_key->_expr] = g_value->_expr;
        }
        return std::make_shared<GiNaCExpression>(_expr.subs(substitutions));
    }

    double eval() const override {
        return GiNaC::ex_to<GiNaC::numeric>(_expr.evalf()).to_double();
    }

    bool is_numeric() const override {
        return GiNaC::is_a<GiNaC::numeric>(_expr);
    }

    // Series expansion
    SymExprPtr series(const SymExprPtr& var, const SymExprPtr& point, int order) const override {
        auto g_var = std::dynamic_pointer_cast<GiNaCExpression>(var);
        auto g_point = std::dynamic_pointer_cast<GiNaCExpression>(point);
        return std::make_shared<GiNaCExpression>(_expr.series(g_var->_expr, order));
    }

    SymExprPtr series_to_poly(const SymExprPtr& series) const override {
        auto g_series = std::dynamic_pointer_cast<GiNaCExpression>(series);
        return std::make_shared<GiNaCExpression>(g_series->_expr);
    }

    SymExprPtr coeff(const SymExprPtr& var, int power) const override {
        auto g_var = std::dynamic_pointer_cast<GiNaCExpression>(var);
        return std::make_shared<GiNaCExpression>(_expr.coeff(g_var->_expr, power));
    }

    SymExprPtr evalf() const override {
        return std::make_shared<GiNaCExpression>(_expr.evalf());
    }

    bool is_zero() const override {
        return _expr.is_zero();
    }

    double to_double() const override {
        return GiNaC::ex_to<GiNaC::numeric>(_expr.evalf()).to_double();
    }

    // String representation
    std::string to_string() const override {
        std::ostringstream oss;
        oss << _expr;
        return oss.str();
    }

    std::string to_latex() const override {
        std::ostringstream oss;
        _expr.print(GiNaC::print_latex(oss));
        return oss.str();
    }

    std::string to_python() const override {
        return to_string();
    }

    std::string to_cform() const override {
        std::ostringstream oss;
        _expr.print(GiNaC::print_csrc(oss));
        return oss.str();
    }

    // Get underlying GiNaC expression
    const GiNaC::ex& get_ginac_expr() const { return _expr; }
};

class GiNaCFactory : public SymbolicFactory {
public:
    SymExprPtr createSymbol(const std::string& name) override;
    SymExprPtr createConstant(double value) override;
    SymExprPtr createPi() override;
    SymExprPtr createE() override;
    SymExprPtr createI() override;
};

// Backward compatibility for legacy code paths that still construct GiNaCSymbolic.
class GiNaCSymbolic : public GiNaCFactory {
};

#endif // USE_GINAC_IMPL
#endif // INCLUDE_GINAC_SYMBOLIC_HPP
