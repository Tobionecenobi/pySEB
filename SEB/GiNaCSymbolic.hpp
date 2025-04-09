//===========================================================================
// Included guards
#ifndef INCLUDE_GUARD_GINACSYMBOLIC
#define INCLUDE_GUARD_GINACSYMBOLIC

//===========================================================================
// included dependencies
#include "SymbolicInterface.hpp"
#include <ginac/ginac.h>
#include <sstream>
#include "SpecialFunctions.hpp"

//===========================================================================
// GiNaC implementation of symbolic expressions
class GiNaCExpression : public SymbolicExpression {
public:
    GiNaCExpression(const GiNaC::ex& expr) : _expr(expr) {}
    virtual ~GiNaCExpression() = default;

    // Get the underlying GiNaC expression
    const GiNaC::ex& expr() const { return _expr; }

    // Basic operations
    SymExprPtr add(const SymExprPtr& other) const override {
        auto ginac_other = std::dynamic_pointer_cast<GiNaCExpression>(other);
        return std::make_shared<GiNaCExpression>(_expr + ginac_other->expr());
    }

    SymExprPtr subtract(const SymExprPtr& other) const override {
        auto ginac_other = std::dynamic_pointer_cast<GiNaCExpression>(other);
        return std::make_shared<GiNaCExpression>(_expr - ginac_other->expr());
    }

    SymExprPtr multiply(const SymExprPtr& other) const override {
        auto ginac_other = std::dynamic_pointer_cast<GiNaCExpression>(other);
        return std::make_shared<GiNaCExpression>(_expr * ginac_other->expr());
    }

    SymExprPtr divide(const SymExprPtr& other) const override {
        auto ginac_other = std::dynamic_pointer_cast<GiNaCExpression>(other);
        return std::make_shared<GiNaCExpression>(_expr / ginac_other->expr());
    }

    SymExprPtr negate() const override {
        return std::make_shared<GiNaCExpression>(-_expr);
    }

    SymExprPtr pow(const SymExprPtr& exponent) const override {
        auto ginac_exponent = std::dynamic_pointer_cast<GiNaCExpression>(exponent);
        return std::make_shared<GiNaCExpression>(GiNaC::pow(_expr, ginac_exponent->expr()));
    }

    // Functions
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

    // Special functions needed for scattering
    SymExprPtr bessel_j0() const override {
        return std::make_shared<GiNaCExpression>(BesselJ0(_expr));
    }

    SymExprPtr bessel_j1() const override {
        return std::make_shared<GiNaCExpression>(BesselJ1(_expr));
    }

    SymExprPtr dawson() const override {
        return std::make_shared<GiNaCExpression>(DawsonF(_expr));
    }

    SymExprPtr erf() const override {
        return std::make_shared<GiNaCExpression>(Erf(_expr));
    }

    SymExprPtr erfc() const override {
        return std::make_shared<GiNaCExpression>(Erfc(_expr));
    }

    // Substitution and evaluation
    SymExprPtr subs(const std::string& symbol, const SymExprPtr& value) const override {
        std::cout << "Substituting " << symbol << " with value" << std::endl;
        auto ginac_value = std::dynamic_pointer_cast<GiNaCExpression>(value);
        if (!ginac_value) {
            std::cerr << "Error: Failed to cast value to GiNaCExpression" << std::endl;
            return std::make_shared<GiNaCExpression>(_expr);
        }
        GiNaC::symbol sym(symbol);
        return std::make_shared<GiNaCExpression>(_expr.subs(sym == ginac_value->expr()));
    }

    SymExprPtr subs(const ParameterMap& params) const override {
        GiNaC::exmap substitutions;
        for (const auto& param : params) {
            GiNaC::symbol sym(param.first);
            substitutions[sym] = GiNaC::numeric(param.second);
        }
        return std::make_shared<GiNaCExpression>(_expr.subs(substitutions));
    }

    SymExprPtr subs(const std::map<SymExprPtr, SymExprPtr>& exprMap) const override {
        GiNaC::exmap substitutions;
        for (const auto& pair : exprMap) {
            auto ginacExpr1 = std::dynamic_pointer_cast<GiNaCExpression>(pair.first);
            auto ginacExpr2 = std::dynamic_pointer_cast<GiNaCExpression>(pair.second);
            if (ginacExpr1 && ginacExpr2) {
                substitutions[ginacExpr1->expr()] = ginacExpr2->expr();
            }
        }
        return std::make_shared<GiNaCExpression>(_expr.subs(substitutions));
    }

    double eval() const override {
        if (!is_numeric()) {
            throw std::runtime_error("Cannot evaluate non-numeric expression");
        }
        return GiNaC::ex_to<GiNaC::numeric>(_expr).to_double();
    }

    bool is_numeric() const override {
        return GiNaC::is_a<GiNaC::numeric>(_expr);
    }

    // Series expansion and coefficient extraction
    SymExprPtr series(const SymExprPtr& var, const SymExprPtr& point, int order) const override {
        auto ginac_var = std::dynamic_pointer_cast<GiNaCExpression>(var);
        auto ginac_point = std::dynamic_pointer_cast<GiNaCExpression>(point);
        return std::make_shared<GiNaCExpression>(_expr.series(ginac_var->expr() == ginac_point->expr(), order));
    }

    SymExprPtr series_to_poly(const SymExprPtr& series) const override {
        auto ginac_series = std::dynamic_pointer_cast<GiNaCExpression>(series);
        return std::make_shared<GiNaCExpression>(GiNaC::series_to_poly(ginac_series->expr()));
    }

    SymExprPtr coeff(const SymExprPtr& var, int power) const override {
        auto ginac_var = std::dynamic_pointer_cast<GiNaCExpression>(var);
        return std::make_shared<GiNaCExpression>(_expr.coeff(ginac_var->expr(), power));
    }

    SymExprPtr evalf() const override {
        return std::make_shared<GiNaCExpression>(_expr.evalf());
    }

    bool is_zero() const override {
        return GiNaC::is_zero(_expr);
    }

    double to_double() const override {
        try {
            if (is_numeric()) {
                return GiNaC::ex_to<GiNaC::numeric>(_expr).to_double();
            } else {
                // Try to evaluate the expression first
                GiNaC::ex evaluated = _expr.evalf();
                if (GiNaC::is_a<GiNaC::numeric>(evaluated)) {
                    return GiNaC::ex_to<GiNaC::numeric>(evaluated).to_double();
                }
            }
            throw std::runtime_error("Cannot convert expression to double: " + to_string());
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Error converting to double: ") + e.what());
        }
    }

    // String representation
    std::string to_string() const override {
        std::ostringstream oss;
        oss << _expr;
        return oss.str();
    }

    std::string to_latex() const override {
        std::ostringstream oss;
        oss << GiNaC::latex << _expr;
        return oss.str();
    }

    std::string to_python() const override {
        std::ostringstream oss;
        oss << GiNaC::python << _expr;
        return oss.str();
    }

    std::string to_cform() const override {
        std::ostringstream oss;
        oss << GiNaC::csrc_double << _expr;
        return oss.str();
    }

    // Get the underlying GiNaC expression
    const GiNaC::ex& getGiNaCExpression() const {
        return _expr;
    }

private:
    GiNaC::ex _expr;
};

//===========================================================================
// GiNaC implementation of symbolic factory
class GiNaCSymbolic : public SymbolicFactory {
public:
    GiNaCSymbolic() = default;
    virtual ~GiNaCSymbolic() = default;

    SymExprPtr createSymbol(const std::string& name) override {
        return std::make_shared<GiNaCExpression>(GiNaC::symbol(name));
    }

    SymExprPtr createConstant(double value) override {
        return std::make_shared<GiNaCExpression>(GiNaC::numeric(value));
    }

    SymExprPtr createPi() override {
        return std::make_shared<GiNaCExpression>(GiNaC::Pi);
    }

    SymExprPtr createE() override {
        return std::make_shared<GiNaCExpression>(GiNaC::exp(1));
    }

    SymExprPtr createI() override {
        return std::make_shared<GiNaCExpression>(GiNaC::I);
    }
};

#endif // INCLUDE_GUARD_GINACSYMBOLIC
