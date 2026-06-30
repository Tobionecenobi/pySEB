#ifndef INCLUDE_GUARD_SYMPYSYMBOLIC
#define INCLUDE_GUARD_SYMPYSYMBOLIC

#include "SymbolicInterface.hpp"
#include <string>
#include <sstream>
#include <memory>
#include <stdexcept>

//===========================================================================
// SymPy implementation of symbolic expression
class SymPyExpression : public SymbolicExpression {
public:
    SymPyExpression(const std::string& expr_str) : _expr_str(expr_str) {}
    virtual ~SymPyExpression() = default;

    // Basic operations
    SymExprPtr add(const SymExprPtr& other) const override {
        auto sympy_other = std::dynamic_pointer_cast<SymPyExpression>(other);
        if (!sympy_other) {
            throw std::runtime_error("Cannot add non-SymPyExpression");
        }
        // In C++, we just store the expression as a string
        // The actual computation will be done in Python
        return std::make_shared<SymPyExpression>(
            "(" + _expr_str + " + " + sympy_other->_expr_str + ")"
        );
    }

    SymExprPtr subtract(const SymExprPtr& other) const override {
        auto sympy_other = std::dynamic_pointer_cast<SymPyExpression>(other);
        if (!sympy_other) {
            throw std::runtime_error("Cannot subtract non-SymPyExpression");
        }
        return std::make_shared<SymPyExpression>(
            "(" + _expr_str + " - " + sympy_other->_expr_str + ")"
        );
    }

    SymExprPtr multiply(const SymExprPtr& other) const override {
        auto sympy_other = std::dynamic_pointer_cast<SymPyExpression>(other);
        if (!sympy_other) {
            throw std::runtime_error("Cannot multiply non-SymPyExpression");
        }
        return std::make_shared<SymPyExpression>(
            "(" + _expr_str + " * " + sympy_other->_expr_str + ")"
        );
    }

    SymExprPtr divide(const SymExprPtr& other) const override {
        auto sympy_other = std::dynamic_pointer_cast<SymPyExpression>(other);
        if (!sympy_other) {
            throw std::runtime_error("Cannot divide non-SymPyExpression");
        }
        return std::make_shared<SymPyExpression>(
            "(" + _expr_str + " / " + sympy_other->_expr_str + ")"
        );
    }

    SymExprPtr negate() const override {
        return std::make_shared<SymPyExpression>("(-" + _expr_str + ")");
    }

    SymExprPtr pow(const SymExprPtr& exponent) const override {
        auto sympy_exponent = std::dynamic_pointer_cast<SymPyExpression>(exponent);
        if (!sympy_exponent) {
            throw std::runtime_error("Cannot raise to non-SymPyExpression power");
        }
        return std::make_shared<SymPyExpression>(
            "(" + _expr_str + "**" + sympy_exponent->_expr_str + ")"
        );
    }

    // Non-virtual method to handle double exponents
    SymExprPtr pow(double exponent) const {
        std::ostringstream oss;
        oss << exponent;
        return std::make_shared<SymPyExpression>(
            "(" + _expr_str + "**" + oss.str() + ")"
        );
    }

    // Basic functions
    SymExprPtr exp() const override {
        return std::make_shared<SymPyExpression>("sympy.exp(" + _expr_str + ")");
    }

    SymExprPtr log() const override {
        return std::make_shared<SymPyExpression>("sympy.log(" + _expr_str + ")");
    }

    SymExprPtr sin() const override {
        return std::make_shared<SymPyExpression>("sympy.sin(" + _expr_str + ")");
    }

    SymExprPtr cos() const override {
        return std::make_shared<SymPyExpression>("sympy.cos(" + _expr_str + ")");
    }

    SymExprPtr tan() const override {
        return std::make_shared<SymPyExpression>("sympy.tan(" + _expr_str + ")");
    }

    SymExprPtr asin() const override {
        return std::make_shared<SymPyExpression>("sympy.asin(" + _expr_str + ")");
    }

    SymExprPtr acos() const override {
        return std::make_shared<SymPyExpression>("sympy.acos(" + _expr_str + ")");
    }

    SymExprPtr atan() const override {
        return std::make_shared<SymPyExpression>("sympy.atan(" + _expr_str + ")");
    }

    SymExprPtr sinh() const override {
        return std::make_shared<SymPyExpression>("sympy.sinh(" + _expr_str + ")");
    }

    SymExprPtr cosh() const override {
        return std::make_shared<SymPyExpression>("sympy.cosh(" + _expr_str + ")");
    }

    SymExprPtr tanh() const override {
        return std::make_shared<SymPyExpression>("sympy.tanh(" + _expr_str + ")");
    }

    SymExprPtr sqrt() const override {
        return std::make_shared<SymPyExpression>("sympy.sqrt(" + _expr_str + ")");
    }

    SymExprPtr abs() const override {
        return std::make_shared<SymPyExpression>("sympy.Abs(" + _expr_str + ")");
    }

    // Special functions for scattering
    SymExprPtr bessel_j0() const override {
        return std::make_shared<SymPyExpression>("sympy.besselj(0, " + _expr_str + ")");
    }

    SymExprPtr bessel_j1() const override {
        return std::make_shared<SymPyExpression>("sympy.besselj(1, " + _expr_str + ")");
    }

    SymExprPtr dawson() const override {
        // Dawson function F(x) = exp(-x²)∫₀ˣexp(t²)dt
        return std::make_shared<SymPyExpression>(
            "(sympy.exp(-(" + _expr_str + ")**2) * sympy.sqrt(sympy.pi)/2 * sympy.erfi(" + _expr_str + "))"
        );
    }

    SymExprPtr erf() const override {
        return std::make_shared<SymPyExpression>("sympy.erf(" + _expr_str + ")");
    }

    SymExprPtr erfc() const override {
        return std::make_shared<SymPyExpression>("sympy.erfc(" + _expr_str + ")");
    }

    // Substitution and evaluation
    SymExprPtr subs(const std::string& symbol, const SymExprPtr& value) const override {
        auto sympy_value = std::dynamic_pointer_cast<SymPyExpression>(value);
        if (!sympy_value) {
            throw std::runtime_error("Cannot substitute non-SymPyExpression value");
        }
        // In C++, we just build the expression string
        // The actual substitution will be done in Python
        return std::make_shared<SymPyExpression>(
            _expr_str + ".subs('" + symbol + "', " + sympy_value->_expr_str + ")"
        );
    }

    // Implementation of ParameterMap substitution
    SymExprPtr subs(const ParameterMap& params) const override {
        std::string result = _expr_str;
        for (const auto& param : params) {
            std::ostringstream oss;
            oss << param.second;
            result += ".subs('" + param.first + "', " + oss.str() + ")";
        }
        return std::make_shared<SymPyExpression>(result);
    }

    // Implementation of expression map substitution
    SymExprPtr subs(const std::map<SymExprPtr, SymExprPtr>& exprMap) const override {
        std::string result = _expr_str;
        for (const auto& pair : exprMap) {
            auto sympy_key = std::dynamic_pointer_cast<SymPyExpression>(pair.first);
            auto sympy_value = std::dynamic_pointer_cast<SymPyExpression>(pair.second);
            if (!sympy_key || !sympy_value) {
                throw std::runtime_error("Cannot substitute non-SymPyExpression in expression map");
            }
            result += ".subs(" + sympy_key->_expr_str + ", " + sympy_value->_expr_str + ")";
        }
        return std::make_shared<SymPyExpression>(result);
    }

    // Non-virtual method to handle double value substitutions
    SymExprPtr subs(const std::string& symbol, double value) const {
        std::ostringstream oss;
        oss << value;
        return std::make_shared<SymPyExpression>(
            _expr_str + ".subs('" + symbol + "', " + oss.str() + ")"
        );
    }

    double eval() const override {
        throw std::runtime_error("Cannot evaluate SymPyExpression in C++, use Python");
    }

    bool is_numeric() const override {
        // This is a placeholder - in C++ we can't determine if a SymPy expression is numeric
        return false;
    }

    // Series expansion
    SymExprPtr series(const SymExprPtr& var, const SymExprPtr& point, int order) const override {
        auto sympy_var = std::dynamic_pointer_cast<SymPyExpression>(var);
        auto sympy_point = std::dynamic_pointer_cast<SymPyExpression>(point);
        if (!sympy_var || !sympy_point) {
            throw std::runtime_error("Cannot expand series with non-SymPyExpression arguments");
        }
        std::ostringstream oss;
        oss << order;
        return std::make_shared<SymPyExpression>(
            _expr_str + ".series(" + sympy_var->_expr_str + ", " + sympy_point->_expr_str + ", " + oss.str() + ").removeO()"
        );
    }

    SymExprPtr series_to_poly(const SymExprPtr& series) const override {
        auto sympy_series = std::dynamic_pointer_cast<SymPyExpression>(series);
        if (!sympy_series) {
            throw std::runtime_error("Cannot convert non-SymPyExpression series to polynomial");
        }
        return std::make_shared<SymPyExpression>(
            "sympy.Poly(" + sympy_series->_expr_str + ")"
        );
    }

    SymExprPtr coeff(const SymExprPtr& var, int power) const override {
        auto sympy_var = std::dynamic_pointer_cast<SymPyExpression>(var);
        if (!sympy_var) {
            throw std::runtime_error("Cannot get coefficient with non-SymPyExpression variable");
        }
        std::ostringstream oss;
        oss << power;
        return std::make_shared<SymPyExpression>(
            _expr_str + ".coeff(" + sympy_var->_expr_str + ", " + oss.str() + ")"
        );
    }

    SymExprPtr evalf() const override {
        return std::make_shared<SymPyExpression>(_expr_str + ".evalf()");
    }

    bool is_zero() const override {
        // This is a placeholder - in C++ we can't determine if a SymPy expression is zero
        return false;
    }

    double to_double() const override {
        throw std::runtime_error("Cannot convert SymPyExpression to double in C++, use Python");
    }

    // String representations
    std::string to_string() const override {
        return _expr_str;
    }

    std::string to_latex() const override {
        return "sympy.latex(" + _expr_str + ")";
    }

    std::string to_python() const override {
        return _expr_str;
    }

    std::string to_cform() const override {
        return "sympy.ccode(" + _expr_str + ")";
    }

    // Get the underlying expression string
    const std::string& getExpressionString() const {
        return _expr_str;
    }

private:
    std::string _expr_str;
};

//===========================================================================
// SymPy implementation of symbolic factory
class SymPySymbolic : public SymbolicFactory {
public:
    SymPySymbolic() = default;
    virtual ~SymPySymbolic() = default;

    std::string backendName() const override { return "sympy-string"; }

    SymbolicCapabilities capabilities() const override {
        SymbolicCapabilities caps;
        caps.python_output = true;
        return caps;
    }

    SymExprPtr createSymbol(const std::string& name) override {
        return std::make_shared<SymPyExpression>("sympy.Symbol('" + name + "')");
    }

    SymExprPtr createConstant(double value) override {
        std::ostringstream oss;
        oss << value;
        return std::make_shared<SymPyExpression>(oss.str());
    }

    SymExprPtr createPi() override {
        return std::make_shared<SymPyExpression>("sympy.pi");
    }

    SymExprPtr createE() override {
        return std::make_shared<SymPyExpression>("sympy.E");
    }

    SymExprPtr createI() override {
        return std::make_shared<SymPyExpression>("sympy.I");
    }
};

#endif // INCLUDE_GUARD_SYMPYSYMBOLIC
