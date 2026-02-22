//===========================================================================
// Included guards
#ifndef INCLUDE_GUARD_EXPRESSION
#define INCLUDE_GUARD_EXPRESSION

//===========================================================================
// included dependencies
#include <string>
#include <map>
#include <stdexcept>
#include "SymbolicInterface.hpp"

// Forward declarations
class GiNaCSymbolic;
class GiNaCExpression;

//===========================================================================
// used namespaces
using namespace std;

// Forward declarations
class Expression;

// Define parameter map type
typedef std::map<std::string, double> ParameterMap;

/**
 * Expression class that wraps SymExprPtr and provides operator overloading
 * to make it compatible with the existing GiNaC-based code.
 */
class Expression {
private:
    SymExprPtr expr;

public:
    // Constructors
    Expression() : expr(nullptr) {}
    Expression(const SymExprPtr& e) : expr(e) {}
    Expression(const Expression& other) : expr(other.expr) {}

    // Assignment operators
    Expression& operator=(const Expression& other) {
        expr = other.expr;
        return *this;
    }

    Expression& operator=(int value) {
        expr = SymbolicExpression::constant(value);
        return *this;
    }

    Expression& operator=(double value) {
        expr = SymbolicExpression::constant(value);
        return *this;
    }

    // Conversion to SymExprPtr
    operator SymExprPtr() const { return expr; }
    SymExprPtr get() const { return expr; }

    // Basic operators
    Expression operator+(const Expression& other) const {
        if (!expr || !other.expr) return Expression();
        return Expression(expr->add(other.expr));
    }

    Expression operator+(int value) const {
        if (!expr) return Expression();
        return Expression(expr->add(SymbolicExpression::constant(value)));
    }

    Expression operator+(double value) const {
        if (!expr) return Expression();
        return Expression(expr->add(SymbolicExpression::constant(value)));
    }

    friend Expression operator+(int value, const Expression& e) {
        if (!e.expr) return Expression();
        return Expression(SymbolicExpression::constant(value)->add(e.expr));
    }

    friend Expression operator+(double value, const Expression& e) {
        if (!e.expr) return Expression();
        return Expression(SymbolicExpression::constant(value)->add(e.expr));
    }

    Expression operator-(const Expression& other) const {
        if (!expr || !other.expr) return Expression();
        return Expression(expr->subtract(other.expr));
    }

    Expression operator-(int value) const {
        if (!expr) return Expression();
        return Expression(expr->subtract(SymbolicExpression::constant(value)));
    }

    Expression operator-(double value) const {
        if (!expr) return Expression();
        return Expression(expr->subtract(SymbolicExpression::constant(value)));
    }

    friend Expression operator-(int value, const Expression& e) {
        if (!e.expr) return Expression();
        return Expression(SymbolicExpression::constant(value)->subtract(e.expr));
    }

    friend Expression operator-(double value, const Expression& e) {
        if (!e.expr) return Expression();
        return Expression(SymbolicExpression::constant(value)->subtract(e.expr));
    }

    Expression operator*(const Expression& other) const {
        if (!expr || !other.expr) return Expression();
        return Expression(expr->multiply(other.expr));
    }

    Expression operator*(int value) const {
        if (!expr) return Expression();
        return Expression(expr->multiply(SymbolicExpression::constant(value)));
    }

    Expression operator*(double value) const {
        if (!expr) return Expression();
        return Expression(expr->multiply(SymbolicExpression::constant(value)));
    }

    friend Expression operator*(int value, const Expression& e) {
        if (!e.expr) return Expression();
        return Expression(SymbolicExpression::constant(value)->multiply(e.expr));
    }

    friend Expression operator*(double value, const Expression& e) {
        if (!e.expr) return Expression();
        return Expression(SymbolicExpression::constant(value)->multiply(e.expr));
    }

    Expression operator/(const Expression& other) const {
        if (!expr || !other.expr) return Expression();
        return Expression(expr->divide(other.expr));
    }

    Expression operator/(int value) const {
        if (!expr) return Expression();
        return Expression(expr->divide(SymbolicExpression::constant(value)));
    }

    Expression operator/(double value) const {
        if (!expr) return Expression();
        return Expression(expr->divide(SymbolicExpression::constant(value)));
    }

    friend Expression operator/(int value, const Expression& e) {
        if (!e.expr) return Expression();
        return Expression(SymbolicExpression::constant(value)->divide(e.expr));
    }

    friend Expression operator/(double value, const Expression& e) {
        if (!e.expr) return Expression();
        return Expression(SymbolicExpression::constant(value)->divide(e.expr));
    }

    Expression operator-() const {
        if (!expr) return Expression();
        return Expression(expr->negate());
    }



    // Functions
    Expression exp() const {
        if (!expr) return Expression();
        return Expression(expr->exp());
    }

    Expression log() const {
        if (!expr) return Expression();
        return Expression(expr->log());
    }

    Expression sin() const {
        if (!expr) return Expression();
        return Expression(expr->sin());
    }

    Expression cos() const {
        if (!expr) return Expression();
        return Expression(expr->cos());
    }

    Expression tan() const {
        if (!expr) return Expression();
        return Expression(expr->tan());
    }

    Expression asin() const {
        if (!expr) return Expression();
        return Expression(expr->asin());
    }

    Expression acos() const {
        if (!expr) return Expression();
        return Expression(expr->acos());
    }

    Expression atan() const {
        if (!expr) return Expression();
        return Expression(expr->atan());
    }

    Expression sinh() const {
        if (!expr) return Expression();
        return Expression(expr->sinh());
    }

    Expression cosh() const {
        if (!expr) return Expression();
        return Expression(expr->cosh());
    }

    Expression tanh() const {
        if (!expr) return Expression();
        return Expression(expr->tanh());
    }

    Expression sqrt() const {
        if (!expr) return Expression();
        return Expression(expr->sqrt());
    }

    Expression abs() const {
        if (!expr) return Expression();
        return Expression(expr->abs());
    }

    // Special functions needed for scattering
    Expression bessel_j0() const {
        if (!expr) return Expression();
        return Expression(expr->bessel_j0());
    }

    Expression bessel_j1() const {
        if (!expr) return Expression();
        return Expression(expr->bessel_j1());
    }

    Expression dawson() const {
        if (!expr) return Expression();
        return Expression(expr->dawson());
    }

    Expression erf() const {
        if (!expr) return Expression();
        return Expression(expr->erf());
    }

    Expression erfc() const {
        if (!expr) return Expression();
        return Expression(expr->erfc());
    }

    // Substitution and evaluation
    Expression subs(const std::string& symbol, const Expression& value) const {
        if (!expr) return Expression();
        return Expression(expr->subs(symbol, value.expr));
    }

    Expression subs(const ParameterMap& params) const {
        if (!expr) return Expression();
        return Expression(expr->subs(params));
    }

    Expression subs(const std::map<Expression, Expression>& exprMap) const {
        if (!expr) return Expression();
        std::map<SymExprPtr, SymExprPtr> symExprMap;
        for (const auto& pair : exprMap) {
            if (!pair.first.expr || !pair.second.expr) {
                continue;
            }
            symExprMap[pair.first.expr] = pair.second.expr;
        }
        return Expression(expr->subs(symExprMap));
    }

    // Handle substitution with a comparison expression like q==0.1
    Expression subs(bool comparisonExpr) const {
        if (!expr) return Expression();
        // This is a placeholder implementation
        // In a real implementation, you would parse the comparison expression
        // and perform the substitution
        return *this;
    }

    Expression evalf() const {
        if (!expr) return Expression();
        return Expression(expr->evalf());
    }

    double eval() const {
        if (!expr) throw std::runtime_error("Cannot evaluate null expression");
        return expr->eval();
    }

    bool is_numeric() const {
        if (!expr) return false;
        return expr->is_numeric();
    }

    // GiNaC compatibility methods
    bool is_a_numeric() const {
        return is_numeric();
    }

    double to_numeric() const {
        if (!is_numeric()) {
            throw std::runtime_error("Expression is not numeric");
        }
        return to_double();
    }

    // Series expansion and coefficient extraction
    Expression series(const Expression& var, const Expression& point, int order) const {
        if (!expr) return Expression();
        return Expression(expr->series(var.expr, point.expr, order));
    }

    Expression series_to_poly(const Expression& series) const {
        if (!expr) return Expression();
        return Expression(expr->series_to_poly(series.expr));
    }

    Expression coeff(const Expression& var, int power) const {
        if (!expr) return Expression();
        return Expression(expr->coeff(var.expr, power));
    }

    bool is_zero() const {
        if (!expr) return true;
        return expr->is_zero();
    }

    double to_double() const {
        if (!expr) throw std::runtime_error("Cannot convert null expression to double");
        return expr->to_double();
    }

    // String representation
    std::string to_string() const {
        if (!expr) return "null";
        return expr->to_string();
    }

    // Conversion to string for use as a key in a ParameterList
    operator std::string() const {
        return to_string();
    }

    std::string to_latex() const {
        if (!expr) return "null";
        return expr->to_latex();
    }

    std::string to_python() const {
        if (!expr) return "null";
        return expr->to_python();
    }

    std::string to_cform() const {
        if (!expr) return "null";
        return expr->to_cform();
    }

    // Comparison operators
    bool operator==(const Expression& other) const {
        if (!expr && !other.expr) return true;
        if (!expr || !other.expr) return false;
        return expr->to_string() == other.expr->to_string(); // Simple string comparison
    }

    bool operator<(const Expression& other) const {
        if (!expr && !other.expr) return false;
        if (!expr) return true;
        if (!other.expr) return false;
        // Use string representation for comparison
        return expr->to_string() < other.expr->to_string();
    }

    bool operator==(double value) const {
        if (!expr || !is_numeric()) {
            return false;
        }
        return to_double() == value;
    }

    bool operator!=(const Expression& other) const {
        return !(*this == other);
    }

    friend bool operator==(double value, const Expression& expr) {
        return expr == value;
    }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, const Expression& e) {
        // Just use the string representation for now
        os << e.to_string();
        return os;
    }
};

// Helper functions to create expressions
inline Expression symbol(const std::string& name) {
    return Expression(SymbolicExpression::symbol(name));
}

inline Expression constant(double value) {
    return Expression(SymbolicExpression::constant(value));
}

inline Expression pi() {
    return Expression(SymbolicExpression::pi());
}

inline Expression e() {
    return Expression(SymbolicExpression::e());
}

inline Expression i() {
    return Expression(SymbolicExpression::i());
}

// Power function
inline Expression pow(const Expression& base, const Expression& exponent) {
    if (!base.get()) return Expression();
    return Expression(base.get()->pow(exponent.get()));
}

inline Expression pow(const Expression& base, int exponent) {
    if (!base.get()) return Expression();
    return Expression(base.get()->pow(SymbolicExpression::constant(exponent)));
}

inline Expression pow(const Expression& base, double exponent) {
    if (!base.get()) return Expression();
    return Expression(base.get()->pow(SymbolicExpression::constant(exponent)));
}

// Convenience functions for numeric conversion
inline double to_double(const Expression& e) {
    return e.eval();
}

#endif // INCLUDE_GUARD_EXPRESSION
