//===========================================================================
// Included guards
#ifndef INCLUDE_GUARD_SYMBOLICINTERFACE
#define INCLUDE_GUARD_SYMBOLICINTERFACE

//===========================================================================
// included dependencies
#include <string>
#include <map>
#include <vector>
#include <memory>

//===========================================================================
// Forward declarations
class SymbolicExpression;
typedef std::shared_ptr<SymbolicExpression> SymExprPtr;
typedef std::map<std::string, double> ParameterMap;
typedef std::vector<double> DoubleVector;

//===========================================================================
// Abstract interface for symbolic expressions
class SymbolicExpression {
public:
    virtual ~SymbolicExpression() = default;

    // Basic operations
    virtual SymExprPtr add(const SymExprPtr& other) const = 0;
    virtual SymExprPtr subtract(const SymExprPtr& other) const = 0;
    virtual SymExprPtr multiply(const SymExprPtr& other) const = 0;
    virtual SymExprPtr divide(const SymExprPtr& other) const = 0;
    virtual SymExprPtr negate() const = 0;
    virtual SymExprPtr pow(const SymExprPtr& exponent) const = 0;

    // Functions
    virtual SymExprPtr exp() const = 0;
    virtual SymExprPtr log() const = 0;
    virtual SymExprPtr sin() const = 0;
    virtual SymExprPtr cos() const = 0;
    virtual SymExprPtr tan() const = 0;
    virtual SymExprPtr asin() const = 0;
    virtual SymExprPtr acos() const = 0;
    virtual SymExprPtr atan() const = 0;
    virtual SymExprPtr sinh() const = 0;
    virtual SymExprPtr cosh() const = 0;
    virtual SymExprPtr tanh() const = 0;
    virtual SymExprPtr sqrt() const = 0;
    virtual SymExprPtr abs() const = 0;

    // Special functions needed for scattering
    virtual SymExprPtr bessel_j0() const = 0;
    virtual SymExprPtr bessel_j1() const = 0;
    virtual SymExprPtr dawson() const = 0;
    virtual SymExprPtr erf() const = 0;
    virtual SymExprPtr erfc() const = 0;

    // Substitution and evaluation
    virtual SymExprPtr subs(const std::string& symbol, const SymExprPtr& value) const = 0;
    virtual SymExprPtr subs(const ParameterMap& params) const = 0;
    virtual double eval() const = 0;
    virtual bool is_numeric() const = 0;

    // Series expansion and coefficient extraction
    virtual SymExprPtr series(const SymExprPtr& var, const SymExprPtr& point, int order) const = 0;
    virtual SymExprPtr series_to_poly(const SymExprPtr& series) const = 0;
    virtual SymExprPtr coeff(const SymExprPtr& var, int power) const = 0;
    virtual SymExprPtr evalf() const = 0;
    virtual bool is_zero() const = 0;
    virtual double to_double() const = 0;

    // String representation
    virtual std::string to_string() const = 0;
    virtual std::string to_latex() const = 0;
    virtual std::string to_python() const = 0;
    virtual std::string to_cform() const = 0;

    // Symbol creation and manipulation
    static SymExprPtr symbol(const std::string& name);
    static SymExprPtr constant(double value);
    static SymExprPtr pi();
    static SymExprPtr e();
    static SymExprPtr i();
};

//===========================================================================
// Factory for creating symbolic expressions
class SymbolicFactory {
public:
    virtual ~SymbolicFactory() = default;

    virtual SymExprPtr createSymbol(const std::string& name) = 0;
    virtual SymExprPtr createConstant(double value) = 0;
    virtual SymExprPtr createPi() = 0;
    virtual SymExprPtr createE() = 0;
    virtual SymExprPtr createI() = 0;

    // Get the singleton instance
    static SymbolicFactory* instance();

    // Set the factory implementation
    static void setInstance(SymbolicFactory* factory);

private:
    static SymbolicFactory* _instance;
};

// Inline implementation of static methods
inline SymExprPtr SymbolicExpression::symbol(const std::string& name) {
    return SymbolicFactory::instance()->createSymbol(name);
}

inline SymExprPtr SymbolicExpression::constant(double value) {
    return SymbolicFactory::instance()->createConstant(value);
}

inline SymExprPtr SymbolicExpression::pi() {
    return SymbolicFactory::instance()->createPi();
}

inline SymExprPtr SymbolicExpression::e() {
    return SymbolicFactory::instance()->createE();
}

inline SymExprPtr SymbolicExpression::i() {
    return SymbolicFactory::instance()->createI();
}

#endif // INCLUDE_GUARD_SYMBOLICINTERFACE
