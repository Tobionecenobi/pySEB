#ifndef INCLUDE_GUARD_SYMBOLICPORTABLE
#define INCLUDE_GUARD_SYMBOLICPORTABLE

#include "SymbolicInterface.hpp"
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace sebsym {

class PortableExpression : public SymbolicExpression {
public:
    enum class Kind {
        Constant,
        Symbol,
        Unary,
        Binary,
        Function
    };

    PortableExpression(double value);
    PortableExpression(const std::string& symbol);
    PortableExpression(Kind kind,
                       const std::string& op,
                       const std::vector<SymExprPtr>& args);

    SymExprPtr add(const SymExprPtr& other) const override;
    SymExprPtr subtract(const SymExprPtr& other) const override;
    SymExprPtr multiply(const SymExprPtr& other) const override;
    SymExprPtr divide(const SymExprPtr& other) const override;
    SymExprPtr negate() const override;
    SymExprPtr pow(const SymExprPtr& exponent) const override;

    SymExprPtr exp() const override;
    SymExprPtr log() const override;
    SymExprPtr sin() const override;
    SymExprPtr cos() const override;
    SymExprPtr tan() const override;
    SymExprPtr asin() const override;
    SymExprPtr acos() const override;
    SymExprPtr atan() const override;
    SymExprPtr sinh() const override;
    SymExprPtr cosh() const override;
    SymExprPtr tanh() const override;
    SymExprPtr sqrt() const override;
    SymExprPtr abs() const override;

    SymExprPtr bessel_j0() const override;
    SymExprPtr bessel_j1() const override;
    SymExprPtr dawson() const override;
    SymExprPtr erf() const override;
    SymExprPtr erfc() const override;

    SymExprPtr subs(const std::string& symbol, const SymExprPtr& value) const override;
    SymExprPtr subs(const ParameterMap& params) const override;
    SymExprPtr subs(const std::map<SymExprPtr, SymExprPtr>& exprMap) const override;
    double eval() const override;
    bool is_numeric() const override;

    SymExprPtr series(const SymExprPtr& var, const SymExprPtr& point, int order) const override;
    SymExprPtr series_to_poly(const SymExprPtr& series) const override;
    SymExprPtr coeff(const SymExprPtr& var, int power) const override;
    SymExprPtr evalf() const override;
    bool is_zero() const override;
    double to_double() const override;

    std::string to_string() const override;
    std::string to_latex() const override;
    std::string to_python() const override;
    std::string to_cform() const override;

    Kind kind() const { return _kind; }
    const std::string& op() const { return _op; }
    const std::vector<SymExprPtr>& args() const { return _args; }
    double value() const { return _value; }

private:
    Kind _kind;
    std::string _op;
    std::vector<SymExprPtr> _args;
    double _value = 0.0;

    SymExprPtr clone() const;
    SymExprPtr unary(const std::string& op) const;
    SymExprPtr binary(const std::string& op, const SymExprPtr& other) const;
    SymExprPtr function(const std::string& name) const;

    std::string render_python() const;
    std::string render_plain() const;
    bool has_symbols() const;
};

class PortableFactory : public SymbolicFactory {
public:
    std::string backendName() const override { return "portable"; }
    SymbolicCapabilities capabilities() const override;

    SymExprPtr createSymbol(const std::string& name) override;
    SymExprPtr createConstant(double value) override;
    SymExprPtr createPi() override;
    SymExprPtr createE() override;
    SymExprPtr createI() override;
};

SymbolicFactory* portableFactory();
void registerPortableBackend();

} // namespace sebsym

#endif // INCLUDE_GUARD_SYMBOLICPORTABLE
