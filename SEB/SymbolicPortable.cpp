#include "SymbolicPortable.hpp"

#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dawson.h>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace sebsym {
namespace {

std::shared_ptr<PortableExpression> as_portable(const SymExprPtr& expr)
{
    auto portable = std::dynamic_pointer_cast<PortableExpression>(expr);
    if (!portable) {
        throw std::runtime_error("portable backend cannot combine expressions from another backend");
    }
    return portable;
}

std::string double_to_string(double value)
{
    std::ostringstream oss;
    oss << std::setprecision(17) << value;
    return oss.str();
}

std::string sympy_function_name(const std::string& name)
{
    if (name == "abs") return "Abs";
    if (name == "bessel_j0") return "besselj";
    if (name == "bessel_j1") return "besselj";
    if (name == "dawson") return "dawson";
    return name;
}

double evaluate_function(const std::string& name, double x)
{
    if (name == "exp") return std::exp(x);
    if (name == "log") return std::log(x);
    if (name == "sin") return std::sin(x);
    if (name == "cos") return std::cos(x);
    if (name == "tan") return std::tan(x);
    if (name == "asin") return std::asin(x);
    if (name == "acos") return std::acos(x);
    if (name == "atan") return std::atan(x);
    if (name == "sinh") return std::sinh(x);
    if (name == "cosh") return std::cosh(x);
    if (name == "tanh") return std::tanh(x);
    if (name == "sqrt") return std::sqrt(x);
    if (name == "abs") return std::fabs(x);
    if (name == "bessel_j0") return gsl_sf_bessel_J0(x);
    if (name == "bessel_j1") return gsl_sf_bessel_J1(x);
    if (name == "dawson") return gsl_sf_dawson(x);
    if (name == "erf") return std::erf(x);
    if (name == "erfc") return std::erfc(x);
    throw std::runtime_error("portable backend cannot numerically evaluate function: " + name);
}

} // namespace

PortableExpression::PortableExpression(double value)
    : _kind(Kind::Constant), _value(value)
{
}

PortableExpression::PortableExpression(const std::string& symbol)
    : _kind(Kind::Symbol), _op(symbol)
{
}

PortableExpression::PortableExpression(Kind kind,
                                       const std::string& op,
                                       const std::vector<SymExprPtr>& args)
    : _kind(kind), _op(op), _args(args)
{
}

SymExprPtr PortableExpression::clone() const
{
    if (_kind == Kind::Constant) return std::make_shared<PortableExpression>(_value);
    if (_kind == Kind::Symbol) return std::make_shared<PortableExpression>(_op);
    return std::make_shared<PortableExpression>(_kind, _op, _args);
}

SymExprPtr PortableExpression::unary(const std::string& op) const
{
    return std::make_shared<PortableExpression>(Kind::Unary, op, std::vector<SymExprPtr>{clone()});
}

SymExprPtr PortableExpression::binary(const std::string& op, const SymExprPtr& other) const
{
    as_portable(other);
    return std::make_shared<PortableExpression>(Kind::Binary, op, std::vector<SymExprPtr>{clone(), other});
}

SymExprPtr PortableExpression::function(const std::string& name) const
{
    return std::make_shared<PortableExpression>(Kind::Function, name, std::vector<SymExprPtr>{clone()});
}

SymExprPtr PortableExpression::add(const SymExprPtr& other) const { return binary("+", other); }
SymExprPtr PortableExpression::subtract(const SymExprPtr& other) const { return binary("-", other); }
SymExprPtr PortableExpression::multiply(const SymExprPtr& other) const { return binary("*", other); }
SymExprPtr PortableExpression::divide(const SymExprPtr& other) const { return binary("/", other); }
SymExprPtr PortableExpression::negate() const { return unary("-"); }
SymExprPtr PortableExpression::pow(const SymExprPtr& exponent) const { return binary("**", exponent); }

SymExprPtr PortableExpression::exp() const { return function("exp"); }
SymExprPtr PortableExpression::log() const { return function("log"); }
SymExprPtr PortableExpression::sin() const { return function("sin"); }
SymExprPtr PortableExpression::cos() const { return function("cos"); }
SymExprPtr PortableExpression::tan() const { return function("tan"); }
SymExprPtr PortableExpression::asin() const { return function("asin"); }
SymExprPtr PortableExpression::acos() const { return function("acos"); }
SymExprPtr PortableExpression::atan() const { return function("atan"); }
SymExprPtr PortableExpression::sinh() const { return function("sinh"); }
SymExprPtr PortableExpression::cosh() const { return function("cosh"); }
SymExprPtr PortableExpression::tanh() const { return function("tanh"); }
SymExprPtr PortableExpression::sqrt() const { return function("sqrt"); }
SymExprPtr PortableExpression::abs() const { return function("abs"); }
SymExprPtr PortableExpression::bessel_j0() const { return function("bessel_j0"); }
SymExprPtr PortableExpression::bessel_j1() const { return function("bessel_j1"); }
SymExprPtr PortableExpression::dawson() const { return function("dawson"); }
SymExprPtr PortableExpression::erf() const { return function("erf"); }
SymExprPtr PortableExpression::erfc() const { return function("erfc"); }

SymExprPtr PortableExpression::subs(const std::string& symbol, const SymExprPtr& value) const
{
    as_portable(value);
    if (_kind == Kind::Symbol && _op == symbol) {
        return value;
    }
    if (_kind == Kind::Constant || _kind == Kind::Symbol) {
        return clone();
    }

    std::vector<SymExprPtr> new_args;
    for (const auto& arg : _args) {
        new_args.push_back(arg->subs(symbol, value));
    }
    return std::make_shared<PortableExpression>(_kind, _op, new_args);
}

SymExprPtr PortableExpression::subs(const ParameterMap& params) const
{
    SymExprPtr result = clone();
    for (const auto& param : params) {
        result = result->subs(param.first, std::make_shared<PortableExpression>(param.second));
    }
    return result;
}

SymExprPtr PortableExpression::subs(const std::map<SymExprPtr, SymExprPtr>& exprMap) const
{
    for (const auto& pair : exprMap) {
        if (to_string() == pair.first->to_string()) {
            as_portable(pair.second);
            return pair.second;
        }
    }
    if (_kind == Kind::Constant || _kind == Kind::Symbol) {
        return clone();
    }

    std::vector<SymExprPtr> new_args;
    for (const auto& arg : _args) {
        new_args.push_back(arg->subs(exprMap));
    }
    return std::make_shared<PortableExpression>(_kind, _op, new_args);
}

double PortableExpression::eval() const
{
    if (_kind == Kind::Constant) return _value;
    if (_kind == Kind::Symbol) {
        throw std::runtime_error("portable backend cannot evaluate unresolved symbol: " + _op);
    }
    if (_kind == Kind::Unary) {
        double x = _args[0]->eval();
        if (_op == "-") return -x;
    }
    if (_kind == Kind::Binary) {
        double a = _args[0]->eval();
        double b = _args[1]->eval();
        if (_op == "+") return a + b;
        if (_op == "-") return a - b;
        if (_op == "*") return a * b;
        if (_op == "/") return a / b;
        if (_op == "**") return std::pow(a, b);
    }
    if (_kind == Kind::Function) {
        return evaluate_function(_op, _args[0]->eval());
    }
    throw std::runtime_error("portable backend cannot evaluate expression");
}

bool PortableExpression::has_symbols() const
{
    if (_kind == Kind::Symbol) return true;
    for (const auto& arg : _args) {
        auto portable = as_portable(arg);
        if (portable->has_symbols()) return true;
    }
    return false;
}

bool PortableExpression::is_numeric() const
{
    return !has_symbols();
}

SymExprPtr PortableExpression::series(const SymExprPtr&, const SymExprPtr&, int) const
{
    throw std::runtime_error("portable backend does not support series expansion");
}

SymExprPtr PortableExpression::series_to_poly(const SymExprPtr&) const
{
    throw std::runtime_error("portable backend does not support series-to-polynomial conversion");
}

SymExprPtr PortableExpression::coeff(const SymExprPtr&, int) const
{
    throw std::runtime_error("portable backend does not support coefficient extraction");
}

SymExprPtr PortableExpression::evalf() const
{
    if (is_numeric()) {
        return std::make_shared<PortableExpression>(eval());
    }
    return clone();
}

bool PortableExpression::is_zero() const
{
    return is_numeric() && std::fabs(eval()) < 1e-14;
}

double PortableExpression::to_double() const
{
    return eval();
}

std::string PortableExpression::render_python() const
{
    if (_kind == Kind::Constant) return double_to_string(_value);
    if (_kind == Kind::Symbol) return _op;
    if (_kind == Kind::Unary) return "(-" + as_portable(_args[0])->render_python() + ")";
    if (_kind == Kind::Binary) {
        return "(" + as_portable(_args[0])->render_python() + " " + _op + " " +
               as_portable(_args[1])->render_python() + ")";
    }
    if (_kind == Kind::Function) {
        std::string arg = as_portable(_args[0])->render_python();
        if (_op == "bessel_j0") return "besselj(0, " + arg + ")";
        if (_op == "bessel_j1") return "besselj(1, " + arg + ")";
        return sympy_function_name(_op) + "(" + arg + ")";
    }
    return "";
}

std::string PortableExpression::render_plain() const
{
    return render_python();
}

std::string PortableExpression::to_string() const { return render_plain(); }
std::string PortableExpression::to_latex() const { return render_python(); }
std::string PortableExpression::to_python() const { return render_python(); }
std::string PortableExpression::to_cform() const { return render_python(); }

SymbolicCapabilities PortableFactory::capabilities() const
{
    SymbolicCapabilities caps;
    caps.numeric_evaluation = true;
    caps.python_output = true;
    return caps;
}

SymExprPtr PortableFactory::createSymbol(const std::string& name)
{
    return std::make_shared<PortableExpression>(name);
}

SymExprPtr PortableFactory::createConstant(double value)
{
    return std::make_shared<PortableExpression>(value);
}

SymExprPtr PortableFactory::createPi()
{
    return std::make_shared<PortableExpression>(3.141592653589793238462643383279502884);
}

SymExprPtr PortableFactory::createE()
{
    return std::make_shared<PortableExpression>(2.718281828459045235360287471352662498);
}

SymExprPtr PortableFactory::createI()
{
    throw std::runtime_error("portable backend does not support complex unit I");
}

SymbolicFactory* portableFactory()
{
    static PortableFactory factory;
    return &factory;
}

void registerPortableBackend()
{
    SymbolicFactory::registerBackend("portable", portableFactory());
}

} // namespace sebsym
