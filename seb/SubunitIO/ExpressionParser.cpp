#include "SubunitIO/ExpressionParser.hpp"

#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_erf.h>

#include "Exceptions.hpp"
#include "ExpressionFunctions.hpp"
#include "SpecialFunctions.hpp"

namespace pyseb {
namespace {

using Node = ParsedExpressionNode;
using NodePtr = std::shared_ptr<const Node>;

NodePtr numberNode(double value) {
    auto node = std::make_shared<Node>();
    node->kind = Node::Kind::Number;
    node->number = value;
    return node;
}

NodePtr textNode(Node::Kind kind, const std::string& text, std::vector<NodePtr> children = {}) {
    auto node = std::make_shared<Node>();
    node->kind = kind;
    node->text = text;
    node->children = std::move(children);
    return node;
}

class Parser {
public:
    Parser(const std::string& input, const std::string& location)
        : input_(input), location_(location) {}

    ParsedExpression parse() {
        skipWhitespace();
        if (atEnd()) fail("expression is empty");
        NodePtr root = parseAdditive();
        skipWhitespace();
        if (!atEnd()) fail(std::string("unexpected character '") + input_[position_] + "'");
        return ParsedExpression(input_, root, identifiers_);
    }

private:
    const std::string& input_;
    std::string location_;
    std::size_t position_ = 0;
    std::set<std::string> identifiers_;

    bool atEnd() const { return position_ >= input_.size(); }

    void skipWhitespace() {
        while (!atEnd() && std::isspace(static_cast<unsigned char>(input_[position_]))) {
            ++position_;
        }
    }

    bool consume(char token) {
        skipWhitespace();
        if (!atEnd() && input_[position_] == token) {
            ++position_;
            return true;
        }
        return false;
    }

    [[noreturn]] void fail(const std::string& message) const {
        std::ostringstream full;
        full << location_ << ": " << message << " at character " << (position_ + 1);
        throw SEBException(full.str(), "ParseSubunitExpression()");
    }

    NodePtr parseAdditive() {
        NodePtr lhs = parseMultiplicative();
        while (true) {
            if (consume('+')) {
                lhs = textNode(Node::Kind::Binary, "+", {lhs, parseMultiplicative()});
            } else if (consume('-')) {
                lhs = textNode(Node::Kind::Binary, "-", {lhs, parseMultiplicative()});
            } else {
                return lhs;
            }
        }
    }

    NodePtr parseMultiplicative() {
        NodePtr lhs = parseUnary();
        while (true) {
            if (consume('*')) {
                if (consume('*')) fail("use '^' for exponentiation");
                lhs = textNode(Node::Kind::Binary, "*", {lhs, parseUnary()});
            } else if (consume('/')) {
                lhs = textNode(Node::Kind::Binary, "/", {lhs, parseUnary()});
            } else {
                return lhs;
            }
        }
    }

    NodePtr parseUnary() {
        if (consume('+')) return parseUnary();
        if (consume('-')) return textNode(Node::Kind::Unary, "-", {parseUnary()});
        return parsePower();
    }

    NodePtr parsePower() {
        NodePtr base = parsePrimary();
        if (consume('^')) {
            return textNode(Node::Kind::Binary, "^", {base, parseUnary()});
        }
        return base;
    }

    NodePtr parsePrimary() {
        skipWhitespace();
        if (atEnd()) fail("expected a number, identifier, or '('");

        if (consume('(')) {
            NodePtr nested = parseAdditive();
            if (!consume(')')) fail("expected ')'");
            return nested;
        }

        const unsigned char current = static_cast<unsigned char>(input_[position_]);
        if (std::isdigit(current) || input_[position_] == '.') return parseNumber();
        if (std::isalpha(current) || input_[position_] == '_') return parseIdentifierOrFunction();
        fail(std::string("unexpected character '") + input_[position_] + "'");
    }

    NodePtr parseNumber() {
        const char* begin = input_.c_str() + position_;
        char* end = nullptr;
        errno = 0;
        const double value = std::strtod(begin, &end);
        if (end == begin || errno == ERANGE) fail("invalid numeric literal");
        position_ += static_cast<std::size_t>(end - begin);
        return numberNode(value);
    }

    std::string parseIdentifier() {
        const std::size_t begin = position_;
        ++position_;
        while (!atEnd()) {
            const unsigned char ch = static_cast<unsigned char>(input_[position_]);
            if (!std::isalnum(ch) && input_[position_] != '_') break;
            ++position_;
        }
        return input_.substr(begin, position_ - begin);
    }

    NodePtr parseIdentifierOrFunction() {
        const std::string name = parseIdentifier();
        skipWhitespace();
        if (!consume('(')) {
            identifiers_.insert(name);
            return textNode(Node::Kind::Identifier, name);
        }

        if (!IsSubunitExpressionFunction(name)) {
            fail("unsupported function '" + name + "'");
        }

        std::vector<NodePtr> arguments;
        if (!consume(')')) {
            do {
                arguments.push_back(parseAdditive());
            } while (consume(','));
            if (!consume(')')) fail("expected ')' after function arguments");
        }

        const std::size_t expected = name == "pow" ? 2u : 1u;
        if (arguments.size() != expected) {
            fail("function '" + name + "' expects " + std::to_string(expected) + " argument(s)");
        }
        return textNode(Node::Kind::Function, name, std::move(arguments));
    }
};

sebsym::Expression materialize(
    const NodePtr& node,
    const std::function<sebsym::Expression(const std::string&)>& resolve) {
    if (!node) throw SEBException("Cannot materialize an empty expression");

    switch (node->kind) {
        case Node::Kind::Number:
            return sebsym::constant(node->number);
        case Node::Kind::Identifier:
            if (node->text == "pi") return sebsym::pi();
            if (node->text == "e") return sebsym::e();
            return resolve(node->text);
        case Node::Kind::Unary:
            return -materialize(node->children.at(0), resolve);
        case Node::Kind::Binary: {
            const sebsym::Expression lhs = materialize(node->children.at(0), resolve);
            const sebsym::Expression rhs = materialize(node->children.at(1), resolve);
            if (node->text == "+") return lhs + rhs;
            if (node->text == "-") return lhs - rhs;
            if (node->text == "*") return lhs * rhs;
            if (node->text == "/") return lhs / rhs;
            if (node->text == "^") return sebsym::pow(lhs, rhs);
            break;
        }
        case Node::Kind::Function: {
            const sebsym::Expression arg = materialize(node->children.at(0), resolve);
            if (node->text == "pow") {
                return sebsym::pow(arg, materialize(node->children.at(1), resolve));
            }
            if (node->text == "exp") return arg.exp();
            if (node->text == "log") return arg.log();
            if (node->text == "sqrt") return arg.sqrt();
            if (node->text == "abs") return arg.abs();
            if (node->text == "sin") return arg.sin();
            if (node->text == "cos") return arg.cos();
            if (node->text == "tan") return arg.tan();
            if (node->text == "asin") return arg.asin();
            if (node->text == "acos") return arg.acos();
            if (node->text == "atan") return arg.atan();
            if (node->text == "sinh") return arg.sinh();
            if (node->text == "cosh") return arg.cosh();
            if (node->text == "tanh") return arg.tanh();
            if (node->text == "erf") return arg.erf();
            if (node->text == "erfc") return arg.erfc();
            if (node->text == "bessel_j0") return BesselJ0(arg);
            if (node->text == "bessel_j1") return BesselJ1(arg);
            if (node->text == "dawson") return DawsonF(arg);
            if (node->text == "six") return Six(arg);
            if (node->text == "sinc") return Sinc(arg);
            if (node->text == "jinc") return Jinc(arg);
            if (node->text == "struve_h0") return StruveH0(arg);
            if (node->text == "struve_h1") return StruveH1(arg);
            break;
        }
    }
    throw SEBException("Unsupported parsed expression node", "MaterializeSubunitExpression()");
}

double evaluate(
    const NodePtr& node,
    const std::function<double(const std::string&)>& resolve) {
    if (!node) throw SEBException("Cannot evaluate an empty expression");

    switch (node->kind) {
        case Node::Kind::Number:
            return node->number;
        case Node::Kind::Identifier:
            if (node->text == "pi") return std::acos(-1.0);
            if (node->text == "e") return std::exp(1.0);
            return resolve(node->text);
        case Node::Kind::Unary:
            return -evaluate(node->children.at(0), resolve);
        case Node::Kind::Binary: {
            const double lhs = evaluate(node->children.at(0), resolve);
            const double rhs = evaluate(node->children.at(1), resolve);
            if (node->text == "+") return lhs + rhs;
            if (node->text == "-") return lhs - rhs;
            if (node->text == "*") return lhs * rhs;
            if (node->text == "/") return lhs / rhs;
            if (node->text == "^") return std::pow(lhs, rhs);
            break;
        }
        case Node::Kind::Function: {
            const double arg = evaluate(node->children.at(0), resolve);
            if (node->text == "pow") {
                return std::pow(arg, evaluate(node->children.at(1), resolve));
            }
            if (node->text == "exp") return std::exp(arg);
            if (node->text == "log") return std::log(arg);
            if (node->text == "sqrt") return std::sqrt(arg);
            if (node->text == "abs") return std::abs(arg);
            if (node->text == "sin") return std::sin(arg);
            if (node->text == "cos") return std::cos(arg);
            if (node->text == "tan") return std::tan(arg);
            if (node->text == "asin") return std::asin(arg);
            if (node->text == "acos") return std::acos(arg);
            if (node->text == "atan") return std::atan(arg);
            if (node->text == "sinh") return std::sinh(arg);
            if (node->text == "cosh") return std::cosh(arg);
            if (node->text == "tanh") return std::tanh(arg);
            if (node->text == "erf") return gsl_sf_erf(arg);
            if (node->text == "erfc") return gsl_sf_erfc(arg);
            if (node->text == "bessel_j0") return gsl_sf_bessel_J0(arg);
            if (node->text == "bessel_j1") return gsl_sf_bessel_J1(arg);
            if (node->text == "dawson") return gsl_sf_dawson(arg);
            if (node->text == "six") {
                const double squared = arg * arg;
                if (squared < 1e-8) {
                    return 1.0 - squared / 18.0 + squared * squared / 600.0;
                }
                return gsl_sf_Si(arg) / arg;
            }
            if (node->text == "sinc") {
                const double squared = arg * arg;
                if (squared < 1e-8) {
                    return 1.0 - squared / 6.0 + squared * squared / 120.0;
                }
                return std::sin(arg) / arg;
            }
            if (node->text == "jinc") {
                const double squared = arg * arg;
                if (squared < 1e-8) {
                    return 1.0 - squared / 8.0 + squared * squared / 192.0;
                }
                return 2.0 * gsl_sf_bessel_J1(arg) / arg;
            }
            if (node->text == "struve_h0") return StruveH0Value(arg);
            if (node->text == "struve_h1") return StruveH1Value(arg);
            break;
        }
    }
    throw SEBException(
        "Unsupported parsed expression node during numerical evaluation",
        "EvaluateSubunitExpression()");
}

} // namespace

ParsedExpression::ParsedExpression(
    std::string source,
    std::shared_ptr<const ParsedExpressionNode> root,
    std::set<std::string> identifiers)
    : source_(std::move(source)), root_(std::move(root)), identifiers_(std::move(identifiers)) {}

ParsedExpression ParseSubunitExpression(const std::string& expression, const std::string& location) {
    return Parser(expression, location).parse();
}

sebsym::Expression MaterializeSubunitExpression(
    const ParsedExpression& expression,
    const std::function<sebsym::Expression(const std::string&)>& resolveIdentifier) {
    return materialize(expression.root(), resolveIdentifier);
}

double EvaluateSubunitExpression(
    const ParsedExpression& expression,
    const std::function<double(const std::string&)>& resolveIdentifier) {
    return evaluate(expression.root(), resolveIdentifier);
}

bool IsSubunitExpressionFunction(const std::string& name) {
    static const std::set<std::string> functions = {
        "abs", "acos", "asin", "atan", "bessel_j0", "bessel_j1", "cos", "cosh",
        "dawson", "erf", "erfc", "exp", "jinc", "log", "pow", "sin", "sinc", "sinh", "six", "sqrt",
        "struve_h0", "struve_h1", "tan", "tanh"
    };
    return functions.find(name) != functions.end();
}

} // namespace pyseb
