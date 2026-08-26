#ifndef PYSEB_SUBUNIT_IO_EXPRESSION_PARSER_HPP
#define PYSEB_SUBUNIT_IO_EXPRESSION_PARSER_HPP

#include <functional>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "Expression.hpp"

namespace pyseb {

struct ParsedExpressionNode {
    enum class Kind { Number, Identifier, Unary, Binary, Function };

    Kind kind = Kind::Number;
    double number = 0.0;
    std::string text;
    std::vector<std::shared_ptr<const ParsedExpressionNode>> children;
};

class ParsedExpression {
public:
    ParsedExpression() = default;
    ParsedExpression(
        std::string source,
        std::shared_ptr<const ParsedExpressionNode> root,
        std::set<std::string> identifiers);

    const std::string& source() const { return source_; }
    const std::shared_ptr<const ParsedExpressionNode>& root() const { return root_; }
    const std::set<std::string>& identifiers() const { return identifiers_; }
    bool empty() const { return !root_; }

private:
    std::string source_;
    std::shared_ptr<const ParsedExpressionNode> root_;
    std::set<std::string> identifiers_;
};

ParsedExpression ParseSubunitExpression(
    const std::string& expression,
    const std::string& location = "expression");

sebsym::Expression MaterializeSubunitExpression(
    const ParsedExpression& expression,
    const std::function<sebsym::Expression(const std::string&)>& resolveIdentifier);

bool IsSubunitExpressionFunction(const std::string& name);

} // namespace pyseb

#endif
