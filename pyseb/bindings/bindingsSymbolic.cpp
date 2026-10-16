#include "bindingsSymbolic.hpp"

#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "SEB.hpp"
#include "Subunit.hpp"
#include "Subunits/Subunits.hpp"
#include "SymbolicPortable.hpp"

#include <unordered_map>

namespace {

// Transfer the portable expression tree directly. Parsing its printed form
// duplicates shared subexpressions and hits Python's parenthesis-depth limit.
// Keep memoization local to this conversion so Python objects never outlive
// the interpreter through a static cache.
class PortableSympyConverter {
    py::object sympy = py::module::import("sympy");
    py::dict functions = py::module::import("pyseb.symbolic").attr("_SYMPY_LOCALS");
    std::unordered_map<const SymbolicExpression*, py::object> converted;
    std::unordered_map<double, py::object> numbers;
    std::unordered_map<std::string, py::object> symbols;

    void collect(const SymExprPtr& expression, const std::string& operation,
                 std::vector<py::object>& arguments)
    {
        const auto found = converted.find(expression.get());
        if (found != converted.end()) {
            arguments.push_back(found->second);
            return;
        }
        const auto node = std::dynamic_pointer_cast<sebsym::PortableExpression>(expression);
        if (node && node->kind() == sebsym::PortableExpression::Kind::Binary && node->op() == operation) {
            for (const auto& child : node->args()) collect(child, operation, arguments);
        } else arguments.push_back(convert(expression));
    }

public:
    py::object convert(const SymExprPtr& expression)
    {
        const auto found = converted.find(expression.get());
        if (found != converted.end()) return found->second;
        const auto node = std::dynamic_pointer_cast<sebsym::PortableExpression>(expression);
        if (!node) throw std::runtime_error("Expected a portable expression tree");
        using Kind = sebsym::PortableExpression::Kind;
        py::object result;
        if (node->kind() == Kind::Constant) {
            // Preserve the existing printer's integer/17-digit float semantics.
            // Parsing one number is independent of expression size or nesting.
            const auto number = numbers.find(node->value());
            if (number != numbers.end()) result = number->second;
            else {
                result = sympy.attr("sympify")(node->to_python());
                numbers.emplace(node->value(), result);
            }
        } else if (node->kind() == Kind::Symbol) {
            const auto symbol = symbols.find(node->op());
            if (symbol != symbols.end()) result = symbol->second;
            else {
                result = sympy.attr("Symbol")(node->op());
                symbols.emplace(node->op(), result);
            }
        } else if (node->kind() == Kind::Binary && (node->op() == "+" || node->op() == "*")) {
            // Construct an associative operation once, rather than repeatedly
            // asking SymPy to flatten and canonicalize growing binary results.
            std::vector<py::object> arguments;
            for (const auto& child : node->args()) collect(child, node->op(), arguments);
            py::tuple values(arguments.size());
            for (std::size_t i=0; i<arguments.size(); ++i) values[i] = arguments[i];
            result = sympy.attr(node->op() == "+" ? "Add" : "Mul")(*values);
        } else {
            std::vector<py::object> args;
            for (const auto& child : node->args()) args.push_back(convert(child));
            const auto& op = node->op();
            if (node->kind() == Kind::Unary) {
                result = args.at(0).attr("__neg__")();
            } else if (node->kind() == Kind::Binary) {
                static const std::map<std::string, std::string> methods = {
                    {"+", "__add__"}, {"-", "__sub__"}, {"*", "__mul__"},
                    {"/", "__truediv__"}, {"**", "__pow__"}};
                result = args.at(0).attr(methods.at(op).c_str())(args.at(1));
            } else if (op == "Integral") {
                result = sympy.attr("Integral")(args.at(3),
                    py::make_tuple(args.at(0), args.at(1), args.at(2)));
            } else if (op == "hyper0f1_regularized") {
                result = sympy.attr("hyper")(py::tuple(), py::make_tuple(args.at(0)), args.at(1))
                    .attr("__truediv__")(sympy.attr("gamma")(args.at(0)));
            } else if (op == "struve_h0" || op == "struve_h1") {
                result = functions["struve"](op == "struve_h0" ? 0 : 1, args.at(0));
            } else if (op == "bessel_j0" || op == "bessel_j1") {
                result = sympy.attr("besselj")(op == "bessel_j0" ? 0 : 1, args.at(0));
            } else {
                if (args.size() != 1) throw std::runtime_error("Unsupported portable function arity: " + op);
                const py::str name(op == "abs" ? "Abs" : op);
                py::object function = functions.contains(name)
                    ? py::reinterpret_borrow<py::object>(functions[name])
                    : sympy.attr("Function")(name);
                result = function(args.at(0));
            }
        }
        converted.emplace(expression.get(), result);
        return result;
    }
};

} // namespace

py::object pyseb_expression_to_sympy(const Expression& expr)
{
    if (std::dynamic_pointer_cast<sebsym::PortableExpression>(expr.get())) {
        PortableSympyConverter converter;
        return py::module::import("pyseb.symbolic").attr("SymPyExpression")(converter.convert(expr.get()));
    }
    return py::module::import("pyseb.symbolic").attr("SymPyExpression")(expr.to_python());
}

py::object pyseb_zero_sympy_expression()
{
    return py::module::import("pyseb.symbolic").attr("SymPyExpression")("0");
}

void register_symbolic_world_bindings(py::class_<World>& world)
{
    world.def("FormFactor", [](World& self, const std::string& myself, int depth, int varForm) {
            return pyseb_expression_to_sympy(self.FormFactor(myself, depth, varForm));
        }, py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = static_cast<int>(QVAR))
        .def("FormFactor_Unnormalized", [](World& self, const std::string& myself, int depth, int varForm) {
            return pyseb_expression_to_sympy(self.FormFactor_Unnormalized(myself, depth, varForm));
        }, py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = static_cast<int>(QVAR))
        .def("FormFactorGeneric", [](World& self, const std::string& myself, int depth) {
            return pyseb_expression_to_sympy(self.FormFactorGeneric(myself, depth));
        }, py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH)
        .def("FormFactorAmplitude", [](World& self, refPoint& ref, int depth, int varForm) {
            return pyseb_expression_to_sympy(self.FormFactorAmplitude(ref, depth, varForm));
        }, py::arg("ref"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = static_cast<int>(QVAR))
        .def("PhaseFactor", [](World& self, refPoint& ref1, refPoint& ref2, int depth) {
            return pyseb_expression_to_sympy(self.PhaseFactor(ref1, ref2, depth));
        }, py::arg("ref1"), py::arg("ref2"), py::arg("depth") = WORLDMAXDEPTH)
        .def("RadiusOfGyration2", [](World& self, const std::string& myself) {
            return pyseb_expression_to_sympy(self.RadiusOfGyration2(myself));
        }, py::arg("myself"))
        .def("SMSD_ref2ref", [](World& self, refPoint& ref1, refPoint& ref2) {
            return pyseb_expression_to_sympy(self.SMSD_ref2ref(ref1, ref2));
        }, py::arg("ref1"), py::arg("ref2"))
        .def("SMSD_ref2scat", [](World& self, refPoint& ref) {
            return pyseb_expression_to_sympy(self.SMSD_ref2scat(ref));
        }, py::arg("ref"))
        .def("SMSD_scat2scat", [](World& self) {
            (void)self;
            return pyseb_zero_sympy_expression();
        });
}
