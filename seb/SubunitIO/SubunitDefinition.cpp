#include "SubunitIO/SubunitDefinition.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>

#include "Exceptions.hpp"
#include "fkYAML/node.hpp"

namespace pyseb {
namespace {

using Node = fkyaml::node;
constexpr std::size_t kMaximumModelFileSize = 1024u * 1024u;

void rejectUnknownKeys(
    const Node& node,
    const std::set<std::string>& allowed,
    const std::string& source,
    const std::string& path);

[[noreturn]] void schemaError(
    const std::string& source,
    const std::string& path,
    const std::string& message) {
    throw SEBException(source + ": " + path + ": " + message, "LoadSubunitDefinitionYaml()");
}

void requireMapping(const Node& node, const std::string& source, const std::string& path) {
    if (!node.is_mapping()) schemaError(source, path, "expected a mapping");
}

void requireSequence(const Node& node, const std::string& source, const std::string& path) {
    if (!node.is_sequence()) schemaError(source, path, "expected a sequence");
}

std::string scalarString(const Node& node, const std::string& source, const std::string& path) {
    if (!node.is_string()) schemaError(source, path, "expected a string");
    return node.get_value<std::string>();
}

double scalarDouble(const Node& node, const std::string& source, const std::string& path) {
    if (!node.is_integer() && !node.is_float_number()) {
        schemaError(source, path, "expected a finite number");
    }
    const double value = node.get_value<double>();
    if (!std::isfinite(value)) schemaError(source, path, "expected a finite number");
    return value;
}

std::size_t scalarSize(
    const Node& node,
    const std::string& source,
    const std::string& path) {
    if (!node.is_integer()) schemaError(source, path, "expected a positive integer");
    const long long value = node.get_value<long long>();
    if (value <= 0) schemaError(source, path, "expected a positive integer");
    return static_cast<std::size_t>(value);
}

int qagRuleKey(int points, const std::string& source, const std::string& path) {
    switch (points) {
        case 15: return GSL_INTEG_GAUSS15;
        case 21: return GSL_INTEG_GAUSS21;
        case 31: return GSL_INTEG_GAUSS31;
        case 41: return GSL_INTEG_GAUSS41;
        case 51: return GSL_INTEG_GAUSS51;
        case 61: return GSL_INTEG_GAUSS61;
        default:
            schemaError(source, path, "expected one of 15, 21, 31, 41, 51, or 61");
    }
}

IntegrationOptions integrationOptions(
    const Node& node,
    const IntegrationOptions& inherited,
    const std::string& source,
    const std::string& path) {
    rejectUnknownKeys(
        node,
        {"method", "absolute_tolerance", "relative_tolerance", "workspace_size", "qag_rule"},
        source,
        path);
    IntegrationOptions result = inherited;
    if (node.contains("method")) {
        const std::string method = scalarString(node.at("method"), source, path + ".method");
        if (method == "qag") result.method = IntegrationMethod::QAG;
        else if (method == "cquad") result.method = IntegrationMethod::CQUAD;
        else schemaError(source, path + ".method", "expected 'qag' or 'cquad'");
    }
    if (node.contains("absolute_tolerance")) {
        result.absoluteTolerance = scalarDouble(
            node.at("absolute_tolerance"), source, path + ".absolute_tolerance");
    }
    if (node.contains("relative_tolerance")) {
        result.relativeTolerance = scalarDouble(
            node.at("relative_tolerance"), source, path + ".relative_tolerance");
    }
    if (result.absoluteTolerance < 0.0 || result.relativeTolerance < 0.0 ||
        (result.absoluteTolerance == 0.0 && result.relativeTolerance == 0.0)) {
        schemaError(
            source,
            path,
            "integration tolerances must be non-negative and at least one must be positive");
    }
    if (node.contains("workspace_size")) {
        result.workspaceSize = scalarSize(
            node.at("workspace_size"), source, path + ".workspace_size");
    }
    if (result.workspaceSize < 3) {
        schemaError(source, path + ".workspace_size", "must be at least 3");
    }
    if (node.contains("qag_rule")) {
        if (!node.at("qag_rule").is_integer()) {
            schemaError(source, path + ".qag_rule", "expected an integer");
        }
        result.qagRule = qagRuleKey(
            node.at("qag_rule").get_value<int>(), source, path + ".qag_rule");
    }
    return result;
}

void rejectUnknownKeys(
    const Node& node,
    const std::set<std::string>& allowed,
    const std::string& source,
    const std::string& path) {
    requireMapping(node, source, path);
    for (const auto& item : node.map_items()) {
        if (!item.key().is_string()) schemaError(source, path, "mapping keys must be strings");
        const std::string key = item.key().get_value<std::string>();
        if (allowed.find(key) == allowed.end()) {
            schemaError(source, path + "." + key, "unknown field");
        }
    }
}

const Node& requiredNode(
    const Node& node,
    const char* key,
    const std::string& source,
    const std::string& path) {
    if (!node.contains(key)) schemaError(source, path + "." + key, "required field is missing");
    return node.at(key);
}

bool validSymbolName(const std::string& name) {
    if (name.empty() || !std::isalpha(static_cast<unsigned char>(name.front()))) return false;
    return std::all_of(name.begin() + 1, name.end(), [](char ch) {
        return std::isalnum(static_cast<unsigned char>(ch)) != 0;
    });
}

bool validModelId(const std::string& id) {
    const std::size_t slash = id.find('/');
    if (slash == std::string::npos || slash == 0 || slash + 1 == id.size()) return false;
    if (id.find('/', slash + 1) != std::string::npos) return false;
    for (char ch : id) {
        const unsigned char value = static_cast<unsigned char>(ch);
        if (!std::isalnum(value) && ch != '/' && ch != '-' && ch != '_' && ch != '.') return false;
    }
    return true;
}

void validateNodeRestrictions(
    const Node& node,
    const std::string& source,
    const std::string& path,
    std::size_t depth = 0) {
    if (depth > 64) schemaError(source, path, "YAML nesting exceeds 64 levels");
    if (node.is_anchor() || node.is_alias()) schemaError(source, path, "YAML anchors and aliases are not supported");
    if (node.has_tag_name()) schemaError(source, path, "YAML tags are not supported");
    if (node.is_sequence()) {
        std::size_t index = 0;
        for (const auto& child : node) {
            validateNodeRestrictions(child, source, path + "[" + std::to_string(index++) + "]", depth + 1);
        }
    } else if (node.is_mapping()) {
        for (const auto& item : node.map_items()) {
            validateNodeRestrictions(item.key(), source, path + ".<key>", depth + 1);
            const std::string key = item.key().is_string() ? item.key().get_value<std::string>() : "<value>";
            validateNodeRestrictions(item.value(), source, path + "." + key, depth + 1);
        }
    }
}

bool isDocumentIndicator(const std::string& line, const char* indicator) {
    if (line.compare(0, 3, indicator) != 0) return false;
    return line.size() == 3 || std::isspace(static_cast<unsigned char>(line[3])) || line[3] == '#';
}

void rejectMultipleDocuments(const std::string& yaml, const std::string& source) {
    std::istringstream input(yaml);
    std::string line;
    std::size_t lineNumber = 0;
    bool seenContent = false;
    bool ended = false;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::size_t first = line.find_first_not_of(" \t");
        if (first == std::string::npos || line[first] == '#') continue;

        // YAML document indicators are only recognized at column zero. This
        // deliberately leaves indented "---" text in block scalars alone.
        if (first == 0 && isDocumentIndicator(line, "---")) {
            if (seenContent || ended) {
                schemaError(source, "$", "multiple YAML documents are not supported (line " +
                    std::to_string(lineNumber) + ", column 1)");
            }
            const std::size_t inlineContent = line.find_first_not_of(" \t", 3);
            if (inlineContent != std::string::npos && line[inlineContent] != '#') seenContent = true;
            continue;
        }
        if (first == 0 && isDocumentIndicator(line, "...")) {
            ended = true;
            continue;
        }
        if (ended) {
            schemaError(source, "$", "multiple YAML documents are not supported (line " +
                std::to_string(lineNumber) + ", column 1)");
        }
        seenContent = true;
    }
}

std::vector<std::string> stringSequence(
    const Node& node,
    const std::string& source,
    const std::string& path) {
    requireSequence(node, source, path);
    std::vector<std::string> result;
    std::size_t index = 0;
    for (const auto& value : node) {
        result.push_back(scalarString(value, source, path + "[" + std::to_string(index++) + "]"));
    }
    return result;
}

ParsedExpression parsedExpression(
    const Node& node,
    const std::string& source,
    const std::string& path) {
    return ParseSubunitExpression(scalarString(node, source, path), source + ": " + path);
}

std::map<std::string, ParsedExpression> expressionMap(
    const Node& node,
    const std::string& source,
    const std::string& path) {
    requireMapping(node, source, path);
    std::map<std::string, ParsedExpression> result;
    for (const auto& item : node.map_items()) {
        const std::string key = scalarString(item.key(), source, path + ".<key>");
        if (!validSymbolName(key)) schemaError(source, path + "." + key, "invalid identifier");
        result.emplace(key, parsedExpression(item.value(), source, path + "." + key));
    }
    return result;
}

std::map<std::string, IntegralDefinition> integralMap(
    const Node& node,
    const IntegrationOptions& inherited,
    const std::string& source,
    const std::string& path) {
    requireMapping(node, source, path);
    std::map<std::string, IntegralDefinition> result;
    for (const auto& item : node.map_items()) {
        IntegralDefinition integral;
        integral.name = scalarString(item.key(), source, path + ".<key>");
        const std::string itemPath = path + "." + integral.name;
        if (!validSymbolName(integral.name)) {
            schemaError(source, itemPath, "invalid identifier");
        }
        rejectUnknownKeys(
            item.value(),
            {"variable", "lower", "upper", "variables", "integrand", "integration"},
            source,
            itemPath);
        const bool hasLegacyDimension = item.value().contains("variable") ||
            item.value().contains("lower") || item.value().contains("upper");
        const bool hasDimensions = item.value().contains("variables");
        if (hasLegacyDimension && hasDimensions) {
            schemaError(source, itemPath, "use either the scalar integral fields or variables, not both");
        }
        if (hasLegacyDimension) {
            if (!item.value().contains("variable") || !item.value().contains("lower") ||
                !item.value().contains("upper")) {
                schemaError(source, itemPath, "variable, lower, and upper are required together");
            }
            IntegralDefinition::Dimension dimension;
            dimension.variable = scalarString(
                item.value().at("variable"), source, itemPath + ".variable");
            if (!validSymbolName(dimension.variable)) {
                schemaError(source, itemPath + ".variable", "invalid identifier");
            }
            dimension.lower = parsedExpression(
                item.value().at("lower"), source, itemPath + ".lower");
            dimension.upper = parsedExpression(
                item.value().at("upper"), source, itemPath + ".upper");
            integral.dimensions.push_back(std::move(dimension));
            integral.variable = integral.dimensions.front().variable;
            integral.lower = integral.dimensions.front().lower;
            integral.upper = integral.dimensions.front().upper;
        } else if (hasDimensions) {
            const Node& dimensions = item.value().at("variables");
            requireSequence(dimensions, source, itemPath + ".variables");
            std::size_t index = 0;
            for (const auto& value : dimensions) {
                const std::string dimensionPath =
                    itemPath + ".variables[" + std::to_string(index++) + "]";
                rejectUnknownKeys(value, {"variable", "lower", "upper"}, source, dimensionPath);
                IntegralDefinition::Dimension dimension;
                dimension.variable = scalarString(
                    requiredNode(value, "variable", source, dimensionPath),
                    source,
                    dimensionPath + ".variable");
                if (!validSymbolName(dimension.variable)) {
                    schemaError(source, dimensionPath + ".variable", "invalid identifier");
                }
                dimension.lower = parsedExpression(
                    requiredNode(value, "lower", source, dimensionPath),
                    source,
                    dimensionPath + ".lower");
                dimension.upper = parsedExpression(
                    requiredNode(value, "upper", source, dimensionPath),
                    source,
                    dimensionPath + ".upper");
                integral.dimensions.push_back(std::move(dimension));
            }
            if (integral.dimensions.empty()) {
                schemaError(source, itemPath + ".variables", "expected at least one dimension");
            }
        } else {
            schemaError(source, itemPath, "variable/lower/upper or variables is required");
        }
        integral.integrand = parsedExpression(
            requiredNode(item.value(), "integrand", source, itemPath),
            source,
            itemPath + ".integrand");
        integral.integration = item.value().contains("integration")
            ? integrationOptions(
                item.value().at("integration"), inherited, source, itemPath + ".integration")
            : inherited;
        result.emplace(integral.name, std::move(integral));
    }
    return result;
}

ReferencePair parseReferencePair(
    const Node& node,
    const std::string& source,
    const std::string& path) {
    const auto values = stringSequence(node, source, path);
    if (values.size() != 2) schemaError(source, path, "expected exactly two reference names");
    return CanonicalReferencePair(values[0], values[1]);
}

std::map<ReferencePair, ParsedExpression> pairExpressionSequence(
    const Node& node,
    const std::string& source,
    const std::string& path) {
    requireSequence(node, source, path);
    std::map<ReferencePair, ParsedExpression> result;
    std::size_t index = 0;
    for (const auto& item : node) {
        const std::string itemPath = path + "[" + std::to_string(index++) + "]";
        rejectUnknownKeys(item, {"references", "expression"}, source, itemPath);
        const ReferencePair pair = parseReferencePair(
            requiredNode(item, "references", source, itemPath), source, itemPath + ".references");
        const ParsedExpression expression = parsedExpression(
            requiredNode(item, "expression", source, itemPath), source, itemPath + ".expression");
        if (!result.emplace(pair, expression).second) {
            schemaError(source, itemPath + ".references", "duplicate unordered reference pair");
        }
    }
    return result;
}

void validateUniqueNames(SubunitDefinition& definition) {
    const std::string& source = definition.source;
    const std::set<std::string> reserved = {
        "q", "pi", "e", "abs", "acos", "asin", "atan", "bessel_j0", "bessel_j1",
        "cos", "cosh", "dawson", "erf", "erfc", "exp", "jinc", "log", "pow", "sin",
        "sinc", "sinh", "six", "sqrt", "struve_h0", "struve_h1", "tan", "tanh"
    };
    std::set<std::string> names = reserved;
    const auto insertName = [&](const std::string& name, const std::string& path) {
        if (!validSymbolName(name)) schemaError(source, path, "names must start with a letter and contain only letters and digits");
        if (!names.insert(name).second) schemaError(source, path, "name is duplicated or reserved");
    };

    for (const auto& parameter : definition.parameters) insertName(parameter.first, "parameters." + parameter.first);
    for (const auto& variable : definition.variables) insertName(variable.first, "variables." + variable.first);
    for (const auto& named : definition.definitions) insertName(named.first, "definitions." + named.first);
    for (const auto& integral : definition.integrals) insertName(integral.first, "integrals." + integral.first);
    for (const auto& integral : definition.integrals) {
        std::set<std::string> dimensionNames;
        std::size_t index = 0;
        for (const auto& dimension : integral.second.dimensions) {
            const std::string path = "integrals." + integral.first + ".variables[" +
                std::to_string(index++) + "].variable";
            if (!validSymbolName(dimension.variable)) {
                schemaError(source, path, "names must start with a letter and contain only letters and digits");
            }
            if (names.count(dimension.variable)) {
                schemaError(source, path, "bound variable collides with a model-level or reserved name");
            }
            if (!dimensionNames.insert(dimension.variable).second) {
                schemaError(source, path, "bound variable is duplicated");
            }
        }
    }

    std::set<std::string> references;
    for (const auto& ref : definition.specificReferences) {
        if (!validSymbolName(ref)) schemaError(source, "references.specific", "invalid reference name '" + ref + "'");
        if (!references.insert(ref).second) schemaError(source, "references", "duplicate reference '" + ref + "'");
    }
    for (const auto& ref : definition.distributedReferences) {
        if (!validSymbolName(ref)) schemaError(source, "references.distributed", "invalid reference name '" + ref + "'");
        if (!references.insert(ref).second) schemaError(source, "references", "duplicate reference '" + ref + "'");
    }
    if (references.empty()) schemaError(source, "references", "at least one reference is required");
}

void validateExpressionDependencies(const SubunitDefinition& definition) {
    const std::string& source = definition.source;
    std::set<std::string> parameters;
    for (const auto& value : definition.parameters) parameters.insert(value.first);
    std::set<std::string> variables;
    for (const auto& value : definition.variables) variables.insert(value.first);
    std::set<std::string> namedDefinitions;
    for (const auto& value : definition.definitions) namedDefinitions.insert(value.first);
    std::set<std::string> integrals;
    for (const auto& value : definition.integrals) integrals.insert(value.first);

    const auto knownBase = [&](const std::string& id) {
        return id == "q" || id == "pi" || id == "e" || parameters.count(id) ||
            variables.count(id) || namedDefinitions.count(id) || integrals.count(id);
    };
    const auto checkKnown = [&](const ParsedExpression& expression, const std::string& path) {
        for (const auto& id : expression.identifiers()) {
            if (!knownBase(id)) schemaError(source, path, "unknown identifier '" + id + "'");
        }
    };

    for (const auto& variable : definition.variables) {
        for (const auto& id : variable.second.identifiers()) {
            if (namedDefinitions.count(id)) {
                schemaError(source, "variables." + variable.first, "variables cannot depend on reusable definitions");
            }
            if (id != "q" && id != "pi" && id != "e" && !parameters.count(id) && !variables.count(id)) {
                schemaError(source, "variables." + variable.first, "unknown identifier '" + id + "'");
            }
        }
    }
    for (const auto& named : definition.definitions) checkKnown(named.second, "definitions." + named.first);
    checkKnown(definition.formFactor, "expressions.form_factor");
    for (const auto& item : definition.amplitudes) checkKnown(item.second, "expressions.amplitudes." + item.first);
    for (const auto& item : definition.phases) checkKnown(item.second, "expressions.phases");
    checkKnown(definition.radiusOfGyrationSquared, "sizes.radius_of_gyration_squared");
    for (const auto& item : definition.referenceToScatterer) checkKnown(item.second, "sizes.reference_to_scatterer." + item.first);
    for (const auto& item : definition.referenceToReference) checkKnown(item.second, "sizes.reference_to_reference");

    const auto checkIntegralExpression = [&](const ParsedExpression& expression,
                                             const std::set<std::string>& localNames,
                                             const std::string& path) {
        for (const auto& id : expression.identifiers()) {
            if (localNames.count(id)) continue;
            if (!knownBase(id)) schemaError(source, path, "unknown identifier '" + id + "'");
        }
    };
    for (const auto& item : definition.integrals) {
        const IntegralDefinition& integral = item.second;
        const std::string path = "integrals." + item.first;
        std::set<std::string> localNames;
        std::size_t index = 0;
        for (const auto& dimension : integral.dimensions) {
            const std::string dimensionPath = path + ".variables[" +
                std::to_string(index++) + "]";
            checkIntegralExpression(dimension.lower, localNames, dimensionPath + ".lower");
            checkIntegralExpression(dimension.upper, localNames, dimensionPath + ".upper");
            for (const auto& bound : {std::cref(dimension.lower), std::cref(dimension.upper)}) {
                for (const auto& id : bound.get().identifiers()) {
                    if (id == "q") {
                        schemaError(source, dimensionPath, "integration bounds cannot depend on q");
                    }
                    if (id == dimension.variable) {
                        schemaError(source, dimensionPath, "integration bounds cannot depend on the bound variable");
                    }
                    if (integrals.count(id)) {
                        schemaError(source, dimensionPath, "integration bounds cannot depend on another integral");
                    }
                }
            }
            localNames.insert(dimension.variable);
        }
        checkIntegralExpression(integral.integrand, localNames, path + ".integrand");
        for (const auto& id : integral.integrand.identifiers()) {
            if (integrals.count(id)) {
                schemaError(source, path + ".integrand", "nested integrals are not supported");
            }
        }
    }

    enum class Visit { Unvisited, Visiting, Done };
    std::map<std::string, Visit> visits;
    std::function<void(const std::string&)> visit = [&](const std::string& name) {
        if (visits[name] == Visit::Done) return;
        if (visits[name] == Visit::Visiting) schemaError(source, name, "cyclic expression dependency");
        visits[name] = Visit::Visiting;
        const auto variable = definition.variables.find(name);
        const ParsedExpression& expression = variable != definition.variables.end()
            ? variable->second : definition.definitions.at(name);
        for (const auto& id : expression.identifiers()) {
            if (variables.count(id) || namedDefinitions.count(id)) visit(id);
        }
        visits[name] = Visit::Done;
    };
    for (const auto& value : definition.variables) visit(value.first);
    for (const auto& value : definition.definitions) visit(value.first);

    std::map<std::string, bool> integralDependencies;
    std::function<bool(const std::string&)> dependsOnIntegral = [&](const std::string& name) {
        const auto cached = integralDependencies.find(name);
        if (cached != integralDependencies.end()) return cached->second;
        const auto variable = definition.variables.find(name);
        const ParsedExpression& expression = variable != definition.variables.end()
            ? variable->second : definition.definitions.at(name);
        bool depends = false;
        for (const auto& id : expression.identifiers()) {
            if (integrals.count(id)) depends = true;
            if (variables.count(id) || namedDefinitions.count(id)) {
                depends = depends || dependsOnIntegral(id);
            }
        }
        integralDependencies[name] = depends;
        return depends;
    };
    for (const auto& value : definition.variables) dependsOnIntegral(value.first);
    for (const auto& value : definition.definitions) dependsOnIntegral(value.first);

    for (const auto& item : definition.integrals) {
        const IntegralDefinition& integral = item.second;
        const auto requireIntegralFree = [&](const ParsedExpression& expression,
                                             const std::string& path) {
            for (const auto& id : expression.identifiers()) {
                if ((variables.count(id) || namedDefinitions.count(id)) && dependsOnIntegral(id)) {
                    schemaError(source, path, "nested integral dependency through '" + id + "'");
                }
            }
        };
        for (std::size_t index = 0; index < integral.dimensions.size(); ++index) {
            const auto& dimension = integral.dimensions[index];
            const std::string path = "integrals." + item.first + ".variables[" +
                std::to_string(index) + "]";
            requireIntegralFree(dimension.lower, path + ".lower");
            requireIntegralFree(dimension.upper, path + ".upper");
        }
        requireIntegralFree(integral.integrand, "integrals." + item.first + ".integrand");
    }

    std::map<std::string, bool> qDependencies;
    std::function<bool(const std::string&)> dependsOnQ = [&](const std::string& name) {
        const auto cached = qDependencies.find(name);
        if (cached != qDependencies.end()) return cached->second;
        const auto variable = definition.variables.find(name);
        const ParsedExpression& expression = variable != definition.variables.end()
            ? variable->second : definition.definitions.at(name);
        bool depends = expression.identifiers().count("q") != 0;
        for (const auto& id : expression.identifiers()) {
            if (variables.count(id) || namedDefinitions.count(id)) depends = depends || dependsOnQ(id);
        }
        qDependencies[name] = depends;
        return depends;
    };
    for (const auto& item : definition.integrals) {
        const auto requireBoundIndependent = [&](const ParsedExpression& expression,
                                                 const std::string& path) {
            for (const auto& id : expression.identifiers()) {
                if ((variables.count(id) || namedDefinitions.count(id)) && dependsOnQ(id)) {
                    schemaError(source, path, "integration bounds cannot depend on q through '" + id + "'");
                }
            }
        };
        for (std::size_t index = 0; index < item.second.dimensions.size(); ++index) {
            const auto& dimension = item.second.dimensions[index];
            const std::string path = "integrals." + item.first + ".variables[" +
                std::to_string(index) + "]";
            requireBoundIndependent(dimension.lower, path + ".lower");
            requireBoundIndependent(dimension.upper, path + ".upper");
        }
    }
    const auto requireSizeIndependent = [&](const ParsedExpression& expression, const std::string& path) {
        if (expression.identifiers().count("q")) schemaError(source, path, "size expressions cannot depend on q");
        for (const auto& id : expression.identifiers()) {
            if (integrals.count(id)) {
                schemaError(source, path, "size expressions cannot depend on integrals");
            }
            if ((variables.count(id) || namedDefinitions.count(id)) && dependsOnQ(id)) {
                schemaError(source, path, "size expressions cannot depend on q through '" + id + "'");
            }
            if ((variables.count(id) || namedDefinitions.count(id)) && dependsOnIntegral(id)) {
                schemaError(source, path, "size expressions cannot depend on integrals through '" + id + "'");
            }
        }
    };
    requireSizeIndependent(definition.radiusOfGyrationSquared, "sizes.radius_of_gyration_squared");
    for (const auto& item : definition.referenceToScatterer) {
        requireSizeIndependent(item.second, "sizes.reference_to_scatterer." + item.first);
    }
    for (const auto& item : definition.referenceToReference) {
        requireSizeIndependent(item.second, "sizes.reference_to_reference");
    }
}

void validateReferenceCompleteness(const SubunitDefinition& definition) {
    const std::string& source = definition.source;
    std::set<std::string> specific(definition.specificReferences.begin(), definition.specificReferences.end());
    std::set<std::string> references = specific;
    references.insert(definition.distributedReferences.begin(), definition.distributedReferences.end());

    const auto validateSingleMap = [&](const auto& values, const std::string& path) {
        for (const auto& ref : references) {
            if (!values.count(ref)) schemaError(source, path + "." + ref, "required expression is missing");
        }
        for (const auto& value : values) {
            if (!references.count(value.first)) schemaError(source, path + "." + value.first, "unknown reference");
        }
    };
    validateSingleMap(definition.amplitudes, "expressions.amplitudes");
    validateSingleMap(definition.referenceToScatterer, "sizes.reference_to_scatterer");

    std::set<ReferencePair> requiredPairs;
    for (auto first = references.begin(); first != references.end(); ++first) {
        for (auto second = first; second != references.end(); ++second) {
            if (*first == *second && specific.count(*first)) continue;
            requiredPairs.emplace(*first, *second);
        }
    }
    const auto validatePairMap = [&](const auto& values, const std::string& path) {
        for (const auto& pair : requiredPairs) {
            if (!values.count(pair)) {
                schemaError(source, path, "required pair ['" + pair.first + "', '" + pair.second + "'] is missing");
            }
        }
        for (const auto& value : values) {
            if (!references.count(value.first.first) || !references.count(value.first.second)) {
                schemaError(source, path, "pair contains an unknown reference");
            }
            if (!requiredPairs.count(value.first)) {
                schemaError(source, path, "pair is redundant for an identical specific reference");
            }
        }
    };
    validatePairMap(definition.phases, "expressions.phases");
    validatePairMap(definition.referenceToReference, "sizes.reference_to_reference");
}

SubunitValidationCase parseValidationCase(
    const Node& node,
    const std::string& source,
    const std::string& path) {
    rejectUnknownKeys(
        node,
        {"name", "quantity", "parameters", "references", "q", "beta", "expected",
         "absolute_tolerance", "relative_tolerance"},
        source,
        path);
    SubunitValidationCase result;
    result.name = node.contains("name") ? scalarString(node.at("name"), source, path + ".name") : path;
    result.quantity = scalarString(requiredNode(node, "quantity", source, path), source, path + ".quantity");
    result.expected = scalarDouble(requiredNode(node, "expected", source, path), source, path + ".expected");
    if (node.contains("q")) {
        result.q = scalarDouble(node.at("q"), source, path + ".q");
        result.hasQ = true;
    }
    if (node.contains("beta")) result.beta = scalarDouble(node.at("beta"), source, path + ".beta");
    if (node.contains("references")) result.references = stringSequence(node.at("references"), source, path + ".references");
    if (node.contains("absolute_tolerance")) {
        result.absoluteTolerance = scalarDouble(node.at("absolute_tolerance"), source, path + ".absolute_tolerance");
    }
    if (node.contains("relative_tolerance")) {
        result.relativeTolerance = scalarDouble(node.at("relative_tolerance"), source, path + ".relative_tolerance");
    }
    if (node.contains("parameters")) {
        const Node& parameters = node.at("parameters");
        requireMapping(parameters, source, path + ".parameters");
        for (const auto& item : parameters.map_items()) {
            const std::string name = scalarString(item.key(), source, path + ".parameters.<key>");
            result.parameters.emplace(name, scalarDouble(item.value(), source, path + ".parameters." + name));
        }
    }
    static const std::set<std::string> quantities = {
        "form_factor", "form_factor_unnormalized", "amplitude", "amplitude_unnormalized",
        "phase", "radius_of_gyration_squared", "reference_to_scatterer", "reference_to_reference"
    };
    if (!quantities.count(result.quantity)) schemaError(source, path + ".quantity", "unsupported validation quantity");
    const bool qRequired = result.quantity == "form_factor" || result.quantity == "form_factor_unnormalized" ||
        result.quantity == "amplitude" || result.quantity == "amplitude_unnormalized" || result.quantity == "phase";
    if (qRequired != result.hasQ) {
        schemaError(source, path + ".q", qRequired ? "q is required for this quantity" : "q is not valid for this quantity");
    }
    const std::size_t requiredReferences =
        result.quantity == "amplitude" || result.quantity == "amplitude_unnormalized" ||
        result.quantity == "reference_to_scatterer" ? 1u :
        result.quantity == "phase" || result.quantity == "reference_to_reference" ? 2u : 0u;
    if (result.references.size() != requiredReferences) {
        schemaError(source, path + ".references", "quantity expects " + std::to_string(requiredReferences) + " reference(s)");
    }
    if ((result.absoluteTolerance != -1.0 && result.absoluteTolerance <= 0.0) ||
        (result.relativeTolerance != -1.0 && result.relativeTolerance <= 0.0)) {
        schemaError(source, path, "validation tolerances must be positive");
    }
    return result;
}

} // namespace

ReferencePair CanonicalReferencePair(std::string first, std::string second) {
    if (second < first) std::swap(first, second);
    return {std::move(first), std::move(second)};
}

SubunitDefinition LoadSubunitDefinitionFile(const std::string& path) {
    static const std::string extension = ".pyseb.yaml";
    if (path.size() < extension.size() ||
        path.compare(path.size() - extension.size(), extension.size(), extension) != 0) {
        throw SEBException(path + ": model files must use the .pyseb.yaml extension", "LoadSubunitDefinitionFile()");
    }
    std::ifstream input(path, std::ios::binary);
    if (!input) throw SEBException("Unable to open subunit model file " + path, "LoadSubunitDefinitionFile()");
    input.seekg(0, std::ios::end);
    const std::streamoff size = input.tellg();
    if (size < 0 || static_cast<std::size_t>(size) > kMaximumModelFileSize) {
        throw SEBException(path + ": model file exceeds the 1 MiB limit", "LoadSubunitDefinitionFile()");
    }
    input.seekg(0, std::ios::beg);
    const std::string yaml((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    return LoadSubunitDefinitionYaml(yaml, path);
}

SubunitDefinition LoadSubunitDefinitionYaml(const std::string& yaml, const std::string& source) {
    if (yaml.size() > kMaximumModelFileSize) {
        throw SEBException(source + ": model file exceeds the 1 MiB limit", "LoadSubunitDefinitionYaml()");
    }
    rejectMultipleDocuments(yaml, source);

    Node root;
    try {
        root = Node::deserialize(yaml);
    } catch (const std::exception& error) {
        throw SEBException(source + ": " + error.what(), "LoadSubunitDefinitionYaml()");
    }
    validateNodeRestrictions(root, source, "$");
    rejectUnknownKeys(
        root,
        {"format", "schema_version", "id", "api_name", "model_version", "metadata", "behavior", "parameters",
         "variables", "definitions", "integration", "integrals", "references", "expressions", "sizes", "validation"},
        source,
        "$");

    const std::string format = scalarString(requiredNode(root, "format", source, "$"), source, "format");
    if (format != "pyseb-subunit") schemaError(source, "format", "expected 'pyseb-subunit'");
    const Node& schemaVersion = requiredNode(root, "schema_version", source, "$");
    if (!schemaVersion.is_integer() || schemaVersion.get_value<int>() != 1) {
        schemaError(source, "schema_version", "only schema version 1 is supported");
    }

    SubunitDefinition definition;
    definition.source = source;
    definition.id = scalarString(requiredNode(root, "id", source, "$"), source, "id");
    definition.apiName = scalarString(requiredNode(root, "api_name", source, "$"), source, "api_name");
    definition.modelVersion = scalarString(requiredNode(root, "model_version", source, "$"), source, "model_version");
    if (!validModelId(definition.id)) schemaError(source, "id", "expected a namespaced ID such as 'my-lab/MyModel'");
    if (!validSymbolName(definition.apiName)) {
        schemaError(source, "api_name", "must start with a letter and contain only letters and digits");
    }
    if (definition.modelVersion.empty()) schemaError(source, "model_version", "must not be empty");

    if (root.contains("metadata")) {
        const Node& metadata = root.at("metadata");
        rejectUnknownKeys(metadata, {"title", "description", "authors", "citations", "license"}, source, "metadata");
        if (metadata.contains("title")) definition.metadata.title = scalarString(metadata.at("title"), source, "metadata.title");
        if (metadata.contains("description")) definition.metadata.description = scalarString(metadata.at("description"), source, "metadata.description");
        if (metadata.contains("authors")) definition.metadata.authors = stringSequence(metadata.at("authors"), source, "metadata.authors");
        if (metadata.contains("citations")) definition.metadata.citations = stringSequence(metadata.at("citations"), source, "metadata.citations");
        if (metadata.contains("license")) definition.metadata.license = scalarString(metadata.at("license"), source, "metadata.license");
    }

    if (root.contains("behavior")) {
        const Node& behavior = root.at("behavior");
        rejectUnknownKeys(behavior, {"scattering"}, source, "behavior");
        if (behavior.contains("scattering")) {
            const std::string scattering = scalarString(behavior.at("scattering"), source, "behavior.scattering");
            if (scattering == "invisible") definition.invisible = true;
            else if (scattering != "normal") schemaError(source, "behavior.scattering", "expected 'normal' or 'invisible'");
        }
    }

    const Node& parameters = requiredNode(root, "parameters", source, "$");
    requireMapping(parameters, source, "parameters");
    for (const auto& item : parameters.map_items()) {
        ParameterDefinition parameter;
        parameter.name = scalarString(item.key(), source, "parameters.<key>");
        const std::string path = "parameters." + parameter.name;
        rejectUnknownKeys(item.value(), {"unit", "description"}, source, path);
        if (item.value().contains("unit")) parameter.unit = scalarString(item.value().at("unit"), source, path + ".unit");
        if (item.value().contains("description")) parameter.description = scalarString(item.value().at("description"), source, path + ".description");
        definition.parameters.emplace(parameter.name, parameter);
    }

    if (root.contains("variables")) definition.variables = expressionMap(root.at("variables"), source, "variables");
    if (root.contains("definitions")) definition.definitions = expressionMap(root.at("definitions"), source, "definitions");
    if (root.contains("integration")) {
        definition.integration = integrationOptions(
            root.at("integration"), definition.integration, source, "integration");
    }
    if (root.contains("integrals")) {
        definition.integrals = integralMap(
            root.at("integrals"), definition.integration, source, "integrals");
    }

    const Node& references = requiredNode(root, "references", source, "$");
    rejectUnknownKeys(references, {"specific", "distributed"}, source, "references");
    definition.specificReferences = stringSequence(requiredNode(references, "specific", source, "references"), source, "references.specific");
    definition.distributedReferences = stringSequence(requiredNode(references, "distributed", source, "references"), source, "references.distributed");

    const Node& expressions = requiredNode(root, "expressions", source, "$");
    rejectUnknownKeys(expressions, {"form_factor", "amplitudes", "phases"}, source, "expressions");
    definition.formFactor = parsedExpression(requiredNode(expressions, "form_factor", source, "expressions"), source, "expressions.form_factor");
    definition.amplitudes = expressionMap(requiredNode(expressions, "amplitudes", source, "expressions"), source, "expressions.amplitudes");
    definition.phases = pairExpressionSequence(requiredNode(expressions, "phases", source, "expressions"), source, "expressions.phases");

    const Node& sizes = requiredNode(root, "sizes", source, "$");
    rejectUnknownKeys(sizes, {"radius_of_gyration_squared", "reference_to_scatterer", "reference_to_reference"}, source, "sizes");
    definition.radiusOfGyrationSquared = parsedExpression(
        requiredNode(sizes, "radius_of_gyration_squared", source, "sizes"), source, "sizes.radius_of_gyration_squared");
    definition.referenceToScatterer = expressionMap(
        requiredNode(sizes, "reference_to_scatterer", source, "sizes"), source, "sizes.reference_to_scatterer");
    definition.referenceToReference = pairExpressionSequence(
        requiredNode(sizes, "reference_to_reference", source, "sizes"), source, "sizes.reference_to_reference");

    if (root.contains("validation")) {
        const Node& validation = root.at("validation");
        rejectUnknownKeys(validation, {"absolute_tolerance", "relative_tolerance", "cases"}, source, "validation");
        if (validation.contains("absolute_tolerance")) {
            definition.validationAbsoluteTolerance = scalarDouble(validation.at("absolute_tolerance"), source, "validation.absolute_tolerance");
        }
        if (validation.contains("relative_tolerance")) {
            definition.validationRelativeTolerance = scalarDouble(validation.at("relative_tolerance"), source, "validation.relative_tolerance");
        }
        if (definition.validationAbsoluteTolerance <= 0.0 || definition.validationRelativeTolerance <= 0.0) {
            schemaError(source, "validation", "tolerances must be positive");
        }
        if (validation.contains("cases")) {
            const Node& cases = validation.at("cases");
            requireSequence(cases, source, "validation.cases");
            std::size_t index = 0;
            for (const auto& item : cases) {
                definition.validationCases.push_back(parseValidationCase(
                    item, source, "validation.cases[" + std::to_string(index++) + "]"));
            }
        }
    }

    validateUniqueNames(definition);
    validateExpressionDependencies(definition);
    validateReferenceCompleteness(definition);
    std::set<std::string> validationReferences(
        definition.specificReferences.begin(), definition.specificReferences.end());
    validationReferences.insert(
        definition.distributedReferences.begin(), definition.distributedReferences.end());
    for (const auto& testCase : definition.validationCases) {
        for (const auto& parameter : testCase.parameters) {
            if (!definition.parameters.count(parameter.first)) {
                schemaError(source, "validation.cases." + testCase.name, "unknown parameter '" + parameter.first + "'");
            }
        }
        for (const auto& reference : testCase.references) {
            if (!validationReferences.count(reference)) {
                schemaError(source, "validation.cases." + testCase.name, "unknown reference '" + reference + "'");
            }
        }
    }
    if (definition.invisible) {
        const auto isLiteralZero = [](const ParsedExpression& expression) {
            return expression.root() &&
                expression.root()->kind == ParsedExpressionNode::Kind::Number &&
                expression.root()->number == 0.0;
        };
        if (!isLiteralZero(definition.formFactor)) {
            schemaError(source, "expressions.form_factor", "invisible models must use a literal zero expression");
        }
        for (const auto& amplitude : definition.amplitudes) {
            if (!isLiteralZero(amplitude.second)) {
                schemaError(source, "expressions.amplitudes." + amplitude.first, "invisible models must use a literal zero expression");
            }
        }
    }
    return definition;
}

} // namespace pyseb
