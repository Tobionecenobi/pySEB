#include "Subunits/FileDefined.hpp"

#include <functional>
#include <map>
#include <set>

#include "Exceptions.hpp"
#include "SubunitIO/ExpressionParser.hpp"

FileDefinedSubunit::FileDefinedSubunit(const pyseb::SubunitDefinition& definition)
    : SubUnit(), definition_(definition) {
    type = SUBUNITCHILD;
    stype = FILEDEFINEDSUBUNIT;
    setName(definition.id);
}

void FileDefinedSubunit::Init(string name, string tag, SymbolInterface* GEX) {
    SubUnit::Init(name, tag, GEX);
    const string modelTag = getTag();

    for (const auto& reference : definition_.specificReferences) setReferencePointName(reference);
    for (const auto& reference : definition_.distributedReferences) setDistReferencePointType(reference);

    std::map<std::string, Expression> parameterSymbols;
    for (const auto& parameter : definition_.parameters) {
        Expression symbol = GLEX->getSymbol(parameter.first, modelTag);
        parameterSymbols.emplace(parameter.first, symbol);
        parameters.insert(symbol);
    }

    std::map<std::string, Expression> variableSymbols;
    for (const auto& variable : definition_.variables) {
        Expression symbol = GLEX->getSymbol(variable.first, modelTag);
        variableSymbols.emplace(variable.first, symbol);
        xparameters.insert(symbol);
    }

    const auto baseResolve = [&](const std::string& identifier) -> Expression {
        if (identifier == "q") return GLEX->getSymbol("q");
        const auto parameter = parameterSymbols.find(identifier);
        if (parameter != parameterSymbols.end()) return parameter->second;
        throw SEBException(
            definition_.source + ": unresolved identifier '" + identifier + "'",
            "FileDefinedSubunit::Init()");
    };

    std::map<std::string, Expression> expandedVariables;
    std::set<std::string> expandingVariables;
    std::function<Expression(const std::string&)> expandVariable;
    expandVariable = [&](const std::string& name) -> Expression {
        const auto cached = expandedVariables.find(name);
        if (cached != expandedVariables.end()) return cached->second;
        if (!expandingVariables.insert(name).second) {
            throw SEBException(definition_.source + ": cyclic variable dependency at '" + name + "'");
        }
        const Expression value = pyseb::MaterializeSubunitExpression(
            definition_.variables.at(name),
            [&](const std::string& identifier) -> Expression {
                if (definition_.variables.count(identifier)) return expandVariable(identifier);
                return baseResolve(identifier);
            });
        expandingVariables.erase(name);
        expandedVariables.emplace(name, value);
        return value;
    };
    for (const auto& variable : definition_.variables) {
        expand[variableSymbols.at(variable.first)] = expandVariable(variable.first);
    }

    std::map<std::string, Expression> namedExpressions;
    std::set<std::string> materializingDefinitions;
    std::function<Expression(const std::string&)> materializeDefinition;
    std::function<Expression(const std::string&)> symbolicResolve;
    symbolicResolve = [&](const std::string& identifier) -> Expression {
        const auto variable = variableSymbols.find(identifier);
        if (variable != variableSymbols.end()) return variable->second;
        if (definition_.definitions.count(identifier)) return materializeDefinition(identifier);
        return baseResolve(identifier);
    };
    materializeDefinition = [&](const std::string& name) -> Expression {
        const auto cached = namedExpressions.find(name);
        if (cached != namedExpressions.end()) return cached->second;
        if (!materializingDefinitions.insert(name).second) {
            throw SEBException(definition_.source + ": cyclic reusable definition at '" + name + "'");
        }
        const Expression value = pyseb::MaterializeSubunitExpression(
            definition_.definitions.at(name), symbolicResolve);
        materializingDefinitions.erase(name);
        namedExpressions.emplace(name, value);
        return value;
    };

    std::map<std::string, Expression> expandedNamedExpressions;
    std::set<std::string> expandingDefinitions;
    std::function<Expression(const std::string&)> materializeExpandedDefinition;
    std::function<Expression(const std::string&)> expandedResolve;
    expandedResolve = [&](const std::string& identifier) -> Expression {
        if (definition_.variables.count(identifier)) return expandVariable(identifier);
        if (definition_.definitions.count(identifier)) return materializeExpandedDefinition(identifier);
        return baseResolve(identifier);
    };
    materializeExpandedDefinition = [&](const std::string& name) -> Expression {
        const auto cached = expandedNamedExpressions.find(name);
        if (cached != expandedNamedExpressions.end()) return cached->second;
        if (!expandingDefinitions.insert(name).second) {
            throw SEBException(definition_.source + ": cyclic reusable definition at '" + name + "'");
        }
        const Expression value = pyseb::MaterializeSubunitExpression(
            definition_.definitions.at(name), expandedResolve);
        expandingDefinitions.erase(name);
        expandedNamedExpressions.emplace(name, value);
        return value;
    };

    FormFactorExpression = pyseb::MaterializeSubunitExpression(definition_.formFactor, symbolicResolve);
    for (const auto& amplitude : definition_.amplitudes) {
        FormFactorAmplitudeExpressions[amplitude.first] =
            pyseb::MaterializeSubunitExpression(amplitude.second, symbolicResolve);
    }
    for (const auto& phase : definition_.phases) {
        PhaseFactorExpressions[phase.first.first][phase.first.second] =
            pyseb::MaterializeSubunitExpression(phase.second, symbolicResolve);
    }

    RadiusOfGyration2 = pyseb::MaterializeSubunitExpression(
        definition_.radiusOfGyrationSquared, expandedResolve);
    for (const auto& size : definition_.referenceToScatterer) {
        sigmaMSDref2scat[size.first] = pyseb::MaterializeSubunitExpression(size.second, expandedResolve);
    }
    for (const auto& size : definition_.referenceToReference) {
        sigmaMSDref2ref[size.first.first][size.first.second] =
            pyseb::MaterializeSubunitExpression(size.second, expandedResolve);
    }
}

Expression FileDefinedSubunit::FormFactor(
    ParameterList& betas,
    ParameterList& params,
    int varForm) {
    if (definition_.invisible) return constant(0);
    return SubUnit::FormFactor(betas, params, varForm);
}

Expression FileDefinedSubunit::FormFactorAmplitude(
    refPoint reference,
    ParameterList& betas,
    ParameterList& params,
    int varForm) {
    if (definition_.invisible) {
        if (!testReference(reference)) throw SEBException("Invalid reference point " + reference);
        return constant(0);
    }
    return SubUnit::FormFactorAmplitude(reference, betas, params, varForm);
}

double FileDefinedSubunit::NumericTotalBeta(const ParameterList& values) {
    if (definition_.invisible) return 0.0;
    return SubUnit::NumericTotalBeta(values);
}
