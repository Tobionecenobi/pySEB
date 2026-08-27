#include "Subunits/FileDefined.hpp"

#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <sstream>

#include "Exceptions.hpp"
#include "ExpressionFunctions.hpp"
#include "SubunitIO/ExpressionParser.hpp"

FileDefinedSubunit::FileDefinedSubunit(const pyseb::SubunitDefinition& definition)
    : IntegratedSubunit(), definition_(definition) {
    type = SUBUNITCHILD;
    stype = FILEDEFINEDSUBUNIT;
    setName(definition.id);
    for (const auto& integral : definition_.integrals) {
        integrators_.emplace(
            integral.first,
            std::unique_ptr<NumericalIntegrator>(
                new NumericalIntegrator(integral.second.integration)));
    }
}

void FileDefinedSubunit::Init(string name, string tag, SymbolInterface* GEX) {
    IntegratedSubunit::Init(name, tag, GEX);
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
    std::map<std::string, Expression> namedIntegrals;
    std::set<std::string> materializingDefinitions;
    std::set<std::string> materializingIntegrals;
    std::function<Expression(const std::string&)> materializeDefinition;
    std::function<Expression(const std::string&)> materializeIntegral;
    std::function<Expression(const std::string&)> symbolicResolve;
    symbolicResolve = [&](const std::string& identifier) -> Expression {
        const auto variable = variableSymbols.find(identifier);
        if (variable != variableSymbols.end()) return variable->second;
        if (definition_.definitions.count(identifier)) return materializeDefinition(identifier);
        if (definition_.integrals.count(identifier)) return materializeIntegral(identifier);
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
    materializeIntegral = [&](const std::string& name) -> Expression {
        const auto cached = namedIntegrals.find(name);
        if (cached != namedIntegrals.end()) return cached->second;
        if (!materializingIntegrals.insert(name).second) {
            throw SEBException(definition_.source + ": cyclic integral dependency at '" + name + "'");
        }
        const pyseb::IntegralDefinition& integral = definition_.integrals.at(name);
        const Expression variable = GLEX->getSymbol(integral.variable, modelTag);
        const auto integralResolve = [&](const std::string& identifier) -> Expression {
            if (identifier == integral.variable) return variable;
            return symbolicResolve(identifier);
        };
        const Expression lower = pyseb::MaterializeSubunitExpression(
            integral.lower, symbolicResolve);
        const Expression upper = pyseb::MaterializeSubunitExpression(
            integral.upper, symbolicResolve);
        const Expression integrand = pyseb::MaterializeSubunitExpression(
            integral.integrand, integralResolve);
        const Expression value = integrate(variable, lower, upper, integrand);
        materializingIntegrals.erase(name);
        namedIntegrals.emplace(name, value);
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

    if (requiresParsedNumericalEvaluation(definition_.formFactor)) {
        setNumericalFormFactorFunction(
            [this](double q, const ParameterList& values) {
                return evaluateNumerically(definition_.formFactor, q, values);
            });
    }
    for (const auto& amplitude : definition_.amplitudes) {
        if (!requiresParsedNumericalEvaluation(amplitude.second)) continue;
        const std::string reference = amplitude.first;
        setNumericalFormFactorAmplitudeFunction(
            reference,
            [this, reference](double q, const ParameterList& values) {
                return evaluateNumerically(definition_.amplitudes.at(reference), q, values);
            });
    }
    for (const auto& phase : definition_.phases) {
        if (!requiresParsedNumericalEvaluation(phase.second)) continue;
        const pyseb::ReferencePair references = phase.first;
        setNumericalPhaseFactorFunction(
            references.first,
            references.second,
            [this, references](double q, const ParameterList& values) {
                return evaluateNumerically(definition_.phases.at(references), q, values);
            });
    }
}

bool FileDefinedSubunit::requiresParsedNumericalEvaluation(
    const pyseb::ParsedExpression& expression) const {
    std::set<std::string> visited;
    std::function<bool(const pyseb::ParsedExpression&)> inspect;
    inspect = [&](const pyseb::ParsedExpression& current) {
        std::function<bool(const std::shared_ptr<const pyseb::ParsedExpressionNode>&)> usesAstFunction;
        usesAstFunction = [&](const std::shared_ptr<const pyseb::ParsedExpressionNode>& node) {
            if (!node) return false;
            if (node->kind == pyseb::ParsedExpressionNode::Kind::Function &&
                (node->text == "struve_h0" || node->text == "struve_h1")) {
                return true;
            }
            for (const auto& child : node->children) {
                if (usesAstFunction(child)) return true;
            }
            return false;
        };
        if (usesAstFunction(current.root())) return true;
        for (const auto& identifier : current.identifiers()) {
            if (definition_.integrals.count(identifier)) return true;
            const auto named = definition_.definitions.find(identifier);
            if (named != definition_.definitions.end() && visited.insert(identifier).second &&
                inspect(named->second)) {
                return true;
            }
        }
        return false;
    };
    return inspect(expression);
}

double FileDefinedSubunit::evaluateNumerically(
    const pyseb::ParsedExpression& expression,
    double q,
    const ParameterList& values) {
    std::map<std::string, double> variableValues;
    std::map<std::string, double> definitionValues;
    std::map<std::string, double> integralValues;
    std::set<std::string> evaluating;

    std::function<double(const pyseb::ParsedExpression&, const std::string*, double)> evaluate;
    std::function<double(const std::string&)> resolve;
    std::function<double(const std::string&)> evaluateVariable;
    std::function<double(const std::string&)> evaluateDefinition;
    std::function<double(const std::string&)> evaluateIntegral;

    evaluate = [&](const pyseb::ParsedExpression& current,
                   const std::string* localName,
                   double localValue) {
        return pyseb::EvaluateSubunitExpression(
            current,
            [&](const std::string& identifier) {
                if (localName && identifier == *localName) return localValue;
                return resolve(identifier);
            });
    };
    evaluateVariable = [&](const std::string& name) {
        const auto cached = variableValues.find(name);
        if (cached != variableValues.end()) return cached->second;
        if (!evaluating.insert("variable:" + name).second) {
            throw SEBException(definition_.source + ": cyclic variable dependency at '" + name + "'");
        }
        const double value = evaluate(definition_.variables.at(name), nullptr, 0.0);
        evaluating.erase("variable:" + name);
        variableValues.emplace(name, value);
        return value;
    };
    evaluateDefinition = [&](const std::string& name) {
        const auto cached = definitionValues.find(name);
        if (cached != definitionValues.end()) return cached->second;
        if (!evaluating.insert("definition:" + name).second) {
            throw SEBException(definition_.source + ": cyclic definition dependency at '" + name + "'");
        }
        const double value = evaluate(definition_.definitions.at(name), nullptr, 0.0);
        evaluating.erase("definition:" + name);
        definitionValues.emplace(name, value);
        return value;
    };
    evaluateIntegral = [&](const std::string& name) {
        const auto cached = integralValues.find(name);
        if (cached != integralValues.end()) return cached->second;
        const pyseb::IntegralDefinition& integral = definition_.integrals.at(name);
        const double lower = evaluate(integral.lower, nullptr, 0.0);
        const double upper = evaluate(integral.upper, nullptr, 0.0);
        if (!std::isfinite(lower) || !std::isfinite(upper)) {
            throw SEBException(
                definition_.source + ": integrals." + name +
                    ": evaluated bounds must be finite",
                "FileDefinedSubunit numerical integration");
        }
        try {
            const double value = integrators_.at(name)->integrate(
                [&](double integrationVariable) {
                    return evaluate(
                        integral.integrand,
                        &integral.variable,
                        integrationVariable);
                },
                lower,
                upper).value;
            integralValues.emplace(name, value);
            return value;
        } catch (const std::exception& error) {
            std::ostringstream message;
            message << definition_.source << ": integrals." << name
                    << " on [" << lower << ", " << upper << "]: " << error.what();
            throw SEBException(message.str(), "FileDefinedSubunit numerical integration");
        }
    };
    resolve = [&](const std::string& identifier) {
        if (identifier == "q") return q;
        if (definition_.variables.count(identifier)) return evaluateVariable(identifier);
        if (definition_.definitions.count(identifier)) return evaluateDefinition(identifier);
        if (definition_.integrals.count(identifier)) return evaluateIntegral(identifier);
        if (definition_.parameters.count(identifier)) {
            const std::string tagged = identifier + "_" + getTag();
            const auto value = values.find(tagged);
            if (value == values.end()) {
                throw SEBException(
                    definition_.source + ": missing numerical parameter " + tagged,
                    "FileDefinedSubunit numerical evaluation");
            }
            if (!std::isfinite(value->second)) {
                throw SEBException(
                    definition_.source + ": numerical parameter " + tagged + " must be finite",
                    "FileDefinedSubunit numerical evaluation");
            }
            return value->second;
        }
        throw SEBException(
            definition_.source + ": unresolved numerical identifier '" + identifier + "'",
            "FileDefinedSubunit numerical evaluation");
    };

    return evaluate(expression, nullptr, 0.0);
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
