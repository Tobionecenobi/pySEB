#include "SubunitIO/SubunitDefinition.hpp"

#include <cmath>
#include <sstream>

#include "Exceptions.hpp"
#include "Subunits/FileDefined.hpp"
#include "SymbolInterface.hpp"

namespace pyseb {
namespace {

double evaluateCase(
    FileDefinedSubunit& subunit,
    const SubunitDefinition& definition,
    const SubunitValidationCase& testCase) {
    ParameterList values;
    for (const auto& parameter : testCase.parameters) {
        values[parameter.first + "_validation"] = parameter.second;
    }
    values["beta_validation"] = testCase.beta;

    if (testCase.quantity == "form_factor" || testCase.quantity == "form_factor_unnormalized") {
        const double raw = subunit.NumericFormFactorUnnormalized(testCase.q, values);
        if (testCase.quantity == "form_factor_unnormalized") return raw;
        const double beta = subunit.NumericTotalBeta(values);
        if (std::abs(beta) <= 1e-14) throw SEBException("cannot normalize a zero-beta model");
        return raw / (beta * beta);
    }
    if (testCase.quantity == "amplitude" || testCase.quantity == "amplitude_unnormalized") {
        const double raw = subunit.NumericFormFactorAmplitudeUnnormalized(
            testCase.references.at(0), testCase.q, values);
        if (testCase.quantity == "amplitude_unnormalized") return raw;
        const double beta = subunit.NumericTotalBeta(values);
        if (std::abs(beta) <= 1e-14) throw SEBException("cannot normalize a zero-beta model");
        return raw / beta;
    }
    if (testCase.quantity == "phase") {
        return subunit.NumericPhaseFactor(
            testCase.references.at(0), testCase.references.at(1), testCase.q, values);
    }
    if (testCase.quantity == "radius_of_gyration_squared") {
        return subunit.NumericRadiusOfGyration2(values);
    }
    if (testCase.quantity == "reference_to_scatterer") {
        return subunit.NumericSigmaMSDRef2Scat(testCase.references.at(0), values);
    }
    if (testCase.quantity == "reference_to_reference") {
        return subunit.NumericSigmaMSDRef2Ref(
            testCase.references.at(0), testCase.references.at(1), values);
    }
    throw SEBException("Unsupported validation quantity " + testCase.quantity);
}

} // namespace

SubunitValidationReport ValidateSubunitDefinition(const SubunitDefinition& definition) {
    SubunitValidationReport report;
    report.modelId = definition.id;
    report.caseCount = definition.validationCases.size();
    if (definition.validationCases.empty()) {
        report.warnings.push_back("schema is valid but the model contains no reference validation cases");
        return report;
    }

    SymbolInterface symbols;
    FileDefinedSubunit subunit(definition);
    subunit.Init("validation", "validation", &symbols);

    for (const auto& testCase : definition.validationCases) {
        try {
            const double actual = evaluateCase(subunit, definition, testCase);
            const double absoluteTolerance = testCase.absoluteTolerance > 0.0
                ? testCase.absoluteTolerance : definition.validationAbsoluteTolerance;
            const double relativeTolerance = testCase.relativeTolerance > 0.0
                ? testCase.relativeTolerance : definition.validationRelativeTolerance;
            const double allowed = absoluteTolerance + relativeTolerance * std::abs(testCase.expected);
            const double difference = std::abs(actual - testCase.expected);
            if (!std::isfinite(actual) || difference > allowed) {
                std::ostringstream message;
                message.precision(17);
                message << "expected " << testCase.expected << ", got " << actual
                        << " (difference " << difference << ", tolerance " << allowed << ")";
                report.failures.push_back({testCase.name, message.str()});
            }
        } catch (const std::exception& error) {
            report.failures.push_back({testCase.name, error.what()});
        }
    }
    return report;
}

SubunitValidationReport ValidateSubunitFile(const std::string& path) {
    return ValidateSubunitDefinition(LoadSubunitDefinitionFile(path));
}

} // namespace pyseb
