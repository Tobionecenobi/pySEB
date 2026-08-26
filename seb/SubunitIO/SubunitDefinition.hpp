#ifndef PYSEB_SUBUNIT_IO_SUBUNIT_DEFINITION_HPP
#define PYSEB_SUBUNIT_IO_SUBUNIT_DEFINITION_HPP

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "SubunitIO/ExpressionParser.hpp"

namespace pyseb {

using ReferencePair = std::pair<std::string, std::string>;

struct ParameterDefinition {
    std::string name;
    std::string unit;
    std::string description;
};

struct SubunitMetadata {
    std::string title;
    std::string description;
    std::vector<std::string> authors;
    std::vector<std::string> citations;
    std::string license;
};

struct SubunitValidationCase {
    std::string name;
    std::string quantity;
    std::map<std::string, double> parameters;
    std::vector<std::string> references;
    double q = 0.0;
    bool hasQ = false;
    double beta = 1.0;
    double expected = 0.0;
    double absoluteTolerance = -1.0;
    double relativeTolerance = -1.0;
};

struct SubunitDefinition {
    int schemaVersion = 1;
    std::string id;
    std::string apiName;
    std::string modelVersion;
    std::string source;
    SubunitMetadata metadata;
    bool invisible = false;

    std::map<std::string, ParameterDefinition> parameters;
    std::map<std::string, ParsedExpression> variables;
    std::map<std::string, ParsedExpression> definitions;
    std::vector<std::string> specificReferences;
    std::vector<std::string> distributedReferences;

    ParsedExpression formFactor;
    std::map<std::string, ParsedExpression> amplitudes;
    std::map<ReferencePair, ParsedExpression> phases;
    ParsedExpression radiusOfGyrationSquared;
    std::map<std::string, ParsedExpression> referenceToScatterer;
    std::map<ReferencePair, ParsedExpression> referenceToReference;

    double validationAbsoluteTolerance = 1e-10;
    double validationRelativeTolerance = 1e-8;
    std::vector<SubunitValidationCase> validationCases;
};

struct SubunitModelInfo {
    std::string id;
    std::string apiName;
    std::string modelVersion;
    std::string title;
    std::string source;
    bool bundled = false;
};

struct SubunitValidationFailure {
    std::string caseName;
    std::string message;
};

struct SubunitValidationReport {
    std::string modelId;
    std::size_t caseCount = 0;
    std::vector<SubunitValidationFailure> failures;
    std::vector<std::string> warnings;

    bool ok() const { return failures.empty(); }
};

SubunitDefinition LoadSubunitDefinitionFile(const std::string& path);
SubunitDefinition LoadSubunitDefinitionYaml(
    const std::string& yaml,
    const std::string& source = "<string>");

SubunitValidationReport ValidateSubunitDefinition(const SubunitDefinition& definition);
SubunitValidationReport ValidateSubunitFile(const std::string& path);

ReferencePair CanonicalReferencePair(std::string first, std::string second);

} // namespace pyseb

#endif
