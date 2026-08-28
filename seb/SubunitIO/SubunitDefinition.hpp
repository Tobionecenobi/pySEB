#ifndef PYSEB_SUBUNIT_IO_SUBUNIT_DEFINITION_HPP
#define PYSEB_SUBUNIT_IO_SUBUNIT_DEFINITION_HPP

#include <map>
#include <string>
#include <utility>
#include <vector>
#include <array>

#include "NumericalIntegration.hpp"
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

struct IntegralDefinition {
    std::string name;
    struct Dimension {
        std::string variable;
        ParsedExpression lower;
        ParsedExpression upper;
    };
    std::vector<Dimension> dimensions;
    std::string variable;
    ParsedExpression lower;
    ParsedExpression upper;
    ParsedExpression integrand;
    IntegrationOptions integration;
};

/* Version-1, renderer-neutral visualization description.  Expressions are
   intentionally kept as ParsedExpression so exporters can evaluate them
   without embedding a renderer or the source YAML. */
enum class VisualizationGeometryKind { Curve, Surface, RandomWalk };
struct VisualizationGeometry {
    std::string name;
    VisualizationGeometryKind kind = VisualizationGeometryKind::Curve;
    ParsedExpression uLower, uUpper, vLower, vUpper;
    std::array<ParsedExpression,3> coordinates;
    std::size_t samples = 0;
    std::string distribution = "gaussian";
    std::string closure = "open";
    ParsedExpression targetRg;
    double width = 0.01;
    bool periodic = false;
};
struct VisualizationReferenceVariant {
    std::string geometry;
    std::string sampling = "uniform_parameter";
};
struct VisualizationReference {
    std::string name;
    std::string kind = "fixed";       // fixed, curve_fraction, distributed
    std::string geometry;
    ParsedExpression fraction;
    std::array<ParsedExpression,3> position;
    std::string sampling = "uniform_parameter";
    std::vector<std::string> patches;
    std::vector<double> weights;
    std::map<std::string, VisualizationReferenceVariant> variants;
};
struct VisualizationDefinition {
    bool present = false;
    std::map<std::string, ParsedExpression> definitions;
    std::map<std::string, VisualizationGeometry> geometry;
    std::map<std::string, VisualizationReference> references;
    std::array<double,3> color{{0.7,0.7,0.7}};
    double opacity = 1.0;
    bool doubleSided = true;
    double curveWidth = 0.01;
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
    IntegrationOptions integration;
    std::map<std::string, IntegralDefinition> integrals;
    std::vector<std::string> specificReferences;
    std::vector<std::string> distributedReferences;

    ParsedExpression formFactor;
    std::map<std::string, ParsedExpression> amplitudes;
    std::map<ReferencePair, ParsedExpression> phases;
    ParsedExpression radiusOfGyrationSquared;
    std::map<std::string, ParsedExpression> referenceToScatterer;
    std::map<ReferencePair, ParsedExpression> referenceToReference;

    VisualizationDefinition visualization;

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
