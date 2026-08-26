#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <string>

#include "SEB.hpp"
#include "SubunitIO/BundledModels.hpp"
#include "SubunitIO/ExpressionParser.hpp"
#include "SubunitIO/SubunitDefinition.hpp"

namespace {

const char kMinimalModel[] = R"YAML(
format: pyseb-subunit
schema_version: 1
id: test/Minimal
model_version: "1"
parameters: {}
references:
  specific: [center]
  distributed: []
expressions:
  form_factor: "1"
  amplitudes: {center: "1"}
  phases: []
sizes:
  radius_of_gyration_squared: "0"
  reference_to_scatterer: {center: "0"}
  reference_to_reference: []
)YAML";

TEST(FileDefinedExpression, UsesMathematicalPowerPrecedence) {
    SymbolicFactory::instance();
    const auto parsed = pyseb::ParseSubunitExpression("-2^2 + 2^-2");
    const auto expression = pyseb::MaterializeSubunitExpression(
        parsed,
        [](const std::string& name) -> Expression {
            throw SEBException("unexpected identifier " + name);
        });
    EXPECT_NEAR(expression.eval(), -3.75, 1e-12);

    const auto powers = pyseb::ParseSubunitExpression("2^3^2 + 1e-3 + 2E-3");
    const auto powerExpression = pyseb::MaterializeSubunitExpression(
        powers,
        [](const std::string& name) -> Expression {
            throw SEBException("unexpected identifier " + name);
        });
    EXPECT_NEAR(powerExpression.eval(), 512.003, 1e-12);
}

TEST(FileDefinedExpression, RejectsUnsupportedFunctionsAndTrailingInput) {
    EXPECT_THROW(pyseb::ParseSubunitExpression("system(1)"), SEBException);
    EXPECT_THROW(pyseb::ParseSubunitExpression("2 + 3 garbage"), SEBException);
    EXPECT_THROW(pyseb::ParseSubunitExpression("pow(2)"), SEBException);
}

TEST(FileDefinedSchema, LoadsAndValidatesBundledModels) {
    const auto point = pyseb::LoadSubunitDefinitionFile(
        std::string(PYSEB_MODEL_DIR) + "/Point.pyseb.yaml");
    const auto polymer = pyseb::LoadSubunitDefinitionFile(
        std::string(PYSEB_MODEL_DIR) + "/GaussianPolymer.pyseb.yaml");
    EXPECT_EQ(point.id, "pyseb/Point");
    EXPECT_TRUE(point.invisible);
    EXPECT_EQ(polymer.amplitudes.size(), 4u);
    EXPECT_TRUE(pyseb::ValidateSubunitDefinition(point).ok());
    EXPECT_TRUE(pyseb::ValidateSubunitDefinition(polymer).ok());
}

TEST(FileDefinedSchema, RejectsUnknownFieldsAndDependencyCycles) {
    std::string unknown = kMinimalModel;
    unknown += "unknown_field: true\n";
    EXPECT_THROW(pyseb::LoadSubunitDefinitionYaml(unknown, "unknown.yaml"), SEBException);

    std::string cyclic = kMinimalModel;
    const std::string marker = "parameters: {}";
    cyclic.replace(cyclic.find(marker), marker.size(),
        "parameters: {}\nvariables: {x: \"y\", y: \"x\"}");
    EXPECT_THROW(pyseb::LoadSubunitDefinitionYaml(cyclic, "cycle.yaml"), SEBException);
}

TEST(FileDefinedSchema, RejectsUnsupportedYamlFeatures) {
    std::string tagged = kMinimalModel;
    tagged.replace(tagged.find("id: test/Minimal"), 16,
                   "id: !model test/Minimal");
    EXPECT_THROW(pyseb::LoadSubunitDefinitionYaml(tagged, "tagged.yaml"), SEBException);

    std::string anchored = kMinimalModel;
    anchored.replace(anchored.find("parameters: {}"), 14,
                     "parameters: &params {}");
    EXPECT_THROW(pyseb::LoadSubunitDefinitionYaml(anchored, "anchored.yaml"), SEBException);

    std::string multiple = kMinimalModel;
    multiple += "---\nformat: pyseb-subunit\n";
    EXPECT_THROW(pyseb::LoadSubunitDefinitionYaml(multiple, "multiple.yaml"), SEBException);
}

TEST(FileDefinedWorld, DirectAndBundledModelsMatchLegacyName) {
    const std::string model = std::string(PYSEB_MODEL_DIR) + "/GaussianPolymer.pyseb.yaml";
    ParameterList values{{"beta_poly", 2.0}, {"Rg_poly", 2.0}};

    World legacy("legacy");
    legacy.Add("GaussianPolymer", "poly");
    World canonical("canonical");
    canonical.Add("pyseb/GaussianPolymer", "poly");
    World direct("direct");
    direct.AddFile(model, "poly");

    for (double q : {0.0, 1e-6, 0.1, 0.5, 1.0}) {
        const double expected = legacy.EvaluateFormFactor("poly", values, q);
        EXPECT_NEAR(canonical.EvaluateFormFactor("poly", values, q), expected, 1e-12);
        EXPECT_NEAR(direct.EvaluateFormFactor("poly", values, q), expected, 1e-12);
    }
    EXPECT_NEAR(direct.EvaluateRadiusOfGyration2("poly", values), 4.0, 1e-12);
}

TEST(FileDefinedWorld, InvisiblePointPreservesCountingAndPhaseBehavior) {
    World world("points");
    world.Add("Point", "legacy");
    world.Add("pyseb/Point", "canonical");
    EXPECT_DOUBLE_EQ(world.Count("legacy.point").eval(), 0.0);
    EXPECT_DOUBLE_EQ(world.Count("canonical.point").eval(), 0.0);
    EXPECT_DOUBLE_EQ(world.PhaseFactor("canonical.point", "canonical.point").eval(), 1.0);
}

TEST(FileDefinedWorld, RegistrationRejectsConflictingIds) {
    World world("registry");
    const std::string model = std::string(PYSEB_MODEL_DIR) + "/GaussianPolymer.pyseb.yaml";
    EXPECT_THROW(world.RegisterSubunitFile(model), SEBException);
    const auto models = world.ListSubunitModels();
    EXPECT_FALSE(models.empty());
}

} // namespace
