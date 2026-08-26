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
api_name: Minimal
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
    EXPECT_EQ(point.apiName, "Point");
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

    std::string missingApiName = kMinimalModel;
    missingApiName.erase(
        missingApiName.find("api_name: Minimal\n"),
        std::string("api_name: Minimal\n").size());
    EXPECT_THROW(pyseb::LoadSubunitDefinitionYaml(missingApiName, "missing-api.yaml"), SEBException);

    std::string invalidApiName = kMinimalModel;
    invalidApiName.replace(
        invalidApiName.find("api_name: Minimal"),
        std::string("api_name: Minimal").size(),
        "api_name: invalid_name");
    EXPECT_THROW(pyseb::LoadSubunitDefinitionYaml(invalidApiName, "invalid-api.yaml"), SEBException);

    std::string duplicateApiName = kMinimalModel;
    duplicateApiName.replace(
        duplicateApiName.find("api_name: Minimal"),
        std::string("api_name: Minimal").size(),
        "api_name: Minimal\napi_name: Duplicate");
    EXPECT_THROW(pyseb::LoadSubunitDefinitionYaml(duplicateApiName, "duplicate-api.yaml"), SEBException);
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

TEST(FileDefinedWorld, BundledCatalogueContainsEveryConfiguredModel) {
    const auto models = pyseb::BundledSubunitModels();
    ASSERT_EQ(models.size(), static_cast<std::size_t>(PYSEB_BUNDLED_MODEL_COUNT));
    World world("catalogue");
    std::size_t index = 0;
    for (const auto& info : models) {
        EXPECT_FALSE(info.id.empty());
        EXPECT_FALSE(info.apiName.empty());
        EXPECT_TRUE(info.bundled);
        EXPECT_TRUE(pyseb::ValidateSubunitDefinition(
            pyseb::BundledSubunitDefinition(info.id)).ok());
        const std::string suffix = std::to_string(index++);
        world.Add(info.id, "canonical" + suffix);
        world.Add(info.apiName, "alias" + suffix);
        EXPECT_EQ(world.getSubunit("canonical" + suffix)->getSubunitType(), FILEDEFINEDSUBUNIT);
        EXPECT_EQ(world.getSubunit("alias" + suffix)->getSubunitType(), FILEDEFINEDSUBUNIT);
    }
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

// SubUnit::NumericSigmaMSDRef2Ref must only short-circuit self-pairs of a
// *specific* reference to zero; a distributed reference paired with itself
// (e.g. GaussianLoop's "contour") has a genuine, nonzero sigma<R^2> that the
// bundled model supplies and must be evaluated, not skipped.
TEST(FileDefinedWorld, DistributedSelfReferenceSigmaIsEvaluatedNotZeroed) {
    World world("loop");
    world.Add("GaussianLoop", "loop");
    const ParameterList values{{"Rg_loop", 2.0}};
    SubUnit* loop = world.getSubunit("loop");
    EXPECT_NEAR(loop->NumericSigmaMSDRef2Ref("contour", "contour", values), 8.0, 1e-12);

    world.Add("ThinRod", "rod");
    const ParameterList rodValues{{"L_rod", 10.0}};
    SubUnit* rod = world.getSubunit("rod");
    EXPECT_DOUBLE_EQ(rod->NumericSigmaMSDRef2Ref("end1", "end1", rodValues), 0.0);
}

static_assert(POINT == 2, "POINT must retain its legacy enum slot");
static_assert(GAUSSIANPOLYMER == 3, "GAUSSIANPOLYMER must retain its legacy enum slot");

} // namespace
