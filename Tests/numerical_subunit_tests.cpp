#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "SEB.hpp"

namespace {

double sinc(double value)
{
    if (std::abs(value) < 1e-8) {
        return 1.0 - value * value / 6.0;
    }
    return std::sin(value) / value;
}

}

// Verify that the two callback contracts are equivalent after World applies
// the beta weighting: normalized callbacks return F and A, whereas
// unnormalized callbacks return beta^2 F and beta A directly.
TEST(NumericalSubunitTest, NormalizedAndUnnormalizedCallbacksAgree)
{
    NumericalSubunit* normalized =
        new NumericalSubunit(NormalizationMode::Normalized);
    normalized->addReferencePoint("center");
    normalized->setTotalBeta(2.0);
    normalized->setFormFactorFunction(
        [](double q, const ParameterList&) { return std::exp(-q * q); }
    );
    normalized->setFormFactorAmplitudeFunction(
        "center",
        [](double q, const ParameterList&) { return std::exp(-0.5 * q * q); }
    );
    normalized->setRadiusOfGyration2(3.0);
    normalized->setSigmaMSDRef2Scat("center", 3.0);

    NumericalSubunit* unnormalized =
        new NumericalSubunit(NormalizationMode::Unnormalized);
    unnormalized->addReferencePoint("center");
    unnormalized->setTotalBeta(2.0);
    unnormalized->setFormFactorFunction(
        [](double q, const ParameterList&) { return 4.0 * std::exp(-q * q); }
    );
    unnormalized->setFormFactorAmplitudeFunction(
        "center",
        [](double q, const ParameterList&) { return 2.0 * std::exp(-0.5 * q * q); }
    );
    unnormalized->setRadiusOfGyration2(3.0);
    unnormalized->setSigmaMSDRef2Scat("center", 3.0);

    EXPECT_NO_THROW(normalized->ValidateNumerically());
    EXPECT_NO_THROW(unnormalized->ValidateNumerically());

    World world;
    world.Add(normalized, "normalized");
    world.Add(unnormalized, "unnormalized");

    const ParameterList parameters;
    const double q = 0.3;
    EXPECT_NEAR(
        world.EvaluateFormFactor("normalized", parameters, q),
        world.EvaluateFormFactor("unnormalized", parameters, q),
        1e-12
    );
    EXPECT_NEAR(
        world.EvaluateFormFactorAmplitude("normalized.center", parameters, q),
        world.EvaluateFormFactorAmplitude("unnormalized.center", parameters, q),
        1e-12
    );
}

// Verify that a distributed reference type may define a non-trivial phase
// factor between two independently selected points from the same distribution.
// Unlike a specific reference point compared with itself, this term is required.
TEST(NumericalSubunitTest, DistributedReferenceCanDefineSelfPhaseFactor)
{
    NumericalSubunit unit(NormalizationMode::Normalized);
    unit.addDistributedReferencePointType("surface");
    unit.setTotalBeta(1.0);
    unit.setFormFactorFunction(
        [](double, const ParameterList&) { return 1.0; }
    );
    unit.setFormFactorAmplitudeFunction(
        "surface",
        [](double, const ParameterList&) { return 1.0; }
    );
    unit.setPhaseFactorFunction(
        "surface",
        "surface",
        [](double, const ParameterList&) { return 1.0; }
    );
    unit.setRadiusOfGyration2(0.0);
    unit.setSigmaMSDRef2Scat("surface", 0.0);
    unit.setSigmaMSDRef2Ref("surface", "surface", 0.0);

    EXPECT_NO_THROW(unit.ValidateNumerically());
}

// Verify that a Debye cloud containing one finite sphere reproduces the
// existing analytic SolidSphere form factor, center amplitude, and Rg^2 over
// several q values, including the exact q=0 limit.
TEST(DebyeSphereCloudTest, SingleSphereMatchesRayleighSphere)
{
    std::vector<SphereScatterer> scatterers{
        SphereScatterer(0.0, 0.0, 0.0, 2.0, 3.0)
    };

    World world;
    world.Add(new DebyeSphereCloud(scatterers), "cloud");
    world.Add(new SolidSphere(), "sphere");

    ParameterList parameters{
        {"beta_sphere", 3.0},
        {"R_sphere", 2.0}
    };

    for (const double q : std::vector<double>{0.0, 0.05, 0.2, 0.7}) {
        EXPECT_NEAR(
            world.EvaluateFormFactor("cloud", parameters, q),
            world.EvaluateFormFactor("sphere", parameters, q),
            1e-11
        );
        EXPECT_NEAR(
            world.EvaluateFormFactorAmplitude("cloud.center", parameters, q),
            world.EvaluateFormFactorAmplitude("sphere.center", parameters, q),
            1e-11
        );
    }

    EXPECT_NEAR(
        world.EvaluateRadiusOfGyration2("cloud", parameters),
        12.0 / 5.0,
        1e-12
    );
}

// Verify the three Debye distance calculations for two point scatterers:
// scatterer-to-scatterer for F, reference-to-scatterer for A, and
// reference-to-reference for Psi. Also check their associated size measures.
TEST(DebyeSphereCloudTest, TwoPointScatterersAndReferencesUseDebyeDistances)
{
    DebyeSphereCloud* cloud = new DebyeSphereCloud({
        SphereScatterer(0.0, 0.0, 0.0, 0.0, 1.0),
        SphereScatterer(2.0, 0.0, 0.0, 0.0, 1.0)
    });
    cloud->addReferencePoint("left", 0.0, 0.0, 0.0);
    cloud->addReferencePoint("right", 2.0, 0.0, 0.0);

    World world;
    world.Add(cloud, "cloud");
    const ParameterList parameters;
    const double q = 0.4;

    EXPECT_NEAR(
        world.EvaluateFormFactor("cloud", parameters, q),
        (2.0 + 2.0 * sinc(2.0 * q)) / 4.0,
        1e-12
    );
    EXPECT_NEAR(
        world.EvaluateFormFactorAmplitude("cloud.left", parameters, q),
        (1.0 + sinc(2.0 * q)) / 2.0,
        1e-12
    );
    EXPECT_NEAR(
        world.EvaluatePhaseFactor("cloud.left", "cloud.right", parameters, q),
        sinc(2.0 * q),
        1e-12
    );
    EXPECT_NEAR(
        world.EvaluateSMSDRef2Ref("cloud.left", "cloud.right", parameters),
        4.0,
        1e-12
    );
    EXPECT_NEAR(
        world.EvaluateRadiusOfGyration2("cloud", parameters),
        1.0,
        1e-12
    );
}

// Verify the zero-contrast edge case. The raw Debye intensity remains defined
// when the particle betas cancel, but normalized scattering and Rg^2 must fail
// because they require division by the total beta.
TEST(DebyeSphereCloudTest, ZeroTotalBetaKeepsUnnormalizedScattering)
{
    World world;
    world.Add(
        new DebyeSphereCloud({
            SphereScatterer(0.0, 0.0, 0.0, 0.0, 1.0),
            SphereScatterer(1.0, 0.0, 0.0, 0.0, -1.0)
        }),
        "cloud"
    );
    const ParameterList parameters;

    EXPECT_TRUE(std::isfinite(
        world.EvaluateFormFactorUnnormalized("cloud", parameters, 0.2)
    ));
    EXPECT_THROW(
        world.EvaluateFormFactor("cloud", parameters, 0.2),
        SEBException
    );
    EXPECT_THROW(
        world.EvaluateRadiusOfGyration2("cloud", parameters),
        SEBException
    );
}

// Verify that World composes an analytic SolidSphere and a numerical Debye
// cloud using the standard interference expression
// I_total = I_1 + I_2 + 2 A_1 Psi A_2, followed by total-beta normalization.
TEST(NumericalWorldTest, MixedAnalyticAndNumericalStructureComposes)
{
    DebyeSphereCloud* cloud = new DebyeSphereCloud({
        SphereScatterer(0.0, 0.0, 0.0, 0.0, 1.0),
        SphereScatterer(2.0, 0.0, 0.0, 0.0, 1.0)
    });
    cloud->addReferencePoint("left", 0.0, 0.0, 0.0);

    World world;
    const GraphID graph = world.Add(cloud, "cloud");
    world.Link(new SolidSphere(), "sphere.center", "cloud.left");
    world.Add(graph, "mixed");

    const ParameterList parameters{
        {"beta_sphere", 3.0},
        {"R_sphere", 1.0}
    };
    const double q = 0.2;

    const double cloudI =
        world.EvaluateFormFactorUnnormalized("cloud", parameters, q);
    const double sphereI =
        world.EvaluateFormFactorUnnormalized("sphere", parameters, q);
    const double cloudA =
        world.EvaluateFormFactorAmplitudeUnnormalized("cloud.left", parameters, q);
    const double sphereA =
        world.EvaluateFormFactorAmplitudeUnnormalized("sphere.center", parameters, q);
    const double expected = (cloudI + sphereI + 2.0 * cloudA * sphereA) / 25.0;

    EXPECT_NEAR(
        world.EvaluateFormFactor("mixed", parameters, q),
        expected,
        1e-11
    );
}
