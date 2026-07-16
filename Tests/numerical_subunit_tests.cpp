#include <cmath>
#include <stdexcept>
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

double radicalInverse(std::size_t index, unsigned int base)
{
    double result = 0.0;
    double factor = 1.0 / static_cast<double>(base);
    while (index > 0) {
        result += factor * static_cast<double>(index % base);
        index /= base;
        factor /= static_cast<double>(base);
    }
    return result;
}

std::vector<SphereScatterer> sampleSolidSphere(
    double radius,
    double totalBeta,
    std::size_t scattererCount)
{
    const double pi = 3.14159265358979323846;
    const std::size_t pairCount = scattererCount / 2;
    const double betaPerScatterer =
        totalBeta / static_cast<double>(2 * pairCount);

    std::vector<SphereScatterer> scatterers;
    scatterers.reserve(2 * pairCount);
    for (std::size_t index = 1; index <= pairCount; ++index) {
        // Halton coordinates give a deterministic, approximately uniform
        // volume sampling. Antipodal pairs keep the cloud center exactly zero.
        const double radialFraction =
            (static_cast<double>(index) - 0.5) / static_cast<double>(pairCount);
        const double cosTheta = 1.0 - 2.0 * radicalInverse(index, 2);
        const double phi = 2.0 * pi * radicalInverse(index, 3);
        const double sampledRadius = radius * std::cbrt(radialFraction);
        const double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);

        const double x = sampledRadius * sinTheta * std::cos(phi);
        const double y = sampledRadius * sinTheta * std::sin(phi);
        const double z = sampledRadius * cosTheta;

        scatterers.emplace_back(x, y, z, 0.0, betaPerScatterer);
        scatterers.emplace_back(-x, -y, -z, 0.0, betaPerScatterer);
    }
    return scatterers;
}

std::vector<SphereScatterer> sampleThinRod(
    double length,
    double totalBeta,
    std::size_t scattererCount)
{
    const double betaPerScatterer =
        totalBeta / static_cast<double>(scattererCount);

    std::vector<SphereScatterer> scatterers;
    scatterers.reserve(scattererCount);
    for (std::size_t index = 0; index < scattererCount; ++index) {
        // Midpoint sampling gives equal scattering length per contour segment
        // and avoids placing a finite point mass exactly at either endpoint.
        const double contour =
            (static_cast<double>(index) + 0.5) /
            static_cast<double>(scattererCount);
        scatterers.emplace_back(
            length * contour,
            0.0,
            0.0,
            0.0,
            betaPerScatterer
        );
    }
    return scatterers;
}

DebyeSphereCloud* sampledSphereCloud(
    double radius,
    double beta)
{
    DebyeSphereCloud* cloud =
        new DebyeSphereCloud(sampleSolidSphere(radius, beta, 1200));
    cloud->addReferencePoint("surface", radius, 0.0, 0.0);
    return cloud;
}

DebyeSphereCloud* sampledRodCloud(
    double length,
    double beta,
    std::size_t scattererCount = 600)
{
    DebyeSphereCloud* cloud =
        new DebyeSphereCloud(
            sampleThinRod(length, beta, scattererCount)
        );
    cloud->addReferencePoint("end1", 0.0, 0.0, 0.0);
    cloud->addReferencePoint("middle", length / 2.0, 0.0, 0.0);
    cloud->addReferencePoint("end2", length, 0.0, 0.0);
    return cloud;
}

enum class RodDendrimerLayerMode {
    AllCloud,
    AllSymbolic,
    AlternatingWithSymbolicRoot
};

bool useCloudRodAtLevel(
    RodDendrimerLayerMode mode,
    int level,
    int levels)
{
    if (mode == RodDendrimerLayerMode::AllCloud) {
        return true;
    }
    if (mode == RodDendrimerLayerMode::AllSymbolic) {
        return false;
    }

    const int distanceFromRoot = levels - 1 - level;
    return distanceFromRoot % 2 == 1;
}

void buildBinaryRodDendrimer(
    World& world,
    RodDendrimerLayerMode mode,
    int levels,
    double length,
    double beta,
    std::size_t cloudScattererCount)
{
    if (levels < 1) {
        throw std::invalid_argument("Dendrimer must contain at least one level");
    }

    string currentRoot = "level0Rod";
    GraphID currentGraph = useCloudRodAtLevel(mode, 0, levels)
        ? world.Add(
            sampledRodCloud(length, beta, cloudScattererCount),
            currentRoot,
            "rod"
        )
        : world.Add(new ThinRod(), currentRoot, "rod");

    for (int level = 1; level < levels; ++level) {
        const string newRoot = "level" + std::to_string(level) + "Rod";
        const GraphID nextGraph = useCloudRodAtLevel(mode, level, levels)
            ? world.Add(
                sampledRodCloud(length, beta, cloudScattererCount),
                newRoot,
                "rod"
            )
            : world.Add(new ThinRod(), newRoot, "rod");

        const string leftBranch = "level" + std::to_string(level) + "Left";
        const string rightBranch = "level" + std::to_string(level) + "Right";
        world.Link(
            currentGraph,
            leftBranch + ":" + currentRoot + ".end1",
            newRoot + ".end2"
        );
        world.Link(
            currentGraph,
            rightBranch + ":" + currentRoot + ".end1",
            newRoot + ".end2"
        );

        currentGraph = nextGraph;
        currentRoot = newRoot;
    }

    world.Add(currentGraph, "dendrimer");
}

void expectRodSphereStructureMatchesAnalytic(
    World& candidate,
    const ParameterList& parameters,
    double formFactorTolerance,
    double amplitudeTolerance,
    double rg2Tolerance)
{
    World analytic;
    const GraphID graph = analytic.Add(new ThinRod(), "rod");
    analytic.Link(
        new SolidSphere(),
        "sphere.surface#join",
        "rod.end1"
    );
    analytic.Add(graph, "pair");

    for (const double q : std::vector<double>{0.0, 0.05, 0.15, 0.3}) {
        EXPECT_NEAR(
            candidate.EvaluateFormFactor("pair", parameters, q),
            analytic.EvaluateFormFactor("pair", parameters, q),
            formFactorTolerance
        );
        EXPECT_NEAR(
            candidate.EvaluateFormFactorAmplitude(
                "pair:rod.end1",
                parameters,
                q
            ),
            analytic.EvaluateFormFactorAmplitude(
                "pair:rod.end1",
                parameters,
                q
            ),
            amplitudeTolerance
        );
    }

    EXPECT_NEAR(
        candidate.EvaluateRadiusOfGyration2("pair", parameters),
        analytic.EvaluateRadiusOfGyration2("pair", parameters),
        rg2Tolerance
    );
}

void expectRodRodStructureMatchesAnalytic(
    World& candidate,
    const ParameterList& parameters,
    double formFactorTolerance,
    double amplitudeTolerance,
    double rg2Tolerance)
{
    World analytic;
    const GraphID graph = analytic.Add(new ThinRod(), "rod1");
    analytic.Link(
        new ThinRod(),
        "rod2.end1",
        "rod1.end2"
    );
    analytic.Add(graph, "pair");

    for (const double q : std::vector<double>{0.0, 0.05, 0.15, 0.3}) {
        EXPECT_NEAR(
            candidate.EvaluateFormFactor("pair", parameters, q),
            analytic.EvaluateFormFactor("pair", parameters, q),
            formFactorTolerance
        );
        EXPECT_NEAR(
            candidate.EvaluateFormFactorAmplitude(
                "pair:rod1.end1",
                parameters,
                q
            ),
            analytic.EvaluateFormFactorAmplitude(
                "pair:rod1.end1",
                parameters,
                q
            ),
            amplitudeTolerance
        );
    }

    EXPECT_NEAR(
        candidate.EvaluateRadiusOfGyration2("pair", parameters),
        analytic.EvaluateRadiusOfGyration2("pair", parameters),
        rg2Tolerance
    );
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

// Verify that a deterministic point cloud filling a solid spherical volume
// converges to the analytic SolidSphere form factor, center amplitude, and
// Rg^2. The tolerance accounts for the finite spatial discretization.
TEST(DebyeSphereCloudTest, SampledSolidSphereCloudApproximatesRayleighSphere)
{
    const double radius = 2.0;
    const double beta = 3.0;
    const std::vector<SphereScatterer> scatterers =
        sampleSolidSphere(radius, beta, 1200);

    World world;
    world.Add(new DebyeSphereCloud(scatterers), "cloud");
    world.Add(new SolidSphere(), "sphere");

    ParameterList parameters{
        {"beta_sphere", beta},
        {"R_sphere", radius}
    };

    for (const double q : std::vector<double>{0.0, 0.1, 0.3, 0.6}) {
        EXPECT_NEAR(
            world.EvaluateFormFactor("cloud", parameters, q),
            world.EvaluateFormFactor("sphere", parameters, q),
            1e-3
        );
        EXPECT_NEAR(
            world.EvaluateFormFactorAmplitude("cloud.center", parameters, q),
            world.EvaluateFormFactorAmplitude("sphere.center", parameters, q),
            6e-4
        );
    }

    EXPECT_NEAR(
        world.EvaluateRadiusOfGyration2("cloud", parameters),
        12.0 / 5.0,
        1e-2
    );
}

// Verify the basic rod discretization independently of structure composition.
// A line-sampled point cloud must reproduce the analytic ThinRod form factor,
// amplitudes from its named references, end-to-end phase factor, and Rg^2.
TEST(DebyeSphereCloudTest, SampledThinRodCloudApproximatesAnalyticThinRod)
{
    const double length = 5.0;
    const double beta = 2.0;

    World world;
    world.Add(sampledRodCloud(length, beta), "cloud");
    world.Add(new ThinRod(), "rod");

    const ParameterList parameters{
        {"beta_rod", beta},
        {"L_rod", length}
    };

    for (const double q : std::vector<double>{0.0, 0.05, 0.15, 0.3, 0.6}) {
        EXPECT_NEAR(
            world.EvaluateFormFactor("cloud", parameters, q),
            world.EvaluateFormFactor("rod", parameters, q),
            1e-4
        );
        for (const string& reference :
             std::vector<string>{"end1", "middle", "end2"}) {
            EXPECT_NEAR(
                world.EvaluateFormFactorAmplitude(
                    "cloud." + reference,
                    parameters,
                    q
                ),
                world.EvaluateFormFactorAmplitude(
                    "rod." + reference,
                    parameters,
                    q
                ),
                1e-4
            );
        }
        EXPECT_NEAR(
            world.EvaluatePhaseFactor(
                "cloud.end1",
                "cloud.end2",
                parameters,
                q
            ),
            world.EvaluatePhaseFactor(
                "rod.end1",
                "rod.end2",
                parameters,
                q
            ),
            1e-12
        );
    }

    EXPECT_NEAR(
        world.EvaluateRadiusOfGyration2("cloud", parameters),
        world.EvaluateRadiusOfGyration2("rod", parameters),
        1e-5
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

// Verify that numerical clouds and analytic subunits compose identically at
// structure level. Two sampled solid spheres are linked at surface reference
// points and compared with two analytic SolidSphere subunits linked through
// their corresponding distributed surface references.
TEST(NumericalWorldTest, TwoConnectedSphereCloudsApproximateTwoAnalyticSpheres)
{
    const double radius = 2.0;
    const double beta = 3.0;

    DebyeSphereCloud* cloud1 =
        new DebyeSphereCloud(sampleSolidSphere(radius, beta, 1200));
    DebyeSphereCloud* cloud2 =
        new DebyeSphereCloud(sampleSolidSphere(radius, beta, 1200));
    cloud1->addReferencePoint("surface", radius, 0.0, 0.0);
    cloud2->addReferencePoint("surface", radius, 0.0, 0.0);

    World cloudWorld;
    const GraphID cloudGraph = cloudWorld.Add(cloud1, "cloud1");
    cloudWorld.Link(cloud2, "cloud2.surface", "cloud1.surface");
    cloudWorld.Add(cloudGraph, "cloudPair");

    World analyticWorld;
    const GraphID analyticGraph =
        analyticWorld.Add(new SolidSphere(), "sphere1");
    analyticWorld.Link(
        new SolidSphere(),
        "sphere2.surface#join",
        "sphere1.surface#join"
    );
    analyticWorld.Add(analyticGraph, "analyticPair");

    const ParameterList parameters{
        {"beta_sphere1", beta},
        {"beta_sphere2", beta},
        {"R_sphere1", radius},
        {"R_sphere2", radius}
    };

    for (const double q : std::vector<double>{0.0, 0.1, 0.3, 0.6}) {
        EXPECT_NEAR(
            cloudWorld.EvaluateFormFactor("cloudPair", parameters, q),
            analyticWorld.EvaluateFormFactor("analyticPair", parameters, q),
            2e-3
        );
        EXPECT_NEAR(
            cloudWorld.EvaluateFormFactorAmplitude(
                "cloudPair:cloud1.surface",
                parameters,
                q
            ),
            analyticWorld.EvaluateFormFactorAmplitude(
                "analyticPair:sphere1.surface#join",
                parameters,
                q
            ),
            1e-3
        );
    }

    EXPECT_NEAR(
        cloudWorld.EvaluateRadiusOfGyration2("cloudPair", parameters),
        analyticWorld.EvaluateRadiusOfGyration2("analyticPair", parameters),
        2e-2
    );
}

// Verify hybrid composition: replacing only one sphere in an analytic
// two-sphere structure by a sampled cloud should preserve the total form
// factor, connection-point amplitude, and Rg^2 within sampling accuracy.
TEST(NumericalWorldTest, ConnectedCloudAndAnalyticSphereApproximateTwoAnalyticSpheres)
{
    const double radius = 2.0;
    const double beta = 3.0;

    DebyeSphereCloud* cloud =
        new DebyeSphereCloud(sampleSolidSphere(radius, beta, 1200));
    cloud->addReferencePoint("surface", radius, 0.0, 0.0);

    World mixedWorld;
    const GraphID mixedGraph = mixedWorld.Add(cloud, "cloud");
    mixedWorld.Link(
        new SolidSphere(),
        "sphere2.surface#join",
        "cloud.surface"
    );
    mixedWorld.Add(mixedGraph, "mixedPair");

    World analyticWorld;
    const GraphID analyticGraph =
        analyticWorld.Add(new SolidSphere(), "sphere1");
    analyticWorld.Link(
        new SolidSphere(),
        "sphere2.surface#join",
        "sphere1.surface#join"
    );
    analyticWorld.Add(analyticGraph, "analyticPair");

    const ParameterList parameters{
        {"beta_sphere1", beta},
        {"beta_sphere2", beta},
        {"R_sphere1", radius},
        {"R_sphere2", radius}
    };

    for (const double q : std::vector<double>{0.0, 0.1, 0.3, 0.6}) {
        EXPECT_NEAR(
            mixedWorld.EvaluateFormFactor("mixedPair", parameters, q),
            analyticWorld.EvaluateFormFactor("analyticPair", parameters, q),
            1e-3
        );
        EXPECT_NEAR(
            mixedWorld.EvaluateFormFactorAmplitude(
                "mixedPair:cloud.surface",
                parameters,
                q
            ),
            analyticWorld.EvaluateFormFactorAmplitude(
                "analyticPair:sphere1.surface#join",
                parameters,
                q
            ),
            6e-4
        );
    }

    EXPECT_NEAR(
        mixedWorld.EvaluateRadiusOfGyration2("mixedPair", parameters),
        analyticWorld.EvaluateRadiusOfGyration2("analyticPair", parameters),
        1e-2
    );
}

// Verify a fully numerical rod-sphere structure. A line-sampled rod cloud and
// a volume-sampled sphere cloud are connected at rod.end1 and sphere.surface,
// then compared with the corresponding ThinRod-SolidSphere structure.
TEST(NumericalWorldTest, ConnectedRodCloudAndSphereCloudApproximateAnalyticPair)
{
    const double rodLength = 5.0;
    const double rodBeta = 2.0;
    const double sphereRadius = 2.0;
    const double sphereBeta = 3.0;

    World candidate;
    const GraphID graph = candidate.Add(
        sampledRodCloud(rodLength, rodBeta),
        "rod"
    );
    candidate.Link(
        sampledSphereCloud(sphereRadius, sphereBeta),
        "sphere.surface",
        "rod.end1"
    );
    candidate.Add(graph, "pair");

    const ParameterList parameters{
        {"beta_rod", rodBeta},
        {"L_rod", rodLength},
        {"beta_sphere", sphereBeta},
        {"R_sphere", sphereRadius}
    };
    expectRodSphereStructureMatchesAnalytic(
        candidate,
        parameters,
        3e-3,
        2e-3,
        3e-2
    );
}

// Verify the cloud-symbolic ordering. Replacing only the rod by its sampled
// cloud must give the same rod-sphere structure scattering as ThinRod linked
// to SolidSphere.
TEST(NumericalWorldTest, ConnectedRodCloudAndAnalyticSphereApproximateAnalyticPair)
{
    const double rodLength = 5.0;
    const double rodBeta = 2.0;
    const double sphereRadius = 2.0;
    const double sphereBeta = 3.0;

    World candidate;
    const GraphID graph = candidate.Add(
        sampledRodCloud(rodLength, rodBeta),
        "rod"
    );
    candidate.Link(
        new SolidSphere(),
        "sphere.surface#join",
        "rod.end1"
    );
    candidate.Add(graph, "pair");

    const ParameterList parameters{
        {"beta_rod", rodBeta},
        {"L_rod", rodLength},
        {"beta_sphere", sphereBeta},
        {"R_sphere", sphereRadius}
    };
    expectRodSphereStructureMatchesAnalytic(
        candidate,
        parameters,
        2e-3,
        1e-3,
        2e-2
    );
}

// Verify the symbolic-cloud ordering. Replacing only the sphere by its sampled
// cloud must give the same rod-sphere structure scattering as ThinRod linked
// to SolidSphere.
TEST(NumericalWorldTest, ConnectedAnalyticRodAndSphereCloudApproximateAnalyticPair)
{
    const double rodLength = 5.0;
    const double rodBeta = 2.0;
    const double sphereRadius = 2.0;
    const double sphereBeta = 3.0;

    World candidate;
    const GraphID graph = candidate.Add(new ThinRod(), "rod");
    candidate.Link(
        sampledSphereCloud(sphereRadius, sphereBeta),
        "sphere.surface",
        "rod.end1"
    );
    candidate.Add(graph, "pair");

    const ParameterList parameters{
        {"beta_rod", rodBeta},
        {"L_rod", rodLength},
        {"beta_sphere", sphereBeta},
        {"R_sphere", sphereRadius}
    };
    expectRodSphereStructureMatchesAnalytic(
        candidate,
        parameters,
        2e-3,
        1e-3,
        2e-2
    );
}

// Verify a fully numerical two-rod structure under World's factorized,
// orientationally averaged composition. Two sampled rod clouds linked
// end-to-end must reproduce two linked analytic ThinRod subunits.
TEST(NumericalWorldTest, ConnectedRodCloudsApproximateConnectedAnalyticRods)
{
    const double rod1Length = 4.0;
    const double rod1Beta = 2.0;
    const double rod2Length = 6.0;
    const double rod2Beta = 1.5;

    World candidate;
    const GraphID graph = candidate.Add(
        sampledRodCloud(rod1Length, rod1Beta),
        "rod1"
    );
    candidate.Link(
        sampledRodCloud(rod2Length, rod2Beta),
        "rod2.end1",
        "rod1.end2"
    );
    candidate.Add(graph, "pair");

    const ParameterList parameters{
        {"beta_rod1", rod1Beta},
        {"L_rod1", rod1Length},
        {"beta_rod2", rod2Beta},
        {"L_rod2", rod2Length}
    };
    expectRodRodStructureMatchesAnalytic(
        candidate,
        parameters,
        3e-3,
        2e-3,
        2e-2
    );
}

// Verify the cloud-symbolic ordering for two rods. Replacing rod1 by a
// line-sampled cloud must preserve the scattering of the analytic
// ThinRod-ThinRod structure.
TEST(NumericalWorldTest, ConnectedRodCloudAndAnalyticRodApproximateAnalyticRods)
{
    const double rod1Length = 4.0;
    const double rod1Beta = 2.0;
    const double rod2Length = 6.0;
    const double rod2Beta = 1.5;

    World candidate;
    const GraphID graph = candidate.Add(
        sampledRodCloud(rod1Length, rod1Beta),
        "rod1"
    );
    candidate.Link(
        new ThinRod(),
        "rod2.end1",
        "rod1.end2"
    );
    candidate.Add(graph, "pair");

    const ParameterList parameters{
        {"beta_rod1", rod1Beta},
        {"L_rod1", rod1Length},
        {"beta_rod2", rod2Beta},
        {"L_rod2", rod2Length}
    };
    expectRodRodStructureMatchesAnalytic(
        candidate,
        parameters,
        2e-3,
        1e-3,
        1e-2
    );
}

// Verify the symbolic-cloud ordering for two rods. Replacing rod2 by a
// line-sampled cloud must preserve the scattering of the analytic
// ThinRod-ThinRod structure.
TEST(NumericalWorldTest, ConnectedAnalyticRodAndRodCloudApproximateAnalyticRods)
{
    const double rod1Length = 4.0;
    const double rod1Beta = 2.0;
    const double rod2Length = 6.0;
    const double rod2Beta = 1.5;

    World candidate;
    const GraphID graph = candidate.Add(new ThinRod(), "rod1");
    candidate.Link(
        sampledRodCloud(rod2Length, rod2Beta),
        "rod2.end1",
        "rod1.end2"
    );
    candidate.Add(graph, "pair");

    const ParameterList parameters{
        {"beta_rod1", rod1Beta},
        {"L_rod1", rod1Length},
        {"beta_rod2", rod2Beta},
        {"L_rod2", rod2Length}
    };
    expectRodRodStructureMatchesAnalytic(
        candidate,
        parameters,
        2e-3,
        1e-3,
        1e-2
    );
}

// Verify hierarchical composition for a binary rod dendrimer at the configured
// depth, with each rod producing two child rods from end2. A dendrimer built
// from sampled rod clouds must reproduce the same structure built from
// analytic ThinRod subunits.
TEST(NumericalWorldTest, BinaryRodCloudDendrimerApproximatesAnalyticDendrimer)
{
    const int levels = 7;
    const double rodLength = 3.0;
    const double rodBeta = 1.0;
    const std::size_t effectiveRodCount =
        (static_cast<std::size_t>(1) << levels) - 1;
    const string rootReference =
        "dendrimer:level" + std::to_string(levels - 1) + "Rod.end1";

    World cloudWorld;
    buildBinaryRodDendrimer(
        cloudWorld,
        RodDendrimerLayerMode::AllCloud,
        levels,
        rodLength,
        rodBeta,
        300
    );

    World analyticWorld;
    buildBinaryRodDendrimer(
        analyticWorld,
        RodDendrimerLayerMode::AllSymbolic,
        levels,
        rodLength,
        rodBeta,
        0
    );

    const ParameterList parameters{
        {"beta_rod", rodBeta},
        {"L_rod", rodLength}
    };

    EXPECT_NEAR(
        cloudWorld.EvaluateFormFactorUnnormalized(
            "dendrimer",
            parameters,
            0.0
        ),
        static_cast<double>(effectiveRodCount * effectiveRodCount),
        1e-10
    );

    for (const double q : std::vector<double>{0.0, 0.05, 0.15, 0.3}) {
        EXPECT_NEAR(
            cloudWorld.EvaluateFormFactor("dendrimer", parameters, q),
            analyticWorld.EvaluateFormFactor("dendrimer", parameters, q),
            4e-3
        );
        EXPECT_NEAR(
            cloudWorld.EvaluateFormFactorAmplitude(
                rootReference,
                parameters,
                q
            ),
            analyticWorld.EvaluateFormFactorAmplitude(
                rootReference,
                parameters,
                q
            ),
            3e-3
        );
    }

    EXPECT_NEAR(
        cloudWorld.EvaluateRadiusOfGyration2("dendrimer", parameters),
        analyticWorld.EvaluateRadiusOfGyration2("dendrimer", parameters),
        3e-2
    );
}

// Verify that numerical and analytic subunits can alternate through every
// recursive layer of a dendrimer. The root layer is symbolic, the next layer
// uses rod clouds, and the pattern repeats down to the leaves.
TEST(NumericalWorldTest, AlternatingRodCloudAndSymbolicDendrimerApproximatesAnalyticDendrimer)
{
    const int levels = 7;
    const double rodLength = 3.0;
    const double rodBeta = 1.0;
    const std::size_t effectiveRodCount =
        (static_cast<std::size_t>(1) << levels) - 1;
    const string rootReference =
        "dendrimer:level" + std::to_string(levels - 1) + "Rod.end1";

    World alternatingWorld;
    buildBinaryRodDendrimer(
        alternatingWorld,
        RodDendrimerLayerMode::AlternatingWithSymbolicRoot,
        levels,
        rodLength,
        rodBeta,
        300
    );

    World analyticWorld;
    buildBinaryRodDendrimer(
        analyticWorld,
        RodDendrimerLayerMode::AllSymbolic,
        levels,
        rodLength,
        rodBeta,
        0
    );

    const ParameterList parameters{
        {"beta_rod", rodBeta},
        {"L_rod", rodLength}
    };

    EXPECT_NEAR(
        alternatingWorld.EvaluateFormFactorUnnormalized(
            "dendrimer",
            parameters,
            0.0
        ),
        static_cast<double>(effectiveRodCount * effectiveRodCount),
        1e-10
    );

    for (const double q : std::vector<double>{0.0, 0.05, 0.15, 0.3}) {
        EXPECT_NEAR(
            alternatingWorld.EvaluateFormFactor(
                "dendrimer",
                parameters,
                q
            ),
            analyticWorld.EvaluateFormFactor(
                "dendrimer",
                parameters,
                q
            ),
            4e-3
        );
        EXPECT_NEAR(
            alternatingWorld.EvaluateFormFactorAmplitude(
                rootReference,
                parameters,
                q
            ),
            analyticWorld.EvaluateFormFactorAmplitude(
                rootReference,
                parameters,
                q
            ),
            3e-3
        );
    }

    EXPECT_NEAR(
        alternatingWorld.EvaluateRadiusOfGyration2(
            "dendrimer",
            parameters
        ),
        analyticWorld.EvaluateRadiusOfGyration2(
            "dendrimer",
            parameters
        ),
        3e-2
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
