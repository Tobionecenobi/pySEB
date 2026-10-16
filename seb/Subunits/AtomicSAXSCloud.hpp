#ifndef INCLUDE_GUARD_ATOMICSAXSCLOUD
#define INCLUDE_GUARD_ATOMICSAXSCLOUD

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Geometry/Point3D.hpp"
#include "Numerical.hpp"

// Cromer--Mann coefficients for
// f(q) = sum_k a_k exp[-b_k (q / 4 pi)^2] + c.
struct CromerMannCoefficients {
    std::array<double, 4> a{{0.0, 0.0, 0.0, 0.0}};
    std::array<double, 4> b{{0.0, 0.0, 0.0, 0.0}};
    double c = 0.0;

    CromerMannCoefficients() {}
    CromerMannCoefficients(
        const std::array<double, 4>& aValues,
        const std::array<double, 4>& bValues,
        double cValue)
        : a(aValues), b(bValues), c(cValue) {}

    double evaluate(double q) const;
    double valueAtZero() const;
    double secondMoment() const;
};

struct AtomicSAXSScatterer {
    CartesianPoint3D center;
    CromerMannCoefficients vacuum;
    double excludedVolume = 0.0;
    double occupancy = 1.0;
    double temperatureFactor = 0.0;
    double solventAccessibility = 0.0;
    int implicitHydrogens = 0;
    // Radius of the spherically averaged bonded-H distribution; zero
    // preserves the coincident-group approximation for custom scatterers.
    double implicitHydrogenDistance = 0.0;
    // Positive for spatial hydration quadrature points, zero for atoms.
    double hydrationVolume = 0.0;
    // Pepsi shell convention multiplies geometric volume by the current bulk density.
    bool hydrationUsesSolventDensity = false;
    std::string element;

    AtomicSAXSScatterer() {}
    AtomicSAXSScatterer(
        const CartesianPoint3D& centerValue,
        const CromerMannCoefficients& vacuumValue,
        double excludedVolumeValue,
        double occupancyValue = 1.0,
        double temperatureFactorValue = 0.0,
        double solventAccessibilityValue = 0.0,
        int implicitHydrogensValue = 0,
        const std::string& elementValue = "")
        : center(centerValue),
          vacuum(vacuumValue),
          excludedVolume(excludedVolumeValue),
          occupancy(occupancyValue),
          temperatureFactor(temperatureFactorValue),
          solventAccessibility(solventAccessibilityValue),
          implicitHydrogens(implicitHydrogensValue),
          element(elementValue) {}
};

enum class AtomicSAXSHydrationModel { AtomCentered, SpatialShell };
enum class AtomicSAXSHydrationWaterModel {
    NeutralWater, BondedWater, FoXSReciprocal, FoXSFast
};

struct AtomicSAXSOptions {
    double solventElectronDensity = 0.334;
    double excludedVolumeScale = 1.0;
    bool includeExcludedSolvent = true;
    bool useTemperatureFactors = false;

    // Hydration is deliberately opt-in. Builders calculate fractional solvent
    // accessibility by default; direct callers may supply it themselves.
    bool includeHydration = false;
    AtomicSAXSHydrationModel hydrationModel = AtomicSAXSHydrationModel::AtomCentered;
    AtomicSAXSHydrationWaterModel hydrationWaterModel = AtomicSAXSHydrationWaterModel::NeutralWater;
    double hydrationShellContrast = 0.0;
    double hydrationVolumePerExposedAtom = 30.0;
};

class AtomicSAXSCloud : public NumericalSubunit {
public:
    AtomicSAXSCloud();
    explicit AtomicSAXSCloud(
        const std::vector<AtomicSAXSScatterer>& scatterers,
        const AtomicSAXSOptions& options = AtomicSAXSOptions());
    virtual ~AtomicSAXSCloud() {}

    AtomicSAXSCloud(const AtomicSAXSCloud&) = delete;
    AtomicSAXSCloud& operator=(const AtomicSAXSCloud&) = delete;

    void addScatterer(const AtomicSAXSScatterer& scatterer);
    void addReferencePoint(
        const refPoint& name,
        const CartesianPoint3D& coordinate);
    void addReferencePoint(
        const refPoint& name,
        double x,
        double y,
        double z);

    const std::vector<AtomicSAXSScatterer>& getScatterers() const;
    std::map<refPoint, CartesianPoint3D> getReferenceCoordinates() const;
    std::size_t scattererCount() const;

    AtomicSAXSOptions getOptions() const;
    void setOptions(const AtomicSAXSOptions& value);
    void setSolventElectronDensity(double value);
    void setExcludedVolumeScale(double value);
    void setIncludeExcludedSolvent(bool value);
    void setUseTemperatureFactors(bool value);
    void setHydration(
        bool enabled,
        double shellContrast = 0.0,
        double volumePerExposedAtom = 30.0);

    double atomicAmplitude(std::size_t index, double q) const;

    static double sinc(double x);
    static CromerMannCoefficients hydrogenCoefficients();

private:
    struct PairDistance {
        std::size_t first;
        std::size_t second;
        double distance;
        double distance2;
    };

    struct AmplitudeMoment {
        double beta;
        double secondMoment;
    };

    std::vector<AtomicSAXSScatterer> atoms;
    AtomicSAXSOptions options;
    mutable std::map<refPoint, CartesianPoint3D> referenceCoordinates;
    mutable std::vector<PairDistance> atomPairDistances;
    mutable std::map<refPoint, std::vector<double>> referenceAtomDistances;
    mutable std::map<refPoint, std::vector<double>> referenceAtomDistances2;
    mutable std::map<std::pair<refPoint, refPoint>, double> referencePairDistances;
    mutable std::map<std::pair<refPoint, refPoint>, double> referencePairDistances2;
    mutable bool geometryDirty = true;

    void configureCoreCallbacks();
    void configureReferenceCallbacks(const refPoint& name);
    void validate() const;
    void ensureGeometryPrepared() const;
    AmplitudeMoment amplitudeMoment(const AtomicSAXSScatterer& atom) const;
    AmplitudeMoment waterMoment() const;
    double normalizedWaterAmplitude(double q) const;
    double totalBeta() const;
    double nonzeroTotalBeta(const std::string& operation) const;
    double formFactorUnnormalized(double q) const;
    double amplitudeUnnormalized(const refPoint& reference, double q) const;
    double phaseFactor(
        const refPoint& first,
        const refPoint& second,
        double q) const;
    double radiusOfGyration2() const;
    double sigmaReferenceToScatterers(const refPoint& reference) const;
    double sigmaReferenceToReference(
        const refPoint& first,
        const refPoint& second) const;
};

#endif
