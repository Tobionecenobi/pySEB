#include "AtomicSAXSCloud.hpp"

#include <cmath>
#include <limits>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void validatePoint(const CartesianPoint3D& point, const std::string& description)
{
    if (!std::isfinite(point.x) ||
        !std::isfinite(point.y) ||
        !std::isfinite(point.z)) {
        throw SEBException(
            "Non-finite coordinate in " + description,
            "AtomicSAXSCloud input validation"
        );
    }
}

double squaredDistance(
    const CartesianPoint3D& left,
    const CartesianPoint3D& right)
{
    const double dx = left.x - right.x;
    const double dy = left.y - right.y;
    const double dz = left.z - right.z;
    return dx * dx + dy * dy + dz * dz;
}

void validateOptions(const AtomicSAXSOptions& options)
{
    if (options.hydrationWaterModel != AtomicSAXSHydrationWaterModel::NeutralWater &&
        options.hydrationWaterModel != AtomicSAXSHydrationWaterModel::BondedWater &&
        options.hydrationWaterModel != AtomicSAXSHydrationWaterModel::FoXSReciprocal &&
        options.hydrationWaterModel != AtomicSAXSHydrationWaterModel::FoXSFast) {
        throw SEBException("Unknown hydration water model", "AtomicSAXSCloud options");
    }
    if (options.hydrationModel != AtomicSAXSHydrationModel::AtomCentered &&
        options.hydrationModel != AtomicSAXSHydrationModel::SpatialShell) {
        throw SEBException("Unknown hydration model", "AtomicSAXSCloud options");
    }
    if (!std::isfinite(options.solventElectronDensity) ||
        options.solventElectronDensity < 0.0) {
        throw SEBException(
            "Solvent electron density must be finite and non-negative",
            "AtomicSAXSCloud options"
        );
    }
    if (!std::isfinite(options.excludedVolumeScale) ||
        options.excludedVolumeScale <= 0.0) {
        throw SEBException(
            "Excluded-volume scale must be finite and positive",
            "AtomicSAXSCloud options"
        );
    }
    if (!std::isfinite(options.hydrationShellContrast)) {
        throw SEBException(
            "Hydration-shell contrast must be finite",
            "AtomicSAXSCloud options"
        );
    }
    if (!std::isfinite(options.hydrationVolumePerExposedAtom) ||
        options.hydrationVolumePerExposedAtom < 0.0) {
        throw SEBException(
            "Hydration volume per exposed atom must be finite and non-negative",
            "AtomicSAXSCloud options"
        );
    }
}

} // namespace

double CromerMannCoefficients::evaluate(double q) const
{
    if (!std::isfinite(q)) {
        throw SEBException("q must be finite", "CromerMannCoefficients::evaluate()");
    }
    const double s2 = q * q / (16.0 * kPi * kPi);
    double result = c;
    for (std::size_t index = 0; index < a.size(); ++index) {
        result += a[index] * std::exp(-b[index] * s2);
    }
    return result;
}

double CromerMannCoefficients::valueAtZero() const
{
    double result = c;
    for (double coefficient : a) {
        result += coefficient;
    }
    return result;
}

double CromerMannCoefficients::secondMoment() const
{
    double weightedWidth = 0.0;
    for (std::size_t index = 0; index < a.size(); ++index) {
        weightedWidth += a[index] * b[index];
    }
    return 3.0 * weightedWidth / (8.0 * kPi * kPi);
}

AtomicSAXSCloud::AtomicSAXSCloud()
    : NumericalSubunit(NormalizationMode::Unnormalized)
{
    stype = ATOMICSAXSCLOUD;
    setReferencePointName("center");
    referenceCoordinates["center"] = CartesianPoint3D();
    configureCoreCallbacks();
    configureReferenceCallbacks("center");
}

AtomicSAXSCloud::AtomicSAXSCloud(
    const std::vector<AtomicSAXSScatterer>& scatterers,
    const AtomicSAXSOptions& optionValues)
    : AtomicSAXSCloud()
{
    atoms = scatterers;
    options = optionValues;
    geometryDirty = true;
}

void AtomicSAXSCloud::addScatterer(const AtomicSAXSScatterer& scatterer)
{
    atoms.push_back(scatterer);
    geometryDirty = true;
}

void AtomicSAXSCloud::addReferencePoint(
    const refPoint& name,
    const CartesianPoint3D& coordinate)
{
    if (hasReference(name)) {
        throw SEBException(
            "Reference point " + name + " already exists",
            "AtomicSAXSCloud::addReferencePoint()"
        );
    }
    validatePoint(coordinate, "reference point " + name);
    setReferencePointName(name);
    referenceCoordinates[name] = coordinate;
    configureReferenceCallbacks(name);
    geometryDirty = true;
}

void AtomicSAXSCloud::addReferencePoint(
    const refPoint& name,
    double x,
    double y,
    double z)
{
    addReferencePoint(name, CartesianPoint3D(x, y, z));
}

const std::vector<AtomicSAXSScatterer>& AtomicSAXSCloud::getScatterers() const
{
    return atoms;
}

std::map<refPoint, CartesianPoint3D>
AtomicSAXSCloud::getReferenceCoordinates() const
{
    ensureGeometryPrepared();
    return referenceCoordinates;
}

std::size_t AtomicSAXSCloud::scattererCount() const
{
    return atoms.size();
}

AtomicSAXSOptions AtomicSAXSCloud::getOptions() const
{
    return options;
}

void AtomicSAXSCloud::setOptions(const AtomicSAXSOptions& value)
{
    validateOptions(value);
    if (value.includeHydration &&
        value.hydrationModel == AtomicSAXSHydrationModel::SpatialShell) {
        bool found = false;
        for (const auto& atom : atoms) found = found || atom.hydrationVolume > 0.0;
        if (!found) throw SEBException(
            "Spatial hydration requires shell points; rebuild with SpatialShell",
            "AtomicSAXSCloud::setOptions()");
    }
    options = value;
}

void AtomicSAXSCloud::setSolventElectronDensity(double value)
{
    AtomicSAXSOptions updated = options;
    updated.solventElectronDensity = value;
    setOptions(updated);
}

void AtomicSAXSCloud::setExcludedVolumeScale(double value)
{
    AtomicSAXSOptions updated = options;
    updated.excludedVolumeScale = value;
    setOptions(updated);
}

void AtomicSAXSCloud::setIncludeExcludedSolvent(bool value)
{
    options.includeExcludedSolvent = value;
}

void AtomicSAXSCloud::setUseTemperatureFactors(bool value)
{
    options.useTemperatureFactors = value;
}

void AtomicSAXSCloud::setHydration(
    bool enabled,
    double shellContrast,
    double volumePerExposedAtom)
{
    AtomicSAXSOptions updated = options;
    updated.includeHydration = enabled;
    updated.hydrationShellContrast = shellContrast;
    updated.hydrationVolumePerExposedAtom = volumePerExposedAtom;
    setOptions(updated);
}

double AtomicSAXSCloud::sinc(double x)
{
    const double x2 = x * x;
    if (x2 < 1e-8) {
        return 1.0 - x2 / 6.0 + x2 * x2 / 120.0;
    }
    return std::sin(x) / x;
}

CromerMannCoefficients AtomicSAXSCloud::hydrogenCoefficients()
{
    return CromerMannCoefficients(
        {{0.493002, 0.322912, 0.140191, 0.040810}},
        {{10.5109, 26.1257, 3.14236, 57.7997}},
        0.003038
    );
}

void AtomicSAXSCloud::configureCoreCallbacks()
{
    setTotalBetaProvider([this](const ParameterList&) {
        return totalBeta();
    });
    setFormFactorFunction([this](double q, const ParameterList&) {
        return formFactorUnnormalized(q);
    });
    setRadiusOfGyration2Provider([this](const ParameterList&) {
        return radiusOfGyration2();
    });
}

void AtomicSAXSCloud::configureReferenceCallbacks(const refPoint& name)
{
    setFormFactorAmplitudeFunction(
        name,
        [this, name](double q, const ParameterList&) {
            return amplitudeUnnormalized(name, q);
        }
    );
    setSigmaMSDRef2ScatProvider(
        name,
        [this, name](const ParameterList&) {
            return sigmaReferenceToScatterers(name);
        }
    );

    for (const auto& reference : referenceCoordinates) {
        if (reference.first == name) {
            continue;
        }
        const refPoint other = reference.first;
        setPhaseFactorFunction(
            name,
            other,
            [this, name, other](double q, const ParameterList&) {
                return phaseFactor(name, other, q);
            }
        );
        setSigmaMSDRef2RefProvider(
            name,
            other,
            [this, name, other](const ParameterList&) {
                return sigmaReferenceToReference(name, other);
            }
        );
    }
}

void AtomicSAXSCloud::validate() const
{
    validateOptions(options);
    if (atoms.empty()) {
        throw SEBException(
            "An atomic SAXS cloud must contain at least one scatterer",
            "AtomicSAXSCloud input validation"
        );
    }

    for (const auto& atom : atoms) {
        validatePoint(atom.center, "atom center");
        for (std::size_t index = 0; index < atom.vacuum.a.size(); ++index) {
            if (!std::isfinite(atom.vacuum.a[index]) ||
                !std::isfinite(atom.vacuum.b[index]) ||
                atom.vacuum.b[index] < 0.0) {
                throw SEBException(
                    "Cromer--Mann coefficients must be finite and b values non-negative",
                    "AtomicSAXSCloud input validation"
                );
            }
        }
        if (!std::isfinite(atom.vacuum.c) ||
            !std::isfinite(atom.excludedVolume) || atom.excludedVolume < 0.0 ||
            !std::isfinite(atom.occupancy) || atom.occupancy < 0.0 ||
            !std::isfinite(atom.temperatureFactor) || atom.temperatureFactor < 0.0 ||
            !std::isfinite(atom.solventAccessibility) ||
                atom.solventAccessibility < 0.0 || atom.solventAccessibility > 1.0 ||
            atom.implicitHydrogens < 0 ||
            !std::isfinite(atom.implicitHydrogenDistance) ||
            atom.implicitHydrogenDistance < 0.0) {
            throw SEBException(
                "Invalid atomic SAXS scatterer parameters",
                "AtomicSAXSCloud input validation"
            );
        }
        if (!std::isfinite(atom.hydrationVolume) || atom.hydrationVolume < 0.0) {
            throw SEBException("Hydration point volume must be finite and non-negative",
                               "AtomicSAXSCloud input validation");
        }
        if (atom.hydrationVolume > 0.0 &&
            (atom.vacuum.valueAtZero() != 0.0 || atom.vacuum.secondMoment() != 0.0 ||
             atom.excludedVolume != 0.0 || atom.implicitHydrogens != 0 ||
             atom.solventAccessibility != 0.0)) {
            throw SEBException("Shell points must contain only hydration volume",
                               "AtomicSAXSCloud input validation");
        }
    }
    if (options.includeHydration &&
        options.hydrationModel == AtomicSAXSHydrationModel::SpatialShell) {
        bool found = false;
        for (const auto& atom : atoms) found = found || atom.hydrationVolume > 0.0;
        if (!found) throw SEBException("Spatial hydration requires shell points",
                                      "AtomicSAXSCloud input validation");
    }
}

void AtomicSAXSCloud::ensureGeometryPrepared() const
{
    if (!geometryDirty) {
        validateOptions(options);
        return;
    }
    validate();

    CartesianPoint3D centroid;
    std::size_t atomCount = 0;
    for (const auto& atom : atoms) {
        if (atom.hydrationVolume > 0.0) continue;
        ++atomCount;
        centroid.x += atom.center.x;
        centroid.y += atom.center.y;
        centroid.z += atom.center.z;
    }
    const double inverseCount = 1.0 / static_cast<double>(atomCount ? atomCount : atoms.size());
    if (!atomCount) {
        for (const auto& atom : atoms) {
            centroid.x += atom.center.x;
            centroid.y += atom.center.y;
            centroid.z += atom.center.z;
        }
    }
    centroid.x *= inverseCount;
    centroid.y *= inverseCount;
    centroid.z *= inverseCount;
    referenceCoordinates["center"] = centroid;

    atomPairDistances.clear();
    for (std::size_t first = 0; first < atoms.size(); ++first) {
        for (std::size_t second = first + 1; second < atoms.size(); ++second) {
            const double distance2 = squaredDistance(
                atoms[first].center,
                atoms[second].center
            );
            atomPairDistances.push_back(
                PairDistance{first, second, std::sqrt(distance2), distance2}
            );
        }
    }

    referenceAtomDistances.clear();
    referenceAtomDistances2.clear();
    for (const auto& reference : referenceCoordinates) {
        std::vector<double> distances;
        std::vector<double> distances2;
        distances.reserve(atoms.size());
        distances2.reserve(atoms.size());
        for (const auto& atom : atoms) {
            const double distance2 = squaredDistance(reference.second, atom.center);
            distances.push_back(std::sqrt(distance2));
            distances2.push_back(distance2);
        }
        referenceAtomDistances[reference.first] = distances;
        referenceAtomDistances2[reference.first] = distances2;
    }

    referencePairDistances.clear();
    referencePairDistances2.clear();
    for (auto first = referenceCoordinates.begin();
         first != referenceCoordinates.end(); ++first) {
        auto second = first;
        ++second;
        for (; second != referenceCoordinates.end(); ++second) {
            const auto pair = referencePair(first->first, second->first);
            const double distance2 = squaredDistance(first->second, second->second);
            referencePairDistances[pair] = std::sqrt(distance2);
            referencePairDistances2[pair] = distance2;
        }
    }

    geometryDirty = false;
}

double AtomicSAXSCloud::atomicAmplitude(std::size_t index, double q) const
{
    ensureGeometryPrepared();
    if (index >= atoms.size()) {
        throw SEBException(
            "Atomic scatterer index is out of range",
            "AtomicSAXSCloud::atomicAmplitude()"
        );
    }
    const AtomicSAXSScatterer& atom = atoms[index];
    const bool spatial = options.hydrationModel == AtomicSAXSHydrationModel::SpatialShell;
    if (atom.hydrationVolume > 0.0 && (!spatial || !options.includeHydration)) return 0.0;
    const CromerMannCoefficients hydrogen = hydrogenCoefficients();
    double vacuum = atom.vacuum.evaluate(q) +
        static_cast<double>(atom.implicitHydrogens) * hydrogen.evaluate(q) *
            sinc(q * atom.implicitHydrogenDistance);
    if (options.useTemperatureFactors && atom.temperatureFactor > 0.0) {
        vacuum *= std::exp(
            -atom.temperatureFactor * q * q / (16.0 * kPi * kPi)
        );
    }

    double result = vacuum;
    if (options.includeExcludedSolvent && atom.excludedVolume > 0.0) {
        const double volume = options.excludedVolumeScale * atom.excludedVolume;
        result -= options.solventElectronDensity * volume *
            std::exp(-q * q * std::pow(volume, 2.0 / 3.0) / (4.0 * kPi));
    }

    if (options.includeHydration &&
        ((spatial && atom.hydrationVolume > 0.0) ||
         (!spatial && atom.solventAccessibility > 0.0))) {
        const double volume = spatial ? atom.hydrationVolume :
            options.hydrationVolumePerExposedAtom * atom.solventAccessibility;
        const double density = atom.hydrationUsesSolventDensity ? options.solventElectronDensity : 1.0;
        result += options.hydrationShellContrast * volume * density * normalizedWaterAmplitude(q);
    }
    return atom.occupancy * result;
}

AtomicSAXSCloud::AmplitudeMoment AtomicSAXSCloud::amplitudeMoment(
    const AtomicSAXSScatterer& atom) const
{
    const bool spatial = options.hydrationModel == AtomicSAXSHydrationModel::SpatialShell;
    if (atom.hydrationVolume > 0.0 && (!spatial || !options.includeHydration))
        return AmplitudeMoment{0.0, 0.0};
    const CromerMannCoefficients hydrogen = hydrogenCoefficients();
    double vacuumBeta = atom.vacuum.valueAtZero() +
        static_cast<double>(atom.implicitHydrogens) * hydrogen.valueAtZero();
    double moment = atom.vacuum.secondMoment() +
        static_cast<double>(atom.implicitHydrogens) *
            (hydrogen.secondMoment() + hydrogen.valueAtZero() *
             atom.implicitHydrogenDistance * atom.implicitHydrogenDistance);
    if (options.useTemperatureFactors && atom.temperatureFactor > 0.0) {
        moment += 3.0 * vacuumBeta * atom.temperatureFactor /
            (8.0 * kPi * kPi);
    }
    double beta = vacuumBeta;

    if (options.includeExcludedSolvent && atom.excludedVolume > 0.0) {
        const double volume = options.excludedVolumeScale * atom.excludedVolume;
        const double excludedBeta = -options.solventElectronDensity * volume;
        beta += excludedBeta;
        moment += excludedBeta * 3.0 * std::pow(volume, 2.0 / 3.0) /
            (2.0 * kPi);
    }

    if (options.includeHydration &&
        ((spatial && atom.hydrationVolume > 0.0) ||
         (!spatial && atom.solventAccessibility > 0.0))) {
        const AmplitudeMoment water = waterMoment();
        const double volume = spatial ? atom.hydrationVolume :
            options.hydrationVolumePerExposedAtom * atom.solventAccessibility;
        const double density = atom.hydrationUsesSolventDensity ? options.solventElectronDensity : 1.0;
        const double hydrationBeta = options.hydrationShellContrast * volume * density;
        beta += hydrationBeta;
        moment += hydrationBeta * water.secondMoment / water.beta;
    }

    return AmplitudeMoment{
        atom.occupancy * beta,
        atom.occupancy * moment
    };
}

AtomicSAXSCloud::AmplitudeMoment AtomicSAXSCloud::waterMoment() const
{
    const CromerMannCoefficients oxygen(
        {{3.0485, 2.2868, 1.5463, 0.8670}},
        {{13.2771, 5.7011, 0.3239, 32.9089}}, 0.2508);
    const auto hydrogen = hydrogenCoefficients();
    double beta = oxygen.valueAtZero() + 2*hydrogen.valueAtZero();
    double moment = oxygen.secondMoment() + 2*hydrogen.secondMoment();
    if (options.hydrationWaterModel == AtomicSAXSHydrationWaterModel::BondedWater) {
        // O-H distance: Grudinin et al. (2017), Table 1 / equation (8).
        moment += 2*hydrogen.valueAtZero()*0.967*0.967;
    } else if (options.hydrationWaterModel == AtomicSAXSHydrationWaterModel::FoXSReciprocal ||
               options.hydrationWaterModel == AtomicSAXSHydrationWaterModel::FoXSFast) {
        // FoXS water is a solution form factor: O+2H minus their displaced
        // solvent, not a neutral-water form factor. Native units are c2*fw.
        const double vo=9.13, vh=5.15;
        beta -= options.solventElectronDensity*(vo+2*vh);
        moment -= options.solventElectronDensity*3/(8*kPi)*
            (vo*std::pow(vo, 2.0/3.0)+2*vh*std::pow(vh, 2.0/3.0));
        if (options.hydrationWaterModel == AtomicSAXSHydrationWaterModel::FoXSFast)
            moment = 0.69*beta; // amplitude exp(-0.115*q^2)
    }
    if (std::abs(beta) < 1e-12) throw SEBException(
        "Cannot normalize a zero-contrast hydration water form factor",
        "AtomicSAXSCloud water model");
    return AmplitudeMoment{beta, moment};
}

double AtomicSAXSCloud::normalizedWaterAmplitude(double q) const
{
    const AmplitudeMoment water = waterMoment();
    if (options.hydrationWaterModel == AtomicSAXSHydrationWaterModel::FoXSFast)
        return std::exp(-0.115*q*q);
    const CromerMannCoefficients oxygen(
        {{3.0485, 2.2868, 1.5463, 0.8670}},
        {{13.2771, 5.7011, 0.3239, 32.9089}}, 0.2508);
    const auto hydrogen = hydrogenCoefficients();
    const double bond = options.hydrationWaterModel == AtomicSAXSHydrationWaterModel::BondedWater
        ? sinc(q*0.967) : 1.0;
    double value = oxygen.evaluate(q) + 2*hydrogen.evaluate(q)*bond;
    if (options.hydrationWaterModel == AtomicSAXSHydrationWaterModel::FoXSReciprocal) {
        const double vo=9.13, vh=5.15;
        value -= options.solventElectronDensity*(
            vo*std::exp(-q*q*std::pow(vo,2.0/3.0)/(16*kPi)) +
            2*vh*std::exp(-q*q*std::pow(vh,2.0/3.0)/(16*kPi)));
    }
    return value/water.beta;
}

double AtomicSAXSCloud::totalBeta() const
{
    ensureGeometryPrepared();
    double beta = 0.0;
    for (const auto& atom : atoms) {
        beta += amplitudeMoment(atom).beta;
    }
    return beta;
}

double AtomicSAXSCloud::nonzeroTotalBeta(const std::string& operation) const
{
    const double beta = totalBeta();
    double scale = 1.0;
    for (const auto& atom : atoms) {
        scale += std::abs(amplitudeMoment(atom).beta);
    }
    if (std::abs(beta) <=
        std::numeric_limits<double>::epsilon() * scale * 16.0) {
        throw SEBException(
            "Total solution contrast is zero for " + operation,
            "AtomicSAXSCloud"
        );
    }
    return beta;
}

double AtomicSAXSCloud::formFactorUnnormalized(double q) const
{
    ensureGeometryPrepared();
    std::vector<double> amplitudes;
    amplitudes.reserve(atoms.size());
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        amplitudes.push_back(atomicAmplitude(index, q));
    }

    double result = 0.0;
    // Spatial quadrature adds many small terms, often with cancellation
    // against displaced solvent. Keep low-q derivatives stable as well.
    double compensation = 0.0;
    const auto accumulate = [&](double term) {
        const double corrected = term - compensation;
        const double updated = result + corrected;
        compensation = (updated - result) - corrected;
        result = updated;
    };
    for (double amplitude : amplitudes) {
        accumulate(amplitude * amplitude);
    }
    for (const auto& pair : atomPairDistances) {
        accumulate(2.0 * amplitudes[pair.first] * amplitudes[pair.second] *
            sinc(q * pair.distance));
    }
    return result;
}

double AtomicSAXSCloud::amplitudeUnnormalized(
    const refPoint& reference,
    double q) const
{
    ensureGeometryPrepared();
    const auto distances = referenceAtomDistances.find(reference);
    if (distances == referenceAtomDistances.end()) {
        throw SEBException(
            "Unknown reference point " + reference,
            "AtomicSAXSCloud::amplitudeUnnormalized()"
        );
    }

    double result = 0.0;
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        result += atomicAmplitude(index, q) *
            sinc(q * distances->second[index]);
    }
    return result;
}

double AtomicSAXSCloud::phaseFactor(
    const refPoint& first,
    const refPoint& second,
    double q) const
{
    if (first == second) {
        return 1.0;
    }
    ensureGeometryPrepared();
    const auto pair = referencePair(first, second);
    const auto distance = referencePairDistances.find(pair);
    if (distance == referencePairDistances.end()) {
        throw SEBException(
            "Unknown reference-point pair " + first + ", " + second,
            "AtomicSAXSCloud::phaseFactor()"
        );
    }
    return sinc(q * distance->second);
}

double AtomicSAXSCloud::radiusOfGyration2() const
{
    ensureGeometryPrepared();
    const double beta = nonzeroTotalBeta("radius of gyration");
    double internal = 0.0;
    // Moments depend on the atom and current solvent options, not on its pair.
    // Reuse them within this evaluation without retaining stale mutable state.
    std::vector<double> atomBetas;
    atomBetas.reserve(atoms.size());
    for (const auto& atom : atoms) {
        const AmplitudeMoment moment = amplitudeMoment(atom);
        internal += moment.secondMoment;
        atomBetas.push_back(moment.beta);
    }
    internal /= beta;

    double centers = 0.0;
    for (const auto& pair : atomPairDistances) {
        centers += atomBetas[pair.first] * atomBetas[pair.second] * pair.distance2;
    }
    centers /= beta * beta;
    return internal + centers;
}

double AtomicSAXSCloud::sigmaReferenceToScatterers(
    const refPoint& reference) const
{
    ensureGeometryPrepared();
    const double beta = nonzeroTotalBeta(
        "reference-to-scatterer mean-square distance"
    );
    const auto distances2 = referenceAtomDistances2.find(reference);
    if (distances2 == referenceAtomDistances2.end()) {
        throw SEBException(
            "Unknown reference point " + reference,
            "AtomicSAXSCloud::sigmaReferenceToScatterers()"
        );
    }

    double result = 0.0;
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        const AmplitudeMoment moment = amplitudeMoment(atoms[index]);
        result += moment.secondMoment + moment.beta * distances2->second[index];
    }
    return result / beta;
}

double AtomicSAXSCloud::sigmaReferenceToReference(
    const refPoint& first,
    const refPoint& second) const
{
    if (first == second) {
        return 0.0;
    }
    ensureGeometryPrepared();
    const auto pair = referencePair(first, second);
    const auto distance2 = referencePairDistances2.find(pair);
    if (distance2 == referencePairDistances2.end()) {
        throw SEBException(
            "Unknown reference-point pair " + first + ", " + second,
            "AtomicSAXSCloud::sigmaReferenceToReference()"
        );
    }
    return distance2->second;
}
