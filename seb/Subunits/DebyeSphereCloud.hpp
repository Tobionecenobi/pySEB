#ifndef INCLUDE_GUARD_DEBYESPHERECLOUD
#define INCLUDE_GUARD_DEBYESPHERECLOUD

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Numerical.hpp"

struct CartesianPoint3D {
    double x;
    double y;
    double z;

    CartesianPoint3D(double xValue = 0.0,
                     double yValue = 0.0,
                     double zValue = 0.0)
        : x(xValue), y(yValue), z(zValue) {}
};

struct SphereScatterer {
    CartesianPoint3D center;
    double radius;
    double beta;

    SphereScatterer(
        const CartesianPoint3D& centerValue = CartesianPoint3D(),
        double radiusValue = 0.0,
        double betaValue = 0.0)
        : center(centerValue), radius(radiusValue), beta(betaValue) {}

    SphereScatterer(
        double x,
        double y,
        double z,
        double radiusValue,
        double betaValue)
        : center(x, y, z), radius(radiusValue), beta(betaValue) {}
};

class DebyeSphereCloud : public NumericalSubunit {
public:
    DebyeSphereCloud()
        : NumericalSubunit(NormalizationMode::Unnormalized)
    {
        stype = DEBYESPHERECLOUD;
        setReferencePointName("center");
        referenceCoordinates["center"] = CartesianPoint3D();
        configureCoreCallbacks();
        configureReferenceCallbacks("center");
    }

    explicit DebyeSphereCloud(const std::vector<SphereScatterer>& scatterers)
        : DebyeSphereCloud()
    {
        particles = scatterers;
        dirty = true;
    }

    virtual ~DebyeSphereCloud() {}

    DebyeSphereCloud(const DebyeSphereCloud&) = delete;
    DebyeSphereCloud& operator=(const DebyeSphereCloud&) = delete;

    void addScatterer(const SphereScatterer& scatterer) {
        particles.push_back(scatterer);
        dirty = true;
    }

    void addScatterer(
        double x,
        double y,
        double z,
        double radius,
        double beta)
    {
        addScatterer(SphereScatterer(x, y, z, radius, beta));
    }

    void addReferencePoint(
        const refPoint& name,
        const CartesianPoint3D& coordinate)
    {
        if (hasReference(name)) {
            throw SEBException(
                "Reference point " + name + " already exists",
                "DebyeSphereCloud::addReferencePoint()"
            );
        }
        validatePoint(coordinate, "reference point " + name);
        setReferencePointName(name);
        referenceCoordinates[name] = coordinate;
        configureReferenceCallbacks(name);
        dirty = true;
    }

    void addReferencePoint(
        const refPoint& name,
        double x,
        double y,
        double z)
    {
        addReferencePoint(name, CartesianPoint3D(x, y, z));
    }

    const std::vector<SphereScatterer>& getScatterers() const {
        return particles;
    }

    std::map<refPoint, CartesianPoint3D> getReferenceCoordinates() const {
        ensurePrepared();
        return referenceCoordinates;
    }

    std::size_t scattererCount() const {
        return particles.size();
    }

    static double sinc(double x) {
        const double x2 = x * x;
        if (x2 < 1e-8) {
            return 1.0 - x2 / 6.0 + x2 * x2 / 120.0;
        }
        return std::sin(x) / x;
    }

    static double sphereAmplitude(double x) {
        const double x2 = x * x;
        if (x2 < 1e-6) {
            return 1.0 - x2 / 10.0 + x2 * x2 / 280.0;
        }
        return 3.0 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
    }

private:
    struct PairDistance {
        std::size_t first;
        std::size_t second;
        double distance;
        double distance2;
    };

    std::vector<SphereScatterer> particles;
    mutable std::map<refPoint, CartesianPoint3D> referenceCoordinates;
    mutable std::vector<PairDistance> scattererPairDistances;
    mutable std::map<refPoint, std::vector<double>> referenceScattererDistances;
    mutable std::map<refPoint, std::vector<double>> referenceScattererDistances2;
    mutable std::map<std::pair<refPoint, refPoint>, double> referencePairDistances;
    mutable std::map<std::pair<refPoint, refPoint>, double> referencePairDistances2;
    mutable bool dirty = true;

    void configureCoreCallbacks() {
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

    void configureReferenceCallbacks(const refPoint& name) {
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

    static void validatePoint(
        const CartesianPoint3D& point,
        const string& description)
    {
        if (!std::isfinite(point.x) ||
            !std::isfinite(point.y) ||
            !std::isfinite(point.z)) {
            throw SEBException(
                "Non-finite coordinate in " + description,
                "DebyeSphereCloud input validation"
            );
        }
    }

    void validateParticles() const {
        if (particles.empty()) {
            throw SEBException(
                "A Debye sphere cloud must contain at least one scatterer",
                "DebyeSphereCloud input validation"
            );
        }

        for (const auto& particle : particles) {
            validatePoint(particle.center, "sphere center");
            if (!std::isfinite(particle.radius) || particle.radius < 0.0) {
                throw SEBException(
                    "Sphere radii must be finite and non-negative",
                    "DebyeSphereCloud input validation"
                );
            }
            if (!std::isfinite(particle.beta)) {
                throw SEBException(
                    "Sphere excess scattering lengths must be finite",
                    "DebyeSphereCloud input validation"
                );
            }
        }
    }

    static double squaredDistance(
        const CartesianPoint3D& left,
        const CartesianPoint3D& right)
    {
        const double dx = left.x - right.x;
        const double dy = left.y - right.y;
        const double dz = left.z - right.z;
        return dx * dx + dy * dy + dz * dz;
    }

    void ensurePrepared() const {
        if (!dirty) {
            return;
        }

        validateParticles();

        CartesianPoint3D centroid;
        for (const auto& particle : particles) {
            centroid.x += particle.center.x;
            centroid.y += particle.center.y;
            centroid.z += particle.center.z;
        }
        const double inverseCount = 1.0 / static_cast<double>(particles.size());
        centroid.x *= inverseCount;
        centroid.y *= inverseCount;
        centroid.z *= inverseCount;
        referenceCoordinates["center"] = centroid;

        scattererPairDistances.clear();
        for (std::size_t first = 0; first < particles.size(); ++first) {
            for (std::size_t second = first + 1; second < particles.size(); ++second) {
                const double distance2 = squaredDistance(
                    particles[first].center,
                    particles[second].center
                );
                scattererPairDistances.push_back(
                    PairDistance{
                        first,
                        second,
                        std::sqrt(distance2),
                        distance2
                    }
                );
            }
        }

        referenceScattererDistances.clear();
        referenceScattererDistances2.clear();
        for (const auto& reference : referenceCoordinates) {
            std::vector<double> distances;
            std::vector<double> distances2;
            distances.reserve(particles.size());
            distances2.reserve(particles.size());
            for (const auto& particle : particles) {
                const double distance2 = squaredDistance(
                    reference.second,
                    particle.center
                );
                distances.push_back(std::sqrt(distance2));
                distances2.push_back(distance2);
            }
            referenceScattererDistances[reference.first] = distances;
            referenceScattererDistances2[reference.first] = distances2;
        }

        referencePairDistances.clear();
        referencePairDistances2.clear();
        for (auto first = referenceCoordinates.begin();
             first != referenceCoordinates.end();
             ++first) {
            auto second = first;
            ++second;
            for (; second != referenceCoordinates.end(); ++second) {
                const auto pair = referencePair(first->first, second->first);
                const double distance2 = squaredDistance(first->second, second->second);
                referencePairDistances[pair] = std::sqrt(distance2);
                referencePairDistances2[pair] = distance2;
            }
        }

        dirty = false;
    }

    double totalBeta() const {
        ensurePrepared();
        double beta = 0.0;
        for (const auto& particle : particles) {
            beta += particle.beta;
        }
        return beta;
    }

    double nonzeroTotalBeta(const string& operation) const {
        const double beta = totalBeta();
        double scale = 1.0;
        for (const auto& particle : particles) {
            scale += std::abs(particle.beta);
        }
        if (std::abs(beta) <= std::numeric_limits<double>::epsilon() * scale * 16.0) {
            throw SEBException(
                "Total excess scattering length is zero for " + operation,
                "DebyeSphereCloud"
            );
        }
        return beta;
    }

    std::vector<double> sphereAmplitudes(double q) const {
        std::vector<double> amplitudes;
        amplitudes.reserve(particles.size());
        for (const auto& particle : particles) {
            amplitudes.push_back(sphereAmplitude(q * particle.radius));
        }
        return amplitudes;
    }

    double formFactorUnnormalized(double q) const {
        ensurePrepared();
        const std::vector<double> amplitudes = sphereAmplitudes(q);

        double result = 0.0;
        for (std::size_t index = 0; index < particles.size(); ++index) {
            const double weightedAmplitude =
                particles[index].beta * amplitudes[index];
            result += weightedAmplitude * weightedAmplitude;
        }

        for (const auto& pair : scattererPairDistances) {
            result += 2.0 *
                particles[pair.first].beta *
                particles[pair.second].beta *
                amplitudes[pair.first] *
                amplitudes[pair.second] *
                sinc(q * pair.distance);
        }
        return result;
    }

    double amplitudeUnnormalized(const refPoint& reference, double q) const {
        ensurePrepared();
        const auto distances = referenceScattererDistances.find(reference);
        if (distances == referenceScattererDistances.end()) {
            throw SEBException(
                "Unknown reference point " + reference,
                "DebyeSphereCloud::amplitudeUnnormalized()"
            );
        }

        const std::vector<double> amplitudes = sphereAmplitudes(q);
        double result = 0.0;
        for (std::size_t index = 0; index < particles.size(); ++index) {
            result += particles[index].beta *
                amplitudes[index] *
                sinc(q * distances->second[index]);
        }
        return result;
    }

    double phaseFactor(
        const refPoint& first,
        const refPoint& second,
        double q) const
    {
        if (first == second) {
            return 1.0;
        }
        ensurePrepared();
        const auto pair = referencePair(first, second);
        const auto distance = referencePairDistances.find(pair);
        if (distance == referencePairDistances.end()) {
            throw SEBException(
                "Unknown reference-point pair " + first + ", " + second,
                "DebyeSphereCloud::phaseFactor()"
            );
        }
        return sinc(q * distance->second);
    }

    double radiusOfGyration2() const {
        ensurePrepared();
        const double beta = nonzeroTotalBeta("radius of gyration");

        double internal = 0.0;
        for (const auto& particle : particles) {
            internal += particle.beta * particle.radius * particle.radius;
        }
        internal *= 3.0 / (5.0 * beta);

        double centers = 0.0;
        for (const auto& pair : scattererPairDistances) {
            centers += particles[pair.first].beta *
                particles[pair.second].beta *
                pair.distance2;
        }
        centers /= beta * beta;
        return internal + centers;
    }

    double sigmaReferenceToScatterers(const refPoint& reference) const {
        ensurePrepared();
        const double beta = nonzeroTotalBeta(
            "reference-to-scatterer mean-square distance"
        );
        const auto distances2 = referenceScattererDistances2.find(reference);
        if (distances2 == referenceScattererDistances2.end()) {
            throw SEBException(
                "Unknown reference point " + reference,
                "DebyeSphereCloud::sigmaReferenceToScatterers()"
            );
        }

        double result = 0.0;
        for (std::size_t index = 0; index < particles.size(); ++index) {
            result += particles[index].beta *
                (3.0 * particles[index].radius * particles[index].radius / 5.0 +
                 distances2->second[index]);
        }
        return result / beta;
    }

    double sigmaReferenceToReference(
        const refPoint& first,
        const refPoint& second) const
    {
        if (first == second) {
            return 0.0;
        }
        ensurePrepared();
        const auto pair = referencePair(first, second);
        const auto distance2 = referencePairDistances2.find(pair);
        if (distance2 == referencePairDistances2.end()) {
            throw SEBException(
                "Unknown reference-point pair " + first + ", " + second,
                "DebyeSphereCloud::sigmaReferenceToReference()"
            );
        }
        return distance2->second;
    }
};

#endif
