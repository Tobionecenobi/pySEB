// Geometry realization, reference sampling, and pose mathematics.
#include "USDExportInternal.hpp"
#include "Exceptions.hpp"
#include "SubunitIO/ExpressionParser.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace pyseb {
namespace usd_detail {

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

double length(const Vec3& value) {
    return std::sqrt(dot(value, value));
}

Vec3 rotate(const Quaternion& q, const Vec3& value) {
    const Vec3 u{q.x, q.y, q.z};
    return value
        + cross(u, value) * (2.0 * q.w)
        + cross(u, cross(u, value)) * 2.0;
}

Quaternion multiply(const Quaternion& first, const Quaternion& second) {
    return {
        first.w * second.w - first.x * second.x
            - first.y * second.y - first.z * second.z,
        first.w * second.x + first.x * second.w
            + first.y * second.z - first.z * second.y,
        first.w * second.y - first.x * second.z
            + first.y * second.w + first.z * second.x,
        first.w * second.z + first.x * second.y
            - first.y * second.x + first.z * second.w
    };
}

unsigned long long hashString(unsigned long long hash, const std::string& value) {
    for (unsigned char ch : value) hash = (hash ^ ch) * 1099511628211ULL;
    return hash;
}

unsigned long long derivedSeed(unsigned long long seed, const std::string& key) {
    return hashString(seed ^ 1469598103934665603ULL, key);
}

double unitRandom(std::mt19937_64& generator) {
    return std::generate_canonical<double, 53>(generator);
}

Quaternion uniformRotation(unsigned long long seed, const std::string& key) {
    std::mt19937_64 generator(derivedSeed(seed, key));
    const double u1 = unitRandom(generator);
    const double u2 = unitRandom(generator);
    const double u3 = unitRandom(generator);
    const double a = 2.0 * M_PI * u2;
    const double b = 2.0 * M_PI * u3;
    Quaternion result;
    result.x = std::sqrt(1.0-u1) * std::sin(a);
    result.y = std::sqrt(1.0-u1) * std::cos(a);
    result.z = std::sqrt(u1) * std::sin(b);
    result.w = std::sqrt(u1) * std::cos(b);
    return result;
}

VisualizationReferenceName parseVisualizationReference(const std::string& reference) {
    VisualizationReferenceName result;
    const std::size_t first=reference.find('#');
    if (first == std::string::npos) {
        result.base=reference;
        return result;
    }
    result.base=reference.substr(0,first);
    const std::size_t second=reference.find('#',first+1);
    if (second == std::string::npos) {
        result.label=reference.substr(first+1);
        return result;
    }
    if (reference.find('#',second+1) != std::string::npos) throw SEBException("visualization reference contains too many '#' separators", "World::ExportUSD");
    result.variant=reference.substr(first+1,second-first-1);
    result.label=reference.substr(second+1);
    if (result.variant.empty() || result.label.empty()) throw SEBException("visualization variant and label must be non-empty", "World::ExportUSD");
    return result;
}

Vec3 interpolate(const Vec3& a, const Vec3& b, double fraction) {
    return a * (1.0-fraction) + b * fraction;
}

Vec3 sampleCurveParameter(const GeometryPatch& patch, double parameter) {
    if (patch.points.empty()) throw SEBException("cannot sample an empty curve", "World::ExportUSD");
    if (patch.points.size() == 1) return patch.points.front();
    parameter = std::max(0.0, std::min(1.0, parameter));
    const double index = parameter * static_cast<double>(patch.points.size()-1);
    const std::size_t first = static_cast<std::size_t>(std::floor(index));
    const std::size_t second = std::min(first + 1, patch.points.size()-1);
    return interpolate(patch.points[first], patch.points[second], index-first);
}

Vec3 sampleCurveLength(const GeometryPatch& patch, double draw) {
    if (patch.points.size() < 2) return sampleCurveParameter(patch, draw);
    std::vector<double> lengths;
    double total = 0.0;
    for (std::size_t i=1; i<patch.points.size(); ++i) {
        total += length(patch.points[i]-patch.points[i-1]);
        lengths.push_back(total);
    }
    if (total <= 0.0) return patch.points.front();
    const double target = draw * total;
    const auto found = std::lower_bound(lengths.begin(), lengths.end(), target);
    const std::size_t segment = static_cast<std::size_t>(found-lengths.begin());
    const double before = segment == 0 ? 0.0 : lengths[segment-1];
    const double span = lengths[segment]-before;
    return interpolate(patch.points[segment], patch.points[segment+1], span > 0.0 ? (target-before)/span : 0.0);
}

Vec3 sampleSurfaceArea(const GeometryPatch& patch, std::mt19937_64& generator) {
    if (patch.uCount < 2 || patch.vCount < 2) throw SEBException("surface requires at least two samples per dimension", "World::ExportUSD");
    struct Triangle { Vec3 a,b,c; double cumulative; };
    std::vector<Triangle> triangles;
    double total = 0.0;
    for (std::size_t u=0; u+1<patch.uCount; ++u) {
        for (std::size_t v=0; v+1<patch.vCount; ++v) {
            const std::size_t i=u*patch.vCount+v;
            const Vec3 p00=patch.points[i], p01=patch.points[i+1];
            const Vec3 p10=patch.points[i+patch.vCount], p11=patch.points[i+patch.vCount+1];
            for (const auto& vertices : {std::array<Vec3,3>{{p00,p01,p11}}, std::array<Vec3,3>{{p00,p11,p10}}}) {
                const double area=0.5*length(cross(vertices[1]-vertices[0],vertices[2]-vertices[0]));
                if (area > 0.0) { total += area; triangles.push_back({vertices[0],vertices[1],vertices[2],total}); }
            }
        }
    }
    if (triangles.empty() || total <= 0.0) throw SEBException("surface has zero sampled area", "World::ExportUSD");
    const double target=unitRandom(generator)*total;
    const auto found=std::lower_bound(triangles.begin(),triangles.end(),target,[](const Triangle& t,double value){return t.cumulative<value;});
    double r1=unitRandom(generator), r2=unitRandom(generator);
    const double root=std::sqrt(r1);
    return found->a*(1.0-root)+found->b*(root*(1.0-r2))+found->c*(root*r2);
}

Vec3 sampleSurfaceBoundary(const GeometryPatch& patch, double draw) {
    if (patch.uCount == 0 || patch.vCount == 0) throw SEBException("cannot sample an empty surface", "World::ExportUSD");
    GeometryPatch boundary;
    boundary.kind=pyseb::VisualizationGeometryKind::Curve;
    const std::size_t offset=(patch.uCount-1)*patch.vCount;
    boundary.points.assign(patch.points.begin()+offset,patch.points.begin()+offset+patch.vCount);
    return sampleCurveLength(boundary,draw);
}

GeometryPatch realizeGeometry(
    const pyseb::VisualizationGeometry& definition,
    std::size_t sampleCount,
    const ExpressionEvaluator& evaluate,
    unsigned long long seed,
    const std::string& instanceName,
    const std::string& geometryName) {
    GeometryPatch patch;
    patch.kind = definition.kind;

    if (definition.kind == pyseb::VisualizationGeometryKind::RandomWalk) {
        std::mt19937_64 generator(
            derivedSeed(seed, instanceName + "/geometry/" + geometryName));
        std::normal_distribution<double> normal(0.0, 1.0);
        patch.points.resize(sampleCount);

        for (std::size_t index = 1; index < sampleCount; ++index) {
            Vec3 step{normal(generator), normal(generator), normal(generator)};
            if (definition.distribution == "fixed_length") {
                const double size = length(step);
                if (size > 0.0) step = step * (1.0 / size);
            }
            patch.points[index] = patch.points[index - 1] + step;
        }

        if (definition.closure == "bridge") {
            const Vec3 end = patch.points.back();
            for (std::size_t index = 0; index < sampleCount; ++index) {
                patch.points[index] = patch.points[index]
                    - end * (static_cast<double>(index) / (sampleCount - 1));
            }
        }

        Vec3 center;
        for (const Vec3& point : patch.points) center = center + point;
        center = center * (1.0 / sampleCount);
        double radiusOfGyrationSquared = 0.0;
        for (Vec3& point : patch.points) {
            point = point - center;
            radiusOfGyrationSquared += dot(point, point);
        }

        const double radiusOfGyration = std::sqrt(
            radiusOfGyrationSquared / sampleCount);
        const double target = definition.targetRg.root()
            ? evaluate(definition.targetRg, {})
            : 1.0;
        if (radiusOfGyration > 0.0) {
            for (Vec3& point : patch.points) {
                point = point * (target / radiusOfGyration);
            }
        }
        return patch;
    }

    if (definition.kind == pyseb::VisualizationGeometryKind::Curve) {
        patch.uCount = sampleCount;
        const double lower = evaluate(definition.uLower, {});
        const double upper = evaluate(definition.uUpper, {});
        for (std::size_t index = 0; index < sampleCount; ++index) {
            const double u = lower + (upper - lower)
                * static_cast<double>(index) / (sampleCount - 1);
            const std::map<std::string, double> local{{"u", u}, {"t", u}};
            patch.points.push_back({
                evaluate(definition.coordinates[0], local),
                evaluate(definition.coordinates[1], local),
                evaluate(definition.coordinates[2], local)
            });
        }
        return patch;
    }

    patch.uCount = sampleCount;
    patch.vCount = sampleCount;
    const double uLower = evaluate(definition.uLower, {});
    const double uUpper = evaluate(definition.uUpper, {});
    const double vLower = evaluate(definition.vLower, {});
    const double vUpper = evaluate(definition.vUpper, {});
    for (std::size_t uIndex = 0; uIndex < sampleCount; ++uIndex) {
        const double u = uLower + (uUpper - uLower)
            * static_cast<double>(uIndex) / (sampleCount - 1);
        for (std::size_t vIndex = 0; vIndex < sampleCount; ++vIndex) {
            const double v = vLower + (vUpper - vLower)
                * static_cast<double>(vIndex) / (sampleCount - 1);
            const std::map<std::string, double> local{{"u", u}, {"v", v}};
            patch.points.push_back({
                evaluate(definition.coordinates[0], local),
                evaluate(definition.coordinates[1], local),
                evaluate(definition.coordinates[2], local)
            });
        }
    }
    return patch;
}

void resolveSpecificReferences(
    InstanceScene& instance,
    const pyseb::VisualizationDefinition& visualization,
    const ExpressionEvaluator& evaluate) {
    for (const auto& item : visualization.references) {
        const auto& reference = item.second;
        if (reference.kind == "fixed") {
            instance.resolvedReferences[item.first] = {
                evaluate(reference.position[0], {}),
                evaluate(reference.position[1], {}),
                evaluate(reference.position[2], {})
            };
            continue;
        }
        if (reference.kind == "curve_fraction") {
            const auto patch = instance.patches.find(reference.geometry);
            if (patch == instance.patches.end()) {
                throw SEBException(
                    "reference uses unknown geometry",
                    "World::ExportUSD");
            }
            instance.resolvedReferences[item.first] = sampleCurveParameter(
                patch->second,
                evaluate(reference.fraction, {}));
        }
    }
}

double atomicDisplayedRadius(
    const AtomicSAXSScatterer& scatterer,
    double zeroRadiusMarkerSize)
{
    if (scatterer.excludedVolume <= 0.0) {
        return zeroRadiusMarkerSize;
    }
    constexpr double piValue = 3.141592653589793238462643383279502884;
    return std::cbrt(3.0 * scatterer.excludedVolume / (4.0 * piValue));
}

} // namespace usd_detail
} // namespace pyseb
