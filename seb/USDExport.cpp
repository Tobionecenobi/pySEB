#include "USDExport.hpp"
#include "World.hpp"
#include "Subunits/FileDefined.hpp"
#include "Subunits/DebyeSphereCloud.hpp"
#include "Exceptions.hpp"
#include "SubunitIO/ExpressionParser.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <locale>
#include <map>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <vector>

double MetersPerUnit(LengthUnit unit) {
    switch (unit) {
        case LengthUnit::Meter: return 1.0;
        case LengthUnit::Millimeter: return 1e-3;
        case LengthUnit::Micrometer: return 1e-6;
        case LengthUnit::Nanometer: return 1e-9;
        case LengthUnit::Angstrom: return 1e-10;
    }
    return 1.0;
}

USDExportOptions::USDExportOptions(LengthUnit unitValue)
    : unit(unitValue), metersPerUnit(MetersPerUnit(unitValue)), customUnit(false) {}

USDExportOptions::USDExportOptions(double customMetersPerUnit)
    : unit(LengthUnit::Meter), metersPerUnit(customMetersPerUnit), customUnit(true) {}

namespace {

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
    Vec3 operator+(const Vec3& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }
    Vec3 operator-(const Vec3& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }
    Vec3 operator*(double scale) const {
        return {x * scale, y * scale, z * scale};
    }
};

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

struct Quaternion {
    double w = 1.0, x = 0.0, y = 0.0, z = 0.0;
};

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

struct Pose {
    Quaternion rotation;
    Vec3 translation;
    bool resolved = false;
};

struct GeometryPatch {
    pyseb::VisualizationGeometryKind kind = pyseb::VisualizationGeometryKind::Curve;
    std::vector<Vec3> points;
    std::size_t uCount = 0;
    std::size_t vCount = 0;
};

struct InstanceScene {
    const FileDefinedSubunit* subunit = nullptr;
    std::string name;
    std::string tag;
    std::map<std::string, double> parameters;
    std::map<std::string, GeometryPatch> patches;
    std::map<std::string, Vec3> resolvedReferences;
    std::array<double,3> color{{0.7,0.7,0.7}};
    double opacity = 1.0;
};

using ExpressionEvaluator = std::function<double(
    const pyseb::ParsedExpression&,
    const std::map<std::string, double>&)>;

// Scene realization: turn schema descriptions into sampled, local geometry.
struct LinkEndpoint {
    std::string instance;
    std::string reference;
};

struct SceneLink {
    LinkEndpoint first;
    LinkEndpoint second;
    std::string key;
};

struct RootedLink {
    LinkEndpoint parent;
    LinkEndpoint child;
    std::string key;
};

struct ProxySegment {
    Vec3 first;
    Vec3 second;
};

struct GeometryProxy {
    std::vector<Vec3> points;
    std::vector<ProxySegment> segments;
    Vec3 centroid;
    double radius = 0.0;
};

std::string quote(const std::string& value) {
    std::string result = "\"";
    for (char ch : value) {
        if (ch == '\\' || ch == '\"') result += '\\';
        result += ch;
    }
    return result + "\"";
}

std::string number(double value) {
    if (!std::isfinite(value)) {
        throw SEBException("non-finite visualization value", "World::ExportUSD");
    }
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << std::setprecision(17) << value;
    return stream.str();
}

std::string usdName(const std::string& value) {
    std::string result;
    for (char ch : value) {
        result += (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') ? ch : '_';
    }
    if (result.empty() || std::isdigit(static_cast<unsigned char>(result.front()))) result = "p_" + result;
    return result;
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

bool parseEndpoint(const std::string& value, LinkEndpoint& endpoint) {
    const std::size_t dot = value.find_last_of('.');
    if (dot == std::string::npos || dot + 1 == value.size()) return false;
    std::string path = value.substr(0, dot);
    const std::size_t colon = path.find_last_of(':');
    endpoint.instance = colon == std::string::npos ? path : path.substr(colon + 1);
    endpoint.reference = value.substr(dot + 1);
    return true;
}

struct VisualizationReferenceName {
    std::string base;
    std::string variant;
    std::string label;
};

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

double segmentDistanceSquared(const ProxySegment& first, const ProxySegment& second) {
    const Vec3 u=first.second-first.first;
    const Vec3 v=second.second-second.first;
    const Vec3 w=first.first-second.first;
    const double a=dot(u,u), b=dot(u,v), c=dot(v,v), d=dot(u,w), e=dot(v,w);
    const double denominator=a*c-b*b;
    double numeratorS, denominatorS=denominator;
    double numeratorT, denominatorT=denominator;
    const double epsilon=1e-15;
    if (denominator < epsilon) {
        numeratorS=0.0;
        denominatorS=1.0;
        numeratorT=e;
        denominatorT=c;
    } else {
        numeratorS=b*e-c*d;
        numeratorT=a*e-b*d;
        if (numeratorS < 0.0) {
            numeratorS=0.0;
            numeratorT=e;
            denominatorT=c;
        } else if (numeratorS > denominatorS) {
            numeratorS=denominatorS;
            numeratorT=e+b;
            denominatorT=c;
        }
    }
    if (numeratorT < 0.0) {
        numeratorT=0.0;
        if (-d < 0.0) numeratorS=0.0;
        else if (-d > a) numeratorS=denominatorS;
        else { numeratorS=-d; denominatorS=a; }
    } else if (numeratorT > denominatorT) {
        numeratorT=denominatorT;
        if (-d+b < 0.0) numeratorS=0.0;
        else if (-d+b > a) numeratorS=denominatorS;
        else { numeratorS=-d+b; denominatorS=a; }
    }
    const double s=std::abs(numeratorS)<epsilon ? 0.0 : numeratorS/denominatorS;
    const double t=std::abs(numeratorT)<epsilon ? 0.0 : numeratorT/denominatorT;
    const Vec3 separation=w+u*s-v*t;
    return dot(separation,separation);
}

GeometryProxy localProxy(const InstanceScene& instance) {
    GeometryProxy proxy;
    constexpr std::size_t maxPatchPoints=32;
    constexpr std::size_t maxPatchSegments=32;
    for (const auto& item : instance.patches) {
        const GeometryPatch& patch=item.second;
        if (patch.points.empty()) continue;
        const std::size_t pointCount=std::min(maxPatchPoints,patch.points.size());
        for (std::size_t i=0; i<pointCount; ++i) {
            const std::size_t index=pointCount==1 ? 0 : i*(patch.points.size()-1)/(pointCount-1);
            proxy.points.push_back(patch.points[index]);
        }
        if (patch.kind != pyseb::VisualizationGeometryKind::Surface && patch.points.size()>1) {
            const std::size_t segmentCount=std::min(maxPatchSegments,patch.points.size()-1);
            for (std::size_t i=0; i<segmentCount; ++i) {
                const std::size_t first=i*(patch.points.size()-1)/segmentCount;
                const std::size_t second=(i+1)*(patch.points.size()-1)/segmentCount;
                proxy.segments.push_back({patch.points[first],patch.points[second]});
            }
        }
    }
    if (!proxy.points.empty()) {
        for (const Vec3& point : proxy.points) proxy.centroid=proxy.centroid+point;
        proxy.centroid=proxy.centroid*(1.0/static_cast<double>(proxy.points.size()));
        for (const Vec3& point : proxy.points) proxy.radius=std::max(proxy.radius,length(point-proxy.centroid));
    }
    return proxy;
}

GeometryProxy localProxy(const DebyeSphereCloud& cloud) {
    GeometryProxy proxy;
    const auto& scatterers=cloud.getScatterers();
    constexpr std::size_t maxPoints=64;
    const std::size_t pointCount=std::min(maxPoints,scatterers.size());
    for (std::size_t index=0; index<pointCount; ++index) {
        const std::size_t source=pointCount==1 ? 0 : index*(scatterers.size()-1)/(pointCount-1);
        const auto& center=scatterers[source].center;
        proxy.points.push_back({center.x,center.y,center.z});
    }
    if (!proxy.points.empty()) {
        for (const Vec3& point : proxy.points) proxy.centroid=proxy.centroid+point;
        proxy.centroid=proxy.centroid*(1.0/static_cast<double>(proxy.points.size()));
        for (std::size_t index=0; index<scatterers.size(); ++index) {
            const auto& scatterer=scatterers[index];
            const Vec3 center{scatterer.center.x,scatterer.center.y,scatterer.center.z};
            proxy.radius=std::max(proxy.radius,length(center-proxy.centroid)+scatterer.radius);
        }
    }
    return proxy;
}

struct TriangleMesh {
    std::vector<Vec3> points;
    std::vector<std::array<std::size_t,3>> triangles;
};

double displayedRadius(const SphereScatterer& scatterer,double zeroRadiusMarkerSize) {
    return scatterer.radius>0.0 ? scatterer.radius : zeroRadiusMarkerSize;
}

double automaticEnvelopePadding(const std::vector<SphereScatterer>& scatterers,
                                double zeroRadiusMarkerSize) {
    constexpr std::size_t maxSamples=256;
    const std::size_t sampleCount=std::min(maxSamples,scatterers.size());
    if (sampleCount==0) return zeroRadiusMarkerSize;
    std::vector<std::size_t> samples;
    samples.reserve(sampleCount);
    for (std::size_t index=0; index<sampleCount; ++index) {
        samples.push_back(sampleCount==1 ? 0 : index*(scatterers.size()-1)/(sampleCount-1));
    }
    double meanRadius=0.0;
    for (std::size_t index : samples) meanRadius+=displayedRadius(scatterers[index],zeroRadiusMarkerSize);
    meanRadius/=static_cast<double>(sampleCount);
    if (sampleCount==1) return std::max(0.1*meanRadius,zeroRadiusMarkerSize);

    std::vector<double> best(sampleCount,std::numeric_limits<double>::infinity());
    std::vector<bool> used(sampleCount,false);
    best[0]=0.0;
    double largestConnectionGap=0.0;
    for (std::size_t step=0; step<sampleCount; ++step) {
        std::size_t selected=sampleCount;
        for (std::size_t index=0; index<sampleCount; ++index) {
            if (!used[index] && (selected==sampleCount || best[index]<best[selected])) selected=index;
        }
        if (selected==sampleCount) break;
        used[selected]=true;
        largestConnectionGap=std::max(largestConnectionGap,best[selected]);
        const auto& first=scatterers[samples[selected]];
        const Vec3 firstCenter{first.center.x,first.center.y,first.center.z};
        for (std::size_t other=0; other<sampleCount; ++other) {
            if (used[other]) continue;
            const auto& second=scatterers[samples[other]];
            const Vec3 secondCenter{second.center.x,second.center.y,second.center.z};
            const double gap=std::max(0.0,length(firstCenter-secondCenter)
                -displayedRadius(first,zeroRadiusMarkerSize)
                -displayedRadius(second,zeroRadiusMarkerSize));
            best[other]=std::min(best[other],gap);
        }
    }
    return std::max(0.1*meanRadius,0.55*largestConnectionGap);
}

TriangleMesh debyeEnvelopeMesh(const DebyeSphereCloud& cloud,
                               std::size_t resolution,
                               double requestedPadding,
                               double zeroRadiusMarkerSize) {
    const auto& scatterers=cloud.getScatterers();
    TriangleMesh mesh;
    if (scatterers.empty()) return mesh;
    const double padding=requestedPadding>=0.0
        ? requestedPadding : automaticEnvelopePadding(scatterers,zeroRadiusMarkerSize);
    Vec3 minimum{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    Vec3 maximum{-minimum.x,-minimum.y,-minimum.z};
    for (const auto& scatterer : scatterers) {
        const double radius=displayedRadius(scatterer,zeroRadiusMarkerSize)+padding;
        minimum.x=std::min(minimum.x,scatterer.center.x-radius);
        minimum.y=std::min(minimum.y,scatterer.center.y-radius);
        minimum.z=std::min(minimum.z,scatterer.center.z-radius);
        maximum.x=std::max(maximum.x,scatterer.center.x+radius);
        maximum.y=std::max(maximum.y,scatterer.center.y+radius);
        maximum.z=std::max(maximum.z,scatterer.center.z+radius);
    }
    Vec3 span=maximum-minimum;
    const double longest=std::max({span.x,span.y,span.z,1e-12});
    const double margin=0.06*longest;
    minimum=minimum-Vec3{margin,margin,margin};
    maximum=maximum+Vec3{margin,margin,margin};
    span=maximum-minimum;
    const std::size_t nx=std::max<std::size_t>(8,std::lround(resolution*span.x/std::max({span.x,span.y,span.z})));
    const std::size_t ny=std::max<std::size_t>(8,std::lround(resolution*span.y/std::max({span.x,span.y,span.z})));
    const std::size_t nz=std::max<std::size_t>(8,std::lround(resolution*span.z/std::max({span.x,span.y,span.z})));
    const Vec3 step{span.x/(nx-1),span.y/(ny-1),span.z/(nz-1)};
    const auto gridIndex=[&](std::size_t x,std::size_t y,std::size_t z) {
        return (x*ny+y)*nz+z;
    };
    const auto gridPoint=[&](std::size_t x,std::size_t y,std::size_t z) {
        return Vec3{minimum.x+step.x*x,minimum.y+step.y*y,minimum.z+step.z*z};
    };
    std::vector<double> field(nx*ny*nz,std::numeric_limits<double>::infinity());
    for (std::size_t x=0; x<nx; ++x) for (std::size_t y=0; y<ny; ++y) for (std::size_t z=0; z<nz; ++z) {
        const Vec3 point=gridPoint(x,y,z);
        double value=std::numeric_limits<double>::infinity();
        for (const auto& scatterer : scatterers) {
            const Vec3 center{scatterer.center.x,scatterer.center.y,scatterer.center.z};
            const double radius=displayedRadius(scatterer,zeroRadiusMarkerSize)+padding;
            value=std::min(value,length(point-center)-radius);
        }
        field[gridIndex(x,y,z)]=value;
    }
    const auto interpolateIso=[](const Vec3& first,const Vec3& second,double firstValue,double secondValue) {
        const double denominator=firstValue-secondValue;
        const double fraction=std::abs(denominator)<1e-15 ? 0.5 : std::max(0.0,std::min(1.0,firstValue/denominator));
        return first+(second-first)*fraction;
    };
    std::map<std::array<long long,3>,std::size_t> vertexIndices;
    const auto vertexIndex=[&](const Vec3& point) {
        const double quantization=1e12/longest;
        const std::array<long long,3> key{{
            std::llround((point.x-minimum.x)*quantization),
            std::llround((point.y-minimum.y)*quantization),
            std::llround((point.z-minimum.z)*quantization)}};
        const auto existing=vertexIndices.find(key);
        if (existing!=vertexIndices.end()) return existing->second;
        const std::size_t index=mesh.points.size();
        mesh.points.push_back(point);
        vertexIndices.emplace(key,index);
        return index;
    };
    const auto addTriangle=[&](const Vec3& first,const Vec3& second,const Vec3& third) {
        mesh.triangles.push_back({vertexIndex(first),vertexIndex(second),vertexIndex(third)});
    };
    const std::array<std::array<int,4>,6> tetrahedra{{
        {{0,5,1,6}},{{0,1,2,6}},{{0,2,3,6}},
        {{0,3,7,6}},{{0,7,4,6}},{{0,4,5,6}}
    }};
    const std::array<std::array<int,3>,8> corners{{
        {{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},
        {{0,0,1}},{{1,0,1}},{{1,1,1}},{{0,1,1}}
    }};
    for (std::size_t x=0; x+1<nx; ++x) for (std::size_t y=0; y+1<ny; ++y) for (std::size_t z=0; z+1<nz; ++z) {
        std::array<Vec3,8> points;
        std::array<double,8> values;
        for (std::size_t corner=0; corner<corners.size(); ++corner) {
            const std::size_t cx=x+corners[corner][0],cy=y+corners[corner][1],cz=z+corners[corner][2];
            points[corner]=gridPoint(cx,cy,cz);
            values[corner]=field[gridIndex(cx,cy,cz)];
        }
        for (const auto& tetrahedron : tetrahedra) {
            std::vector<int> inside,outside;
            for (int corner : tetrahedron) (values[corner]<=0.0 ? inside : outside).push_back(corner);
            if (inside.empty() || outside.empty()) continue;
            if (inside.size()==1 || inside.size()==3) {
                const bool reverse=inside.size()==3;
                const int isolated=reverse ? outside.front() : inside.front();
                const auto& others=reverse ? inside : outside;
                const Vec3 first=interpolateIso(points[isolated],points[others[0]],values[isolated],values[others[0]]);
                const Vec3 second=interpolateIso(points[isolated],points[others[1]],values[isolated],values[others[1]]);
                const Vec3 third=interpolateIso(points[isolated],points[others[2]],values[isolated],values[others[2]]);
                if (reverse) addTriangle(first,third,second); else addTriangle(first,second,third);
            } else {
                const Vec3 ac=interpolateIso(points[inside[0]],points[outside[0]],values[inside[0]],values[outside[0]]);
                const Vec3 ad=interpolateIso(points[inside[0]],points[outside[1]],values[inside[0]],values[outside[1]]);
                const Vec3 bc=interpolateIso(points[inside[1]],points[outside[0]],values[inside[1]],values[outside[0]]);
                const Vec3 bd=interpolateIso(points[inside[1]],points[outside[1]],values[inside[1]],values[outside[1]]);
                addTriangle(ac,ad,bd); addTriangle(ac,bd,bc);
            }
        }
    }
    return mesh;
}

GeometryProxy transformProxy(const GeometryProxy& local, const Pose& pose) {
    GeometryProxy world;
    world.radius=local.radius;
    world.centroid=pose.translation+rotate(pose.rotation,local.centroid);
    for (const Vec3& point : local.points) world.points.push_back(pose.translation+rotate(pose.rotation,point));
    for (const ProxySegment& segment : local.segments) {
        world.segments.push_back({
            pose.translation+rotate(pose.rotation,segment.first),
            pose.translation+rotate(pose.rotation,segment.second)
        });
    }
    return world;
}

GeometryProxy maskJointNeighborhood(const GeometryProxy& proxy,
                                    const Vec3& joint,
                                    double radius) {
    if (radius <= 0.0) return proxy;
    GeometryProxy masked;
    masked.centroid=proxy.centroid;
    masked.radius=proxy.radius;
    const double radiusSquared=radius*radius;
    for (const Vec3& point : proxy.points) {
        if (dot(point-joint,point-joint)>radiusSquared) masked.points.push_back(point);
    }
    for (const ProxySegment& segment : proxy.segments) {
        const Vec3 direction=segment.second-segment.first;
        const Vec3 offset=segment.first-joint;
        const double a=dot(direction,direction);
        if (a<=1e-15) continue;
        const double b=2.0*dot(offset,direction);
        const double c=dot(offset,offset)-radiusSquared;
        std::vector<double> cuts{0.0,1.0};
        const double discriminant=b*b-4.0*a*c;
        if (discriminant>0.0) {
            const double root=std::sqrt(discriminant);
            const double first=(-b-root)/(2.0*a);
            const double second=(-b+root)/(2.0*a);
            if (first>0.0 && first<1.0) cuts.push_back(first);
            if (second>0.0 && second<1.0) cuts.push_back(second);
        }
        std::sort(cuts.begin(),cuts.end());
        for (std::size_t index=0; index+1<cuts.size(); ++index) {
            const double lower=cuts[index], upper=cuts[index+1];
            const Vec3 midpoint=segment.first+direction*((lower+upper)*0.5);
            if (dot(midpoint-joint,midpoint-joint)<=radiusSquared) continue;
            masked.segments.push_back({
                segment.first+direction*lower,
                segment.first+direction*upper
            });
        }
    }
    return masked;
}

double proxyDistanceSquared(const GeometryProxy& first, const GeometryProxy& second) {
    double result=std::numeric_limits<double>::infinity();
    for (const ProxySegment& a : first.segments) {
        for (const ProxySegment& b : second.segments) result=std::min(result,segmentDistanceSquared(a,b));
    }
    for (const Vec3& a : first.points) {
        for (const Vec3& b : second.points) result=std::min(result,dot(a-b,a-b));
    }
    return result;
}

double linkedProxyDistanceSquared(const GeometryProxy& first,
                                  const GeometryProxy& second,
                                  const Vec3& joint,
                                  double maskRadius) {
    return proxyDistanceSquared(
        maskJointNeighborhood(first,joint,maskRadius),
        maskJointNeighborhood(second,joint,maskRadius));
}

double readableScore(const GeometryProxy& candidate,
                     const std::map<std::string,GeometryProxy>& placed,
                     const std::string& parent,
                     const Vec3& joint,
                     double minimumClearance) {
    const double scale=std::max(candidate.radius,1e-9);
    double collisionPenalty=0.0;
    double proximityPenalty=0.0;
    std::size_t comparisonCount=0;
    Vec3 placedCentroid;
    std::size_t centroidCount=0;
    for (const auto& item : placed) {
        if (!item.second.points.empty()) {
            placedCentroid=placedCentroid+item.second.centroid;
            ++centroidCount;
        }
        if (item.second.points.empty() || candidate.points.empty()) continue;
        const double maskRadius=std::max(minimumClearance,0.02*scale);
        const double distanceSquared=item.first == parent
            ? linkedProxyDistanceSquared(candidate,item.second,joint,maskRadius)
            : proxyDistanceSquared(candidate,item.second);
        if (!std::isfinite(distanceSquared)) continue;
        const double distance=std::sqrt(std::max(0.0,distanceSquared));
        if (distance < minimumClearance) {
            const double deficit=(minimumClearance-distance)/scale;
            collisionPenalty+=deficit*deficit;
        }
        const double normalizedSquared=distanceSquared/(scale*scale);
        proximityPenalty+=1.0/(0.05+normalizedSquared);
        ++comparisonCount;
    }
    if (comparisonCount) proximityPenalty/=static_cast<double>(comparisonCount);

    double outward=0.0;
    const auto parentProxy=placed.find(parent);
    if (parentProxy != placed.end() && !parentProxy->second.points.empty()) {
        const Vec3 parentDirection=joint-parentProxy->second.centroid;
        const Vec3 childDirection=candidate.centroid-joint;
        const double denominator=length(parentDirection)*length(childDirection);
        if (denominator>1e-12) outward=dot(parentDirection,childDirection)/denominator;
    }

    double spread=0.0;
    if (centroidCount) {
        placedCentroid=placedCentroid*(1.0/static_cast<double>(centroidCount));
        spread=length(candidate.centroid-placedCentroid)/scale;
    }
    return collisionPenalty*1e6+proximityPenalty-2.0*outward-0.1*spread;
}

double pairLayoutScore(const GeometryProxy& first,
                       const GeometryProxy& second,
                       double minimumClearance,
                       const Vec3* joint) {
    if (first.points.empty() || second.points.empty()) return 0.0;
    const double scale=std::max(std::max(first.radius,second.radius),1e-9);
    const double lowerBound=std::max(
        0.0,length(first.centroid-second.centroid)-first.radius-second.radius);
    double distanceSquared;
    if (!joint && lowerBound>std::max(minimumClearance,2.0*scale)) {
        distanceSquared=lowerBound*lowerBound;
    } else if (joint) {
        const double maskRadius=std::max(minimumClearance,0.02*scale);
        distanceSquared=linkedProxyDistanceSquared(first,second,*joint,maskRadius);
    } else {
        distanceSquared=proxyDistanceSquared(first,second);
    }
    if (!std::isfinite(distanceSquared)) return 0.0;
    const double distance=std::sqrt(std::max(0.0,distanceSquared));
    double score=1.0/(0.05+distanceSquared/(scale*scale));
    if (distance<minimumClearance) {
        const double deficit=(minimumClearance-distance)/scale;
        score+=1e6*deficit*deficit;
    }
    return score;
}

double partitionLayoutScore(const std::set<std::string>& subtree,
                            const std::map<std::string,GeometryProxy>& candidateProxies,
                            const std::map<std::string,GeometryProxy>& currentProxies,
                            const std::string& parent,
                            const std::string& child,
                            const Vec3& joint,
                            double minimumClearance) {
    double score=0.0;
    for (const std::string& inside : subtree) {
        const auto first=candidateProxies.find(inside);
        if (first == candidateProxies.end()) continue;
        for (const auto& outside : currentProxies) {
            if (subtree.count(outside.first)) continue;
            const Vec3* sharedJoint=inside==child && outside.first==parent ? &joint : nullptr;
            score+=pairLayoutScore(
                first->second,outside.second,minimumClearance,sharedJoint);
        }
    }
    return score;
}

Pose rotatePoseAround(const Pose& pose, const Quaternion& rotation, const Vec3& pivot) {
    Pose result=pose;
    result.rotation=multiply(rotation,pose.rotation);
    result.translation=pivot+rotate(rotation,pose.translation-pivot);
    return result;
}

// USDA serialization: keep formatting and metadata in one place.
void writePose(std::ostream& output, const Pose& pose) {
    output << "        double3 xformOp:translate = ("
           << number(pose.translation.x) << ", "
           << number(pose.translation.y) << ", "
           << number(pose.translation.z) << ")\n";
    output << "        quatd xformOp:orient = ("
           << number(pose.rotation.w) << ", ("
           << number(pose.rotation.x) << ", "
           << number(pose.rotation.y) << ", "
           << number(pose.rotation.z) << "))\n";
    output << "        uniform token[] xformOpOrder = "
              "[\"xformOp:translate\", \"xformOp:orient\"]\n";
}

void writeReferencePrim(
    std::ostream& output,
    const std::string& referenceName,
    const Vec3& position) {
    output << "        def Xform \"ref_" << usdName(referenceName)
           << "\" {\n"
           << "            visibility = \"invisible\"\n"
           << "            custom string pyseb:reference = "
           << quote(referenceName) << "\n"
           << "            double3 xformOp:translate = ("
           << number(position.x) << ", "
           << number(position.y) << ", "
           << number(position.z) << ")\n"
           << "            uniform token[] xformOpOrder = "
              "[\"xformOp:translate\"]\n"
           << "        }\n";
}

void writeGeometryPrim(
    std::ostream& output,
    const std::string& name,
    const GeometryPatch& patch,
    const std::array<double, 3>& color,
    double opacity) {
    const bool isSurface =
        patch.kind == pyseb::VisualizationGeometryKind::Surface;
    output << "        def " << (isSurface ? "Mesh" : "BasisCurves")
           << " \"" << usdName(name) << "\" {\n";

    if (isSurface) {
        output << "            int[] faceVertexCounts = [";
        for (std::size_t u = 0; u + 1 < patch.uCount; ++u) {
            for (std::size_t v = 0; v + 1 < patch.vCount; ++v) {
                output << "4, ";
            }
        }
        output << "]\n            int[] faceVertexIndices = [";
        for (std::size_t u = 0; u + 1 < patch.uCount; ++u) {
            for (std::size_t v = 0; v + 1 < patch.vCount; ++v) {
                const std::size_t index = u * patch.vCount + v;
                output << index << ", "
                       << index + 1 << ", "
                       << index + patch.vCount + 1 << ", "
                       << index + patch.vCount << ", ";
            }
        }
        output << "]\n";
    } else {
        output << "            int[] curveVertexCounts = ["
               << patch.points.size() << "]\n";
    }

    output << "            point3f[] points = [";
    for (const Vec3& point : patch.points) {
        output << "(" << number(point.x) << ", "
               << number(point.y) << ", "
               << number(point.z) << "), ";
    }
    output << "]\n"
           << "            color3f[] primvars:displayColor = [("
           << number(color[0]) << ", "
           << number(color[1]) << ", "
           << number(color[2]) << ")]\n"
           << "            float[] primvars:displayOpacity = ["
           << number(opacity) << "]\n"
           << "        }\n";
}

void writeCloudPrim(
    std::ostream& output,
    const std::string& structure,
    const std::string& name,
    const DebyeSphereCloud& cloud,
    const USDExportOptions& options,
    const std::map<std::string, Vec3>& references) {
    std::array<double, 3> color{{0.7, 0.7, 0.7}};
    const auto colorOverride = options.colorOverrides.find(name);
    if (colorOverride != options.colorOverrides.end()) {
        color = colorOverride->second;
    }

    double opacity = 1.0;
    const auto opacityOverride = options.opacityOverrides.find(name);
    if (opacityOverride != options.opacityOverrides.end()) {
        opacity = opacityOverride->second;
    }

    const auto& scatterers = cloud.getScatterers();
    output << "        def PointInstancer \"cloud\" {\n"
           << "            rel prototypes = [</" << usdName(structure)
           << "/" << usdName(name) << "/spherePrototype>]\n"
           << "            point3f[] positions = [";
    for (const auto& scatterer : scatterers) {
        output << "(" << number(scatterer.center.x) << ", "
               << number(scatterer.center.y) << ", "
               << number(scatterer.center.z) << "), ";
    }
    output << "]\n            int[] protoIndices = [";
    for (std::size_t index = 0; index < scatterers.size(); ++index) {
        output << "0, ";
    }
    output << "]\n            float3[] scales = [";
    for (const auto& scatterer : scatterers) {
        const double radius = displayedRadius(
            scatterer, options.zeroRadiusMarkerSize);
        output << "(" << number(radius) << ", "
               << number(radius) << ", "
               << number(radius) << "), ";
    }
    output << "]\n"
           << "            custom float[] pyseb:radius = [";
    for (const auto& scatterer : scatterers) {
        output << number(scatterer.radius) << ", ";
    }
    output << "]\n            custom float[] pyseb:beta = [";
    for (const auto& scatterer : scatterers) {
        output << number(scatterer.beta) << ", ";
    }
    output << "]\n            custom int[] pyseb:index = [";
    for (std::size_t index = 0; index < scatterers.size(); ++index) {
        output << index << ", ";
    }
    output << "]\n            color3f[] primvars:displayColor = [";
    for (std::size_t index = 0; index < scatterers.size(); ++index) {
        output << "(" << number(color[0]) << ", "
               << number(color[1]) << ", "
               << number(color[2]) << "), ";
    }
    output << "]\n            float[] primvars:displayOpacity = [";
    for (std::size_t index = 0; index < scatterers.size(); ++index) {
        output << number(opacity) << ", ";
    }
    output << "]\n        }\n"
           << "        def Sphere \"spherePrototype\" { float radius = 1 }\n";

    if (options.debyeEnvelope) {
        const double padding = options.debyeEnvelopePadding >= 0.0
            ? options.debyeEnvelopePadding
            : automaticEnvelopePadding(
                scatterers, options.zeroRadiusMarkerSize);
        const TriangleMesh envelope = debyeEnvelopeMesh(
            cloud,
            options.debyeEnvelopeResolution,
            padding,
            options.zeroRadiusMarkerSize);
        if (!envelope.triangles.empty()) {
            output << "        def Mesh \"envelope\" {\n"
                   << "            custom bool pyseb:visualOnly = true\n"
                   << "            custom string pyseb:construction = "
                      "\"union_of_inflated_scatterer_spheres\"\n"
                   << "            custom double pyseb:padding = "
                   << number(padding) << "\n"
                   << "            custom uint64 pyseb:resolution = "
                   << options.debyeEnvelopeResolution << "\n"
                   << "            uniform bool doubleSided = true\n"
                   << "            int[] faceVertexCounts = [";
            for (std::size_t index = 0; index < envelope.triangles.size(); ++index) {
                output << "3, ";
            }
            output << "]\n            int[] faceVertexIndices = [";
            for (const auto& triangle : envelope.triangles) {
                output << triangle[0] << ", "
                       << triangle[1] << ", "
                       << triangle[2] << ", ";
            }
            output << "]\n            point3f[] points = [";
            for (const Vec3& point : envelope.points) {
                output << "(" << number(point.x) << ", "
                       << number(point.y) << ", "
                       << number(point.z) << "), ";
            }
            output << "]\n"
                   << "            color3f[] primvars:displayColor = [("
                   << number(color[0]) << ", "
                   << number(color[1]) << ", "
                   << number(color[2]) << ")]\n"
                   << "            float[] primvars:displayOpacity = ["
                   << number(options.debyeEnvelopeOpacity * opacity) << "]\n"
                   << "        }\n";
        }
    }

    for (const auto& reference : references) {
        writeReferencePrim(output, reference.first, reference.second);
    }
    output << "        custom string pyseb:model_id = "
           << quote("pyseb/DebyeSphereCloud") << "\n"
           << "        custom string pyseb:instance_path = " << quote(name)
           << "\n"
           << "        custom uint64 pyseb:seed = " << options.seed << "\n"
           << "    }\n";
}

void writeFileInstanceBody(
    std::ostream& output,
    const InstanceScene& instance,
    unsigned long long seed) {
    for (const auto& reference : instance.resolvedReferences) {
        writeReferencePrim(output, reference.first, reference.second);
    }
    for (const auto& patch : instance.patches) {
        writeGeometryPrim(
            output,
            patch.first,
            patch.second,
            instance.color,
            instance.opacity);
    }
    output << "        custom string pyseb:model_id = "
           << quote(instance.subunit->getDefinition().id) << "\n"
           << "        custom string pyseb:instance_path = "
           << quote(instance.name) << "\n"
           << "        custom string pyseb:tag = " << quote(instance.tag)
           << "\n"
           << "        custom uint64 pyseb:seed = " << seed
           << "\n";
    for (const auto& parameter : instance.parameters) {
        output << "        custom double pyseb:param_"
               << usdName(parameter.first) << " = "
               << number(parameter.second) << "\n";
    }
    output << "    }\n";
}

void writeLayoutMetadata(std::ostream& output, const USDExportOptions& options) {
    if (options.layoutMode == USDLayoutMode::Readable) {
        output << "    custom string pyseb:layoutMode = \"readable\"\n"
               << "    custom string pyseb:orientationSemantics = "
                  "\"readability_optimized_free_rotation\"\n"
               << "    custom string pyseb:overlapPolicy = "
                  "\"best_effort_minimized\"\n"
               << "    custom uint64 pyseb:orientationTrials = "
               << options.orientationTrials << "\n"
               << "    custom uint64 pyseb:relaxationSweeps = "
               << options.relaxationSweeps << "\n"
               << "    custom double pyseb:minimumClearance = "
               << number(options.minimumClearance) << "\n";
        return;
    }
    output << "    custom string pyseb:layoutMode = \"random\"\n"
           << "    custom string pyseb:orientationSemantics = "
              "\"representative_free_rotation\"\n"
           << "    custom string pyseb:overlapPolicy = \"permitted\"\n";
}

void writeLinkRelationships(
    std::ostream& output,
    const std::string& structure,
    const std::vector<SceneLink>& links) {
    for (std::size_t index = 0; index < links.size(); ++index) {
        const auto& link = links[index];
        output << "    rel pyseb:link_" << index << " = [</"
               << usdName(structure) << "/" << usdName(link.first.instance)
               << "/ref_" << usdName(link.first.reference) << ">, </"
               << usdName(structure) << "/" << usdName(link.second.instance)
               << "/ref_" << usdName(link.second.reference) << ">]\n";
    }
}

} // namespace

void World::ExportUSD(const std::string& structure,
                      const std::string& path,
                      const ParameterList& parameters,
                      const USDExportOptions& options) const {
    const bool validExtension = path.size() >= 4 &&
        (path.substr(path.size()-4) == ".usd" ||
         (path.size() >= 5 && path.substr(path.size()-5) == ".usda"));
    if (!validExtension) throw SEBException("USD output must use .usd or .usda", "World::ExportUSD");
    const double metersPerUnit = options.metersPerUnit;
    if (!(metersPerUnit > 0.0) || !std::isfinite(metersPerUnit)) {
        throw SEBException("metersPerUnit must be positive and finite", "World::ExportUSD");
    }
    if (options.orientationTrials == 0) {
        throw SEBException("orientationTrials must be at least one", "World::ExportUSD");
    }
    if (options.minimumClearance < 0.0 || !std::isfinite(options.minimumClearance)) {
        throw SEBException("minimumClearance must be non-negative and finite", "World::ExportUSD");
    }
    if (!(options.zeroRadiusMarkerSize>0.0) || !std::isfinite(options.zeroRadiusMarkerSize)) {
        throw SEBException("zeroRadiusMarkerSize must be positive and finite", "World::ExportUSD");
    }
    if (options.debyeEnvelopeResolution<8) {
        throw SEBException("debyeEnvelopeResolution must be at least eight", "World::ExportUSD");
    }
    if (!std::isfinite(options.debyeEnvelopePadding)
        || (options.debyeEnvelopePadding<0.0 && options.debyeEnvelopePadding!=-1.0)) {
        throw SEBException("debyeEnvelopePadding must be non-negative or -1 for automatic", "World::ExportUSD");
    }
    if (!std::isfinite(options.debyeEnvelopeOpacity)
        || options.debyeEnvelopeOpacity<0.0 || options.debyeEnvelopeOpacity>1.0) {
        throw SEBException("debyeEnvelopeOpacity must be between zero and one", "World::ExportUSD");
    }

    std::set<std::string> selected;
    const auto rootIt=nameCatalog.find(structure);
    if (rootIt != nameCatalog.end()) {
        const Structure* root=dynamic_cast<const Structure*>(rootIt->second);
        if (root) {
            std::set<GraphID> visiting;
            std::function<void(GraphID)> collect=[&](GraphID id) {
                if (!visiting.insert(id).second) return;
                const auto graph=subGraphs.find(id);
                if (graph != subGraphs.end()) for (const auto& childName : graph->second) {
                    const auto child=nameCatalog.find(childName);
                    if (child == nameCatalog.end()) continue;
                    const Structure* nested=dynamic_cast<const Structure*>(child->second);
                    if (nested) collect(nested->getGraphID()); else selected.insert(childName);
                }
                visiting.erase(id);
            };
            collect(root->getGraphID());
        } else selected.insert(structure);
    } else {
        for (const auto& item : nameCatalog) if (dynamic_cast<SubUnit*>(item.second)) selected.insert(item.first);
    }
    if (selected.empty()) throw SEBException("structure contains no visualizable leaf instances", "World::ExportUSD");

    std::map<std::string,InstanceScene> instances;
    for (const auto& name : selected) {
        const auto catalog=nameCatalog.find(name);
        const SubUnit* sub=catalog == nameCatalog.end() ? nullptr : dynamic_cast<const SubUnit*>(catalog->second);
        const FileDefinedSubunit* file=dynamic_cast<const FileDefinedSubunit*>(sub);
        if (!file || !file->getDefinition().visualization.present) continue;
        InstanceScene instance;
        instance.subunit=file;
        instance.name=name;
        instance.tag=const_cast<SubUnit*>(sub)->getTag();
        instance.color=file->getDefinition().visualization.color;
        instance.opacity=file->getDefinition().visualization.opacity;
        const auto colorOverride=options.colorOverrides.find(name);
        if (colorOverride != options.colorOverrides.end()) instance.color=colorOverride->second;
        const auto opacityOverride=options.opacityOverrides.find(name);
        if (opacityOverride != options.opacityOverrides.end()) instance.opacity=opacityOverride->second;
        for (const auto& parameter : file->getDefinition().parameters) {
            auto value=parameters.find(parameter.first+"_"+instance.tag);
            if (value == parameters.end()) value=parameters.find(parameter.first);
            if (value == parameters.end()) throw SEBException("missing numerical parameter "+parameter.first+"_"+instance.tag, "World::ExportUSD");
            if (!std::isfinite(value->second)) throw SEBException("visualization parameter must be finite", "World::ExportUSD");
            instance.parameters[parameter.first]=value->second;
        }

        const auto& visualization=file->getDefinition().visualization;
        std::map<std::string,double> definitions;
        std::set<std::string> active;
        std::function<double(const std::string&)> resolve;
        std::function<double(const pyseb::ParsedExpression&,const std::map<std::string,double>&)> evaluate;
        resolve=[&](const std::string& identifier) -> double {
            const auto parameter=instance.parameters.find(identifier);
            if (parameter != instance.parameters.end()) return parameter->second;
            const auto cached=definitions.find(identifier);
            if (cached != definitions.end()) return cached->second;
            const auto definition=visualization.definitions.find(identifier);
            if (definition != visualization.definitions.end()) {
                if (!active.insert(identifier).second) throw SEBException("cyclic visualization definition", "World::ExportUSD");
                const double value=evaluate(definition->second,{});
                active.erase(identifier);
                definitions[identifier]=value;
                return value;
            }
            if (identifier == "pi") return M_PI;
            if (identifier == "e") return M_E;
            throw SEBException("unknown visualization identifier "+identifier, "World::ExportUSD");
        };
        evaluate=[&](const pyseb::ParsedExpression& expression,const std::map<std::string,double>& local) {
            const double value=pyseb::EvaluateSubunitExpression(expression,[&](const std::string& identifier) {
                const auto found=local.find(identifier);
                return found == local.end() ? resolve(identifier) : found->second;
            });
            if (!std::isfinite(value)) throw SEBException("non-finite visualization expression", "World::ExportUSD");
            return value;
        };

        for (const auto& geometry : visualization.geometry) {
            const auto& definition = geometry.second;
            const std::size_t sampleCount = definition.samples
                ? definition.samples
                : definition.kind == pyseb::VisualizationGeometryKind::Surface
                    ? options.surfaceSamples
                    : options.curveSamples;
            if (sampleCount < 2) {
                throw SEBException(
                    "visualization sampling requires at least two points",
                    "World::ExportUSD");
            }
            instance.patches.emplace(
                geometry.first,
                realizeGeometry(
                    definition,
                    sampleCount,
                    evaluate,
                    options.seed,
                    name,
                    geometry.first));
        }

        resolveSpecificReferences(instance, visualization, evaluate);
        instances.emplace(name,std::move(instance));
    }

    std::map<std::string,std::map<std::string,Vec3>> numericalReferences;
    for (const auto& name : selected) {
        const auto catalog=nameCatalog.find(name);
        const DebyeSphereCloud* cloud=catalog == nameCatalog.end()
            ? nullptr : dynamic_cast<const DebyeSphereCloud*>(catalog->second);
        if (!cloud) continue;
        for (const auto& reference : cloud->getReferenceCoordinates()) {
            numericalReferences[name][reference.first]={
                reference.second.x,reference.second.y,reference.second.z};
        }
    }

    std::map<std::string,Pose> poses;
    for (const auto& name : selected) poses[name]=Pose();
    poses.begin()->second.resolved=true;

    std::map<std::string,GeometryProxy> localProxies;
    for (const auto& instance : instances) localProxies.emplace(instance.first,localProxy(instance.second));
    for (const auto& name : selected) {
        const auto catalog=nameCatalog.find(name);
        const DebyeSphereCloud* cloud=catalog == nameCatalog.end()
            ? nullptr : dynamic_cast<const DebyeSphereCloud*>(catalog->second);
        if (cloud) localProxies.emplace(name,localProxy(*cloud));
    }
    std::map<std::string,GeometryProxy> placedProxies;
    const auto rootProxy=localProxies.find(poses.begin()->first);
    if (rootProxy != localProxies.end()) {
        placedProxies.emplace(rootProxy->first,transformProxy(rootProxy->second,poses.begin()->second));
    }

    auto resolveReference=[&](const LinkEndpoint& endpoint) -> Vec3 {
        auto instance=instances.find(endpoint.instance);
        if (instance == instances.end()) {
            const auto numerical=numericalReferences.find(endpoint.instance);
            if (numerical == numericalReferences.end()) throw SEBException("unknown visualization instance "+endpoint.instance,"World::ExportUSD");
            const auto reference=numerical->second.find(endpoint.reference);
            if (reference == numerical->second.end()) throw SEBException("unknown numerical visualization reference "+endpoint.reference,"World::ExportUSD");
            return reference->second;
        }
        auto cached=instance->second.resolvedReferences.find(endpoint.reference);
        if (cached != instance->second.resolvedReferences.end()) return cached->second;
        const auto parts=parseVisualizationReference(endpoint.reference);
        const auto base=instance->second.subunit->getDefinition().visualization.references.find(parts.base);
        if (base == instance->second.subunit->getDefinition().visualization.references.end()) throw SEBException("unknown visualization reference "+endpoint.reference,"World::ExportUSD");
        const auto& reference=base->second;
        if (reference.kind != "distributed") {
            const auto fixed=instance->second.resolvedReferences.find(parts.base);
            if (fixed == instance->second.resolvedReferences.end()) throw SEBException("unresolved visualization reference "+endpoint.reference,"World::ExportUSD");
            return fixed->second;
        }
        if (parts.label.empty()) throw SEBException("distributed link reference requires a #label", "World::ExportUSD");
        std::mt19937_64 generator(derivedSeed(options.seed,endpoint.instance+"/reference/"+endpoint.reference));
        std::string patchName=reference.geometry;
        std::string sampling=reference.sampling;
        if (!parts.variant.empty()) {
            const auto variant=reference.variants.find(parts.variant);
            if (variant == reference.variants.end()) throw SEBException("unknown visualization reference variant '"+parts.variant+"'", "World::ExportUSD");
            patchName=variant->second.geometry;
            sampling=variant->second.sampling;
        } else if (reference.sampling == "weighted_mixture") {
            if (reference.patches.empty() || reference.patches.size()!=reference.weights.size()) throw SEBException("weighted mixture requires equally sized patches and weights", "World::ExportUSD");
            double total=0.0; for(double weight:reference.weights){if(weight<0.0)throw SEBException("mixture weights must be non-negative","World::ExportUSD");total+=weight;}
            if(total<=0.0)throw SEBException("mixture weights must have positive total","World::ExportUSD");
            double target=unitRandom(generator)*total; std::size_t index=0; for(;index+1<reference.weights.size()&&target>=reference.weights[index];++index)target-=reference.weights[index];
            patchName=reference.patches[index];
        }
        const auto patch=instance->second.patches.find(patchName);
        if (patch == instance->second.patches.end()) throw SEBException("distributed reference uses unknown geometry patch", "World::ExportUSD");
        Vec3 point;
        if (patch->second.kind == pyseb::VisualizationGeometryKind::Surface) {
            point=sampling == "arc_length" ? sampleSurfaceBoundary(patch->second,unitRandom(generator)) : sampleSurfaceArea(patch->second,generator);
        } else {
            point=sampling == "uniform_parameter" ? sampleCurveParameter(patch->second,unitRandom(generator)) : sampleCurveLength(patch->second,unitRandom(generator));
        }
        instance->second.resolvedReferences[endpoint.reference]=point;
        return point;
    };

    std::vector<SceneLink> sceneLinks;
    for (const auto& link : links) {
        SceneLink sceneLink;
        if (!parseEndpoint(link.first,sceneLink.first) || !parseEndpoint(link.second,sceneLink.second)) continue;
        if (!selected.count(sceneLink.first.instance) || !selected.count(sceneLink.second.instance)) continue;
        sceneLink.key=link.first+"<->"+link.second;
        sceneLinks.push_back(sceneLink);
    }
    std::sort(sceneLinks.begin(),sceneLinks.end(),[](const SceneLink& a,const SceneLink& b){return a.key<b.key;});
    std::vector<RootedLink> rootedLinks;
    std::size_t resolvedCount=1;
    while (resolvedCount < poses.size()) {
        bool progress=false;
        for (const auto& link : sceneLinks) {
            Pose& first=poses[link.first.instance]; Pose& second=poses[link.second.instance];
            if (first.resolved == second.resolved) continue;
            const bool firstIsParent=first.resolved;
            Pose& parent=firstIsParent?first:second;
            Pose& child=firstIsParent?second:first;
            const LinkEndpoint& parentEndpoint=firstIsParent?link.first:link.second;
            const LinkEndpoint& childEndpoint=firstIsParent?link.second:link.first;
            const Vec3 parentLocal=resolveReference(parentEndpoint);
            const Vec3 childLocal=resolveReference(childEndpoint);
            const Vec3 parentWorld=parent.translation+rotate(parent.rotation,parentLocal);
            const std::string orientationKey=childEndpoint.instance+"/orientation/"+link.key;
            if (options.layoutMode == USDLayoutMode::Random) {
                child.rotation=uniformRotation(options.seed,orientationKey);
                child.translation=parentWorld-rotate(child.rotation,childLocal);
            } else {
                Pose best;
                double bestScore=std::numeric_limits<double>::infinity();
                const auto localProxy=localProxies.find(childEndpoint.instance);
                for (std::size_t trial=0; trial<options.orientationTrials; ++trial) {
                    Pose candidate;
                    candidate.rotation=uniformRotation(options.seed,orientationKey+"/candidate/"+std::to_string(trial));
                    candidate.translation=parentWorld-rotate(candidate.rotation,childLocal);
                    const GeometryProxy candidateProxy=localProxy == localProxies.end()
                        ? GeometryProxy() : transformProxy(localProxy->second,candidate);
                    const double score=readableScore(
                        candidateProxy,placedProxies,parentEndpoint.instance,parentWorld,options.minimumClearance);
                    if (score < bestScore) {
                        bestScore=score;
                        best=candidate;
                    }
                }
                child.rotation=best.rotation;
                child.translation=best.translation;
            }
            child.resolved=true; ++resolvedCount; progress=true;
            const auto childProxy=localProxies.find(childEndpoint.instance);
            if (childProxy != localProxies.end()) {
                placedProxies[childEndpoint.instance]=transformProxy(childProxy->second,child);
            }
            rootedLinks.push_back({parentEndpoint,childEndpoint,link.key});
        }
        if (!progress) {
            for (auto& pose : poses) if (!pose.second.resolved) {
                pose.second.resolved=true;
                const auto proxy=localProxies.find(pose.first);
                if (proxy != localProxies.end()) placedProxies[pose.first]=transformProxy(proxy->second,pose.second);
                ++resolvedCount;
                break;
            }
        }
    }

    if (options.layoutMode == USDLayoutMode::Readable && options.relaxationSweeps>0) {
        std::map<std::string,std::vector<std::string>> children;
        for (const RootedLink& link : rootedLinks) {
            children[link.parent.instance].push_back(link.child.instance);
        }
        for (auto& item : children) std::sort(item.second.begin(),item.second.end());

        std::vector<std::set<std::string>> subtrees;
        for (const RootedLink& link : rootedLinks) {
            std::set<std::string> subtree;
            std::function<void(const std::string&)> collect=[&](const std::string& name) {
                if (!subtree.insert(name).second) return;
                const auto found=children.find(name);
                if (found != children.end()) for (const std::string& child : found->second) collect(child);
            };
            collect(link.child.instance);
            subtrees.push_back(std::move(subtree));
        }

        for (std::size_t sweep=0; sweep<options.relaxationSweeps; ++sweep) {
            for (std::size_t order=0; order<rootedLinks.size(); ++order) {
                const std::size_t index=sweep%2==0 ? order : rootedLinks.size()-1-order;
                const RootedLink& link=rootedLinks[index];
                const std::set<std::string>& subtree=subtrees[index];
                const Vec3 parentLocal=resolveReference(link.parent);
                const Pose& parentPose=poses[link.parent.instance];
                const Vec3 joint=parentPose.translation+rotate(parentPose.rotation,parentLocal);

                Quaternion bestRotation;
                double bestScore=std::numeric_limits<double>::infinity();
                for (std::size_t trial=0; trial<options.orientationTrials; ++trial) {
                    const Quaternion candidateRotation=trial==0 ? Quaternion() : uniformRotation(
                        options.seed,link.key+"/relaxation/"+std::to_string(sweep)+"/"+std::to_string(trial));
                    std::map<std::string,GeometryProxy> candidateProxies;
                    for (const std::string& name : subtree) {
                        const auto local=localProxies.find(name);
                        if (local == localProxies.end()) continue;
                        const Pose candidatePose=rotatePoseAround(poses[name],candidateRotation,joint);
                        candidateProxies.emplace(name,transformProxy(local->second,candidatePose));
                    }
                    const double score=partitionLayoutScore(
                        subtree,candidateProxies,placedProxies,
                        link.parent.instance,link.child.instance,joint,options.minimumClearance);
                    if (score<bestScore) {
                        bestScore=score;
                        bestRotation=candidateRotation;
                    }
                }
                for (const std::string& name : subtree) {
                    poses[name]=rotatePoseAround(poses[name],bestRotation,joint);
                    const auto local=localProxies.find(name);
                    if (local != localProxies.end()) {
                        placedProxies[name]=transformProxy(local->second,poses[name]);
                    }
                }
            }
        }
    }

    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << "#usda 1.0\n(\n    upAxis = \"Z\"\n    metersPerUnit = " << number(metersPerUnit) << "\n)\n\n";
    output << "def Xform \"" << usdName(structure) << "\" {\n";
    writeLayoutMetadata(output, options);
    for (const auto& name : selected) {
        const auto catalog=nameCatalog.find(name);
        const SubUnit* sub=catalog == nameCatalog.end()?nullptr:dynamic_cast<const SubUnit*>(catalog->second);
        if (!sub) continue;
        const Pose& pose=poses[name];
        output << "    def Xform \"" << usdName(name) << "\" {\n";
        output << "        double3 xformOp:translate = ("<<number(pose.translation.x)<<", "<<number(pose.translation.y)<<", "<<number(pose.translation.z)<<")\n";
        output << "        quatd xformOp:orient = ("<<number(pose.rotation.w)<<", ("<<number(pose.rotation.x)<<", "<<number(pose.rotation.y)<<", "<<number(pose.rotation.z)<<"))\n";
        output << "        uniform token[] xformOpOrder = [\"xformOp:translate\", \"xformOp:orient\"]\n";
        const DebyeSphereCloud* cloud=dynamic_cast<const DebyeSphereCloud*>(sub);
        if (cloud) {
            const auto references = numericalReferences.find(name);
            const std::map<std::string, Vec3> emptyReferences;
            writeCloudPrim(
                output,
                structure,
                name,
                *cloud,
                options,
                references == numericalReferences.end()
                    ? emptyReferences
                    : references->second);
            continue;
        }
        const auto instance=instances.find(name);
        if (instance == instances.end()) { output<<"    }\n"; continue; }
        writeFileInstanceBody(output, instance->second, options.seed);
    }
    writeLinkRelationships(output, structure, sceneLinks);
    output<<"}\n";

    const std::string temporary=path+".tmp";
    {
        std::ofstream file(temporary,std::ios::binary);
        if(!file)throw SEBException("unable to open USD output","World::ExportUSD");
        file<<output.str();
        if(!file)throw SEBException("failed writing USD output","World::ExportUSD");
    }
    if(std::rename(temporary.c_str(),path.c_str())!=0){std::remove(temporary.c_str());throw SEBException("unable to commit USD output","World::ExportUSD");}
}
