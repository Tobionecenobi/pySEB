// Implementation-only types and interfaces shared by the USD exporter files.
#ifndef PYSEB_USD_EXPORT_INTERNAL_HPP
#define PYSEB_USD_EXPORT_INTERNAL_HPP

#include "USDExport.hpp"
#include "World.hpp"
#include "Subunits/FileDefined.hpp"
#include "Subunits/DebyeSphereCloud.hpp"
#include "Subunits/AtomicSAXSCloud.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <list>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace pyseb {
namespace usd_detail {

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

struct Quaternion {
    double w = 1.0, x = 0.0, y = 0.0, z = 0.0;
};

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
    std::string sourceName;
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

struct OccurrenceNode {
    std::string path;
    std::string sourceName;
    bool structure = false;
    std::vector<OccurrenceNode> children;
};

struct OccurrenceTree {
    std::vector<OccurrenceNode> roots;
    std::vector<SceneLink> links;
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

struct VisualizationReferenceName {
    std::string base;
    std::string variant;
    std::string label;
};

// Helpers shared across implementation files.
double dot(const Vec3& a, const Vec3& b);

double length(const Vec3& value);

Vec3 rotate(const Quaternion& q, const Vec3& value);

Quaternion multiply(const Quaternion& first, const Quaternion& second);

std::string quote(const std::string& value);

double number(double value);

std::string usdName(const std::string& value);

unsigned long long derivedSeed(unsigned long long seed, const std::string& key);

double unitRandom(std::mt19937_64& generator);

Quaternion uniformRotation(unsigned long long seed, const std::string& key);

std::string usdPrimPath(
    const std::string& root,
    const std::string& occurrence,
    const std::string& child = std::string());

VisualizationReferenceName parseVisualizationReference(const std::string& reference);

Vec3 sampleCurveParameter(const GeometryPatch& patch, double parameter);

Vec3 sampleCurveLength(const GeometryPatch& patch, double draw);

Vec3 sampleSurfaceArea(const GeometryPatch& patch, std::mt19937_64& generator);

Vec3 sampleSurfaceBoundary(const GeometryPatch& patch, double draw);

GeometryPatch realizeGeometry(
    const pyseb::VisualizationGeometry& definition,
    std::size_t sampleCount,
    const ExpressionEvaluator& evaluate,
    unsigned long long seed,
    const std::string& instanceName,
    const std::string& geometryName);

void resolveSpecificReferences(
    InstanceScene& instance,
    const pyseb::VisualizationDefinition& visualization,
    const ExpressionEvaluator& evaluate);

GeometryProxy localProxy(const InstanceScene& instance);

GeometryProxy localProxy(const DebyeSphereCloud& cloud);

double atomicDisplayedRadius(
    const AtomicSAXSScatterer& scatterer,
    double zeroRadiusMarkerSize);

GeometryProxy localProxy(const AtomicSAXSCloud& cloud);

GeometryProxy transformProxy(const GeometryProxy& local, const Pose& pose);

double readableScore(const GeometryProxy& candidate,
                     const std::map<std::string,GeometryProxy>& placed,
                     const std::string& parent,
                     const Vec3& joint,
                     double minimumClearance);

double partitionLayoutScore(const std::set<std::string>& subtree,
                            const std::map<std::string,GeometryProxy>& candidateProxies,
                            const std::map<std::string,GeometryProxy>& currentProxies,
                            const std::string& parent,
                            const std::string& child,
                            const Vec3& joint,
                            double minimumClearance);

Pose rotatePoseAround(const Pose& pose, const Quaternion& rotation, const Vec3& pivot);

void writeCloudPrim(
    std::ostream& output,
    const std::string& structure,
    const std::string& name,
    const DebyeSphereCloud& cloud,
    const USDExportOptions& options,
    const std::map<std::string, Vec3>& references);

void writeAtomicCloudPrim(
    std::ostream& output,
    const std::string& structure,
    const std::string& name,
    const AtomicSAXSCloud& cloud,
    const USDExportOptions& options,
    const std::map<std::string, Vec3>& references);

void writeFileInstanceBody(
    std::ostream& output,
    const InstanceScene& instance,
    unsigned long long seed);

void writeLayoutMetadata(std::ostream& output, const USDExportOptions& options);

void writeLinkRelationships(
    std::ostream& output,
    const std::string& structure,
    const std::vector<SceneLink>& links);

} // namespace usd_detail
} // namespace pyseb

#endif
