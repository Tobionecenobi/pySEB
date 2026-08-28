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
    Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
    Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
    Vec3 operator*(double scale) const { return {x * scale, y * scale, z * scale}; }
};

double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
double length(const Vec3& value) { return std::sqrt(dot(value, value)); }

struct Quaternion {
    double w = 1.0, x = 0.0, y = 0.0, z = 0.0;
};

Vec3 rotate(const Quaternion& q, const Vec3& value) {
    const Vec3 u{q.x, q.y, q.z};
    return value + cross(u, value) * (2.0*q.w) + cross(u, cross(u, value)) * 2.0;
}

Quaternion multiply(const Quaternion& first, const Quaternion& second) {
    return {
        first.w*second.w-first.x*second.x-first.y*second.y-first.z*second.z,
        first.w*second.x+first.x*second.w+first.y*second.z-first.z*second.y,
        first.w*second.y-first.x*second.z+first.y*second.w+first.z*second.x,
        first.w*second.z+first.x*second.y-first.y*second.x+first.z*second.w
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

std::pair<std::string,std::string> splitReference(const std::string& reference) {
    const std::size_t hash = reference.find('#');
    return hash == std::string::npos
        ? std::make_pair(reference, std::string())
        : std::make_pair(reference.substr(0, hash), reference.substr(hash + 1));
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
            const auto& definition=geometry.second;
            GeometryPatch patch;
            patch.kind=definition.kind;
            std::size_t count=definition.samples ? definition.samples :
                (definition.kind == pyseb::VisualizationGeometryKind::Surface ? options.surfaceSamples : options.curveSamples);
            if (count < 2) throw SEBException("visualization sampling requires at least two points", "World::ExportUSD");
            if (definition.kind == pyseb::VisualizationGeometryKind::RandomWalk) {
                std::mt19937_64 generator(derivedSeed(options.seed,name+"/geometry/"+geometry.first));
                std::normal_distribution<double> normal(0.0,1.0);
                patch.points.resize(count);
                for (std::size_t i=1; i<count; ++i) {
                    Vec3 step{normal(generator),normal(generator),normal(generator)};
                    if (definition.distribution == "fixed_length") { const double size=length(step); if (size>0.0) step=step*(1.0/size); }
                    patch.points[i]=patch.points[i-1]+step;
                }
                if (definition.closure == "bridge") {
                    const Vec3 end=patch.points.back();
                    for (std::size_t i=0; i<count; ++i) patch.points[i]=patch.points[i]-end*(static_cast<double>(i)/(count-1));
                }
                Vec3 center;
                for (const auto& point : patch.points) center=center+point;
                center=center*(1.0/count);
                double rg2=0.0;
                for (auto& point : patch.points) { point=point-center; rg2+=dot(point,point); }
                const double rg=std::sqrt(rg2/count);
                const double target=definition.targetRg.root() ? evaluate(definition.targetRg,{}) : 1.0;
                if (rg>0.0) for (auto& point : patch.points) point=point*(target/rg);
            } else if (definition.kind == pyseb::VisualizationGeometryKind::Curve) {
                patch.uCount=count;
                const double lower=evaluate(definition.uLower,{}), upper=evaluate(definition.uUpper,{});
                for (std::size_t i=0; i<count; ++i) {
                    const double u=lower+(upper-lower)*static_cast<double>(i)/(count-1);
                    const std::map<std::string,double> local{{"u",u},{"t",u}};
                    patch.points.push_back({evaluate(definition.coordinates[0],local),evaluate(definition.coordinates[1],local),evaluate(definition.coordinates[2],local)});
                }
            } else {
                patch.uCount=count; patch.vCount=count;
                const double u0=evaluate(definition.uLower,{}), u1=evaluate(definition.uUpper,{});
                const double v0=evaluate(definition.vLower,{}), v1=evaluate(definition.vUpper,{});
                for (std::size_t uIndex=0; uIndex<count; ++uIndex) {
                    const double u=u0+(u1-u0)*static_cast<double>(uIndex)/(count-1);
                    for (std::size_t vIndex=0; vIndex<count; ++vIndex) {
                        const double v=v0+(v1-v0)*static_cast<double>(vIndex)/(count-1);
                        const std::map<std::string,double> local{{"u",u},{"v",v}};
                        patch.points.push_back({evaluate(definition.coordinates[0],local),evaluate(definition.coordinates[1],local),evaluate(definition.coordinates[2],local)});
                    }
                }
            }
            instance.patches.emplace(geometry.first,std::move(patch));
        }

        for (const auto& reference : visualization.references) {
            if (reference.second.kind == "fixed") {
                instance.resolvedReferences[reference.first]={evaluate(reference.second.position[0],{}),evaluate(reference.second.position[1],{}),evaluate(reference.second.position[2],{})};
            } else if (reference.second.kind == "curve_fraction") {
                const auto patch=instance.patches.find(reference.second.geometry);
                if (patch == instance.patches.end()) throw SEBException("reference uses unknown geometry", "World::ExportUSD");
                instance.resolvedReferences[reference.first]=sampleCurveParameter(patch->second,evaluate(reference.second.fraction,{}));
            }
        }
        instances.emplace(name,std::move(instance));
    }

    std::map<std::string,Pose> poses;
    for (const auto& name : selected) poses[name]=Pose();
    poses.begin()->second.resolved=true;

    std::map<std::string,GeometryProxy> localProxies;
    for (const auto& instance : instances) localProxies.emplace(instance.first,localProxy(instance.second));
    std::map<std::string,GeometryProxy> placedProxies;
    const auto rootProxy=localProxies.find(poses.begin()->first);
    if (rootProxy != localProxies.end()) {
        placedProxies.emplace(rootProxy->first,transformProxy(rootProxy->second,poses.begin()->second));
    }

    auto resolveReference=[&](const LinkEndpoint& endpoint) -> Vec3 {
        auto instance=instances.find(endpoint.instance);
        if (instance == instances.end()) return Vec3();
        auto cached=instance->second.resolvedReferences.find(endpoint.reference);
        if (cached != instance->second.resolvedReferences.end()) return cached->second;
        const auto parts=splitReference(endpoint.reference);
        const auto base=instance->second.subunit->getDefinition().visualization.references.find(parts.first);
        if (base == instance->second.subunit->getDefinition().visualization.references.end()) throw SEBException("unknown visualization reference "+endpoint.reference,"World::ExportUSD");
        const auto& reference=base->second;
        if (reference.kind != "distributed") {
            const auto fixed=instance->second.resolvedReferences.find(parts.first);
            if (fixed == instance->second.resolvedReferences.end()) throw SEBException("unresolved visualization reference "+endpoint.reference,"World::ExportUSD");
            return fixed->second;
        }
        if (parts.second.empty()) throw SEBException("distributed link reference requires a #label", "World::ExportUSD");
        std::mt19937_64 generator(derivedSeed(options.seed,endpoint.instance+"/reference/"+endpoint.reference));
        std::string patchName=reference.geometry;
        if (reference.sampling == "weighted_mixture") {
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
            point=reference.sampling == "arc_length" ? sampleSurfaceBoundary(patch->second,unitRandom(generator)) : sampleSurfaceArea(patch->second,generator);
        } else {
            point=reference.sampling == "uniform_parameter" ? sampleCurveParameter(patch->second,unitRandom(generator)) : sampleCurveLength(patch->second,unitRandom(generator));
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
    if (options.layoutMode == USDLayoutMode::Readable) {
        output << "    custom string pyseb:layoutMode = \"readable\"\n";
        output << "    custom string pyseb:orientationSemantics = \"readability_optimized_free_rotation\"\n";
        output << "    custom string pyseb:overlapPolicy = \"best_effort_minimized\"\n";
        output << "    custom uint64 pyseb:orientationTrials = " << options.orientationTrials << "\n";
        output << "    custom uint64 pyseb:relaxationSweeps = " << options.relaxationSweeps << "\n";
        output << "    custom double pyseb:minimumClearance = " << number(options.minimumClearance) << "\n";
    } else {
        output << "    custom string pyseb:layoutMode = \"random\"\n";
        output << "    custom string pyseb:orientationSemantics = \"representative_free_rotation\"\n";
        output << "    custom string pyseb:overlapPolicy = \"permitted\"\n";
    }
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
            output << "        def PointInstancer \"cloud\" {\n            rel prototypes = [</"<<usdName(structure)<<"/"<<usdName(name)<<"/spherePrototype>]\n            point3f[] positions = [";
            for(const auto& scatterer:cloud->getScatterers())output<<"("<<number(scatterer.center.x)<<", "<<number(scatterer.center.y)<<", "<<number(scatterer.center.z)<<"), ";
            output<<"]\n            int[] protoIndices = [";for(std::size_t i=0;i<cloud->getScatterers().size();++i)output<<"0, ";output<<"]\n            float3[] scales = [";
            for(const auto& scatterer:cloud->getScatterers()){double radius=scatterer.radius>0.0?scatterer.radius:options.zeroRadiusMarkerSize;output<<"("<<number(radius)<<", "<<number(radius)<<", "<<number(radius)<<"), ";}
            output<<"]\n            custom float[] pyseb:radius = [";for(const auto& scatterer:cloud->getScatterers())output<<number(scatterer.radius)<<", ";output<<"]\n            custom float[] pyseb:beta = [";for(const auto& scatterer:cloud->getScatterers())output<<number(scatterer.beta)<<", ";output<<"]\n        }\n        def Sphere \"spherePrototype\" { float radius = 1 }\n    }\n";
            continue;
        }
        const auto instance=instances.find(name);
        if (instance == instances.end()) { output<<"    }\n"; continue; }
        for (const auto& reference : instance->second.resolvedReferences) {
            output<<"        def Xform \"ref_"<<usdName(reference.first)<<"\" {\n            visibility = \"invisible\"\n            custom string pyseb:reference = "<<quote(reference.first)<<"\n            double3 xformOp:translate = ("<<number(reference.second.x)<<", "<<number(reference.second.y)<<", "<<number(reference.second.z)<<")\n            uniform token[] xformOpOrder = [\"xformOp:translate\"]\n        }\n";
        }
        for (const auto& patch : instance->second.patches) {
            const bool surface=patch.second.kind==pyseb::VisualizationGeometryKind::Surface;
            output<<"        def "<<(surface?"Mesh":"BasisCurves")<<" \""<<usdName(patch.first)<<"\" {\n";
            if(surface){output<<"            int[] faceVertexCounts = [";for(std::size_t u=0;u+1<patch.second.uCount;++u)for(std::size_t v=0;v+1<patch.second.vCount;++v)output<<"4, ";output<<"]\n            int[] faceVertexIndices = [";for(std::size_t u=0;u+1<patch.second.uCount;++u)for(std::size_t v=0;v+1<patch.second.vCount;++v){std::size_t i=u*patch.second.vCount+v;output<<i<<", "<<i+1<<", "<<i+patch.second.vCount+1<<", "<<i+patch.second.vCount<<", ";}output<<"]\n";}
            else output<<"            int[] curveVertexCounts = ["<<patch.second.points.size()<<"]\n";
            output<<"            point3f[] points = [";for(const auto& point:patch.second.points)output<<"("<<number(point.x)<<", "<<number(point.y)<<", "<<number(point.z)<<"), ";output<<"]\n";
            output<<"            color3f[] primvars:displayColor = [ ("<<number(instance->second.color[0])<<", "<<number(instance->second.color[1])<<", "<<number(instance->second.color[2])<<") ]\n            float[] primvars:displayOpacity = ["<<number(instance->second.opacity)<<"]\n        }\n";
        }
        output<<"        custom string pyseb:model_id = "<<quote(instance->second.subunit->getDefinition().id)<<"\n        custom string pyseb:instance_path = "<<quote(name)<<"\n        custom string pyseb:tag = "<<quote(instance->second.tag)<<"\n        custom uint64 pyseb:seed = "<<options.seed<<"\n";
        for(const auto& parameter:instance->second.parameters)output<<"        custom double pyseb:param_"<<usdName(parameter.first)<<" = "<<number(parameter.second)<<"\n";
        output<<"    }\n";
    }
    for (std::size_t index=0; index<sceneLinks.size(); ++index) {
        const auto& link=sceneLinks[index];
        output<<"    rel pyseb:link_"<<index<<" = [</"<<usdName(structure)<<"/"<<usdName(link.first.instance)<<"/ref_"<<usdName(link.first.reference)<<">, </"<<usdName(structure)<<"/"<<usdName(link.second.instance)<<"/ref_"<<usdName(link.second.reference)<<">]\n";
    }
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
