// Geometry proxies and readable-layout scoring.
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

GeometryProxy localProxy(const AtomicSAXSCloud& cloud) {
    GeometryProxy proxy;
    const auto& scatterers = cloud.getScatterers();
    constexpr std::size_t maxPoints = 64;
    const std::size_t pointCount = std::min(maxPoints, scatterers.size());
    for (std::size_t index = 0; index < pointCount; ++index) {
        const std::size_t source = pointCount == 1
            ? 0 : index * (scatterers.size() - 1) / (pointCount - 1);
        const auto& center = scatterers[source].center;
        proxy.points.push_back({center.x, center.y, center.z});
    }
    if (!proxy.points.empty()) {
        for (const Vec3& point : proxy.points) proxy.centroid = proxy.centroid + point;
        proxy.centroid = proxy.centroid *
            (1.0 / static_cast<double>(proxy.points.size()));
        for (const auto& scatterer : scatterers) {
            const Vec3 center{
                scatterer.center.x,
                scatterer.center.y,
                scatterer.center.z};
            proxy.radius = std::max(
                proxy.radius,
                length(center - proxy.centroid) +
                    atomicDisplayedRadius(scatterer, 0.1));
        }
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

} // namespace usd_detail
} // namespace pyseb
