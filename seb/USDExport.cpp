// Scene assembly and export orchestration.
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

namespace pyseb {
namespace usd_detail {

bool parseEndpoint(const std::string& value, LinkEndpoint& endpoint) {
    const std::size_t dot = value.find_last_of('.');
    if (dot == std::string::npos || dot + 1 == value.size()) return false;
    std::string path = value.substr(0, dot);
    endpoint.instance = path;
    endpoint.reference = value.substr(dot + 1);
    return true;
}

std::vector<std::string> splitPath(const std::string& value, char separator) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    while (start <= value.size()) {
        const std::size_t end = value.find(separator, start);
        const std::size_t length = end == std::string::npos
            ? value.size() - start : end - start;
        if (length == 0) throw SEBException("empty path component", "World::ExportUSD");
        parts.push_back(value.substr(start, length));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return parts;
}

std::string joinPath(const std::vector<std::string>& parts) {
    std::string result;
    for (const auto& part : parts) {
        if (!result.empty()) result += "/";
        result += part;
    }
    return result;
}

std::string prefixedEndpointPath(
    const std::string& prefix,
    const std::string& endpointPath) {
    std::vector<std::string> parts;
    if (!prefix.empty()) parts = splitPath(prefix, '/');
    const auto endpointParts = splitPath(endpointPath, ':');
    parts.insert(parts.end(), endpointParts.begin(), endpointParts.end());
    return joinPath(parts);
}

OccurrenceTree buildOccurrenceTree(
    const std::string& structure,
    const std::map<GraphID, std::list<std::string>>& subGraphs,
    const std::map<std::string, ABSSubUnit*>& nameCatalog,
    const std::list<LinkPair>& links) {
    OccurrenceTree result;

    const auto endpointObject = [](const std::string& endpoint) {
        const std::size_t colon = endpoint.find(':');
        const std::size_t period = endpoint.find('.');
        std::size_t cut = std::string::npos;
        if (colon != std::string::npos) cut = colon;
        if (period != std::string::npos) cut = std::min(cut, period);
        return cut == std::string::npos ? endpoint : endpoint.substr(0, cut);
    };

    std::map<std::string, GraphID> graphByObject;
    for (const auto& graph : subGraphs) {
        for (const auto& member : graph.second) {
            graphByObject[member] = graph.first;
        }
    }

    std::map<GraphID, std::vector<LinkPair>> linksByGraph;
    for (const auto& link : links) {
        const auto first = graphByObject.find(endpointObject(link.first));
        const auto second = graphByObject.find(endpointObject(link.second));
        if (first != graphByObject.end()
            && second != graphByObject.end()
            && first->second == second->second) {
            linksByGraph[first->second].push_back(link);
        }
    }

    const auto pathWithChild = [](const std::string& prefix,
                                  const std::string& child) {
        return prefix.empty() ? child : prefix + "/" + child;
    };
    const auto expandEndpoint = [&](const std::string& prefix,
                                    const std::string& endpoint) {
        LinkEndpoint parsed;
        if (!parseEndpoint(endpoint, parsed)) return LinkEndpoint{};
        return LinkEndpoint{
            prefixedEndpointPath(prefix, parsed.instance),
            parsed.reference};
    };

    std::function<std::vector<OccurrenceNode>(
        GraphID,
        const std::string&,
        std::set<GraphID>&)> expandGraph;
    expandGraph = [&](GraphID id,
                      const std::string& prefix,
                      std::set<GraphID>& active) {
        if (!active.insert(id).second) {
            throw SEBException(
                "cyclic nested structure graph",
                "World::ExportUSD");
        }

        std::vector<OccurrenceNode> nodes;
        const auto graph = subGraphs.find(id);
        if (graph == subGraphs.end()) {
            active.erase(id);
            return nodes;
        }

        for (const auto& childName : graph->second) {
            const auto child = nameCatalog.find(childName);
            if (child == nameCatalog.end()) continue;

            OccurrenceNode node;
            node.path = pathWithChild(prefix, childName);
            node.sourceName = childName;
            const Structure* nested =
                dynamic_cast<const Structure*>(child->second);
            node.structure = nested != nullptr;
            if (nested != nullptr) {
                node.children = expandGraph(
                    nested->getGraphID(), node.path, active);
            }
            nodes.push_back(std::move(node));
        }

        const auto graphLinks = linksByGraph.find(id);
        if (graphLinks != linksByGraph.end()) {
            for (const auto& link : graphLinks->second) {
                SceneLink expanded;
                expanded.first = expandEndpoint(prefix, link.first);
                expanded.second = expandEndpoint(prefix, link.second);
                expanded.key = expanded.first.instance + "."
                    + expanded.first.reference + "<->"
                    + expanded.second.instance + "."
                    + expanded.second.reference;
                result.links.push_back(std::move(expanded));
            }
        }

        active.erase(id);
        return nodes;
    };

    const auto root = nameCatalog.find(structure);
    if (root != nameCatalog.end()) {
        const Structure* rootStructure =
            dynamic_cast<const Structure*>(root->second);
        if (rootStructure != nullptr) {
            std::set<GraphID> active;
            result.roots = expandGraph(
                rootStructure->getGraphID(), std::string(), active);
        } else if (dynamic_cast<const SubUnit*>(root->second) != nullptr) {
            OccurrenceNode node;
            node.path = structure;
            node.sourceName = structure;
            result.roots.push_back(std::move(node));
        }
        return result;
    }

    for (const auto& item : nameCatalog) {
        if (dynamic_cast<const SubUnit*>(item.second) != nullptr) {
            OccurrenceNode node;
            node.path = item.first;
            node.sourceName = item.first;
            result.roots.push_back(std::move(node));
        }
    }
    for (const auto& link : links) {
        SceneLink expanded;
        if (!parseEndpoint(link.first, expanded.first)
            || !parseEndpoint(link.second, expanded.second)) {
            continue;
        }
        expanded.key = link.first + "<->" + link.second;
        result.links.push_back(std::move(expanded));
    }
    return result;
}

void collectLeafOccurrences(
    const std::vector<OccurrenceNode>& nodes,
    std::set<std::string>& paths,
    std::map<std::string, std::string>& sourceNames) {
    for (const auto& node : nodes) {
        if (node.structure) {
            collectLeafOccurrences(node.children, paths, sourceNames);
            continue;
        }
        paths.insert(node.path);
        sourceNames[node.path] = node.sourceName;
    }
}

std::string usdPrimPath(
    const std::string& root,
    const std::string& occurrence,
    const std::string& child) {
    std::vector<std::string> parts{usdName(root)};
    if (!occurrence.empty()) {
        const auto occurrenceParts = splitPath(occurrence, '/');
        parts.insert(parts.end(), occurrenceParts.begin(), occurrenceParts.end());
    }
    if (!child.empty()) parts.push_back(child);
    std::string result;
    for (const auto& part : parts) result += "/" + usdName(part);
    return result;
}

} // namespace usd_detail
} // namespace pyseb

using namespace pyseb::usd_detail;

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

    const OccurrenceTree occurrenceTree = buildOccurrenceTree(
        structure, subGraphs, nameCatalog, links);
    std::set<std::string> selected;
    std::map<std::string, std::string> sourceNames;
    collectLeafOccurrences(
        occurrenceTree.roots, selected, sourceNames);
    if (selected.empty()) throw SEBException("structure contains no visualizable leaf instances", "World::ExportUSD");

    std::vector<SceneLink> sceneLinks;
    for (const auto& link : occurrenceTree.links) {
        if (selected.count(link.first.instance) && selected.count(link.second.instance)) {
            sceneLinks.push_back(link);
        }
    }

    const auto overrideValue = [&](const auto& overrides,
                                   const std::string& occurrence,
                                   const std::string& source) {
        auto found = overrides.find(occurrence);
        if (found != overrides.end()) return &found->second;
        found = overrides.find(structure + "/" + occurrence);
        if (found != overrides.end()) return &found->second;
        found = overrides.find(source);
        return found == overrides.end() ? nullptr : &found->second;
    };

    std::map<std::string,InstanceScene> instances;
    for (const auto& name : selected) {
        const auto source = sourceNames.find(name);
        if (source == sourceNames.end()) continue;
        const auto catalog=nameCatalog.find(source->second);
        const SubUnit* sub=catalog == nameCatalog.end() ? nullptr : dynamic_cast<const SubUnit*>(catalog->second);
        const FileDefinedSubunit* file=dynamic_cast<const FileDefinedSubunit*>(sub);
        if (!file || !file->getDefinition().visualization.present) continue;
        InstanceScene instance;
        instance.subunit=file;
        instance.name=name;
        instance.sourceName=source->second;
        instance.tag=const_cast<SubUnit*>(sub)->getTag();
        instance.color=file->getDefinition().visualization.color;
        instance.opacity=file->getDefinition().visualization.opacity;
        const auto colorOverride=overrideValue(options.colorOverrides, name, instance.sourceName);
        if (colorOverride != nullptr) instance.color=*colorOverride;
        const auto opacityOverride=overrideValue(options.opacityOverrides, name, instance.sourceName);
        if (opacityOverride != nullptr) instance.opacity=*opacityOverride;
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
        const auto source = sourceNames.find(name);
        if (source == sourceNames.end()) continue;
        const auto catalog=nameCatalog.find(source->second);
        const DebyeSphereCloud* cloud=catalog == nameCatalog.end()
            ? nullptr : dynamic_cast<const DebyeSphereCloud*>(catalog->second);
        const AtomicSAXSCloud* atomicCloud=catalog == nameCatalog.end()
            ? nullptr : dynamic_cast<const AtomicSAXSCloud*>(catalog->second);
        if (!cloud && !atomicCloud) continue;
        const auto references = cloud
            ? cloud->getReferenceCoordinates()
            : atomicCloud->getReferenceCoordinates();
        for (const auto& reference : references) {
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
        const auto source = sourceNames.find(name);
        if (source == sourceNames.end()) continue;
        const auto catalog=nameCatalog.find(source->second);
        const DebyeSphereCloud* cloud=catalog == nameCatalog.end()
            ? nullptr : dynamic_cast<const DebyeSphereCloud*>(catalog->second);
        if (cloud) localProxies.emplace(name,localProxy(*cloud));
        const AtomicSAXSCloud* atomicCloud=catalog == nameCatalog.end()
            ? nullptr : dynamic_cast<const AtomicSAXSCloud*>(catalog->second);
        if (atomicCloud) localProxies.emplace(name,localProxy(*atomicCloud));
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

    std::stringstream output;
    output.imbue(std::locale::classic());
    output << std::setprecision(17);
    output << "#usda 1.0\n(\n    upAxis = \"Z\"\n    metersPerUnit = " << number(metersPerUnit) << "\n)\n\n";
    output << "def Xform \"" << usdName(structure) << "\" {\n";
    writeLayoutMetadata(output, options);
    std::function<void(
        const std::vector<OccurrenceNode>&,
        const std::string&)> writeNodes;
    writeNodes = [&](const std::vector<OccurrenceNode>& nodes,
                     const std::string& indent) {
        for (const auto& node : nodes) {
            if (node.structure) {
                output << indent << "def Xform \""
                       << usdName(node.sourceName) << "\" {\n"
                       << indent << "    custom string pyseb:structure_path = "
                       << quote(node.path) << "\n";
                writeNodes(node.children, indent + "    ");
                output << indent << "}\n";
                continue;
            }
            const auto source = sourceNames.find(node.path);
            if (source == sourceNames.end()) continue;
            const auto catalog = nameCatalog.find(source->second);
            const SubUnit* sub = catalog == nameCatalog.end()
                ? nullptr : dynamic_cast<const SubUnit*>(catalog->second);
            if (!sub) continue;
            const Pose& pose = poses[node.path];
            output << indent << "def Xform \""
                   << usdName(node.sourceName) << "\" {\n";
            output << indent << "    double3 xformOp:translate = ("
                   << number(pose.translation.x) << ", "
                   << number(pose.translation.y) << ", "
                   << number(pose.translation.z) << ")\n";
            output << indent << "    quatd xformOp:orient = ("
                   << number(pose.rotation.w) << ", "
                   << number(pose.rotation.x) << ", "
                   << number(pose.rotation.y) << ", "
                   << number(pose.rotation.z) << ")\n";
            output << indent << "    uniform token[] xformOpOrder = "
                      "[\"xformOp:translate\", \"xformOp:orient\"]\n";
            const DebyeSphereCloud* cloud = dynamic_cast<const DebyeSphereCloud*>(sub);
            if (cloud) {
                const auto references = numericalReferences.find(node.path);
                const std::map<std::string, Vec3> emptyReferences;
                writeCloudPrim(
                    output,
                    structure,
                    node.path,
                    *cloud,
                    options,
                    references == numericalReferences.end()
                        ? emptyReferences
                        : references->second);
            } else {
                const AtomicSAXSCloud* atomicCloud =
                    dynamic_cast<const AtomicSAXSCloud*>(sub);
                if (atomicCloud) {
                    const auto references = numericalReferences.find(node.path);
                    const std::map<std::string, Vec3> emptyReferences;
                    writeAtomicCloudPrim(
                        output,
                        structure,
                        node.path,
                        *atomicCloud,
                        options,
                        references == numericalReferences.end()
                            ? emptyReferences
                            : references->second);
                    continue;
                }
                const auto instance = instances.find(node.path);
                if (instance != instances.end()) {
                    writeFileInstanceBody(output, instance->second, options.seed);
                } else {
                    output << indent << "}\n";
                }
            }
        }
    };
    writeNodes(occurrenceTree.roots, "    ");
    writeLinkRelationships(output, structure, sceneLinks);
    output<<"}\n";

    const std::string temporary=path+".tmp";
    {
        std::ofstream file(temporary,std::ios::binary);
        if(!file)throw SEBException("unable to open USD output","World::ExportUSD");
        // Transfer the existing buffer instead of allocating a second complete
        // USDA string. A bidirectional stringstream exposes its get area here.
        file<<output.rdbuf();
        if(!file)throw SEBException("failed writing USD output","World::ExportUSD");
    }
    if(std::rename(temporary.c_str(),path.c_str())!=0){std::remove(temporary.c_str());throw SEBException("unable to commit USD output","World::ExportUSD");}
}
