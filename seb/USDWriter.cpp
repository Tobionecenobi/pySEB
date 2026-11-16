// USD text serialization and cloud envelope meshes.
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

std::string quote(const std::string& value) {
    std::string result = "\"";
    for (char ch : value) {
        if (ch == '\\' || ch == '\"') result += '\\';
        result += ch;
    }
    return result + "\"";
}

// Validate each floating-point field, then format it directly in the export's
// classic-locale, 17-digit stream. Avoid a temporary stream/string per coordinate.
double number(double value) {
    if (!std::isfinite(value)) {
        throw SEBException("non-finite visualization value", "World::ExportUSD");
    }
    return value;
}

std::string usdName(const std::string& value) {
    std::string result;
    for (char ch : value) {
        result += (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') ? ch : '_';
    }
    if (result.empty() || std::isdigit(static_cast<unsigned char>(result.front()))) result = "p_" + result;
    return result;
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

// USDA serialization: keep formatting and metadata in one place.
void writePose(std::ostream& output, const Pose& pose) {
    output << "        double3 xformOp:translate = ("
           << number(pose.translation.x) << ", "
           << number(pose.translation.y) << ", "
           << number(pose.translation.z) << ")\n";
    output << "        quatd xformOp:orient = ("
           << number(pose.rotation.w) << ", "
           << number(pose.rotation.x) << ", "
           << number(pose.rotation.y) << ", "
           << number(pose.rotation.z) << ")\n";
    output << "        uniform token[] xformOpOrder = "
              "[\"xformOp:translate\", \"xformOp:orient\"]\n";
}

void writeReferencePrim(
    std::ostream& output,
    const std::string& referenceName,
    const Vec3& position) {
    output << "        def Xform \"ref_" << usdName(referenceName)
           << "\" {\n"
           << "            token visibility = \"invisible\"\n"
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
           << "            rel prototypes = [<"
           << usdPrimPath(structure, name, "spherePrototype") << ">]\n"
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
           << "        def Sphere \"spherePrototype\" {\n"
           << "            float radius = 1\n"
           << "        }\n";

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

void writeAtomicCloudPrim(
    std::ostream& output,
    const std::string& structure,
    const std::string& name,
    const AtomicSAXSCloud& cloud,
    const USDExportOptions& options,
    const std::map<std::string, Vec3>& references)
{
    std::array<double, 3> color{{0.45, 0.65, 0.85}};
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
           << "            rel prototypes = [<"
           << usdPrimPath(structure, name, "spherePrototype") << ">]\n"
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
        const double radius = atomicDisplayedRadius(
            scatterer, options.zeroRadiusMarkerSize);
        output << "(" << number(radius) << ", "
               << number(radius) << ", " << number(radius) << "), ";
    }
    output << "]\n            custom float[] pyseb:excluded_volume = [";
    for (const auto& scatterer : scatterers) {
        output << number(scatterer.excludedVolume) << ", ";
    }
    output << "]\n            custom float[] pyseb:occupancy = [";
    for (const auto& scatterer : scatterers) {
        output << number(scatterer.occupancy) << ", ";
    }
    output << "]\n            custom int[] pyseb:implicit_hydrogens = [";
    for (const auto& scatterer : scatterers) {
        output << scatterer.implicitHydrogens << ", ";
    }
    output << "]\n            custom float[] pyseb:implicit_hydrogen_distance = [";
    for (const auto& scatterer : scatterers) {
        output << number(scatterer.implicitHydrogenDistance) << ", ";
    }
    output << "]\n            custom float[] pyseb:hydration_volume = [";
    for (const auto& scatterer : scatterers) {
        output << number(scatterer.hydrationVolume) << ", ";
    }
    output << "]\n            custom string[] pyseb:element = [";
    for (const auto& scatterer : scatterers) {
        output << quote(scatterer.element) << ", ";
    }
    output << "]\n            color3f[] primvars:displayColor = [";
    for (std::size_t index = 0; index < scatterers.size(); ++index) {
        output << "(" << number(color[0]) << ", "
               << number(color[1]) << ", " << number(color[2]) << "), ";
    }
    output << "]\n            float[] primvars:displayOpacity = [";
    for (std::size_t index = 0; index < scatterers.size(); ++index) {
        output << number(opacity) << ", ";
    }
    output << "]\n        }\n"
           << "        def Sphere \"spherePrototype\" {\n"
           << "            float radius = 1\n"
           << "        }\n";
    for (const auto& reference : references) {
        writeReferencePrim(output, reference.first, reference.second);
    }
    const AtomicSAXSOptions saxsOptions = cloud.getOptions();
    const char* waterModel = "neutral_water";
    if (saxsOptions.hydrationWaterModel == AtomicSAXSHydrationWaterModel::BondedWater)
        waterModel = "bonded_water";
    else if (saxsOptions.hydrationWaterModel == AtomicSAXSHydrationWaterModel::FoXSReciprocal)
        waterModel = "foxs_reciprocal";
    else if (saxsOptions.hydrationWaterModel == AtomicSAXSHydrationWaterModel::FoXSFast)
        waterModel = "foxs_fast";
    output << "        custom string pyseb:model_id = "
           << quote("pyseb/AtomicSAXSCloud") << "\n"
           << "        custom string pyseb:instance_path = " << quote(name)
           << "\n"
           << "        custom double pyseb:solvent_electron_density = "
           << number(saxsOptions.solventElectronDensity) << "\n"
           << "        custom double pyseb:excluded_volume_scale = "
           << number(saxsOptions.excludedVolumeScale) << "\n"
           << "        custom string pyseb:hydration_model = "
           << quote(saxsOptions.hydrationModel == AtomicSAXSHydrationModel::SpatialShell
                    ? "spatial_shell" : "atom_centered") << "\n"
           << "        custom double pyseb:hydration_shell_contrast = "
           << number(saxsOptions.hydrationShellContrast) << "\n"
           << "        custom string pyseb:hydration_water_model = " << quote(waterModel) << "\n"
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
           << "        custom string pyseb:source_name = "
           << quote(instance.sourceName) << "\n"
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
        output << "    rel pyseb:link_" << index << " = [<"
               << usdPrimPath(structure, link.first.instance,
                              "ref_" + usdName(link.first.reference))
               << ">, <"
               << usdPrimPath(structure, link.second.instance,
                              "ref_" + usdName(link.second.reference))
               << ">]\n";
    }
}

} // namespace usd_detail
} // namespace pyseb
