#include "AtomicSAXSCloudBuilder.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>

#include "freesasa.h"

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

std::string normalizeIdentifier(const std::string& value)
{
    std::size_t first = 0;
    while (first < value.size() &&
           std::isspace(static_cast<unsigned char>(value[first]))) {
        ++first;
    }
    std::size_t last = value.size();
    while (last > first &&
           std::isspace(static_cast<unsigned char>(value[last - 1]))) {
        --last;
    }
    std::string result = value.substr(first, last - first);
    std::transform(
        result.begin(), result.end(), result.begin(),
        [](unsigned char character) {
            return static_cast<char>(std::toupper(character));
        }
    );
    return result;
}

bool isStandardProteinResidue(const std::string& residue)
{
    static const std::set<std::string> standard = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL"
    };
    return standard.find(normalizeIdentifier(residue)) != standard.end();
}

void validateSASAOptions(const AtomicSAXSCloudBuildOptions& options)
{
    if (!std::isfinite(options.sasaProbeRadius) ||
        options.sasaProbeRadius <= 0.0) {
        throw SEBException(
            "SASA probe radius must be finite and positive",
            "AtomicSAXSCloudBuilder::Build()"
        );
    }
    if (options.sasaPointCount < 1) {
        throw SEBException(
            "SASA point count must be positive",
            "AtomicSAXSCloudBuilder::Build()"
        );
    }
}

void calculateSolventAccessibility(
    const std::vector<AtomRecord>& atoms,
    const std::vector<AtomicSAXSParameters>& resolvedParameters,
    const AtomicSAXSCloudBuildOptions& options,
    std::vector<AtomicSAXSScatterer>& scatterers)
{
    validateSASAOptions(options);

    std::vector<std::size_t> activeIndices;
    std::vector<double> xyz;
    std::vector<double> radii;
    activeIndices.reserve(atoms.size());
    xyz.reserve(3 * atoms.size());
    radii.reserve(atoms.size());

    for (std::size_t index = 0; index < atoms.size(); ++index) {
        if (atoms[index].occupancy <= 0.0) {
            scatterers[index].solventAccessibility = 0.0;
            continue;
        }

        double radius = -1.0;
        if (!options.sasaUseProfileRadii && isStandardProteinResidue(atoms[index].residueName)) {
            const std::string residue =
                normalizeIdentifier(atoms[index].residueName);
            const std::string atomName =
                normalizeIdentifier(atoms[index].atomName);
            radius = freesasa_classifier_radius(
                &freesasa_protor_classifier,
                residue.c_str(),
                atomName.c_str()
            );
        }
        if (!(radius > 0.0)) {
            radius = resolvedParameters[index].sasaRadius;
        }
        if (!std::isfinite(radius) || radius <= 0.0) {
            throw SEBException(
                "No positive SASA radius for atom serial " +
                    std::to_string(atoms[index].serial) + " (" +
                    atoms[index].residueName + ":" + atoms[index].atomName +
                    ", element " + atoms[index].element + ")",
                "AtomicSAXSCloudBuilder::Build()"
            );
        }

        activeIndices.push_back(index);
        xyz.push_back(atoms[index].coordinate.x);
        xyz.push_back(atoms[index].coordinate.y);
        xyz.push_back(atoms[index].coordinate.z);
        radii.push_back(radius);
    }

    if (activeIndices.empty()) {
        return;
    }

    freesasa_parameters parameters = freesasa_default_parameters;
    parameters.alg = FREESASA_SHRAKE_RUPLEY;
    parameters.probe_radius = options.sasaProbeRadius;
    parameters.shrake_rupley_n_points = options.sasaPointCount;
    parameters.n_threads = 1;

    freesasa_result* result = freesasa_calc_coord(
        xyz.data(), radii.data(), static_cast<int>(activeIndices.size()),
        &parameters
    );
    if (result == nullptr) {
        throw SEBException(
            "FreeSASA failed to calculate solvent accessibility",
            "AtomicSAXSCloudBuilder::Build()"
        );
    }

    try {
        for (std::size_t active = 0; active < activeIndices.size(); ++active) {
            const double expandedRadius =
                radii[active] + options.sasaProbeRadius;
            const double isolatedArea =
                4.0 * kPi * expandedRadius * expandedRadius;
            double fraction = result->sasa[active] / isolatedArea;
            const double tolerance = 1e-10;
            if (!std::isfinite(fraction) || fraction < -tolerance ||
                fraction > 1.0 + tolerance) {
                throw SEBException(
                    "FreeSASA returned an invalid exposure fraction for atom serial " +
                        std::to_string(atoms[activeIndices[active]].serial),
                    "AtomicSAXSCloudBuilder::Build()"
                );
            }
            fraction = std::max(0.0, std::min(1.0, fraction));
            scatterers[activeIndices[active]].solventAccessibility = fraction;
        }
    } catch (...) {
        freesasa_result_free(result);
        throw;
    }
    freesasa_result_free(result);
}

void appendHydrationShell(
    const std::vector<AtomRecord>& atoms,
    const std::vector<AtomicSAXSParameters>& parameters,
    const AtomicSAXSCloudBuildOptions& options,
    std::vector<AtomicSAXSScatterer>& scatterers)
{
    const double inner = options.hydrationShellInnerRadius;
    const double thickness = options.hydrationShellThickness;
    const double step = options.hydrationShellSpacing;
    const bool pepsi = options.hydrationShellConvention == AtomicSAXSHydrationShellConvention::Pepsi3;
    if (!pepsi && options.hydrationShellConvention != AtomicSAXSHydrationShellConvention::VolumeDensity)
        throw SEBException("Unknown shell convention", "AtomicSAXSCloudBuilder hydration shell");
    if (!std::isfinite(inner) || inner <= 0.0 ||
        !std::isfinite(thickness) || thickness <= 0.0 ||
        !std::isfinite(step) || step <= 0.0 || !std::isfinite(inner + thickness)) {
        throw SEBException("Shell radii, thickness and spacing must be finite and positive",
                           "AtomicSAXSCloudBuilder hydration shell");
    }
    std::vector<CartesianPoint3D> centers;
    std::vector<double> innerRadii;
    CartesianPoint3D origin;
    for (std::size_t index=0; index<atoms.size(); ++index) {
        const auto& atom = atoms[index];
        if (atom.occupancy <= 0.0) continue;
        centers.push_back(atom.coordinate);
        if (pepsi && !(parameters[index].sasaRadius > 0.0)) throw SEBException(
            "Pepsi shell convention requires positive surface radii",
            "AtomicSAXSCloudBuilder hydration shell");
        innerRadii.push_back(inner + (pepsi ? parameters[index].sasaRadius : 0.0));
        origin.x += atom.coordinate.x;
        origin.y += atom.coordinate.y;
        origin.z += atom.coordinate.z;
    }
    if (centers.empty()) throw SEBException(
        "Spatial hydration requires at least one positive-occupancy atom",
        "AtomicSAXSCloudBuilder hydration shell");
    origin.x /= centers.size(); origin.y /= centers.size(); origin.z /= centers.size();
    const double outer = inner + thickness;
    if (pepsi) {
        // Reference grid is anchored at the minimum atomic coordinates,
        // padded by 12 A, and sampled at cell centers (half a cell offset).
        origin = centers.front();
        for (const auto& center : centers) {
            origin.x = std::min(origin.x, center.x);
            origin.y = std::min(origin.y, center.y);
            origin.z = std::min(origin.z, center.z);
        }
        const double padding = std::ceil(std::max(12.0,outer)/step)*step;
        origin.x += -padding+step/2;
        origin.y += -padding+step/2;
        origin.z += -padding+step/2;
    }
    double lower[3] = {0.0, 0.0, 0.0}, upper[3] = {0.0, 0.0, 0.0};
    for (const auto& center : centers) {
        const double xyz[3] = {center.x-origin.x, center.y-origin.y, center.z-origin.z};
        for (int axis = 0; axis < 3; ++axis) {
            lower[axis] = std::min(lower[axis], xyz[axis]-outer);
            upper[axis] = std::max(upper[axis], xyz[axis]+outer);
        }
    }
    int first[3], last[3];
    double cells = 1.0;
    for (int axis = 0; axis < 3; ++axis) {
        const double lo = std::floor(lower[axis]/step);
        const double hi = std::ceil(upper[axis]/step);
        if (lo < -1000000 || hi > 1000000) throw SEBException(
            "Hydration grid is too large; increase spacing",
            "AtomicSAXSCloudBuilder hydration shell");
        first[axis] = static_cast<int>(lo);
        last[axis] = static_cast<int>(hi);
        cells *= hi-lo+1;
    }
    if (cells > 2000000) throw SEBException(
        "Hydration grid exceeds 2000000 candidate cells; increase spacing",
        "AtomicSAXSCloudBuilder hydration shell");
    const std::size_t initialSize = scatterers.size();
    for (int x=first[0]; x<=last[0]; ++x)
        for (int y=first[1]; y<=last[1]; ++y)
            for (int z=first[2]; z<=last[2]; ++z) {
                const CartesianPoint3D point(origin.x+x*step, origin.y+y*step, origin.z+z*step);
                bool inside = false, nearby = false;
                for (std::size_t index=0; index<centers.size(); ++index) {
                    const auto& center = centers[index];
                    const double dx=point.x-center.x, dy=point.y-center.y, dz=point.z-center.z;
                    const double d2=dx*dx+dy*dy+dz*dz;
                    if (d2 < innerRadii[index]*innerRadii[index]) { inside=true; break; }
                    if (d2 < outer*outer) nearby=true;
                }
                if (inside || !nearby) continue;
                AtomicSAXSScatterer pointScatterer;
                pointScatterer.center = point;
                // Keep geometry independent of solvent so density updates agree
                // with rebuilding. Apply Pepsi's density factor during evaluation.
                pointScatterer.hydrationVolume = step*step*step;
                pointScatterer.hydrationUsesSolventDensity = pepsi;
                pointScatterer.element = "water";
                scatterers.push_back(pointScatterer);
                if (scatterers.size()-initialSize > 10000) throw SEBException(
                    "Shell exceeds 10000 points for the pair-distance cache; increase spacing",
                    "AtomicSAXSCloudBuilder hydration shell");
            }
    if (scatterers.size() == initialSize) throw SEBException(
        "Hydration shell grid is empty; decrease spacing",
        "AtomicSAXSCloudBuilder hydration shell");
}

} // namespace

std::unique_ptr<AtomicSAXSCloud> AtomicSAXSCloudBuilder::Build(
    const std::vector<AtomRecord>& atoms,
    const AtomicSAXSParameterProfile& profile,
    const AtomicSAXSOptions& saxsOptions,
    const AtomicSAXSCloudBuildOptions& buildOptions)
{
    if (atoms.empty()) {
        throw SEBException(
            "Cannot build an atomic SAXS cloud from an empty atom list",
            "AtomicSAXSCloudBuilder::Build()"
        );
    }

    validateSASAOptions(buildOptions);

    std::vector<AtomicSAXSScatterer> scatterers;
    scatterers.reserve(atoms.size());
    std::vector<AtomicSAXSParameters> resolvedParameters;
    resolvedParameters.reserve(atoms.size());
    std::map<int, CartesianPoint3D> coordinatesBySerial;
    std::set<int> duplicateSerials;

    for (const AtomRecord& atom : atoms) {
        if (!std::isfinite(atom.coordinate.x) ||
            !std::isfinite(atom.coordinate.y) ||
            !std::isfinite(atom.coordinate.z)) {
            throw SEBException(
                "Non-finite coordinate for atom serial " +
                    std::to_string(atom.serial),
                "AtomicSAXSCloudBuilder::Build()"
            );
        }
        if (!std::isfinite(atom.occupancy) || atom.occupancy < 0.0) {
            throw SEBException(
                "Invalid occupancy for atom serial " +
                    std::to_string(atom.serial),
                "AtomicSAXSCloudBuilder::Build()"
            );
        }
        const AtomicSAXSParameters parameters = profile.Resolve(atom);
        resolvedParameters.push_back(parameters);
        const double occupancy = buildOptions.scaleByOccupancy
            ? atom.occupancy : 1.0;
        scatterers.emplace_back(
            atom.coordinate,
            parameters.vacuum,
            !buildOptions.useImplicitHydrogens && parameters.atomicExcludedVolume >= 0.0
                ? parameters.atomicExcludedVolume : parameters.excludedVolume,
            occupancy,
            atom.temperatureFactor,
            0.0,
            buildOptions.useImplicitHydrogens
                ? parameters.implicitHydrogens : 0,
            atom.element
        );
        scatterers.back().implicitHydrogenDistance =
            buildOptions.useImplicitHydrogens
                ? parameters.implicitHydrogenDistance : 0.0;

        if (coordinatesBySerial.find(atom.serial) != coordinatesBySerial.end()) {
            duplicateSerials.insert(atom.serial);
        } else {
            coordinatesBySerial[atom.serial] = atom.coordinate;
        }
    }

    if (buildOptions.calculateSolventAccessibility) {
        calculateSolventAccessibility(
            atoms, resolvedParameters, buildOptions, scatterers
        );
    }

    if (saxsOptions.hydrationModel == AtomicSAXSHydrationModel::SpatialShell) {
        appendHydrationShell(atoms, resolvedParameters, buildOptions, scatterers);
    }
    std::unique_ptr<AtomicSAXSCloud> cloud(
        new AtomicSAXSCloud(scatterers, saxsOptions)
    );

    for (const auto& reference : buildOptions.referenceAtomSerials) {
        if (reference.first == "center") {
            throw SEBException(
                "The reference name 'center' is reserved by AtomicSAXSCloud",
                "AtomicSAXSCloudBuilder::Build()"
            );
        }
        if (duplicateSerials.find(reference.second) != duplicateSerials.end()) {
            throw SEBException(
                "Atom serial " + std::to_string(reference.second) +
                    " is ambiguous and cannot define reference " +
                    reference.first,
                "AtomicSAXSCloudBuilder::Build()"
            );
        }
        const auto coordinate = coordinatesBySerial.find(reference.second);
        if (coordinate == coordinatesBySerial.end()) {
            throw SEBException(
                "Atom serial " + std::to_string(reference.second) +
                    " was not found for reference " + reference.first,
                "AtomicSAXSCloudBuilder::Build()"
            );
        }
        cloud->addReferencePoint(reference.first, coordinate->second);
    }

    return cloud;
}

std::unique_ptr<AtomicSAXSCloud> AtomicSAXSStructureLoader::Load(
    const std::string& filename,
    const StructureParser& parser,
    const AtomicSAXSParameterProfile& profile,
    const StructureParseOptions& parseOptions,
    const AtomicSAXSOptions& saxsOptions,
    const AtomicSAXSCloudBuildOptions& buildOptions)
{
    return AtomicSAXSCloudBuilder::Build(
        parser.ParseFile(filename, parseOptions),
        profile,
        saxsOptions,
        buildOptions
    );
}

std::unique_ptr<AtomicSAXSCloud> AtomicSAXSStructureLoader::LoadPDB(
    const std::string& filename,
    const AtomicSAXSParameterProfile& profile,
    const StructureParseOptions& parseOptions,
    const AtomicSAXSOptions& saxsOptions,
    const AtomicSAXSCloudBuildOptions& buildOptions)
{
    const PDBParser parser;
    return Load(
        filename,
        parser,
        profile,
        parseOptions,
        saxsOptions,
        buildOptions
    );
}
