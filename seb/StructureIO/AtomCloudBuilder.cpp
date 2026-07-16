#include "AtomCloudBuilder.hpp"

#include <cmath>
#include <set>

std::unique_ptr<DebyeSphereCloud> AtomCloudBuilder::Build(
    const std::vector<AtomRecord>& atoms,
    const AtomParameterProfile& profile,
    const AtomCloudBuildOptions& options)
{
    if (atoms.empty()) {
        throw SEBException(
            "Cannot build a Debye sphere cloud from an empty atom list",
            "AtomCloudBuilder::Build()"
        );
    }

    std::vector<SphereScatterer> scatterers;
    scatterers.reserve(atoms.size());
    std::map<int, CartesianPoint3D> coordinatesBySerial;
    std::set<int> duplicateSerials;

    for (const AtomRecord& atom : atoms) {
        const AtomScattererParameters parameters = profile.Resolve(atom);
        const double occupancyScale =
            options.scaleBetaByOccupancy ? atom.occupancy : 1.0;
        if (!std::isfinite(occupancyScale) || occupancyScale < 0.0) {
            throw SEBException(
                "Invalid occupancy for atom serial " +
                    std::to_string(atom.serial),
                "AtomCloudBuilder::Build()"
            );
        }

        scatterers.emplace_back(
            atom.coordinate,
            parameters.radius,
            parameters.beta * occupancyScale
        );

        if (coordinatesBySerial.find(atom.serial) !=
            coordinatesBySerial.end()) {
            duplicateSerials.insert(atom.serial);
        } else {
            coordinatesBySerial[atom.serial] = atom.coordinate;
        }
    }

    std::unique_ptr<DebyeSphereCloud> cloud(
        new DebyeSphereCloud(scatterers)
    );

    for (const auto& reference : options.referenceAtomSerials) {
        if (reference.first == "center") {
            throw SEBException(
                "The reference name 'center' is reserved by DebyeSphereCloud",
                "AtomCloudBuilder::Build()"
            );
        }
        if (duplicateSerials.find(reference.second) !=
            duplicateSerials.end()) {
            throw SEBException(
                "Atom serial " + std::to_string(reference.second) +
                    " is ambiguous and cannot define reference " +
                    reference.first,
                "AtomCloudBuilder::Build()"
            );
        }

        const auto coordinate = coordinatesBySerial.find(reference.second);
        if (coordinate == coordinatesBySerial.end()) {
            throw SEBException(
                "Atom serial " + std::to_string(reference.second) +
                    " was not found for reference " + reference.first,
                "AtomCloudBuilder::Build()"
            );
        }
        cloud->addReferencePoint(reference.first, coordinate->second);
    }

    return cloud;
}

std::unique_ptr<DebyeSphereCloud> StructureCloudLoader::Load(
    const std::string& filename,
    const StructureParser& parser,
    const AtomParameterProfile& profile,
    const StructureParseOptions& parseOptions,
    const AtomCloudBuildOptions& buildOptions)
{
    return AtomCloudBuilder::Build(
        parser.ParseFile(filename, parseOptions),
        profile,
        buildOptions
    );
}

std::unique_ptr<DebyeSphereCloud> StructureCloudLoader::LoadPDB(
    const std::string& filename,
    const AtomParameterProfile& profile,
    const StructureParseOptions& parseOptions,
    const AtomCloudBuildOptions& buildOptions)
{
    const PDBParser parser;
    return Load(
        filename,
        parser,
        profile,
        parseOptions,
        buildOptions
    );
}
