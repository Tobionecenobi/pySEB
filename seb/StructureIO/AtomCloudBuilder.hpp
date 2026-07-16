#ifndef INCLUDE_GUARD_ATOMCLOUDBUILDER
#define INCLUDE_GUARD_ATOMCLOUDBUILDER

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "AtomParameterProfile.hpp"
#include "PDBParser.hpp"
#include "Subunits/DebyeSphereCloud.hpp"

struct AtomCloudBuildOptions {
    bool scaleBetaByOccupancy = true;
    std::map<refPoint, int> referenceAtomSerials;
};

class AtomCloudBuilder {
public:
    static std::unique_ptr<DebyeSphereCloud> Build(
        const std::vector<AtomRecord>& atoms,
        const AtomParameterProfile& profile,
        const AtomCloudBuildOptions& options = AtomCloudBuildOptions());
};

class StructureCloudLoader {
public:
    static std::unique_ptr<DebyeSphereCloud> Load(
        const std::string& filename,
        const StructureParser& parser,
        const AtomParameterProfile& profile,
        const StructureParseOptions& parseOptions =
            StructureParseOptions(),
        const AtomCloudBuildOptions& buildOptions =
            AtomCloudBuildOptions());

    static std::unique_ptr<DebyeSphereCloud> LoadPDB(
        const std::string& filename,
        const AtomParameterProfile& profile,
        const StructureParseOptions& parseOptions =
            StructureParseOptions(),
        const AtomCloudBuildOptions& buildOptions =
            AtomCloudBuildOptions());
};

#endif
