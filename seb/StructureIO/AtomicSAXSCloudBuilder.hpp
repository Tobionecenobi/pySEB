#ifndef INCLUDE_GUARD_ATOMICSAXSCLOUDBUILDER
#define INCLUDE_GUARD_ATOMICSAXSCLOUDBUILDER

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "AtomicSAXSParameterProfile.hpp"
#include "PDBParser.hpp"
#include "Subunits/AtomicSAXSCloud.hpp"

enum class AtomicSAXSHydrationShellConvention { VolumeDensity, Pepsi3 };

struct AtomicSAXSCloudBuildOptions {
    bool scaleByOccupancy = true;
    bool useImplicitHydrogens = true;
    bool calculateSolventAccessibility = true;
    double sasaProbeRadius = 1.4;
    int sasaPointCount = 960;
    bool sasaUseProfileRadii = false;
    // Grid midpoints between the union of inner-radius spheres and the
    // union of (inner radius + thickness) spheres around active atoms.
    double hydrationShellInnerRadius = 3.0;
    double hydrationShellThickness = 3.0;
    double hydrationShellSpacing = 3.0;
    AtomicSAXSHydrationShellConvention hydrationShellConvention =
        AtomicSAXSHydrationShellConvention::VolumeDensity;
    std::map<refPoint, int> referenceAtomSerials;
};

class AtomicSAXSCloudBuilder {
public:
    static std::unique_ptr<AtomicSAXSCloud> Build(
        const std::vector<AtomRecord>& atoms,
        const AtomicSAXSParameterProfile& profile,
        const AtomicSAXSOptions& saxsOptions = AtomicSAXSOptions(),
        const AtomicSAXSCloudBuildOptions& buildOptions =
            AtomicSAXSCloudBuildOptions());
};

class AtomicSAXSStructureLoader {
public:
    static std::unique_ptr<AtomicSAXSCloud> Load(
        const std::string& filename,
        const StructureParser& parser,
        const AtomicSAXSParameterProfile& profile,
        const StructureParseOptions& parseOptions = StructureParseOptions(),
        const AtomicSAXSOptions& saxsOptions = AtomicSAXSOptions(),
        const AtomicSAXSCloudBuildOptions& buildOptions =
            AtomicSAXSCloudBuildOptions());

    static std::unique_ptr<AtomicSAXSCloud> LoadPDB(
        const std::string& filename,
        const AtomicSAXSParameterProfile& profile,
        const StructureParseOptions& parseOptions = StructureParseOptions(),
        const AtomicSAXSOptions& saxsOptions = AtomicSAXSOptions(),
        const AtomicSAXSCloudBuildOptions& buildOptions =
            AtomicSAXSCloudBuildOptions());
};

#endif
