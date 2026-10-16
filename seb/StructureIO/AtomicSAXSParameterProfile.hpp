#ifndef INCLUDE_GUARD_ATOMICSAXSPARAMETERPROFILE
#define INCLUDE_GUARD_ATOMICSAXSPARAMETERPROFILE

#include <map>
#include <string>

#include "AtomRecord.hpp"
#include "Exceptions.hpp"
#include "Subunits/AtomicSAXSCloud.hpp"

struct AtomicSAXSParameters {
    CromerMannCoefficients vacuum;
    double excludedVolume = 0.0;
    // Bare-atom volume when grouped hydrogens are disabled; negative means
    // the supplied excludedVolume already describes an individual atom.
    double atomicExcludedVolume = -1.0;
    int implicitHydrogens = 0;
    double implicitHydrogenDistance = 0.0;
    double sasaRadius = 0.0;

    AtomicSAXSParameters() {}
    AtomicSAXSParameters(
        const CromerMannCoefficients& vacuumValue,
        double excludedVolumeValue,
        int implicitHydrogensValue = 0,
        double sasaRadiusValue = 0.0)
        : vacuum(vacuumValue),
          excludedVolume(excludedVolumeValue),
          implicitHydrogens(implicitHydrogensValue),
          sasaRadius(sasaRadiusValue) {}
};

enum class AtomicSAXSVolumeModel {
    BondiElement,
    Fraser1978,
    Pontius1996
};

class AtomicSAXSParameterProfile {
public:
    void SetElement(
        const std::string& element,
        const AtomicSAXSParameters& parameters);
    void SetAtom(
        const std::string& residueName,
        const std::string& atomName,
        const AtomicSAXSParameters& parameters);

    bool HasParameters(const AtomRecord& atom) const;
    AtomicSAXSParameters Resolve(const AtomRecord& atom) const;

    static AtomicSAXSParameterProfile ElementDefaults();
    static AtomicSAXSParameterProfile ProteinDefaults(
        AtomicSAXSVolumeModel volumeModel = AtomicSAXSVolumeModel::Fraser1978);

private:
    std::map<std::string, AtomicSAXSParameters> elementRules;
    std::map<std::string, AtomicSAXSParameters> atomRules;

    static std::string NormalizeIdentifier(const std::string& value);
    static std::string AtomKey(
        const std::string& residueName,
        const std::string& atomName);
    static void ValidateParameters(
        const AtomicSAXSParameters& parameters,
        const std::string& description);
};

#endif
