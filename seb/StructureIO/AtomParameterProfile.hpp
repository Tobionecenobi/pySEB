#ifndef INCLUDE_GUARD_ATOMPARAMETERPROFILE
#define INCLUDE_GUARD_ATOMPARAMETERPROFILE

#include <map>
#include <string>

#include "AtomRecord.hpp"
#include "Exceptions.hpp"

struct AtomScattererParameters {
    double radius = 0.0;
    double beta = 0.0;

    AtomScattererParameters() {}
    AtomScattererParameters(double radiusValue, double betaValue)
        : radius(radiusValue), beta(betaValue) {}
};

class AtomParameterProfile {
public:
    void SetElement(
        const std::string& element,
        const AtomScattererParameters& parameters);

    void SetElement(
        const std::string& element,
        double radius,
        double beta)
    {
        SetElement(element, AtomScattererParameters(radius, beta));
    }

    void SetAtom(
        const std::string& residueName,
        const std::string& atomName,
        const AtomScattererParameters& parameters);

    void SetAtom(
        const std::string& residueName,
        const std::string& atomName,
        double radius,
        double beta)
    {
        SetAtom(
            residueName,
            atomName,
            AtomScattererParameters(radius, beta)
        );
    }

    bool HasParameters(const AtomRecord& atom) const;
    AtomScattererParameters Resolve(const AtomRecord& atom) const;

private:
    std::map<std::string, AtomScattererParameters> elementRules;
    std::map<std::string, AtomScattererParameters> atomRules;

    static std::string NormalizeIdentifier(const std::string& value);
    static std::string AtomKey(
        const std::string& residueName,
        const std::string& atomName);
    static void ValidateParameters(
        const AtomScattererParameters& parameters,
        const std::string& description);
};

#endif
