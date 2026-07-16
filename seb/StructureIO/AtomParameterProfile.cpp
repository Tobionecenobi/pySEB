#include "AtomParameterProfile.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>

std::string AtomParameterProfile::NormalizeIdentifier(
    const std::string& value)
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
        result.begin(),
        result.end(),
        result.begin(),
        [](unsigned char character) {
            return static_cast<char>(std::toupper(character));
        }
    );
    return result;
}

std::string AtomParameterProfile::AtomKey(
    const std::string& residueName,
    const std::string& atomName)
{
    return NormalizeIdentifier(residueName) + ":" +
        NormalizeIdentifier(atomName);
}

void AtomParameterProfile::ValidateParameters(
    const AtomScattererParameters& parameters,
    const std::string& description)
{
    if (!std::isfinite(parameters.radius) || parameters.radius < 0.0) {
        throw SEBException(
            "Radius for " + description +
                " must be finite and non-negative",
            "AtomParameterProfile"
        );
    }
    if (!std::isfinite(parameters.beta)) {
        throw SEBException(
            "Beta for " + description + " must be finite",
            "AtomParameterProfile"
        );
    }
}

void AtomParameterProfile::SetElement(
    const std::string& element,
    const AtomScattererParameters& parameters)
{
    const std::string key = NormalizeIdentifier(element);
    if (key.empty()) {
        throw SEBException(
            "Element rule cannot have an empty name",
            "AtomParameterProfile::SetElement()"
        );
    }
    ValidateParameters(parameters, "element " + key);
    elementRules[key] = parameters;
}

void AtomParameterProfile::SetAtom(
    const std::string& residueName,
    const std::string& atomName,
    const AtomScattererParameters& parameters)
{
    if (NormalizeIdentifier(residueName).empty() ||
        NormalizeIdentifier(atomName).empty()) {
        throw SEBException(
            "Atom rule must specify both a residue and an atom name",
            "AtomParameterProfile::SetAtom()"
        );
    }
    const std::string key = AtomKey(residueName, atomName);
    ValidateParameters(parameters, "atom " + key);
    atomRules[key] = parameters;
}

bool AtomParameterProfile::HasParameters(const AtomRecord& atom) const
{
    return atomRules.find(AtomKey(atom.residueName, atom.atomName)) !=
               atomRules.end() ||
        elementRules.find(NormalizeIdentifier(atom.element)) !=
               elementRules.end();
}

AtomScattererParameters AtomParameterProfile::Resolve(
    const AtomRecord& atom) const
{
    const std::string atomKey = AtomKey(atom.residueName, atom.atomName);
    const auto atomRule = atomRules.find(atomKey);
    if (atomRule != atomRules.end()) {
        return atomRule->second;
    }

    const std::string elementKey = NormalizeIdentifier(atom.element);
    const auto elementRule = elementRules.find(elementKey);
    if (elementRule != elementRules.end()) {
        return elementRule->second;
    }

    throw SEBException(
        "No scattering parameters for atom serial " +
            std::to_string(atom.serial) + " (" + atomKey +
            ", element " + elementKey + ")",
        "AtomParameterProfile::Resolve()"
    );
}
