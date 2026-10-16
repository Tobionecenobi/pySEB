#include "AtomicSAXSParameterProfile.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <initializer_list>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double sphereVolume(double radius)
{
    return 4.0 * kPi * radius * radius * radius / 3.0;
}

CromerMannCoefficients coefficients(
    const std::array<double, 4>& a,
    const std::array<double, 4>& b,
    double c)
{
    return CromerMannCoefficients(a, b, c);
}

void setProteinHydrogens(
    AtomicSAXSParameterProfile& profile,
    const AtomicSAXSParameterProfile& elements,
    const std::string& residue,
    const std::string& atomName,
    const std::string& element,
    int hydrogenCount,
    AtomicSAXSVolumeModel volumeModel)
{
    AtomRecord atom;
    atom.residueName = residue;
    atom.atomName = atomName;
    atom.element = element;
    AtomicSAXSParameters parameters = elements.Resolve(atom);
    parameters.atomicExcludedVolume = parameters.excludedVolume;
    parameters.implicitHydrogens = hydrogenCount;
    parameters.implicitHydrogenDistance = 0.0;
    // Allen et al. distances tabulated in Grudinin et al. (2017), Table 1;
    // use their spherical-amplitude group approximation, equation (8).
    if (hydrogenCount > 0) {
        if (element == "C") {
            const double distances[] = {0.0, 1.099, 1.092, 1.059};
            parameters.implicitHydrogenDistance = distances[hydrogenCount];
            if (hydrogenCount == 1 &&
                (residue == "PHE" || residue == "TYR" ||
                 residue == "HIS" || residue == "TRP") && atomName != "CA") {
                parameters.implicitHydrogenDistance = 1.083;
            }
        } else if (element == "N") {
            parameters.implicitHydrogenDistance =
                residue == "LYS" && atomName == "NZ" ? 1.033 : 1.009;
        } else if (element == "O") {
            parameters.implicitHydrogenDistance = 0.967;
        }
        // No S-H distance is supplied by that table: retain zero rather
        // than introduce an undocumented calibration to Pepsi output.
    }
    if (volumeModel == AtomicSAXSVolumeModel::Fraser1978) {
        parameters.excludedVolume += 5.15 * hydrogenCount;
    } else if (volumeModel == AtomicSAXSVolumeModel::Pontius1996) {
        if (element == "C") {
            const double carbonVolumes[] = {16.44, 11.8, 20.9, 33.9};
            parameters.excludedVolume = carbonVolumes[hydrogenCount];
        } else if (element == "N") {
            // Pontius reports N=8.8 and NH=14.1 A^3. The NH2/NH3
            // entries use the same observed 5.3 A^3 hydrogen increment.
            parameters.excludedVolume = 8.8 + 5.3 * hydrogenCount;
        } else if (element == "O" && hydrogenCount == 1) {
            parameters.excludedVolume = 23.9;
        } else if (element == "S" && hydrogenCount == 1) {
            // Pontius' compact comparison table has no separate thiol
            // value. Retain the Fraser S+H group explicitly rather than
            // substituting a van der Waals sphere.
            parameters.excludedVolume = 19.86 + 5.15;
        }
    }
    profile.SetAtom(residue, atomName, parameters);
}

void setElementVolume(
    AtomicSAXSParameterProfile& profile,
    const AtomicSAXSParameterProfile& source,
    const std::string& element,
    double volume)
{
    AtomRecord atom;
    atom.element = element;
    AtomicSAXSParameters parameters = source.Resolve(atom);
    parameters.excludedVolume = volume;
    profile.SetElement(element, parameters);
}

} // namespace

std::string AtomicSAXSParameterProfile::NormalizeIdentifier(
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

std::string AtomicSAXSParameterProfile::AtomKey(
    const std::string& residueName,
    const std::string& atomName)
{
    return NormalizeIdentifier(residueName) + ":" +
        NormalizeIdentifier(atomName);
}

void AtomicSAXSParameterProfile::ValidateParameters(
    const AtomicSAXSParameters& parameters,
    const std::string& description)
{
    if (!std::isfinite(parameters.excludedVolume) ||
        parameters.excludedVolume < 0.0) {
        throw SEBException(
            "Excluded volume for " + description +
                " must be finite and non-negative",
            "AtomicSAXSParameterProfile"
        );
    }
    if (parameters.implicitHydrogens < 0) {
        throw SEBException(
            "Implicit hydrogen count for " + description +
                " must be non-negative",
            "AtomicSAXSParameterProfile"
        );
    }
    if (!std::isfinite(parameters.atomicExcludedVolume) ||
        (parameters.atomicExcludedVolume < 0.0 && parameters.atomicExcludedVolume != -1.0)) {
        throw SEBException("Atomic excluded volume must be non-negative or -1 (unspecified)",
                           "AtomicSAXSParameterProfile");
    }
    if (!std::isfinite(parameters.implicitHydrogenDistance) ||
        parameters.implicitHydrogenDistance < 0.0) {
        throw SEBException(
            "Implicit hydrogen distance must be finite and non-negative",
            "AtomicSAXSParameterProfile");
    }
    if (!std::isfinite(parameters.sasaRadius) || parameters.sasaRadius < 0.0) {
        throw SEBException(
            "SASA radius for " + description +
                " must be finite and non-negative",
            "AtomicSAXSParameterProfile"
        );
    }
    if (!std::isfinite(parameters.vacuum.c)) {
        throw SEBException(
            "Cromer--Mann coefficients for " + description +
                " must be finite",
            "AtomicSAXSParameterProfile"
        );
    }
    for (std::size_t index = 0; index < parameters.vacuum.a.size(); ++index) {
        if (!std::isfinite(parameters.vacuum.a[index]) ||
            !std::isfinite(parameters.vacuum.b[index]) ||
            parameters.vacuum.b[index] < 0.0) {
            throw SEBException(
                "Cromer--Mann coefficients for " + description +
                    " must be finite and b values non-negative",
                "AtomicSAXSParameterProfile"
            );
        }
    }
}

void AtomicSAXSParameterProfile::SetElement(
    const std::string& element,
    const AtomicSAXSParameters& parameters)
{
    const std::string key = NormalizeIdentifier(element);
    if (key.empty()) {
        throw SEBException(
            "Element rule cannot have an empty name",
            "AtomicSAXSParameterProfile::SetElement()"
        );
    }
    ValidateParameters(parameters, "element " + key);
    elementRules[key] = parameters;
}

void AtomicSAXSParameterProfile::SetAtom(
    const std::string& residueName,
    const std::string& atomName,
    const AtomicSAXSParameters& parameters)
{
    if (NormalizeIdentifier(residueName).empty() ||
        NormalizeIdentifier(atomName).empty()) {
        throw SEBException(
            "Atom rule must specify both a residue and an atom name",
            "AtomicSAXSParameterProfile::SetAtom()"
        );
    }
    const std::string key = AtomKey(residueName, atomName);
    ValidateParameters(parameters, "atom " + key);
    atomRules[key] = parameters;
}

bool AtomicSAXSParameterProfile::HasParameters(const AtomRecord& atom) const
{
    return atomRules.find(AtomKey(atom.residueName, atom.atomName)) !=
               atomRules.end() ||
        elementRules.find(NormalizeIdentifier(atom.element)) !=
               elementRules.end();
}

AtomicSAXSParameters AtomicSAXSParameterProfile::Resolve(
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
        "No atomic SAXS parameters for atom serial " +
            std::to_string(atom.serial) + " (" + atomKey +
            ", element " + elementKey + ")",
        "AtomicSAXSParameterProfile::Resolve()"
    );
}

AtomicSAXSParameterProfile AtomicSAXSParameterProfile::ElementDefaults()
{
    AtomicSAXSParameterProfile profile;
    // Cromer--Mann neutral-atom coefficients use q/(4 pi) in inverse
    // Angstrom. Volumes are isolated Bondi van der Waals spheres and are an
    // explicit, replaceable baseline rather than fitted FoXS/CRYSOL tables.
    profile.SetElement(
        "H",
        AtomicSAXSParameters(
            coefficients(
                {{0.493002, 0.322912, 0.140191, 0.040810}},
                {{10.5109, 26.1257, 3.14236, 57.7997}},
                0.003038),
            sphereVolume(1.20), 0, 1.20));
    profile.SetElement(
        "C",
        AtomicSAXSParameters(
            coefficients(
                {{2.3100, 1.0200, 1.5886, 0.8650}},
                {{20.8439, 10.2075, 0.5687, 51.6512}},
                0.2156),
            sphereVolume(1.70), 0, 1.70));
    profile.SetElement(
        "N",
        AtomicSAXSParameters(
            coefficients(
                {{12.2126, 3.1322, 2.0125, 1.1663}},
                {{0.0057, 9.8933, 28.9975, 0.5826}},
                -11.5290),
            sphereVolume(1.55), 0, 1.55));
    profile.SetElement(
        "O",
        AtomicSAXSParameters(
            coefficients(
                {{3.0485, 2.2868, 1.5463, 0.8670}},
                {{13.2771, 5.7011, 0.3239, 32.9089}},
                0.2508),
            sphereVolume(1.52), 0, 1.52));
    profile.SetElement(
        "P",
        AtomicSAXSParameters(
            coefficients(
                {{6.4345, 4.1791, 1.7800, 1.4908}},
                {{1.9067, 27.1570, 0.5260, 68.1645}},
                1.1149),
            sphereVolume(1.80), 0, 1.80));
    profile.SetElement(
        "S",
        AtomicSAXSParameters(
            coefficients(
                {{6.9053, 5.2034, 1.4379, 1.5863}},
                {{1.4679, 22.2151, 0.2536, 56.1720}},
                0.8669),
            sphereVolume(1.80), 0, 1.80));
    return profile;
}

AtomicSAXSParameterProfile AtomicSAXSParameterProfile::ProteinDefaults(
    AtomicSAXSVolumeModel volumeModel)
{
    const AtomicSAXSParameterProfile elements = ElementDefaults();
    AtomicSAXSParameterProfile profile = elements;

    if (volumeModel == AtomicSAXSVolumeModel::Fraser1978) {
        setElementVolume(profile, elements, "H", 5.15);
        setElementVolume(profile, elements, "C", 16.44);
        setElementVolume(profile, elements, "N", 2.49);
        setElementVolume(profile, elements, "O", 9.13);
        setElementVolume(profile, elements, "P", 5.73);
        setElementVolume(profile, elements, "S", 19.86);
    } else if (volumeModel == AtomicSAXSVolumeModel::Pontius1996) {
        // Protein-core Voronoi means from Pontius et al. The compact
        // published comparison does not give replacements for H, carbonyl
        // C, P, or S, so those retain the explicitly documented Fraser
        // values instead of silently reverting to Bondi spheres.
        setElementVolume(profile, elements, "H", 5.15);
        setElementVolume(profile, elements, "C", 16.44);
        setElementVolume(profile, elements, "N", 8.8);
        setElementVolume(profile, elements, "O", 22.3);
        setElementVolume(profile, elements, "P", 5.73);
        setElementVolume(profile, elements, "S", 19.86);
    } else if (volumeModel != AtomicSAXSVolumeModel::BondiElement) {
        throw SEBException(
            "Unknown atomic SAXS volume model",
            "AtomicSAXSParameterProfile::ProteinDefaults()"
        );
    }

    const std::initializer_list<const char*> residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL"
    };
    for (const char* residue : residues) {
        setProteinHydrogens(
            profile, profile, residue, "N", "N",
            std::string(residue) == "PRO" ? 0 : 1, volumeModel);
        setProteinHydrogens(
            profile, profile, residue, "CA", "C",
            std::string(residue) == "GLY" ? 2 : 1, volumeModel);

        // These standard atoms carry no implicit hydrogen, but explicit
        // atom rules keep every standard backbone assignment tied to the
        // selected volume model.
        setProteinHydrogens(
            profile, profile, residue, "C", "C", 0, volumeModel);
        setProteinHydrogens(
            profile, profile, residue, "O", "O", 0, volumeModel);
        setProteinHydrogens(
            profile, profile, residue, "OXT", "O", 0, volumeModel);
    }

    const auto add = [&](const char* residue, const char* atom,
                         const char* element, int hydrogens) {
        setProteinHydrogens(
            profile, profile, residue, atom, element, hydrogens, volumeModel);
    };
    add("ALA", "CB", "C", 3);
    add("ARG", "CB", "C", 2); add("ARG", "CG", "C", 2);
    add("ARG", "CD", "C", 2); add("ARG", "NE", "N", 1);
    add("ARG", "CZ", "C", 0);
    add("ARG", "NH1", "N", 2); add("ARG", "NH2", "N", 2);
    add("ASN", "CB", "C", 2); add("ASN", "CG", "C", 0);
    add("ASN", "OD1", "O", 0); add("ASN", "ND2", "N", 2);
    add("ASP", "CB", "C", 2); add("ASP", "CG", "C", 0);
    add("ASP", "OD1", "O", 0); add("ASP", "OD2", "O", 0);
    add("CYS", "CB", "C", 2); add("CYS", "SG", "S", 1);
    add("GLN", "CB", "C", 2); add("GLN", "CG", "C", 2);
    add("GLN", "CD", "C", 0); add("GLN", "OE1", "O", 0);
    add("GLN", "NE2", "N", 2);
    add("GLU", "CB", "C", 2); add("GLU", "CG", "C", 2);
    add("GLU", "CD", "C", 0); add("GLU", "OE1", "O", 0);
    add("GLU", "OE2", "O", 0);
    add("HIS", "CB", "C", 2); add("HIS", "CG", "C", 0);
    add("HIS", "ND1", "N", 1);
    add("HIS", "CD2", "C", 1); add("HIS", "CE1", "C", 1);
    add("HIS", "NE2", "N", 0);
    add("ILE", "CB", "C", 1); add("ILE", "CG1", "C", 2);
    add("ILE", "CG2", "C", 3); add("ILE", "CD1", "C", 3);
    add("LEU", "CB", "C", 2); add("LEU", "CG", "C", 1);
    add("LEU", "CD1", "C", 3); add("LEU", "CD2", "C", 3);
    add("LYS", "CB", "C", 2); add("LYS", "CG", "C", 2);
    add("LYS", "CD", "C", 2); add("LYS", "CE", "C", 2);
    add("LYS", "NZ", "N", 3);
    add("MET", "CB", "C", 2); add("MET", "CG", "C", 2);
    add("MET", "SD", "S", 0); add("MET", "CE", "C", 3);
    add("PHE", "CB", "C", 2); add("PHE", "CG", "C", 0);
    add("PHE", "CD1", "C", 1);
    add("PHE", "CD2", "C", 1); add("PHE", "CE1", "C", 1);
    add("PHE", "CE2", "C", 1); add("PHE", "CZ", "C", 1);
    add("PRO", "CB", "C", 2); add("PRO", "CG", "C", 2);
    add("PRO", "CD", "C", 2);
    add("SER", "CB", "C", 2); add("SER", "OG", "O", 1);
    add("THR", "CB", "C", 1); add("THR", "OG1", "O", 1);
    add("THR", "CG2", "C", 3);
    add("TRP", "CB", "C", 2); add("TRP", "CG", "C", 0);
    add("TRP", "CD1", "C", 1); add("TRP", "CD2", "C", 0);
    add("TRP", "NE1", "N", 1); add("TRP", "CE3", "C", 1);
    add("TRP", "CE2", "C", 0);
    add("TRP", "CZ2", "C", 1); add("TRP", "CZ3", "C", 1);
    add("TRP", "CH2", "C", 1);
    add("TYR", "CB", "C", 2); add("TYR", "CG", "C", 0);
    add("TYR", "CD1", "C", 1);
    add("TYR", "CD2", "C", 1); add("TYR", "CE1", "C", 1);
    add("TYR", "CE2", "C", 1); add("TYR", "CZ", "C", 0);
    add("TYR", "OH", "O", 1);
    add("VAL", "CB", "C", 1); add("VAL", "CG1", "C", 3);
    add("VAL", "CG2", "C", 3);

    return profile;
}
