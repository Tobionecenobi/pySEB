#include "PDBParser.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <sstream>

namespace {

std::string trim(const std::string& value)
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
    return value.substr(first, last - first);
}

std::string field(
    const std::string& line,
    std::size_t offset,
    std::size_t width)
{
    if (offset >= line.size()) {
        return "";
    }
    return line.substr(offset, std::min(width, line.size() - offset));
}

std::string upper(const std::string& value)
{
    std::string result = value;
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

std::string location(
    const std::string& sourceName,
    std::size_t lineNumber)
{
    return sourceName + ":" + std::to_string(lineNumber);
}

int parseInteger(
    const std::string& text,
    const std::string& description,
    const std::string& sourceName,
    std::size_t lineNumber)
{
    const std::string value = trim(text);
    if (value.empty()) {
        throw SEBException(
            "Missing " + description + " at " +
                location(sourceName, lineNumber),
            "PDBParser::Parse()"
        );
    }

    try {
        std::size_t consumed = 0;
        const int result = std::stoi(value, &consumed);
        if (consumed != value.size()) {
            throw std::invalid_argument("trailing data");
        }
        return result;
    } catch (const std::exception&) {
        throw SEBException(
            "Invalid " + description + " '" + value + "' at " +
                location(sourceName, lineNumber),
            "PDBParser::Parse()"
        );
    }
}

double parseDouble(
    const std::string& text,
    const std::string& description,
    const std::string& sourceName,
    std::size_t lineNumber,
    bool required,
    double defaultValue)
{
    const std::string value = trim(text);
    if (value.empty()) {
        if (!required) {
            return defaultValue;
        }
        throw SEBException(
            "Missing " + description + " at " +
                location(sourceName, lineNumber),
            "PDBParser::Parse()"
        );
    }

    try {
        std::size_t consumed = 0;
        const double result = std::stod(value, &consumed);
        if (consumed != value.size() || !std::isfinite(result)) {
            throw std::invalid_argument("invalid floating-point value");
        }
        return result;
    } catch (const std::exception&) {
        throw SEBException(
            "Invalid " + description + " '" + value + "' at " +
                location(sourceName, lineNumber),
            "PDBParser::Parse()"
        );
    }
}

std::string inferElement(
    const std::string& elementField,
    const std::string& rawAtomName)
{
    const std::string explicitElement = upper(trim(elementField));
    if (!explicitElement.empty()) {
        return explicitElement;
    }

    if (!rawAtomName.empty() && rawAtomName[0] == ' ') {
        for (const char character : rawAtomName) {
            if (std::isalpha(static_cast<unsigned char>(character))) {
                return std::string(
                    1,
                    static_cast<char>(
                        std::toupper(static_cast<unsigned char>(character))
                    )
                );
            }
        }
    }

    std::string inferred;
    for (const char character : rawAtomName) {
        if (std::isalpha(static_cast<unsigned char>(character))) {
            inferred.push_back(
                static_cast<char>(
                    std::toupper(static_cast<unsigned char>(character))
                )
            );
            if (inferred.size() == 2) {
                break;
            }
        } else if (!inferred.empty()) {
            break;
        }
    }
    return inferred;
}

bool isWater(const std::string& residueName)
{
    const std::string residue = upper(trim(residueName));
    return residue == "HOH" || residue == "WAT" || residue == "DOD";
}

bool isPrimaryAlternateLocation(const std::string& alternateLocation)
{
    const std::string alternate = upper(trim(alternateLocation));
    return alternate.empty() || alternate == "A" || alternate == "1";
}

}

std::vector<AtomRecord> PDBParser::Parse(
    std::istream& input,
    const StructureParseOptions& options,
    const std::string& sourceName) const
{
    if (options.modelNumber < 1) {
        throw SEBException(
            "PDB model number must be positive",
            "PDBParser::Parse()"
        );
    }

    std::vector<AtomRecord> atoms;
    std::set<int> modelsSeen;
    bool hasModelRecords = false;
    int currentModel = 1;
    std::string line;
    std::size_t lineNumber = 0;

    while (std::getline(input, line)) {
        ++lineNumber;
        const std::string record = upper(trim(field(line, 0, 6)));

        if (record == "MODEL") {
            hasModelRecords = true;
            currentModel = parseInteger(
                field(line, 10, 4),
                "MODEL number",
                sourceName,
                lineNumber
            );
            modelsSeen.insert(currentModel);
            continue;
        }

        if (record == "ENDMDL") {
            if (hasModelRecords && currentModel == options.modelNumber) {
                break;
            }
            continue;
        }

        if (record != "ATOM" && record != "HETATM") {
            continue;
        }
        if (hasModelRecords && currentModel != options.modelNumber) {
            continue;
        }
        if (!hasModelRecords && options.modelNumber != 1) {
            continue;
        }
        if (record == "HETATM" && !options.includeHetatm) {
            continue;
        }

        const std::string residueName = trim(field(line, 17, 3));
        if (record == "HETATM" &&
            isWater(residueName) &&
            !options.includeWater) {
            continue;
        }

        const std::string alternateLocation = trim(field(line, 16, 1));
        if (options.alternateLocations == AlternateLocationPolicy::Primary &&
            !isPrimaryAlternateLocation(alternateLocation)) {
            continue;
        }

        const std::string rawAtomName = field(line, 12, 4);
        const std::string element = inferElement(
            field(line, 76, 2),
            rawAtomName
        );
        if (element.empty()) {
            throw SEBException(
                "Could not determine element at " +
                    location(sourceName, lineNumber),
                "PDBParser::Parse()"
            );
        }
        if (!options.includeHydrogen &&
            (element == "H" || element == "D")) {
            continue;
        }

        AtomRecord atom;
        atom.recordType = record;
        atom.serial = parseInteger(
            field(line, 6, 5),
            "atom serial",
            sourceName,
            lineNumber
        );
        atom.atomName = trim(rawAtomName);
        atom.alternateLocation = alternateLocation;
        atom.residueName = residueName;
        atom.chainId = trim(field(line, 21, 1));
        atom.residueNumber = parseInteger(
            field(line, 22, 4),
            "residue number",
            sourceName,
            lineNumber
        );
        atom.insertionCode = trim(field(line, 26, 1));
        atom.coordinate = CartesianPoint3D(
            parseDouble(
                field(line, 30, 8),
                "x coordinate",
                sourceName,
                lineNumber,
                true,
                0.0
            ),
            parseDouble(
                field(line, 38, 8),
                "y coordinate",
                sourceName,
                lineNumber,
                true,
                0.0
            ),
            parseDouble(
                field(line, 46, 8),
                "z coordinate",
                sourceName,
                lineNumber,
                true,
                0.0
            )
        );
        atom.occupancy = parseDouble(
            field(line, 54, 6),
            "occupancy",
            sourceName,
            lineNumber,
            false,
            1.0
        );
        if (atom.occupancy < 0.0) {
            throw SEBException(
                "Occupancy cannot be negative at " +
                    location(sourceName, lineNumber),
                "PDBParser::Parse()"
            );
        }
        atom.temperatureFactor = parseDouble(
            field(line, 60, 6),
            "temperature factor",
            sourceName,
            lineNumber,
            false,
            0.0
        );
        atom.element = element;
        atom.charge = trim(field(line, 78, 2));
        atom.modelNumber = hasModelRecords ? currentModel : 1;
        atoms.push_back(atom);
    }

    if (hasModelRecords &&
        modelsSeen.find(options.modelNumber) == modelsSeen.end()) {
        throw SEBException(
            "PDB model " + std::to_string(options.modelNumber) +
                " was not found in " + sourceName,
            "PDBParser::Parse()"
        );
    }
    if (!hasModelRecords && options.modelNumber != 1) {
        throw SEBException(
            "PDB file " + sourceName +
                " has no MODEL records; only model 1 is available",
            "PDBParser::Parse()"
        );
    }

    return atoms;
}
