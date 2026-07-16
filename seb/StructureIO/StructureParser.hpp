#ifndef INCLUDE_GUARD_STRUCTUREPARSER
#define INCLUDE_GUARD_STRUCTUREPARSER

#include <fstream>
#include <istream>
#include <string>
#include <vector>

#include "AtomRecord.hpp"
#include "Exceptions.hpp"

enum class AlternateLocationPolicy {
    Primary,
    All
};

struct StructureParseOptions {
    int modelNumber = 1;
    bool includeHetatm = false;
    bool includeWater = false;
    bool includeHydrogen = true;
    AlternateLocationPolicy alternateLocations =
        AlternateLocationPolicy::Primary;
};

class StructureParser {
public:
    virtual ~StructureParser() {}

    virtual std::vector<AtomRecord> Parse(
        std::istream& input,
        const StructureParseOptions& options = StructureParseOptions(),
        const std::string& sourceName = "<stream>") const = 0;

    std::vector<AtomRecord> ParseFile(
        const std::string& filename,
        const StructureParseOptions& options = StructureParseOptions()) const
    {
        std::ifstream input(filename.c_str());
        if (!input.is_open()) {
            throw SEBException(
                "Could not open structure file " + filename,
                "StructureParser::ParseFile()"
            );
        }
        return Parse(input, options, filename);
    }
};

#endif
