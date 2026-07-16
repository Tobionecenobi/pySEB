#ifndef INCLUDE_GUARD_PDBPARSER
#define INCLUDE_GUARD_PDBPARSER

#include "StructureParser.hpp"

class PDBParser : public StructureParser {
public:
    std::vector<AtomRecord> Parse(
        std::istream& input,
        const StructureParseOptions& options = StructureParseOptions(),
        const std::string& sourceName = "<stream>") const override;
};

#endif
