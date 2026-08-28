#ifndef PYSEB_WORLD_IO_WORLD_DEFINITION_HPP
#define PYSEB_WORLD_IO_WORLD_DEFINITION_HPP

#include <string>
#include <utility>
#include <vector>

namespace pyseb {

struct WorldSubunitDefinition {
    std::string name;
    std::string model;
    std::string tag;
};

struct WorldStructureDefinition {
    std::string name;
    std::string graph;
};

struct WorldGraphDefinition {
    std::string id;
    std::vector<std::string> members;
    std::vector<std::pair<std::string, std::string>> links;
};

struct WorldDefinition {
    int schemaVersion = 1;
    std::string worldId;
    std::vector<WorldSubunitDefinition> subunits;
    std::vector<WorldStructureDefinition> structures;
    std::vector<WorldGraphDefinition> graphs;
    std::vector<std::string> roots;
};

} // namespace pyseb

#endif
