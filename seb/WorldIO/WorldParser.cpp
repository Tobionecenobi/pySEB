#include "WorldIO/WorldParser.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "Exceptions.hpp"
#include "fkYAML/node.hpp"
#include "World.hpp"

namespace pyseb {
namespace {

using Node = fkyaml::node;
constexpr std::size_t kMaximumWorldFileSize = 4u * 1024u * 1024u;

[[noreturn]] void schemaError(
    const std::string& source,
    const std::string& path,
    const std::string& message) {
    throw SEBException(
        source + ": " + path + ": " + message,
        "LoadWorldDefinitionYaml()");
}

void requireMapping(const Node& node, const std::string& source, const std::string& path) {
    if (!node.is_mapping()) schemaError(source, path, "expected a mapping");
}

void requireSequence(const Node& node, const std::string& source, const std::string& path) {
    if (!node.is_sequence()) schemaError(source, path, "expected a sequence");
}

void rejectUnknownKeys(
    const Node& node,
    const std::set<std::string>& allowed,
    const std::string& source,
    const std::string& path) {
    requireMapping(node, source, path);
    for (const auto& item : node.map_items()) {
        if (!item.key().is_string()) schemaError(source, path, "mapping keys must be strings");
        const std::string key = item.key().get_value<std::string>();
        if (!allowed.count(key)) schemaError(source, path + "." + key, "unknown field");
    }
}

const Node& requiredNode(
    const Node& node,
    const char* key,
    const std::string& source,
    const std::string& path) {
    if (!node.contains(key)) schemaError(source, path + "." + key, "required field is missing");
    return node.at(key);
}

std::string scalarString(const Node& node, const std::string& source, const std::string& path) {
    if (!node.is_string()) schemaError(source, path, "expected a string");
    const std::string value = node.get_value<std::string>();
    if (value.empty()) schemaError(source, path, "must not be empty");
    return value;
}

int scalarInteger(const Node& node, const std::string& source, const std::string& path) {
    if (!node.is_integer()) schemaError(source, path, "expected an integer");
    return node.get_value<int>();
}

std::vector<std::string> stringSequence(
    const Node& node,
    const std::string& source,
    const std::string& path) {
    requireSequence(node, source, path);
    std::vector<std::string> result;
    std::size_t index = 0;
    for (const auto& item : node) {
        result.push_back(scalarString(item, source, path + "[" + std::to_string(index++) + "]"));
    }
    return result;
}

std::pair<std::string, std::string> linkPair(
    const Node& node,
    const std::string& source,
    const std::string& path) {
    requireSequence(node, source, path);
    std::vector<std::string> values;
    for (const auto& item : node) values.push_back(scalarString(item, source, path));
    if (values.size() != 2) schemaError(source, path, "expected exactly two endpoints");
    return {values[0], values[1]};
}

void validateRestrictions(
    const Node& node,
    const std::string& source,
    const std::string& path,
    std::size_t depth = 0) {
    if (depth > 64) schemaError(source, path, "YAML nesting exceeds 64 levels");
    if (node.is_anchor() || node.is_alias()) schemaError(source, path, "YAML anchors and aliases are not supported");
    if (node.has_tag_name()) schemaError(source, path, "YAML tags are not supported");
    if (node.is_sequence()) {
        std::size_t index = 0;
        for (const auto& child : node) {
            validateRestrictions(child, source, path + "[" + std::to_string(index++) + "]", depth + 1);
        }
    } else if (node.is_mapping()) {
        for (const auto& item : node.map_items()) {
            validateRestrictions(item.key(), source, path + ".<key>", depth + 1);
            const std::string key = item.key().is_string() ? item.key().get_value<std::string>() : "<value>";
            validateRestrictions(item.value(), source, path + "." + key, depth + 1);
        }
    }
}

std::string topName(const std::string& endpoint) {
    const std::size_t colon = endpoint.find(':');
    const std::size_t period = endpoint.find('.');
    std::size_t cut = std::string::npos;
    if (colon != std::string::npos) cut = colon;
    if (period != std::string::npos) cut = std::min(cut, period);
    return cut == std::string::npos ? endpoint : endpoint.substr(0, cut);
}

bool validName(const std::string& value) {
    if (value.empty() || !std::isalpha(static_cast<unsigned char>(value.front()))) return false;
    return std::all_of(value.begin() + 1, value.end(), [](char ch) {
        return std::isalnum(static_cast<unsigned char>(ch)) != 0;
    });
}

bool validGraphId(const std::string& value) {
    if (value.empty()) return false;
    return std::all_of(value.begin(), value.end(), [](char ch) {
        return std::isalnum(static_cast<unsigned char>(ch)) || ch == '_' || ch == '-';
    });
}

void rejectMultipleDocuments(const std::string& yaml, const std::string& source) {
    std::istringstream input(yaml);
    std::string line;
    bool seenContent = false;
    bool ended = false;
    std::size_t lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::size_t first = line.find_first_not_of(" \t");
        if (first == std::string::npos || line[first] == '#') continue;
        const bool start = first == 0 && line.compare(first, 3, "---") == 0 &&
            (line.size() == 3 || std::isspace(static_cast<unsigned char>(line[3])) || line[3] == '#');
        const bool end = first == 0 && line.compare(first, 3, "...") == 0 &&
            (line.size() == 3 || std::isspace(static_cast<unsigned char>(line[3])) || line[3] == '#');
        if (start) {
            if (seenContent || ended) {
                schemaError(source, "$", "multiple YAML documents are not supported (line " +
                    std::to_string(lineNumber) + ")");
            }
            const std::size_t inlineContent = line.find_first_not_of(" \t", 3);
            if (inlineContent != std::string::npos && line[inlineContent] != '#') seenContent = true;
            continue;
        }
        if (end) {
            if (!seenContent) schemaError(source, "$", "unexpected YAML document terminator");
            ended = true;
            continue;
        }
        if (ended) schemaError(source, "$", "multiple YAML documents are not supported (line " + std::to_string(lineNumber) + ")");
        seenContent = true;
    }
}

std::string quote(const std::string& value) {
    std::string result = "\"";
    for (char ch : value) {
        const unsigned char byte = static_cast<unsigned char>(ch);
        if (ch == '\\' || ch == '"') {
            result += '\\';
            result += ch;
        } else if (ch == '\n') {
            result += "\\n";
        } else if (ch == '\r') {
            result += "\\r";
        } else if (ch == '\t') {
            result += "\\t";
        } else if (byte < 0x20 || byte == 0x7f) {
            static constexpr char hex[] = "0123456789abcdef";
            result += "\\u00";
            result += hex[(byte >> 4) & 0x0f];
            result += hex[byte & 0x0f];
        } else {
            result += ch;
        }
    }
    result += '"';
    return result;
}

template <typename T, typename Less>
void sortBy(T& values, Less less) {
    std::sort(values.begin(), values.end(), less);
}

} // namespace

WorldDefinition LoadWorldDefinitionYaml(const std::string& yaml, const std::string& source) {
    if (yaml.size() > kMaximumWorldFileSize) {
        throw SEBException(source + ": world file exceeds the 4 MiB limit", "LoadWorldDefinitionYaml()");
    }
    rejectMultipleDocuments(yaml, source);

    Node root;
    try {
        root = Node::deserialize(yaml);
    } catch (const std::exception& error) {
        throw SEBException(source + ": " + error.what(), "LoadWorldDefinitionYaml()");
    }
    validateRestrictions(root, source, "$");
    rejectUnknownKeys(root, {"format", "schema_version", "world", "subunits", "structures", "graphs", "roots"}, source, "$");

    const std::string format = scalarString(requiredNode(root, "format", source, "$"), source, "format");
    if (format != "pyseb-world") schemaError(source, "format", "expected 'pyseb-world'");
    const int version = scalarInteger(requiredNode(root, "schema_version", source, "$"), source, "schema_version");
    if (version != 1) schemaError(source, "schema_version", "only schema version 1 is supported");

    WorldDefinition result;
    result.schemaVersion = version;

    const Node& world = requiredNode(root, "world", source, "$");
    rejectUnknownKeys(world, {"id"}, source, "world");
    result.worldId = scalarString(requiredNode(world, "id", source, "world"), source, "world.id");

    const Node& subunits = requiredNode(root, "subunits", source, "$");
    requireSequence(subunits, source, "subunits");
    std::size_t index = 0;
    for (const auto& item : subunits) {
        const std::string path = "subunits[" + std::to_string(index++) + "]";
        rejectUnknownKeys(item, {"name", "model", "tag"}, source, path);
        WorldSubunitDefinition value;
        value.name = scalarString(requiredNode(item, "name", source, path), source, path + ".name");
        value.model = scalarString(requiredNode(item, "model", source, path), source, path + ".model");
        value.tag = item.contains("tag") ? scalarString(item.at("tag"), source, path + ".tag") : value.name;
        result.subunits.push_back(std::move(value));
    }

    const Node& structures = requiredNode(root, "structures", source, "$");
    requireSequence(structures, source, "structures");
    index = 0;
    for (const auto& item : structures) {
        const std::string path = "structures[" + std::to_string(index++) + "]";
        rejectUnknownKeys(item, {"name", "graph"}, source, path);
        WorldStructureDefinition value;
        value.name = scalarString(requiredNode(item, "name", source, path), source, path + ".name");
        value.graph = scalarString(requiredNode(item, "graph", source, path), source, path + ".graph");
        result.structures.push_back(std::move(value));
    }

    const Node& graphs = requiredNode(root, "graphs", source, "$");
    requireSequence(graphs, source, "graphs");
    index = 0;
    for (const auto& item : graphs) {
        const std::string path = "graphs[" + std::to_string(index++) + "]";
        rejectUnknownKeys(item, {"id", "members", "links"}, source, path);
        WorldGraphDefinition value;
        value.id = scalarString(requiredNode(item, "id", source, path), source, path + ".id");
        value.members = stringSequence(requiredNode(item, "members", source, path), source, path + ".members");
        const Node& links = requiredNode(item, "links", source, path);
        requireSequence(links, source, path + ".links");
        std::size_t linkIndex = 0;
        for (const auto& link : links) {
            value.links.push_back(linkPair(link, source, path + ".links[" + std::to_string(linkIndex++) + "]"));
        }
        result.graphs.push_back(std::move(value));
    }

    result.roots = stringSequence(requiredNode(root, "roots", source, "$"), source, "roots");
    ValidateWorldDefinition(result);
    return result;
}

WorldDefinition LoadWorldDefinitionFile(const std::string& path) {
    const std::string extension = ".pyseb-world.yaml";
    if (path.size() < extension.size() || path.compare(path.size() - extension.size(), extension.size(), extension) != 0) {
        throw SEBException(path + ": world files must use the .pyseb-world.yaml extension", "LoadWorldDefinitionFile()");
    }
    std::ifstream input(path, std::ios::binary);
    if (!input) throw SEBException("Unable to open world file " + path, "LoadWorldDefinitionFile()");
    input.seekg(0, std::ios::end);
    const std::streamoff size = input.tellg();
    if (size < 0 || static_cast<std::size_t>(size) > kMaximumWorldFileSize) {
        throw SEBException(path + ": world file exceeds the 4 MiB limit", "LoadWorldDefinitionFile()");
    }
    input.seekg(0, std::ios::beg);
    return LoadWorldDefinitionYaml(
        std::string((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>()), path);
}

void ValidateWorldDefinition(const WorldDefinition& definition) {
    if (definition.schemaVersion != 1) throw SEBException("only world schema version 1 is supported", "ValidateWorldDefinition()");
    if (definition.worldId.empty()) throw SEBException("world id must not be empty", "ValidateWorldDefinition()");
    if (definition.subunits.empty() && definition.structures.empty()) throw SEBException("world must contain at least one object", "ValidateWorldDefinition()");

    std::set<std::string> names;
    std::set<std::string> subunitNames;
    for (const auto& subunit : definition.subunits) {
        if (!validName(subunit.name)) throw SEBException("invalid subunit name '" + subunit.name + "'", "ValidateWorldDefinition()");
        if (subunit.model.empty()) throw SEBException("subunit '" + subunit.name + "' has an empty model", "ValidateWorldDefinition()");
        if (!validName(subunit.tag)) throw SEBException("invalid tag '" + subunit.tag + "'", "ValidateWorldDefinition()");
        if (!names.insert(subunit.name).second) throw SEBException("duplicate object name '" + subunit.name + "'", "ValidateWorldDefinition()");
        subunitNames.insert(subunit.name);
    }
    std::map<std::string, std::string> structureGraphs;
    for (const auto& structure : definition.structures) {
        if (!validName(structure.name)) throw SEBException("invalid structure name '" + structure.name + "'", "ValidateWorldDefinition()");
        if (!names.insert(structure.name).second) throw SEBException("duplicate object name '" + structure.name + "'", "ValidateWorldDefinition()");
        structureGraphs.emplace(structure.name, structure.graph);
    }

    std::map<std::string, const WorldGraphDefinition*> graphs;
    std::map<std::string, std::string> outerGraph;
    for (const auto& graph : definition.graphs) {
        if (!validGraphId(graph.id)) throw SEBException("invalid graph id '" + graph.id + "'", "ValidateWorldDefinition()");
        if (!graphs.emplace(graph.id, &graph).second) throw SEBException("duplicate graph id '" + graph.id + "'", "ValidateWorldDefinition()");
        if (graph.members.empty()) throw SEBException("graph '" + graph.id + "' must contain at least one member", "ValidateWorldDefinition()");
        std::set<std::string> members;
        bool hasSubunit = false;
        bool hasStructure = false;
        for (const auto& member : graph.members) {
            if (!names.count(member)) throw SEBException("graph '" + graph.id + "' references unknown member '" + member + "'", "ValidateWorldDefinition()");
            if (!members.insert(member).second) throw SEBException("graph '" + graph.id + "' contains duplicate member '" + member + "'", "ValidateWorldDefinition()");
            const auto previousGraph = outerGraph.find(member);
            if (previousGraph != outerGraph.end()) throw SEBException("object '" + member + "' belongs to graphs '" + previousGraph->second + "' and '" + graph.id + "'", "ValidateWorldDefinition()");
            outerGraph.emplace(member, graph.id);
            hasSubunit = hasSubunit || subunitNames.count(member) != 0;
            hasStructure = hasStructure || structureGraphs.count(member) != 0;
        }
        if (hasSubunit && hasStructure) throw SEBException("graph '" + graph.id + "' cannot mix subunits and structures", "ValidateWorldDefinition()");
        std::set<std::pair<std::string, std::string>> linkSet;
        for (const auto& link : graph.links) {
            const std::string first = topName(link.first);
            const std::string second = topName(link.second);
            if (first == second) throw SEBException("graph '" + graph.id + "' links an object to itself", "ValidateWorldDefinition()");
            if (!members.count(first) || !members.count(second)) throw SEBException("graph '" + graph.id + "' link endpoint is not a member", "ValidateWorldDefinition()");
            const auto ordered = first < second ? std::make_pair(first, second) : std::make_pair(second, first);
            if (!linkSet.insert(ordered).second) throw SEBException("graph '" + graph.id + "' contains a duplicate link", "ValidateWorldDefinition()");
        }
        if (graph.links.size() != graph.members.size() - 1) throw SEBException("graph '" + graph.id + "' must contain exactly members - 1 links", "ValidateWorldDefinition()");
        std::set<std::string> visited;
        std::vector<std::string> pending(1, graph.members.front());
        while (!pending.empty()) {
            const std::string current = pending.back();
            pending.pop_back();
            if (!visited.insert(current).second) continue;
            for (const auto& link : graph.links) {
                if (topName(link.first) == current) pending.push_back(topName(link.second));
                if (topName(link.second) == current) pending.push_back(topName(link.first));
            }
        }
        if (visited.size() != graph.members.size()) throw SEBException("graph '" + graph.id + "' is disconnected", "ValidateWorldDefinition()");
    }
    if (outerGraph.size() != names.size()) throw SEBException("every subunit and structure must belong to exactly one graph", "ValidateWorldDefinition()");
    for (const auto& entry : structureGraphs) {
        if (!graphs.count(entry.second)) throw SEBException("structure '" + entry.first + "' references unknown graph '" + entry.second + "'", "ValidateWorldDefinition()");
    }

    struct DependencyFrame {
        std::string id;
        std::size_t nextMember = 0;
    };
    std::map<std::string, int> state;
    for (const auto& graph : graphs) {
        if (state[graph.first] == 2) continue;
        state[graph.first] = 1;
        std::vector<DependencyFrame> pending = {{graph.first, 0}};
        while (!pending.empty()) {
            DependencyFrame& frame = pending.back();
            const auto* current = graphs.at(frame.id);
            if (frame.nextMember == current->members.size()) {
                state[frame.id] = 2;
                pending.pop_back();
                continue;
            }

            const std::string& member = current->members[frame.nextMember++];
            const auto structure = structureGraphs.find(member);
            if (structure == structureGraphs.end()) continue;

            const std::string& dependency = structure->second;
            if (state[dependency] == 1) {
                throw SEBException(
                    "cyclic graph dependency at '" + dependency + "'",
                    "ValidateWorldDefinition()");
            }
            if (state[dependency] == 2) continue;
            state[dependency] = 1;
            pending.push_back({dependency, 0});
        }
    }

    std::set<std::string> roots;
    for (const auto& root : definition.roots) {
        if (!roots.insert(root).second) throw SEBException("duplicate root graph '" + root + "'", "ValidateWorldDefinition()");
        if (!graphs.count(root)) throw SEBException("unknown root graph '" + root + "'", "ValidateWorldDefinition()");
    }
    if (roots.empty()) throw SEBException("world must declare at least one root graph", "ValidateWorldDefinition()");
    std::set<std::string> referenced;
    for (const auto& structure : structureGraphs) referenced.insert(structure.second);
    for (const auto& root : roots) {
        if (referenced.count(root)) throw SEBException("root graph '" + root + "' is also wrapped by a structure", "ValidateWorldDefinition()");
    }
    std::set<std::string> reachable;
    std::vector<std::string> queue(definition.roots.begin(), definition.roots.end());
    while (!queue.empty()) {
        const std::string id = queue.back(); queue.pop_back();
        if (!reachable.insert(id).second) continue;
        for (const auto& member : graphs.at(id)->members) {
            const auto structure = structureGraphs.find(member);
            if (structure != structureGraphs.end()) queue.push_back(structure->second);
        }
    }
    if (reachable.size() != graphs.size()) throw SEBException("world contains an unreachable graph", "ValidateWorldDefinition()");
}

std::string SerializeWorldDefinition(const WorldDefinition& input) {
    WorldDefinition definition = input;
    ValidateWorldDefinition(definition);
    sortBy(definition.subunits, [](const auto& a, const auto& b) { return a.name < b.name; });
    sortBy(definition.structures, [](const auto& a, const auto& b) { return a.name < b.name; });
    sortBy(definition.graphs, [](const auto& a, const auto& b) { return a.id < b.id; });
    for (auto& graph : definition.graphs) {
        std::sort(graph.members.begin(), graph.members.end());
        std::sort(graph.links.begin(), graph.links.end());
    }
    std::sort(definition.roots.begin(), definition.roots.end());

    std::ostringstream output;
    output << "format: pyseb-world\n";
    output << "schema_version: 1\n\n";
    output << "world:\n  id: " << quote(definition.worldId) << "\n\n";
    if (definition.subunits.empty()) {
        output << "subunits: []\n";
    } else {
        output << "subunits:\n";
        for (const auto& subunit : definition.subunits) {
            output << "  - name: " << quote(subunit.name) << "\n"
                   << "    model: " << quote(subunit.model) << "\n"
                   << "    tag: " << quote(subunit.tag) << "\n";
        }
    }
    output << "\n";
    if (definition.structures.empty()) {
        output << "structures: []\n";
    } else {
        output << "structures:\n";
        for (const auto& structure : definition.structures) {
            output << "  - name: " << quote(structure.name) << "\n"
                   << "    graph: " << quote(structure.graph) << "\n";
        }
    }
    output << "\ngraphs:\n";
    for (const auto& graph : definition.graphs) {
        output << "  - id: " << quote(graph.id) << "\n    members: [";
        for (std::size_t i = 0; i < graph.members.size(); ++i) {
            if (i) output << ", ";
            output << quote(graph.members[i]);
        }
        if (graph.links.empty()) {
            output << "]\n    links: []\n";
        } else {
            output << "]\n    links:\n";
            for (const auto& link : graph.links) output << "      - [" << quote(link.first) << ", " << quote(link.second) << "]\n";
        }
    }
    output << "\nroots: [";
    for (std::size_t i = 0; i < definition.roots.size(); ++i) {
        if (i) output << ", ";
        output << quote(definition.roots[i]);
    }
    output << "]\n";
    return output.str();
}

void SaveWorldDefinitionFile(const WorldDefinition& definition, const std::string& path) {
    const std::string extension = ".pyseb-world.yaml";
    if (path.size() < extension.size() || path.compare(path.size() - extension.size(), extension.size(), extension) != 0) {
        throw SEBException(path + ": world files must use the .pyseb-world.yaml extension", "SaveWorldDefinitionFile()");
    }
    // Validate and render before opening the destination. This prevents a failed
    // serialization from truncating an existing world file.
    const std::string serialized = SerializeWorldDefinition(definition);
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) throw SEBException("Unable to open world file " + path, "SaveWorldDefinitionFile()");
    output << serialized;
    if (!output) throw SEBException("Unable to write world file " + path, "SaveWorldDefinitionFile()");
}

} // namespace pyseb

namespace pyseb {
namespace {

std::string endpointObject(const std::string& endpoint) {
    const std::size_t colon = endpoint.find(':');
    const std::size_t period = endpoint.find('.');
    std::size_t cut = std::string::npos;
    if (colon != std::string::npos) cut = colon;
    if (period != std::string::npos) cut = std::min(cut, period);
    return cut == std::string::npos ? endpoint : endpoint.substr(0, cut);
}

const WorldGraphDefinition& graphFor(
    const std::map<std::string, const WorldGraphDefinition*>& graphs,
    const std::string& id) {
    const auto found = graphs.find(id);
    if (found == graphs.end()) throw SEBException("Unknown graph " + id, "BuildWorld()");
    return *found->second;
}

} // namespace

std::unique_ptr<World> BuildWorld(
    const WorldDefinition& definition,
    const std::vector<std::string>& modelFiles) {
    ValidateWorldDefinition(definition);

    std::unique_ptr<World> result(new World(definition.worldId));
    for (const auto& file : modelFiles) result->RegisterSubunitFile(file);

    std::map<std::string, const WorldSubunitDefinition*> subunits;
    for (const auto& subunit : definition.subunits) subunits.emplace(subunit.name, &subunit);
    std::map<std::string, const WorldStructureDefinition*> structures;
    for (const auto& structure : definition.structures) structures.emplace(structure.name, &structure);
    std::map<std::string, const WorldGraphDefinition*> graphs;
    for (const auto& graph : definition.graphs) graphs.emplace(graph.id, &graph);

    const auto ensureSerializable = [&](const std::string& name) {
        if (result->getSubunit(name)->getSubunitType() != FILEDEFINEDSUBUNIT) {
            throw SEBException(
                "Subunit model for '" + name + "' is not representable in world schema version 1",
                "BuildWorld()");
        }
    };

    struct BuildFrame {
        std::string id;
        std::size_t nextMember = 0;
    };
    std::map<std::string, int> buildState;
    std::vector<std::string> buildOrder;
    for (const auto& root : definition.roots) {
        if (buildState[root] == 2) continue;
        buildState[root] = 1;
        std::vector<BuildFrame> pending = {{root, 0}};
        while (!pending.empty()) {
            BuildFrame& frame = pending.back();
            const WorldGraphDefinition& graph = graphFor(graphs, frame.id);
            if (frame.nextMember == graph.members.size()) {
                buildState[frame.id] = 2;
                buildOrder.push_back(frame.id);
                pending.pop_back();
                continue;
            }

            const std::string& member = graph.members[frame.nextMember++];
            const auto structure = structures.find(member);
            if (structure == structures.end()) continue;
            const std::string& dependency = structure->second->graph;
            if (buildState[dependency] == 1) {
                throw SEBException(
                    "cyclic graph dependency at '" + dependency + "'",
                    "BuildWorld()");
            }
            if (buildState[dependency] == 2) continue;
            buildState[dependency] = 1;
            pending.push_back({dependency, 0});
        }
    }

    std::map<std::string, GraphID> built;
    for (const auto& id : buildOrder) {
        const WorldGraphDefinition& graph = graphFor(graphs, id);

        bool containsStructures = false;
        for (const auto& member : graph.members) containsStructures = containsStructures || structures.count(member) != 0;

        std::set<std::string> builtMembers;
        GraphID graphId = 0;
        if (!containsStructures) {
            const auto first = subunits.find(graph.members.front());
            if (first == subunits.end()) throw SEBException("Graph " + id + " does not contain subunits", "BuildWorld()");
            graphId = result->Add(first->second->model, first->second->name, first->second->tag);
            ensureSerializable(first->second->name);
            builtMembers.insert(first->second->name);
            while (builtMembers.size() < graph.members.size()) {
                bool progressed = false;
                for (const auto& link : graph.links) {
                    const std::string left = endpointObject(link.first);
                    const std::string right = endpointObject(link.second);
                    std::string newMember;
                    std::string newReference;
                    std::string oldReference;
                    if (builtMembers.count(left) && !builtMembers.count(right)) {
                        newMember = right; newReference = link.second; oldReference = link.first;
                    } else if (builtMembers.count(right) && !builtMembers.count(left)) {
                        newMember = left; newReference = link.first; oldReference = link.second;
                    } else continue;
                    const auto item = subunits.find(newMember);
                    if (item == subunits.end()) throw SEBException("Graph " + id + " has an invalid subunit link", "BuildWorld()");
                    result->Link(item->second->model, newReference, oldReference, item->second->tag);
                    ensureSerializable(newMember);
                    builtMembers.insert(newMember);
                    progressed = true;
                }
                if (!progressed) throw SEBException("Unable to order graph " + id + " for construction", "BuildWorld()");
            }
        } else {
            const auto first = structures.find(graph.members.front());
            if (first == structures.end()) throw SEBException("Graph " + id + " does not contain structures", "BuildWorld()");
            const GraphID inner = built.at(first->second->graph);
            graphId = result->Add(inner, first->second->name);
            builtMembers.insert(first->second->name);
            while (builtMembers.size() < graph.members.size()) {
                bool progressed = false;
                for (const auto& link : graph.links) {
                    const std::string left = endpointObject(link.first);
                    const std::string right = endpointObject(link.second);
                    std::string newMember;
                    std::string newReference;
                    std::string oldReference;
                    if (builtMembers.count(left) && !builtMembers.count(right)) {
                        newMember = right; newReference = link.second; oldReference = link.first;
                    } else if (builtMembers.count(right) && !builtMembers.count(left)) {
                        newMember = left; newReference = link.first; oldReference = link.second;
                    } else continue;
                    const auto item = structures.find(newMember);
                    if (item == structures.end()) throw SEBException("Graph " + id + " has an invalid structure link", "BuildWorld()");
                    const GraphID memberGraph = built.at(item->second->graph);
                    result->Link(memberGraph, newReference, oldReference);
                    builtMembers.insert(newMember);
                    progressed = true;
                }
                if (!progressed) throw SEBException("Unable to order graph " + id + " for construction", "BuildWorld()");
            }
        }

        result->SetGraphLabel(graphId, id);
        built.emplace(id, graphId);
    }

    return result;
}

std::unique_ptr<World> LoadWorld(const std::string& path, const std::vector<std::string>& modelFiles) {
    return BuildWorld(LoadWorldDefinitionFile(path), modelFiles);
}

std::unique_ptr<World> WorldFromYaml(const std::string& yaml, const std::vector<std::string>& modelFiles) {
    return BuildWorld(LoadWorldDefinitionYaml(yaml), modelFiles);
}

void SaveWorld(const World& world, const std::string& path) {
    SaveWorldDefinitionFile(world.Describe(), path);
}

std::string WorldToYaml(const World& world) {
    return SerializeWorldDefinition(world.Describe());
}

} // namespace pyseb
