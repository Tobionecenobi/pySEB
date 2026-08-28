#ifndef PYSEB_WORLD_IO_WORLD_PARSER_HPP
#define PYSEB_WORLD_IO_WORLD_PARSER_HPP

#include <string>
#include <vector>
#include <memory>

#include "WorldIO/WorldDefinition.hpp"

class World;

namespace pyseb {

WorldDefinition LoadWorldDefinitionFile(const std::string& path);
WorldDefinition LoadWorldDefinitionYaml(
    const std::string& yaml,
    const std::string& source = "<string>");

void ValidateWorldDefinition(const WorldDefinition& definition);

std::string SerializeWorldDefinition(const WorldDefinition& definition);
void SaveWorldDefinitionFile(
    const WorldDefinition& definition,
    const std::string& path);

std::unique_ptr<World> BuildWorld(
    const WorldDefinition& definition,
    const std::vector<std::string>& modelFiles = {});
std::unique_ptr<World> LoadWorld(
    const std::string& path,
    const std::vector<std::string>& modelFiles = {});
std::unique_ptr<World> WorldFromYaml(
    const std::string& yaml,
    const std::vector<std::string>& modelFiles = {});
void SaveWorld(const World& world, const std::string& path);
std::string WorldToYaml(const World& world);

} // namespace pyseb

#endif
