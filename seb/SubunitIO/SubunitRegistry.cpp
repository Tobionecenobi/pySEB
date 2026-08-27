#include "SubunitIO/SubunitRegistry.hpp"

#include <algorithm>
#include <set>

#ifdef _WIN32
#define POINT PYSEB_WIN32_POINT
#include <windows.h>
#undef POINT
#else
#include <dirent.h>
#endif

#include "Exceptions.hpp"
#include "ExpressionFunctions.hpp"
#include "SpecialFunctions.hpp"
#include "Subunit.hpp"
#include "SubunitIO/BundledModels.hpp"
#include "Subunits/Subunits.hpp"
#include "Subunits/CreateSubunit.hpp"

namespace pyseb {
namespace {

bool hasModelExtension(const std::string& name) {
    static const std::string extension = ".pyseb.yaml";
    return name.size() >= extension.size() &&
        name.compare(name.size() - extension.size(), extension.size(), extension) == 0;
}

std::string joinPath(const std::string& directory, const std::string& name) {
    if (directory.empty()) return name;
    const char last = directory.back();
    if (last == '/' || last == '\\') return directory + name;
#ifdef _WIN32
    return directory + "\\" + name;
#else
    return directory + "/" + name;
#endif
}

std::vector<std::string> modelFiles(const std::string& directory) {
    std::vector<std::string> result;
#ifdef _WIN32
    const std::string pattern = joinPath(directory, "*");
    WIN32_FIND_DATAA data;
    HANDLE handle = FindFirstFileA(pattern.c_str(), &data);
    if (handle == INVALID_HANDLE_VALUE) {
        throw SEBException("Unable to open subunit model directory " + directory, "SubunitRegistry::RegisterDirectory()");
    }
    do {
        const std::string name = data.cFileName;
        if (!(data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) && hasModelExtension(name)) {
            result.push_back(joinPath(directory, name));
        }
    } while (FindNextFileA(handle, &data));
    FindClose(handle);
#else
    DIR* handle = opendir(directory.c_str());
    if (!handle) {
        throw SEBException("Unable to open subunit model directory " + directory, "SubunitRegistry::RegisterDirectory()");
    }
    while (dirent* entry = readdir(handle)) {
        const std::string name = entry->d_name;
        if (hasModelExtension(name)) result.push_back(joinPath(directory, name));
    }
    closedir(handle);
#endif
    std::sort(result.begin(), result.end());
    return result;
}

const std::vector<std::string>& builtinNames() {
    // Every bundled analytic model is now a FileDefinedSubunit YAML model.
    static const std::vector<std::string> names = {
        "SymbolicSubunit"
    };
    return names;
}

} // namespace

SubunitRegistry::SubunitRegistry() {
    for (const auto& model : BundledSubunitModels()) {
        const SubunitDefinition& definition = BundledSubunitDefinition(model.id);
        CheckDefinitionNames(definition, {}, "SubunitRegistry()");
        definitions_.emplace(definition.id, definition);
        aliases_.emplace(definition.apiName, definition.id);
    }
}

SubunitModelInfo SubunitRegistry::Info(const SubunitDefinition& definition, bool bundled) {
    SubunitModelInfo info;
    info.id = definition.id;
    info.apiName = definition.apiName;
    info.modelVersion = definition.modelVersion;
    info.title = definition.metadata.title;
    info.source = definition.source;
    info.bundled = bundled;
    return info;
}

bool SubunitRegistry::IsBuiltinName(const std::string& id) {
    const auto& names = builtinNames();
    return std::find(names.begin(), names.end(), id) != names.end();
}

bool SubunitRegistry::HasRegisteredName(const std::string& name) const {
    return definitions_.count(name) != 0 || aliases_.count(name) != 0 || IsBuiltinName(name);
}

void SubunitRegistry::CheckDefinitionNames(
    const SubunitDefinition& definition,
    const std::set<std::string>& pendingNames,
    const std::string& operation) const {
    for (const auto& name : {definition.id, definition.apiName}) {
        if (HasRegisteredName(name) || pendingNames.count(name)) {
            throw SEBException(
                "Subunit model name '" + name + "' is already registered",
                operation);
        }
    }
}

SubunitModelInfo SubunitRegistry::RegisterFile(const std::string& path) {
    SubunitDefinition definition = LoadSubunitDefinitionFile(path);
    CheckDefinitionNames(definition, {}, "SubunitRegistry::RegisterFile()");
    const SubunitModelInfo info = Info(definition, false);
    aliases_.emplace(definition.apiName, definition.id);
    definitions_.emplace(definition.id, std::move(definition));
    return info;
}

std::vector<SubunitModelInfo> SubunitRegistry::RegisterDirectory(const std::string& path) {
    std::vector<SubunitDefinition> pending;
    std::set<std::string> pendingNames;
    for (const auto& file : modelFiles(path)) {
        SubunitDefinition definition = LoadSubunitDefinitionFile(file);
        CheckDefinitionNames(definition, pendingNames, "SubunitRegistry::RegisterDirectory()");
        pendingNames.insert(definition.id);
        pendingNames.insert(definition.apiName);
        pending.push_back(std::move(definition));
    }

    std::vector<SubunitModelInfo> result;
    result.reserve(pending.size());
    for (auto& definition : pending) {
        result.push_back(Info(definition, false));
        aliases_.emplace(definition.apiName, definition.id);
        definitions_.emplace(definition.id, std::move(definition));
    }
    return result;
}

SubUnit* SubunitRegistry::Create(const std::string& id) const {
    const auto alias = aliases_.find(id);
    const std::string& canonicalId = alias == aliases_.end() ? id : alias->second;
    const auto definition = definitions_.find(canonicalId);
    if (definition != definitions_.end()) return new FileDefinedSubunit(definition->second);
    return CreateSubunit(id);
}

bool SubunitRegistry::Has(const std::string& id) const {
    return HasRegisteredName(id);
}

std::vector<SubunitModelInfo> SubunitRegistry::List() const {
    std::vector<SubunitModelInfo> result;
    for (const auto& name : builtinNames()) {
        SubunitModelInfo info;
        info.id = name;
        info.apiName = name;
        info.modelVersion = "builtin";
        info.title = name;
        info.source = "<builtin>";
        info.bundled = true;
        result.push_back(info);
    }
    for (const auto& definition : definitions_) {
        result.push_back(Info(definition.second, definition.second.source.find("<bundled:") == 0));
    }
    std::sort(result.begin(), result.end(), [](const SubunitModelInfo& lhs, const SubunitModelInfo& rhs) {
        return lhs.id < rhs.id;
    });
    return result;
}

} // namespace pyseb
