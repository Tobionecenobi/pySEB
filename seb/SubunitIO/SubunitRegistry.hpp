#ifndef PYSEB_SUBUNIT_IO_SUBUNIT_REGISTRY_HPP
#define PYSEB_SUBUNIT_IO_SUBUNIT_REGISTRY_HPP

#include <map>
#include <set>
#include <string>
#include <vector>

#include "SubunitIO/SubunitDefinition.hpp"

class SubUnit;

namespace pyseb {

class SubunitRegistry {
public:
    SubunitRegistry();

    SubunitModelInfo RegisterFile(const std::string& path);
    std::vector<SubunitModelInfo> RegisterDirectory(const std::string& path);
    SubUnit* Create(const std::string& id) const;
    bool Has(const std::string& id) const;
    std::vector<SubunitModelInfo> List() const;

private:
    std::map<std::string, SubunitDefinition> definitions_;
    std::map<std::string, std::string> aliases_;

    static bool IsBuiltinName(const std::string& id);
    static SubunitModelInfo Info(const SubunitDefinition& definition, bool bundled);
    bool HasRegisteredName(const std::string& name) const;
    void CheckDefinitionNames(
        const SubunitDefinition& definition,
        const std::set<std::string>& pendingNames,
        const std::string& operation) const;
};

} // namespace pyseb

#endif
