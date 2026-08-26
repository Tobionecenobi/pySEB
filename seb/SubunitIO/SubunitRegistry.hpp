#ifndef PYSEB_SUBUNIT_IO_SUBUNIT_REGISTRY_HPP
#define PYSEB_SUBUNIT_IO_SUBUNIT_REGISTRY_HPP

#include <map>
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

    static bool IsBuiltinName(const std::string& id);
    static SubunitModelInfo Info(const SubunitDefinition& definition, bool bundled);
};

} // namespace pyseb

#endif
