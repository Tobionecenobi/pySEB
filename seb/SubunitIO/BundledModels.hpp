#ifndef PYSEB_SUBUNIT_IO_BUNDLED_MODELS_HPP
#define PYSEB_SUBUNIT_IO_BUNDLED_MODELS_HPP

#include <string>
#include <vector>

#include "SubunitIO/SubunitDefinition.hpp"

namespace pyseb {

const SubunitDefinition& BundledSubunitDefinition(const std::string& id);
std::vector<SubunitModelInfo> BundledSubunitModels();

} // namespace pyseb

#endif
