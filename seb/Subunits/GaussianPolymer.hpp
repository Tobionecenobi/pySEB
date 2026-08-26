#ifndef INCLUDE_GUARD_GAUSSIANPOLYMER
#define INCLUDE_GUARD_GAUSSIANPOLYMER

#include "FileDefined.hpp"
#include "SubunitIO/BundledModels.hpp"

class GaussianPolymer : public FileDefinedSubunit {
public:
    GaussianPolymer()
        : FileDefinedSubunit(pyseb::BundledSubunitDefinition("pyseb/GaussianPolymer")) {
        stype = GAUSSIANPOLYMER;
        setName("GaussianPolymer");
    }
    ~GaussianPolymer() override = default;
};

#endif
