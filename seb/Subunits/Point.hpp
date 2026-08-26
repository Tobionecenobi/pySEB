#ifndef INCLUDE_GUARD_POINT
#define INCLUDE_GUARD_POINT

#include "FileDefined.hpp"
#include "SubunitIO/BundledModels.hpp"

class Point : public FileDefinedSubunit {
public:
    Point() : FileDefinedSubunit(pyseb::BundledSubunitDefinition("pyseb/Point")) {
        stype = POINT;
        setName("Point");
    }
    ~Point() override = default;
};

#endif
