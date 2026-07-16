#ifndef INCLUDE_GUARD_ATOMRECORD
#define INCLUDE_GUARD_ATOMRECORD

#include <string>

#include "Geometry/Point3D.hpp"

struct AtomRecord {
    std::string recordType;
    int serial = 0;
    std::string atomName;
    std::string alternateLocation;
    std::string residueName;
    std::string chainId;
    int residueNumber = 0;
    std::string insertionCode;
    CartesianPoint3D coordinate;
    double occupancy = 1.0;
    double temperatureFactor = 0.0;
    std::string element;
    std::string charge;
    int modelNumber = 1;
};

#endif
