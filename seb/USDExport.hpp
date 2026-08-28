#ifndef PYSEB_USD_EXPORT_HPP
#define PYSEB_USD_EXPORT_HPP

#include <map>
#include <string>
#include <array>

enum class LengthUnit { Meter, Millimeter, Micrometer, Nanometer, Angstrom };
enum class USDLayoutMode { Random, Readable };

struct USDExportOptions {
    explicit USDExportOptions(LengthUnit unit);
    explicit USDExportOptions(double customMetersPerUnit);
    LengthUnit unit = LengthUnit::Meter;
    double metersPerUnit = 1.0;
    bool customUnit = false;
    unsigned long long seed = 0;
    USDLayoutMode layoutMode = USDLayoutMode::Random;
    std::size_t orientationTrials = 64;
    std::size_t relaxationSweeps = 4;
    double minimumClearance = 0.0;
    std::size_t curveSamples = 96;
    std::size_t surfaceSamples = 48;
    bool referenceMarkers = false;
    double zeroRadiusMarkerSize = 0.01;
    std::map<std::string, std::array<double,3>> colorOverrides;
    std::map<std::string, double> opacityOverrides;
};

double MetersPerUnit(LengthUnit unit);

#endif
