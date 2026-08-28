#ifndef PYSEB_USD_EXPORT_HPP
#define PYSEB_USD_EXPORT_HPP

#include <array>
#include <cstddef>
#include <map>
#include <string>

/// Units accepted by the text USD exporter.
enum class LengthUnit { Meter, Millimeter, Micrometer, Nanometer, Angstrom };
/// Orientation policy for one representative exported realization.
enum class USDLayoutMode { Random, Readable };

struct USDExportOptions {
    // Units and deterministic realization controls.
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
    // Sampling and optional diagnostic/reference display.
    std::size_t curveSamples = 96;
    std::size_t surfaceSamples = 48;
    bool referenceMarkers = false;
    double zeroRadiusMarkerSize = 0.01;
    // Debye-cloud envelope settings. Padding -1 selects automatic connectivity padding.
    bool debyeEnvelope = true;
    std::size_t debyeEnvelopeResolution = 24;
    double debyeEnvelopePadding = -1.0;
    double debyeEnvelopeOpacity = 0.15;
    // Per-instance display overrides, keyed by the complete instance path.
    std::map<std::string, std::array<double,3>> colorOverrides;
    std::map<std::string, double> opacityOverrides;
};

double MetersPerUnit(LengthUnit unit);

#endif
