#ifndef PYSEB_SPHERE_AMPLITUDE_HPP
#define PYSEB_SPHERE_AMPLITUDE_HPP

#include <cmath>

namespace pyseb {

// Normalized solid-sphere amplitude 3 (sin(x) - x cos(x)) / x^3.
// The direct numerator loses precision near zero. Through x^8, the first
// omitted term at |x| = 0.1 is below 6e-19.
inline double SphereAmplitudeValue(double x) {
    if (std::abs(x) < 0.1) {
        const double x2 = x * x;
        return 1.0 + x2 * (-1.0 / 10.0 + x2 * (1.0 / 280.0
            + x2 * (-1.0 / 15120.0 + x2 / 1330560.0)));
    }
    return 3.0 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
}

} // namespace pyseb

#endif
