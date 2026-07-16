#ifndef INCLUDE_GUARD_NUMERICALSUBUNIT
#define INCLUDE_GUARD_NUMERICALSUBUNIT

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <utility>

#include "Symbolic.hpp"

enum class NormalizationMode {
    Normalized,
    Unnormalized
};

class NumericalSubunit : public SymbolicSubunit {
public:
    using QFunction = std::function<double(double, const ParameterList&)>;
    using ScalarProvider = std::function<double(const ParameterList&)>;

    explicit NumericalSubunit(NormalizationMode mode = NormalizationMode::Normalized)
        : SymbolicSubunit(), normalizationMode(mode)
    {
        stype = NUMERICALSUBUNIT;
    }

    virtual ~NumericalSubunit() {}

    NormalizationMode getNormalizationMode() const {
        return normalizationMode;
    }

    void addReferencePoint(const refPoint& ref) {
        setReferencePointName(ref);
    }

    void addDistributedReferencePointType(const refPoint& ref) {
        setDistReferencePointType(ref);
    }

    void setTotalBetaProvider(ScalarProvider provider) {
        totalBetaProvider = std::move(provider);
    }

    void setTotalBeta(double beta) {
        setTotalBetaProvider([beta](const ParameterList&) { return beta; });
    }

    void setFormFactorFunction(QFunction function) {
        formFactorFunction = std::move(function);
    }

    void setFormFactorAmplitudeFunction(const refPoint& ref, QFunction function) {
        requireKnownReference(ref, "NumericalSubunit::setFormFactorAmplitudeFunction()");
        amplitudeFunctions[canonicalReference(ref)] = std::move(function);
    }

    void setPhaseFactorFunction(refPoint r1, refPoint r2, QFunction function) {
        requireKnownReference(r1, "NumericalSubunit::setPhaseFactorFunction()");
        requireKnownReference(r2, "NumericalSubunit::setPhaseFactorFunction()");
        if (r1 == r2 && hasSpecificReference(r1)) {
            throw SEBException(
                "A phase-factor callback is not needed for an identical specific reference point",
                "NumericalSubunit::setPhaseFactorFunction()"
            );
        }
        phaseFunctions[referencePair(r1, r2)] = std::move(function);
    }

    void setRadiusOfGyration2Provider(ScalarProvider provider) {
        radiusOfGyration2Provider = std::move(provider);
    }

    void setRadiusOfGyration2(double value) {
        setRadiusOfGyration2Provider([value](const ParameterList&) { return value; });
    }

    void setSigmaMSDRef2ScatProvider(const refPoint& ref, ScalarProvider provider) {
        requireKnownReference(ref, "NumericalSubunit::setSigmaMSDRef2ScatProvider()");
        sigmaRefToScattererProviders[canonicalReference(ref)] = std::move(provider);
    }

    void setSigmaMSDRef2Scat(const refPoint& ref, double value) {
        setSigmaMSDRef2ScatProvider(
            ref,
            [value](const ParameterList&) { return value; }
        );
    }

    void setSigmaMSDRef2RefProvider(refPoint r1, refPoint r2, ScalarProvider provider) {
        requireKnownReference(r1, "NumericalSubunit::setSigmaMSDRef2RefProvider()");
        requireKnownReference(r2, "NumericalSubunit::setSigmaMSDRef2RefProvider()");
        if (r1 == r2 && hasSpecificReference(r1)) {
            throw SEBException(
                "The mean-square distance from a specific reference point to itself is zero by definition",
                "NumericalSubunit::setSigmaMSDRef2RefProvider()"
            );
        }
        sigmaRefToRefProviders[referencePair(r1, r2)] = std::move(provider);
    }

    void setSigmaMSDRef2Ref(refPoint r1, refPoint r2, double value) {
        setSigmaMSDRef2RefProvider(
            r1,
            r2,
            [value](const ParameterList&) { return value; }
        );
    }

    double NumericTotalBeta(const ParameterList& values) override {
        requireProvider(totalBetaProvider, "total beta");
        return finiteValue(totalBetaProvider(values), "total beta");
    }

    double NumericFormFactorUnnormalized(
        double q,
        const ParameterList& values) override
    {
        requireQFunction(formFactorFunction, "form factor");
        const double raw = finiteValue(
            formFactorFunction(finiteQ(q), values),
            "form factor"
        );
        if (normalizationMode == NormalizationMode::Unnormalized) {
            return raw;
        }

        const double beta = nonzeroBeta(values, "normalized form factor");
        return beta * beta * raw;
    }

    double NumericFormFactorAmplitudeUnnormalized(
        refPoint r,
        double q,
        const ParameterList& values) override
    {
        requireKnownReference(r, "NumericalSubunit::NumericFormFactorAmplitudeUnnormalized()");
        r = canonicalReference(r);
        const auto function = amplitudeFunctions.find(r);
        if (function == amplitudeFunctions.end() || !function->second) {
            throw SEBException(
                "Missing numerical form-factor amplitude callback for reference point " + r,
                "NumericalSubunit::NumericFormFactorAmplitudeUnnormalized()"
            );
        }

        const double raw = finiteValue(
            function->second(finiteQ(q), values),
            "form-factor amplitude"
        );
        if (normalizationMode == NormalizationMode::Unnormalized) {
            return raw;
        }

        return nonzeroBeta(values, "normalized form-factor amplitude") * raw;
    }

    double NumericPhaseFactor(
        refPoint r1,
        refPoint r2,
        double q,
        const ParameterList& values) override
    {
        requireKnownReference(r1, "NumericalSubunit::NumericPhaseFactor()");
        requireKnownReference(r2, "NumericalSubunit::NumericPhaseFactor()");
        if (r1 == r2 && hasSpecificReference(r1)) {
            return 1.0;
        }

        const auto pair = referencePair(r1, r2);
        const auto function = phaseFunctions.find(pair);
        if (function == phaseFunctions.end() || !function->second) {
            throw SEBException(
                "Missing numerical phase-factor callback for reference points " +
                    pair.first + " and " + pair.second,
                "NumericalSubunit::NumericPhaseFactor()"
            );
        }
        return finiteValue(
            function->second(finiteQ(q), values),
            "phase factor"
        );
    }

    double NumericRadiusOfGyration2(const ParameterList& values) override {
        requireProvider(radiusOfGyration2Provider, "radius of gyration squared");
        return finiteValue(
            radiusOfGyration2Provider(values),
            "radius of gyration squared"
        );
    }

    double NumericSigmaMSDRef2Scat(
        refPoint r,
        const ParameterList& values) override
    {
        requireKnownReference(r, "NumericalSubunit::NumericSigmaMSDRef2Scat()");
        r = canonicalReference(r);
        const auto provider = sigmaRefToScattererProviders.find(r);
        if (provider == sigmaRefToScattererProviders.end() || !provider->second) {
            throw SEBException(
                "Missing sigma<R^2> provider for reference point " + r,
                "NumericalSubunit::NumericSigmaMSDRef2Scat()"
            );
        }
        return finiteValue(provider->second(values), "reference-to-scatterer sigma<R^2>");
    }

    double NumericSigmaMSDRef2Ref(
        refPoint r1,
        refPoint r2,
        const ParameterList& values) override
    {
        requireKnownReference(r1, "NumericalSubunit::NumericSigmaMSDRef2Ref()");
        requireKnownReference(r2, "NumericalSubunit::NumericSigmaMSDRef2Ref()");
        if (r1 == r2 && hasSpecificReference(r1)) {
            return 0.0;
        }

        const auto pair = referencePair(r1, r2);
        const auto provider = sigmaRefToRefProviders.find(pair);
        if (provider == sigmaRefToRefProviders.end() || !provider->second) {
            throw SEBException(
                "Missing reference-to-reference sigma<R^2> provider for " +
                    pair.first + " and " + pair.second,
                "NumericalSubunit::NumericSigmaMSDRef2Ref()"
            );
        }
        return finiteValue(provider->second(values), "reference-to-reference sigma<R^2>");
    }

    void ValidateNumerically(
        const ParameterList& values = ParameterList(),
        double tolerance = 1e-8)
    {
        if (tolerance <= 0.0 || !std::isfinite(tolerance)) {
            throw SEBException(
                "Validation tolerance must be finite and positive",
                "NumericalSubunit::ValidateNumerically()"
            );
        }

        const double beta = NumericTotalBeta(values);
        const double formFactorAtZero = NumericFormFactorUnnormalized(0.0, values);
        validateNear(
            formFactorAtZero,
            beta * beta,
            tolerance,
            "form factor at q=0"
        );
        NumericRadiusOfGyration2(values);

        ReferencePointSet references(refsSpecific);
        references.insert(refsDistributed.begin(), refsDistributed.end());

        for (const auto& r : references) {
            validateNear(
                NumericFormFactorAmplitudeUnnormalized(r, 0.0, values),
                beta,
                tolerance,
                "form-factor amplitude at q=0 for " + r
            );
            NumericSigmaMSDRef2Scat(r, values);
        }

        for (auto r1 = references.begin(); r1 != references.end(); ++r1) {
            for (auto r2 = r1; r2 != references.end(); ++r2) {
                if (*r1 == *r2 && hasSpecificReference(*r1)) {
                    continue;
                }
                validateNear(
                    NumericPhaseFactor(*r1, *r2, 0.0, values),
                    1.0,
                    tolerance,
                    "phase factor at q=0 for " + *r1 + " and " + *r2
                );
                NumericSigmaMSDRef2Ref(*r1, *r2, values);
            }
        }
    }

protected:
    NormalizationMode normalizationMode;
    ScalarProvider totalBetaProvider;
    QFunction formFactorFunction;
    std::map<refPoint, QFunction> amplitudeFunctions;
    std::map<std::pair<refPoint, refPoint>, QFunction> phaseFunctions;
    ScalarProvider radiusOfGyration2Provider;
    std::map<refPoint, ScalarProvider> sigmaRefToScattererProviders;
    std::map<std::pair<refPoint, refPoint>, ScalarProvider> sigmaRefToRefProviders;

    refPoint canonicalReference(refPoint r) const {
        const auto hash = r.find("#");
        if (hash != string::npos) {
            r = r.substr(0, hash);
        }
        return r;
    }

    std::pair<refPoint, refPoint> referencePair(refPoint r1, refPoint r2) const {
        r1 = canonicalReference(r1);
        r2 = canonicalReference(r2);
        if (r2 < r1) {
            std::swap(r1, r2);
        }
        return std::make_pair(r1, r2);
    }

    void requireKnownReference(const refPoint& ref, const string& where) const {
        const refPoint canonical = canonicalReference(ref);
        if (refsSpecific.find(ref) == refsSpecific.end() &&
            refsSpecific.find(canonical) == refsSpecific.end() &&
            refsDistributed.find(canonical) == refsDistributed.end()) {
            throw SEBException(
                "Unknown numerical reference point " + ref,
                where
            );
        }
    }

private:
    static double finiteQ(double q) {
        if (!std::isfinite(q)) {
            throw SEBException(
                "q must be finite",
                "NumericalSubunit numerical evaluation"
            );
        }
        return q;
    }

    static double finiteValue(double value, const string& what) {
        if (!std::isfinite(value)) {
            throw SEBException(
                "Numerical " + what + " produced a non-finite value",
                "NumericalSubunit numerical evaluation"
            );
        }
        return value;
    }

    static void requireProvider(const ScalarProvider& provider, const string& what) {
        if (!provider) {
            throw SEBException(
                "Missing numerical provider for " + what,
                "NumericalSubunit numerical evaluation"
            );
        }
    }

    static void requireQFunction(const QFunction& function, const string& what) {
        if (!function) {
            throw SEBException(
                "Missing numerical callback for " + what,
                "NumericalSubunit numerical evaluation"
            );
        }
    }

    double nonzeroBeta(const ParameterList& values, const string& what) {
        const double beta = NumericTotalBeta(values);
        if (std::abs(beta) <= 1e-14) {
            throw SEBException(
                "Total excess scattering length is zero for " + what,
                "NumericalSubunit numerical evaluation"
            );
        }
        return beta;
    }

    static void validateNear(
        double actual,
        double expected,
        double tolerance,
        const string& what)
    {
        const double scale = std::max(1.0, std::abs(expected));
        if (std::abs(actual - expected) > tolerance * scale) {
            throw SEBException(
                "Numerical validation failed for " + what +
                    ": expected " + std::to_string(expected) +
                    ", got " + std::to_string(actual),
                "NumericalSubunit::ValidateNumerically()"
            );
        }
    }
};

#endif
