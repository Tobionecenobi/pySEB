#ifndef INCLUDE_GUARD_INTEGRATEDSUBUNIT
#define INCLUDE_GUARD_INTEGRATEDSUBUNIT

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <utility>

#include "NumericalIntegration.hpp"
#include "Subunit.hpp"

class IntegratedSubunit : public SubUnit {
public:
    using QFunction = std::function<double(double, const ParameterList&)>;
    using ScalarProvider = std::function<double(const ParameterList&)>;

    IntegratedSubunit() : SubUnit(), integrator_() {}
    virtual ~IntegratedSubunit() {}

    double NumericFormFactorUnnormalized(
        double q,
        const ParameterList& values) override
    {
        if (!formFactorFunction_) {
            return SubUnit::NumericFormFactorUnnormalized(q, values);
        }
        q = finiteQ(q);
        const double beta = NumericTotalBeta(values);
        if (q == 0.0) {
            return beta * beta;
        }
        return beta * beta * finiteValue(
            formFactorFunction_(q, values),
            "form factor"
        );
    }

    double NumericFormFactorAmplitudeUnnormalized(
        refPoint ref,
        double q,
        const ParameterList& values) override
    {
        const refPoint canonical = canonicalReference(ref);
        const auto function = amplitudeFunctions_.find(canonical);
        if (function == amplitudeFunctions_.end() || !function->second) {
            return SubUnit::NumericFormFactorAmplitudeUnnormalized(
                ref,
                q,
                values
            );
        }
        requireKnownReference(ref);
        q = finiteQ(q);
        const double beta = NumericTotalBeta(values);
        if (q == 0.0) {
            return beta;
        }
        return beta * finiteValue(
            function->second(q, values),
            "form-factor amplitude"
        );
    }

    double NumericPhaseFactor(
        refPoint r1,
        refPoint r2,
        double q,
        const ParameterList& values) override
    {
        requireKnownReference(r1);
        requireKnownReference(r2);
        const refPoint c1 = canonicalReference(r1);
        const refPoint c2 = canonicalReference(r2);
        if (r1 == r2 && hasSpecificReference(r1)) {
            return 1.0;
        }

        const auto pair = referencePair(c1, c2);
        const auto function = phaseFunctions_.find(pair);
        if (function == phaseFunctions_.end() || !function->second) {
            return SubUnit::NumericPhaseFactor(r1, r2, q, values);
        }
        q = finiteQ(q);
        if (q == 0.0) {
            return 1.0;
        }
        return finiteValue(
            function->second(q, values),
            "phase factor"
        );
    }

    double NumericRadiusOfGyration2(
        const ParameterList& values) override
    {
        if (!radiusOfGyration2Provider_) {
            return SubUnit::NumericRadiusOfGyration2(values);
        }
        return finiteValue(
            radiusOfGyration2Provider_(values),
            "radius of gyration squared"
        );
    }

    double NumericSigmaMSDRef2Scat(
        refPoint ref,
        const ParameterList& values) override
    {
        const refPoint canonical = canonicalReference(ref);
        const auto provider = sigmaRefToScattererProviders_.find(canonical);
        if (provider == sigmaRefToScattererProviders_.end() ||
            !provider->second) {
            return SubUnit::NumericSigmaMSDRef2Scat(ref, values);
        }
        requireKnownReference(ref);
        return finiteValue(
            provider->second(values),
            "reference-to-scatterer sigma<R^2>"
        );
    }

    double NumericSigmaMSDRef2Ref(
        refPoint r1,
        refPoint r2,
        const ParameterList& values) override
    {
        const refPoint c1 = canonicalReference(r1);
        const refPoint c2 = canonicalReference(r2);
        if (r1 == r2 && hasSpecificReference(r1)) {
            return 0.0;
        }
        const auto pair = referencePair(c1, c2);
        const auto provider = sigmaRefToRefProviders_.find(pair);
        if (provider == sigmaRefToRefProviders_.end() ||
            !provider->second) {
            return SubUnit::NumericSigmaMSDRef2Ref(r1, r2, values);
        }
        requireKnownReference(r1);
        requireKnownReference(r2);
        return finiteValue(
            provider->second(values),
            "reference-to-reference sigma<R^2>"
        );
    }

protected:
    void setIntegrationOptions(const IntegrationOptions& options)
    {
        integrator_.setOptions(options);
    }

    IntegrationResult integrateNumerically(
        const NumericalIntegrator::Function& function,
        double lower,
        double upper)
    {
        return integrator_.integrate(function, lower, upper);
    }

    void setNumericalFormFactorFunction(QFunction function)
    {
        formFactorFunction_ = std::move(function);
    }

    void setNumericalFormFactorAmplitudeFunction(
        const refPoint& ref,
        QFunction function)
    {
        requireKnownReference(ref);
        amplitudeFunctions_[canonicalReference(ref)] = std::move(function);
    }

    void setNumericalPhaseFactorFunction(
        refPoint r1,
        refPoint r2,
        QFunction function)
    {
        requireKnownReference(r1);
        requireKnownReference(r2);
        const refPoint c1 = canonicalReference(r1);
        const refPoint c2 = canonicalReference(r2);
        if (r1 == r2 && hasSpecificReference(r1)) {
            throw SEBException(
                "A phase-factor callback is not needed for an identical "
                "specific reference point",
                "IntegratedSubunit::setNumericalPhaseFactorFunction()"
            );
        }
        phaseFunctions_[referencePair(c1, c2)] = std::move(function);
    }

    void setNumericalRadiusOfGyration2Provider(ScalarProvider provider)
    {
        radiusOfGyration2Provider_ = std::move(provider);
    }

    void setNumericalSigmaMSDRef2ScatProvider(
        const refPoint& ref,
        ScalarProvider provider)
    {
        requireKnownReference(ref);
        sigmaRefToScattererProviders_[canonicalReference(ref)] =
            std::move(provider);
    }

    void setNumericalSigmaMSDRef2RefProvider(
        refPoint r1,
        refPoint r2,
        ScalarProvider provider)
    {
        requireKnownReference(r1);
        requireKnownReference(r2);
        sigmaRefToRefProviders_[referencePair(r1, r2)] =
            std::move(provider);
    }

    double numericParameter(
        const std::string& baseName,
        const ParameterList& values) const
    {
        const std::string key = baseName + "_" + tag;
        const auto value = values.find(key);
        if (value == values.end()) {
            throw SEBException(
                "Missing numerical parameter " + key,
                "IntegratedSubunit::numericParameter()"
            );
        }
        return finiteValue(value->second, "parameter " + key);
    }

    static double requirePositive(
        double value,
        const std::string& parameterName)
    {
        if (value <= 0.0) {
            throw SEBException(
                "Numerical parameter " + parameterName +
                    " must be positive",
                "IntegratedSubunit numerical evaluation"
            );
        }
        return value;
    }

private:
    NumericalIntegrator integrator_;
    QFunction formFactorFunction_;
    std::map<refPoint, QFunction> amplitudeFunctions_;
    std::map<std::pair<refPoint, refPoint>, QFunction> phaseFunctions_;
    ScalarProvider radiusOfGyration2Provider_;
    std::map<refPoint, ScalarProvider> sigmaRefToScattererProviders_;
    std::map<std::pair<refPoint, refPoint>, ScalarProvider>
        sigmaRefToRefProviders_;

    refPoint canonicalReference(refPoint ref) const
    {
        const auto hash = ref.find("#");
        if (hash != std::string::npos) {
            ref = ref.substr(0, hash);
        }
        return ref;
    }

    std::pair<refPoint, refPoint> referencePair(
        refPoint r1,
        refPoint r2) const
    {
        r1 = canonicalReference(r1);
        r2 = canonicalReference(r2);
        if (r2 < r1) {
            std::swap(r1, r2);
        }
        return std::make_pair(r1, r2);
    }

    void requireKnownReference(const refPoint& ref)
    {
        const refPoint canonical = canonicalReference(ref);
        if (!hasSpecificReference(ref) &&
            !hasSpecificReference(canonical) &&
            !hasDistributedReference(canonical)) {
            throw SEBException(
                "Unknown integrated reference point " + ref,
                "IntegratedSubunit numerical evaluation"
            );
        }
    }

    static double finiteQ(double q)
    {
        if (!std::isfinite(q)) {
            throw SEBException(
                "q must be finite",
                "IntegratedSubunit numerical evaluation"
            );
        }
        return q;
    }

    static double finiteValue(
        double value,
        const std::string& description)
    {
        if (!std::isfinite(value)) {
            throw SEBException(
                "Numerical " + description + " produced a non-finite value",
                "IntegratedSubunit numerical evaluation"
            );
        }
        return value;
    }
};

#endif
