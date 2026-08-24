#include "NumericalIntegration.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "Exceptions.hpp"

namespace {

struct FunctionPayload {
    const NumericalIntegrator::Function* function;
    std::exception_ptr exception;
    bool producedNonFiniteValue = false;
};

double gslFunctionAdapter(double x, void* rawPayload) noexcept
{
    FunctionPayload* payload = static_cast<FunctionPayload*>(rawPayload);
    try {
        const double value = (*payload->function)(x);
        if (!std::isfinite(value)) {
            payload->producedNonFiniteValue = true;
        }
        return value;
    } catch (...) {
        payload->exception = std::current_exception();
        return std::numeric_limits<double>::quiet_NaN();
    }
}

const char* methodName(IntegrationMethod method)
{
    switch (method) {
        case IntegrationMethod::QAG:
            return "QAG";
        case IntegrationMethod::CQUAD:
            return "CQUAD";
    }
    return "unknown";
}

void disableGslAbortHandler()
{
    static std::once_flag flag;
    std::call_once(flag, []() { gsl_set_error_handler_off(); });
}

} // namespace

NumericalIntegrator::NumericalIntegrator(const IntegrationOptions& options)
    : options_(options)
{
    disableGslAbortHandler();
    validateOptions(options_);
    allocateWorkspace();
}

NumericalIntegrator::~NumericalIntegrator()
{
    freeWorkspaces();
}

void NumericalIntegrator::validateOptions(const IntegrationOptions& options)
{
    if (!std::isfinite(options.absoluteTolerance) ||
        options.absoluteTolerance < 0.0) {
        throw SEBException(
            "Absolute integration tolerance must be finite and non-negative",
            "NumericalIntegrator::setOptions()"
        );
    }
    if (!std::isfinite(options.relativeTolerance) ||
        options.relativeTolerance < 0.0) {
        throw SEBException(
            "Relative integration tolerance must be finite and non-negative",
            "NumericalIntegrator::setOptions()"
        );
    }
    if (options.absoluteTolerance == 0.0 &&
        options.relativeTolerance == 0.0) {
        throw SEBException(
            "At least one integration tolerance must be positive",
            "NumericalIntegrator::setOptions()"
        );
    }
    if (options.workspaceSize < 3) {
        throw SEBException(
            "Integration workspace size must be at least 3",
            "NumericalIntegrator::setOptions()"
        );
    }
    if (options.qagRule < GSL_INTEG_GAUSS15 ||
        options.qagRule > GSL_INTEG_GAUSS61) {
        throw SEBException(
            "Invalid GSL QAG rule",
            "NumericalIntegrator::setOptions()"
        );
    }
}

void NumericalIntegrator::allocateWorkspace()
{
    if (options_.method == IntegrationMethod::QAG) {
        qagWorkspace_ =
            gsl_integration_workspace_alloc(options_.workspaceSize);
    } else {
        cquadWorkspace_ =
            gsl_integration_cquad_workspace_alloc(options_.workspaceSize);
    }

    if (!qagWorkspace_ && !cquadWorkspace_) {
        throw SEBException(
            "Unable to allocate GSL integration workspace",
            "NumericalIntegrator::allocateWorkspace()"
        );
    }
}

void NumericalIntegrator::freeWorkspaces()
{
    if (qagWorkspace_) {
        gsl_integration_workspace_free(qagWorkspace_);
        qagWorkspace_ = nullptr;
    }
    if (cquadWorkspace_) {
        gsl_integration_cquad_workspace_free(cquadWorkspace_);
        cquadWorkspace_ = nullptr;
    }
}

void NumericalIntegrator::setOptions(const IntegrationOptions& options)
{
    validateOptions(options);
    const bool reallocate =
        options.method != options_.method ||
        options.workspaceSize != options_.workspaceSize;
    if (reallocate) {
        freeWorkspaces();
    }
    options_ = options;
    if (reallocate) {
        allocateWorkspace();
    }
}

IntegrationResult NumericalIntegrator::integrate(
    const Function& function,
    double lower,
    double upper)
{
    if (!function) {
        throw SEBException(
            "Integration function is empty",
            "NumericalIntegrator::integrate()"
        );
    }
    if (!std::isfinite(lower) || !std::isfinite(upper)) {
        throw SEBException(
            "Finite-interval integration requires finite bounds",
            "NumericalIntegrator::integrate()"
        );
    }
    if (lower == upper) {
        return IntegrationResult();
    }

    double sign = 1.0;
    if (upper < lower) {
        std::swap(lower, upper);
        sign = -1.0;
    }

    FunctionPayload payload{&function, nullptr, false};
    gsl_function gslFunction;
    gslFunction.function = &gslFunctionAdapter;
    gslFunction.params = &payload;

    IntegrationResult result;
    int status = GSL_SUCCESS;
    if (options_.method == IntegrationMethod::QAG) {
        status = gsl_integration_qag(
            &gslFunction,
            lower,
            upper,
            options_.absoluteTolerance,
            options_.relativeTolerance,
            options_.workspaceSize,
            options_.qagRule,
            qagWorkspace_,
            &result.value,
            &result.estimatedError
        );
    } else {
        status = gsl_integration_cquad(
            &gslFunction,
            lower,
            upper,
            options_.absoluteTolerance,
            options_.relativeTolerance,
            cquadWorkspace_,
            &result.value,
            &result.estimatedError,
            &result.evaluations
        );
    }

    if (payload.exception) {
        try {
            std::rethrow_exception(payload.exception);
        } catch (const SEBException&) {
            throw;
        } catch (const std::exception& error) {
            throw SEBException(
                std::string("Integration callback failed: ") + error.what(),
                "NumericalIntegrator::integrate()"
            );
        } catch (...) {
            throw SEBException(
                "Integration callback failed with an unknown exception",
                "NumericalIntegrator::integrate()"
            );
        }
    }

    if (payload.producedNonFiniteValue) {
        throw SEBException(
            "Integration callback produced a non-finite value",
            "NumericalIntegrator::integrate()"
        );
    }

    result.value *= sign;
    if (status != GSL_SUCCESS || !std::isfinite(result.value) ||
        !std::isfinite(result.estimatedError)) {
        std::ostringstream message;
        message << methodName(options_.method)
                << " integration failed on [" << lower << ", " << upper << "]"
                << ": " << gsl_strerror(status)
                << ", estimated error=" << result.estimatedError;
        throw SEBException(
            message.str(),
            "NumericalIntegrator::integrate()"
        );
    }

    return result;
}
