#ifndef INCLUDE_GUARD_NUMERICALINTEGRATION
#define INCLUDE_GUARD_NUMERICALINTEGRATION

#include <cstddef>
#include <functional>

#include <gsl/gsl_integration.h>

enum class IntegrationMethod { QAG, CQUAD };

struct IntegrationOptions {
    IntegrationMethod method = IntegrationMethod::CQUAD;
    double absoluteTolerance = 1e-12;
    double relativeTolerance = 1e-8;
    std::size_t workspaceSize = 1000;
    int qagRule = 6; // GSL_INTEG_GAUSS61
};

struct IntegrationResult {
    double value = 0.0;
    double estimatedError = 0.0;
    std::size_t evaluations = 0;
};

class NumericalIntegrator {
public:
    using Function = std::function<double(double)>;

    explicit NumericalIntegrator(
        const IntegrationOptions& options = IntegrationOptions()
    );
    ~NumericalIntegrator();

    NumericalIntegrator(const NumericalIntegrator&) = delete;
    NumericalIntegrator& operator=(const NumericalIntegrator&) = delete;

    void setOptions(const IntegrationOptions& options);
    const IntegrationOptions& getOptions() const { return options_; }
    IntegrationResult integrate(
        const Function& function,
        double lower,
        double upper
    );

private:
    IntegrationOptions options_;
    gsl_integration_workspace* qagWorkspace_ = nullptr;
    gsl_integration_cquad_workspace* cquadWorkspace_ = nullptr;

    void allocateWorkspaces();
    void freeWorkspaces();
    static void validateOptions(const IntegrationOptions& options);
};

#endif
