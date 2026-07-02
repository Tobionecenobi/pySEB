#include "SpecialFunctionsWrapper.hpp"
#include "SymbolicPortable.hpp"

#ifdef USE_GINAC_IMPL
#include "GiNaCSymbolic.hpp"
#endif

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_erf.h>
#include <stdexcept>
#include <vector>

// Implementation of the special functions
namespace {

Expression portable_function(const std::string& name, const std::vector<SymExprPtr>& args) {
    if (SymbolicFactory::activeBackendName() != "portable") {
        throw std::runtime_error("symbolic backend '" + SymbolicFactory::activeBackendName() +
                                 "' does not support " + name);
    }

    return Expression(std::make_shared<sebsym::PortableExpression>(
        sebsym::PortableExpression::Kind::Function,
        name,
        args));
}

#ifdef USE_GINAC_IMPL
std::shared_ptr<GiNaCExpression> as_ginac_expression(const Expression& expression, const std::string& function_name) {
    auto ginac_expression = std::dynamic_pointer_cast<GiNaCExpression>(expression.get());
    if (!ginac_expression) {
        throw std::runtime_error(function_name + " requires a GiNaC expression");
    }
    return ginac_expression;
}

Expression ginac_expression(const GiNaC::ex& expression) {
    return Expression(std::make_shared<GiNaCExpression>(expression));
}
#endif

} // namespace

Expression csc(const Expression& x) {
    return Expression(SymbolicExpression::constant(1)) / x.sin();
}

Expression sec(const Expression& x) {
    return Expression(SymbolicExpression::constant(1)) / x.cos();
}

Expression BesselJ0_wrapper(const Expression& x) {
    if (SymbolicFactory::activeBackendName() != "portable" && x.is_numeric()) {
        return Expression(SymbolicExpression::constant(gsl_sf_bessel_J0(x.eval())));
    }
    return Expression(x.bessel_j0());
}

// Ensure BesselJ0 is available regardless of symbolic backend
Expression BesselJ0(const Expression& x) {
    return BesselJ0_wrapper(x);
}

Expression BesselJ1_wrapper(const Expression& x) {
    if (SymbolicFactory::activeBackendName() != "portable" && x.is_numeric()) {
        return Expression(SymbolicExpression::constant(gsl_sf_bessel_J1(x.eval())));
    }
    return Expression(x.bessel_j1());
}

// Ensure BesselJ1 is available regardless of symbolic backend
Expression BesselJ1(const Expression& x) {
    return BesselJ1_wrapper(x);
}

Expression BesselJ2(const Expression& x) {
    return (Expression(SymbolicExpression::constant(2.0)) * x.bessel_j1() - x.bessel_j0()) / x;
}

Expression DawsonF_wrapper(const Expression& x) {
    if (SymbolicFactory::activeBackendName() != "portable" && x.is_numeric()) {
        return Expression(SymbolicExpression::constant(gsl_sf_dawson(x.eval())));
    }
    return Expression(x.dawson());
}

// Ensure DawsonF is available regardless of symbolic backend
Expression DawsonF(const Expression& x) {
    return DawsonF_wrapper(x);
}

Expression Six(const Expression& x) {
#ifdef USE_GINAC_IMPL
    if (SymbolicFactory::activeBackendName() == "ginac") {
        auto g_x = as_ginac_expression(x, "Six");
        return ginac_expression(::Six(g_x->get_ginac_expr()));
    }
#endif

    Expression si = portable_function("Si", std::vector<SymExprPtr>{x.get()});
    return si / x;
}

Expression Erf_wrapper(const Expression& x) {
    if (SymbolicFactory::activeBackendName() != "portable" && x.is_numeric()) {
        return Expression(SymbolicExpression::constant(gsl_sf_erf(x.eval())));
    }
    return Expression(x.erf());
}

// Ensure Erf is available regardless of symbolic backend
Expression Erf(const Expression& x) {
    return Erf_wrapper(x);
}

Expression Erfc_wrapper(const Expression& x) {
    if (SymbolicFactory::activeBackendName() != "portable" && x.is_numeric()) {
        return Expression(SymbolicExpression::constant(gsl_sf_erfc(x.eval())));
    }
    return Expression(x.erfc());
}

// Ensure Erfc is available regardless of symbolic backend
Expression Erfc(const Expression& x) {
    return Erfc_wrapper(x);
}

Expression Hypergeometric0F1Regularized(const Expression& a, const Expression& x) {
#ifdef USE_GINAC_IMPL
    if (SymbolicFactory::activeBackendName() == "ginac") {
        auto g_a = as_ginac_expression(a, "Hypergeometric0F1Regularized");
        auto g_x = as_ginac_expression(x, "Hypergeometric0F1Regularized");
        return ginac_expression(::Hypergeometric0F1Regularized(g_a->get_ginac_expr(), g_x->get_ginac_expr()));
    }
#endif

    return portable_function("hyper0f1_regularized", std::vector<SymExprPtr>{a.get(), x.get()});
}

Expression StruveH0(const Expression& x) {
#ifdef USE_GINAC_IMPL
    if (SymbolicFactory::activeBackendName() == "ginac") {
        auto g_x = as_ginac_expression(x, "StruveH0");
        return ginac_expression(::StruveH0(g_x->get_ginac_expr()));
    }
#endif

    return portable_function("struve_h0", std::vector<SymExprPtr>{x.get()});
}

Expression StruveH1(const Expression& x) {
#ifdef USE_GINAC_IMPL
    if (SymbolicFactory::activeBackendName() == "ginac") {
        auto g_x = as_ginac_expression(x, "StruveH1");
        return ginac_expression(::StruveH1(g_x->get_ginac_expr()));
    }
#endif

    return portable_function("struve_h1", std::vector<SymExprPtr>{x.get()});
}
