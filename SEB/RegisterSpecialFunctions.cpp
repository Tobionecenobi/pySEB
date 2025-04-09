#include <ginac/ginac.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_hyperg.h>

using namespace GiNaC;

// Serial numbers for our special functions
DECLARE_FUNCTION_1P(BesselJ0)
DECLARE_FUNCTION_1P(BesselJ1)
DECLARE_FUNCTION_1P(DawsonF)
DECLARE_FUNCTION_1P(Erf)
DECLARE_FUNCTION_1P(Erfc)

// Implementation of the special functions
static ex BesselJ0_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        return numeric(gsl_sf_bessel_J0(ex_to<numeric>(x).to_double()));
    }
    return BesselJ0(x).hold();
}

static ex BesselJ1_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        return numeric(gsl_sf_bessel_J1(ex_to<numeric>(x).to_double()));
    }
    return BesselJ1(x).hold();
}

static ex DawsonF_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        return numeric(gsl_sf_dawson(ex_to<numeric>(x).to_double()));
    }
    return DawsonF(x).hold();
}

static ex Erf_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        return numeric(gsl_sf_erf(ex_to<numeric>(x).to_double()));
    }
    return Erf(x).hold();
}

static ex Erfc_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        return numeric(gsl_sf_erfc(ex_to<numeric>(x).to_double()));
    }
    return Erfc(x).hold();
}

// Register the functions with GiNaC
REGISTER_FUNCTION(BesselJ0, eval_func(BesselJ0_eval).
                 latex_name("\\operatorname{BesselJ0}"));

REGISTER_FUNCTION(BesselJ1, eval_func(BesselJ1_eval).
                 latex_name("\\operatorname{BesselJ1}"));

REGISTER_FUNCTION(DawsonF, eval_func(DawsonF_eval).
                 latex_name("\\operatorname{DawsonF}"));

REGISTER_FUNCTION(Erf, eval_func(Erf_eval).
                 latex_name("\\operatorname{Erf}"));

REGISTER_FUNCTION(Erfc, eval_func(Erfc_eval).
                 latex_name("\\operatorname{Erfc}"));
