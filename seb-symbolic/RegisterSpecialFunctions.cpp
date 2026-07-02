#ifdef USE_GINAC

#include <ginac/ginac.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <cmath>

using namespace GiNaC;

// Serial numbers for our special functions
DECLARE_FUNCTION_1P(BesselJ0)
DECLARE_FUNCTION_1P(BesselJ1)
DECLARE_FUNCTION_1P(BesselJ2)
DECLARE_FUNCTION_1P(DawsonF)
DECLARE_FUNCTION_1P(Si)
DECLARE_FUNCTION_1P(Six)
DECLARE_FUNCTION_1P(Erf)
DECLARE_FUNCTION_1P(Erfc)
DECLARE_FUNCTION_2P(Hypergeometric0F1Regularized)
DECLARE_FUNCTION_1P(StruveH0)
DECLARE_FUNCTION_1P(StruveH1)

namespace {

void STVH0(double X, double *SH0) {
    double A0, BY0, P0, Q0, R, S, T, T2, TA0;
    int K, KM;
    constexpr double PI = 3.14159265358979323846264338327950288419716939;

    S = 1.0;
    R = 1.0;
    if (X <= 20.0) {
        A0 = 2.0 * X / PI;
        for (K = 1; K < 61; K++) {
            R = -R * X / (2.0 * K + 1.0) * X / (2.0 * K + 1.0);
            S = S + R;
            if (std::fabs(R) < std::fabs(S) * 1.0e-12) break;
        }
        *SH0 = A0 * S;
    } else {
        KM = int(0.5 * (X + 1.0));
        if (X >= 50.0) KM = 25;
        for (K = 1; K <= KM; K++) {
            R = -R * std::pow((2.0 * K - 1.0) / X, 2.0);
            S = S + R;
            if (std::fabs(R) < std::fabs(S) * 1.0e-12) break;
        }
        T = 4.0 / X;
        T2 = T * T;
        P0 = ((((-.37043e-5 * T2 + .173565e-4) * T2 - .487613e-4) * T2 + .17343e-3) * T2 - 0.1753062e-2) * T2 + .3989422793;
        Q0 = T * (((((.32312e-5 * T2 - 0.142078e-4) * T2 + 0.342468e-4) * T2 - 0.869791e-4) * T2 + 0.4564324e-3) * T2 - 0.0124669441);
        TA0 = X - 0.25 * PI;
        BY0 = 2.0 / std::sqrt(X) * (P0 * std::sin(TA0) + Q0 * std::cos(TA0));
        *SH0 = 2.0 / (PI * X) * S + BY0;
    }
}

void STVH1(double X, double *SH1) {
    double A0, BY1, P1, Q1, R, S, T, T2, TA1;
    int K, KM;
    constexpr double PI = 3.14159265358979323846264338327950288419716939;

    R = 1.0;
    if (X <= 20.0) {
        S = 0.0;
        A0 = -2.0 / PI;
        for (K = 1; K <= 60; K++) {
            R = -R * X * X / (4.0 * K * K - 1.0);
            S = S + R;
            if (std::fabs(R) < std::fabs(S) * 1.0e-12) break;
        }
        *SH1 = A0 * S;
    } else {
        S = 1.0;
        KM = int(0.5 * X);
        if (X > 50.0) KM = 25;
        for (K = 1; K <= KM; K++) {
            R = -R * (4.0 * K * K - 1.0) / (X * X);
            S = S + R;
            if (std::fabs(R) < std::fabs(S) * 1.0e-12) break;
        }
        T = 4.0 / X;
        T2 = T * T;
        P1 = ((((0.42414e-5 * T2 - 0.20092e-4) * T2 + 0.580759e-4) * T2 - 0.223203e-3) * T2 + 0.29218256e-2) * T2 + 0.3989422819;
        Q1 = T * (((((-0.36594e-5 * T2 + 0.1622e-4) * T2 - 0.398708e-4) * T2 + 0.1064741e-3) * T2 - 0.63904e-3) * T2 + 0.0374008364);
        TA1 = X - 0.75 * PI;
        BY1 = 2.0 / std::sqrt(X) * (P1 * std::sin(TA1) + Q1 * std::cos(TA1));
        *SH1 = 2.0 / PI * (1.0 + S / (X * X)) + BY1;
    }
}

} // namespace

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

static ex BesselJ2_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        return numeric(gsl_sf_bessel_Jn(2, ex_to<numeric>(x).to_double()));
    }
    return BesselJ2(x).hold();
}

static ex DawsonF_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        return numeric(gsl_sf_dawson(ex_to<numeric>(x).to_double()));
    }
    return DawsonF(x).hold();
}

static ex Si_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        return numeric(gsl_sf_Si(ex_to<numeric>(x).to_double()));
    }
    return Si(x).hold();
}

static ex Six_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        const double value = ex_to<numeric>(x).to_double();
        if (std::fabs(value) < 1e-8) return numeric(1.0);
        return numeric(gsl_sf_Si(value) / value);
    }
    return Six(x).hold();
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

static ex Hypergeometric0F1Regularized_eval(const ex& a, const ex& x) {
    if (is_a<numeric>(a) && is_a<numeric>(x)) {
        const double a_value = ex_to<numeric>(a).to_double();
        const double x_value = ex_to<numeric>(x).to_double();
        return numeric(gsl_sf_hyperg_0F1(a_value, x_value) / gsl_sf_gamma(a_value));
    }
    return Hypergeometric0F1Regularized(a, x).hold();
}

static ex StruveH0_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        double y = 0.0;
        STVH0(ex_to<numeric>(x).to_double(), &y);
        return numeric(y);
    }
    return StruveH0(x).hold();
}

static ex StruveH1_eval(const ex& x) {
    if (is_a<numeric>(x)) {
        double y = 0.0;
        STVH1(ex_to<numeric>(x).to_double(), &y);
        return numeric(y);
    }
    return StruveH1(x).hold();
}

// Register the functions with GiNaC
REGISTER_FUNCTION(BesselJ0, eval_func(BesselJ0_eval).
                 evalf_func(BesselJ0_eval).
                 latex_name("\\operatorname{BesselJ0}"));

REGISTER_FUNCTION(BesselJ1, eval_func(BesselJ1_eval).
                 evalf_func(BesselJ1_eval).
                 latex_name("\\operatorname{BesselJ1}"));

REGISTER_FUNCTION(BesselJ2, eval_func(BesselJ2_eval).
                 evalf_func(BesselJ2_eval).
                 latex_name("\\operatorname{BesselJ2}"));

REGISTER_FUNCTION(DawsonF, eval_func(DawsonF_eval).
                 evalf_func(DawsonF_eval).
                 latex_name("\\operatorname{DawsonF}"));

REGISTER_FUNCTION(Si, eval_func(Si_eval).
                 evalf_func(Si_eval).
                 latex_name("\\operatorname{Si}"));

REGISTER_FUNCTION(Six, eval_func(Six_eval).
                 evalf_func(Six_eval).
                 latex_name("\\operatorname{Six}"));

REGISTER_FUNCTION(Erf, eval_func(Erf_eval).
                 evalf_func(Erf_eval).
                 latex_name("\\operatorname{Erf}"));

REGISTER_FUNCTION(Erfc, eval_func(Erfc_eval).
                 evalf_func(Erfc_eval).
                 latex_name("\\operatorname{Erfc}"));

REGISTER_FUNCTION(Hypergeometric0F1Regularized, eval_func(Hypergeometric0F1Regularized_eval).
                 evalf_func(Hypergeometric0F1Regularized_eval).
                 latex_name("\\operatorname{Hypergeometric0F1Regularized}"));

REGISTER_FUNCTION(StruveH0, eval_func(StruveH0_eval).
                 evalf_func(StruveH0_eval).
                 latex_name("\\operatorname{StruveH0}"));

REGISTER_FUNCTION(StruveH1, eval_func(StruveH1_eval).
                 evalf_func(StruveH1_eval).
                 latex_name("\\operatorname{StruveH1}"));

#endif // USE_GINAC
