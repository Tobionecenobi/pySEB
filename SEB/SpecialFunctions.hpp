/*
   ginac_extension  introduces various special functions which are required to
   evaluate scattering expressions numerically. We use GNU Scientific LIbrary's
   implementation for these to extend Ginac.

   https://www.gnu.org/software/gsl/doc/html/specfunc.html
*/

//===========================================================================
// Included guards
#ifndef INCLUDE_GINAC_EXTENSION_SPECIALFUNCTIONS
#define INCLUDE_GINAC_EXTENSION_SPECIALFUNCTIONS

#include "Expression.hpp"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_hyperg.h>

Expression csc(const Expression& x);
Expression sec(const Expression& x);
Expression power(const Expression& x, const Expression& a);
Expression BesselJ0(const Expression& x);
Expression BesselJ1(const Expression& x);
Expression BesselJ2(const Expression& x);
Expression DawsonF(const Expression& x);
Expression Six(const Expression& x);
Expression Erf(const Expression& x);
Expression Erfc(const Expression& x);
Expression Hypergeometric0F1Regularized(const Expression& a, const Expression& x);
Expression StruveH0(const Expression& x);
Expression StruveH1(const Expression& x);

#ifdef USE_GINAC
#include <ginac/ginac.h>

// GiNaC function declarations
DECLARE_FUNCTION_1P(csc)
DECLARE_FUNCTION_1P(sec)
DECLARE_FUNCTION_2P(power)
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
#endif

#endif // INCLUDE_GINAC_EXTENSION_SPECIALFUNCTIONS
