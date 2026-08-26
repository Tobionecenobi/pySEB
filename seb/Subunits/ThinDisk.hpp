//===========================================================================
// Included guards
#ifndef INCLUDE_GUARD_THINDISK
#define INCLUDE_GUARD_THINDISK

#include <cmath>
#include <gsl/gsl_sf_bessel.h>

#include "ExpressionFunctions.hpp"
#include "Integrated.hpp"

/*===========================================================================

   This file implements the scattering terms of a infinitesimal thin disk.

   The Shell has two reference points:
          center     The geometric point at the center of the disk.
          surface    A randomly selected point on the surface of the disk.

   We choose R the radius as the characteristic length scales.
   Hence the dimensionless variable is x=qR.

   See SolidSphere.nb / SolidSphere.pdf for derivations.
============================================================================= */

class ThinDisk : public IntegratedSubunit {
    public:
    /*Constructor*/
    public:
        ThinDisk() : IntegratedSubunit(){
            type = SUBUNITCHILD;
           stype = THINDISK;
        }
    virtual ~ThinDisk(){};

    virtual void Init(string name, string tag, SymbolInterface *GEX)
    {
        // Initialize base class
        IntegratedSubunit::Init(name, tag, GEX);
        string n = getTag();

        // ========================================================================================
        // Setup reference points

        // specific reference points for a ThinDisk
        setReferencePointName("center");

        // distributed reference points for a ThinDisk
        setDistReferencePointType("surface");
        setDistReferencePointType("rim");

        // Define symbols
        Expression q    = GLEX->getSymbol("q");
        Expression x    = GLEX->getSymbol("x", n);
        Expression R    = GLEX->getSymbol("R", n);
        Expression t    = GLEX->getSymbol("t", n);

        // ========================================================================================
        // Define mapping between dimensionless variables and scattering expressions using R and q as variables:

        expand[x] = q*R;            // x=qR

        xparameters.insert(x);      // Structural parameter.
        parameters.insert(R);       // Structural parameter.

        // ========================================================================================
        // Scattering expressions

        Expression Ac = integrate(t, 0, pi()/2, 2*BesselJ1(x*sin(t))/x );                                       // Amplitude relative to center of a thin spherical shell

        FormFactorExpression = 2*(x-BesselJ1(2*x))/pow(x,3);

        // Form factor amplitude expression relative to all reference points.

        FormFactorAmplitudeExpressions["center"]   = Ac;
        FormFactorAmplitudeExpressions["surface"]  = FormFactorExpression;
        FormFactorAmplitudeExpressions["rim"]      = integrate(t, 0, Pi()/2, 2*BesselJ0(x*sin(t))*BesselJ1(x*sin(t))/x );

        // ========================================================================================
        // Phase factors
        // We need phase factors for all reference point pairs except a specific reference point and itself, since PhaseFactor[X][X]=1.

        // between the center and surface
        PhaseFactorExpressions["center"]["surface"]    = Ac;
        PhaseFactorExpressions["surface"]["surface"]   = FormFactorExpression;
        PhaseFactorExpressions["center"]["rim"]        = sin(x)/x;
        PhaseFactorExpressions["rim"]["surface"]       = integrate(t, 0, Pi()/2, 2*BesselJ0(x*sin(t))*BesselJ1(x*sin(t))/x );
        PhaseFactorExpressions["rim"]["rim"]           = BesselJ0(2*x)+Pi()/2*(BesselJ1(2*x)*StruveH0(2*x)-BesselJ0(2*x)*StruveH1(2*x));

        // ========================================================================================
        // Sizes

        // Radius of gyration expression for polymer, simple since we express the polymer as function of R:
        RadiusOfGyration2 = R*R/2;

        // sigma <R^2> for all distances between reference points and scatterers.
        // These are exactly the Guinier expansions of the corresponding amplitudes
        sigmaMSDref2scat["center"]    = R*R/2;
        sigmaMSDref2scat["surface"]   = R*R;
        sigmaMSDref2scat["rim"]       = 3*R*R/2;

        // sigma <R^2> for all distances between pairs of reference points.
        // These are exactly the Guinier expansions of the corresponding phasefactors

        // between the center and surface
        sigmaMSDref2ref["center"]["surface"]   = R*R/2;
        sigmaMSDref2ref["surface"]["surface"]  = R*R;
        sigmaMSDref2ref["center"]["rim"]       = R*R;
        sigmaMSDref2ref["rim"]["surface"]      = 3*R*R/2;;
        sigmaMSDref2ref["rim"]["rim"]          = 2*R*R;

        IntegrationOptions integrationOptions;
        integrationOptions.method = IntegrationMethod::QAG;
        integrationOptions.absoluteTolerance = 1e-12;
        integrationOptions.relativeTolerance = 1e-8;
        integrationOptions.workspaceSize = 1000;
        setIntegrationOptions(integrationOptions);

        setNumericalFormFactorAmplitudeFunction(
            "center",
            [this](double qValue, const ParameterList& values) {
                return numericalCenterAmplitude(qValue, values);
            }
        );
        setNumericalFormFactorAmplitudeFunction(
            "rim",
            [this](double qValue, const ParameterList& values) {
                return numericalRimAmplitude(qValue, values);
            }
        );
        setNumericalPhaseFactorFunction(
            "center",
            "surface",
            [this](double qValue, const ParameterList& values) {
                return numericalCenterAmplitude(qValue, values);
            }
        );
        setNumericalPhaseFactorFunction(
            "rim",
            "surface",
            [this](double qValue, const ParameterList& values) {
                return numericalRimAmplitude(qValue, values);
            }
        );
    }

private:
    static double jinc(double value)
    {
        const double value2 = value * value;
        if (value2 < 1e-8) {
            return 1.0 - value2 / 8.0 + value2 * value2 / 192.0;
        }
        return 2.0 * gsl_sf_bessel_J1(value) / value;
    }

    double numericalCenterAmplitude(
        double qValue,
        const ParameterList& values)
    {
        const double radius = requirePositive(
            numericParameter("R", values),
            "R_" + getTag()
        );
        const double halfPi = std::acos(-1.0) / 2.0;
        return integrateNumerically(
            [qValue, radius](double theta) {
                const double sinTheta = std::sin(theta);
                return sinTheta * jinc(qValue * radius * sinTheta);
            },
            0.0,
            halfPi
        ).value;
    }

    double numericalRimAmplitude(
        double qValue,
        const ParameterList& values)
    {
        const double radius = requirePositive(
            numericParameter("R", values),
            "R_" + getTag()
        );
        const double halfPi = std::acos(-1.0) / 2.0;
        return integrateNumerically(
            [qValue, radius](double theta) {
                const double sinTheta = std::sin(theta);
                const double z = qValue * radius * sinTheta;
                return sinTheta * gsl_sf_bessel_J0(z) * jinc(z);
            },
            0.0,
            halfPi
        ).value;
    }

};

#endif // INCLUDE_GUARD_THINDISK
