/*
   ENUM for sub-unit types
*/

#ifndef INCLUDE_GUARD_MASTERHEADERCONSTANTS
#define INCLUDE_GUARD_MASTERHEADERCONSTANTS

// ENUM identifying all the sub-unit types
enum subunittypes{
    UNKNOWN,
    SYMBOLICSUBUNIT,
    POINT,               // Reserved legacy slot; Point is now FILEDEFINEDSUBUNIT.
    GAUSSIANPOLYMER,     // Reserved legacy slot; GaussianPolymer is now FILEDEFINEDSUBUNIT.
    GAUSSIANLOOP,
    THINCIRCLE,
    THINROD,
    THINSPHERICALSHELL,
    THINDISK,
    SOLIDSPHERE,
    SOLIDSPHERICALSHELL,
    SOLIDCYLINDER,
    NUMERICALSUBUNIT,
    DEBYESPHERECLOUD,
    FILEDEFINEDSUBUNIT
};

#endif
