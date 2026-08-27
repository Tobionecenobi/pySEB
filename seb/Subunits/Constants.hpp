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
    GAUSSIANLOOP,        // Reserved legacy slot; GaussianLoop is now FILEDEFINEDSUBUNIT.
    THINCIRCLE,          // Reserved legacy slot; ThinCircle is now FILEDEFINEDSUBUNIT.
    THINROD,             // Reserved legacy slot; ThinRod is now FILEDEFINEDSUBUNIT.
    THINSPHERICALSHELL,  // Reserved legacy slot; ThinSphericalShell is now FILEDEFINEDSUBUNIT.
    THINDISK,            // Reserved legacy slot; ThinDisk is now FILEDEFINEDSUBUNIT.
    SOLIDSPHERE,         // Reserved legacy slot; SolidSphere is now FILEDEFINEDSUBUNIT.
    SOLIDSPHERICALSHELL, // Reserved legacy slot; SolidSphericalShell is now FILEDEFINEDSUBUNIT.
    SOLIDCYLINDER,       // Reserved legacy slot; SolidCylinder is now FILEDEFINEDSUBUNIT.
    NUMERICALSUBUNIT,
    DEBYESPHERECLOUD,
    FILEDEFINEDSUBUNIT
};

#endif
