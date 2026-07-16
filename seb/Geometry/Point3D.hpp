#ifndef INCLUDE_GUARD_POINT3D
#define INCLUDE_GUARD_POINT3D

struct CartesianPoint3D {
    double x;
    double y;
    double z;

    CartesianPoint3D(double xValue = 0.0,
                     double yValue = 0.0,
                     double zValue = 0.0)
        : x(xValue), y(yValue), z(zValue) {}
};

#endif
