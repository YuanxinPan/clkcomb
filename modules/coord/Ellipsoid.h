#ifndef ELLIPSOID_H
#define ELLIPSOID_H

class Ellipsoid
{
public:
    Ellipsoid() {}
    Ellipsoid(double _a, double _b): a(_a), b(_b)
    {
        e2 = ((a*a-b*b)/(a*a));
    }

    double a, b, e2;
};

static const Ellipsoid WGS84(6378137.0, 6356752.314245);

#endif

