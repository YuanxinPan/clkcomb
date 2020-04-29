#ifndef COORDCONVERTER_H
#define COORDCONVERTER_H

#include "CartCoor.h"

#include <string>
#include <vector>

class CoordConverter
{
public:
    struct erp_t {
        double t;       // MJD
        double dx, dy;  // polar motion
        double ut1tai;  // UT1 - TAI
    };

    enum Type {
        ITRS2ICRS,
        ICRS2ITRS,
        TOPO2ECEF,
        ECEF2TOPO
    };

private:
    CoordConverter(const CoordConverter &);
    CoordConverter &operator=(const CoordConverter &);

public:
    CoordConverter(){}
    ~CoordConverter(){}

    //erp_t *erp(double t); // MJD
    static bool read_igserp(const std::string &path);

    // update ERP parameters to MJD::t
    static bool update(double t, double *xpole, double *ypole);

    static CartCoor convert(CartCoor pos, enum Type tp);
    static void *mat_enu2ecef(const double *xyz);
    static CartCoor convert(void *pmat, CartCoor pos, enum Type tp);

private:
    static std::vector<erp_t> erps_;
};

#endif
