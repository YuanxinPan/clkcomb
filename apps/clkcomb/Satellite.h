#ifndef SATELLITE_H
#define SATELLITE_H

#include <string>
#include <pppx/const.h>

class Satellite
{
public:
    Satellite() {}
    explicit Satellite(const std::string &_prn) : prn(_prn) {}

    std::string prn;
    double f[2];
    double waveLen[2];
};

#endif
