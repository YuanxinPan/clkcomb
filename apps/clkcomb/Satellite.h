#ifndef SATELLITE_H
#define SATELLITE_H

#include <string>
#include <const.h>

class Satellite
{
public:
    Satellite() {}
    explicit Satellite(const std::string &_prn) : prn(_prn) {}

    std::string prn;
    std::string svn;
    std::string blk;
    double f[2];
    double waveLen[2];
    std::string obstp[4]; // obs code in BINEX
};

#endif
