#ifndef RINEXCLK_H
#define RINEXCLK_H

#include "../chrono/chrono.h"

#include <stdio.h>
#include <array>
#include <string>
#include <vector>

class RinexClk
{
private:
    struct Coef_t
    {
        Coef_t(): valid(false), a(0.0) {}
        bool valid;
        double a; //[3]; // bias, vel, acc
    };

// private:
//     RinexClk(const RinexClk &);
//     RinexClk &operator=(const RinexClk &);

public:
    RinexClk(): clkFile_(nullptr), interval_(0) {}
    ~RinexClk(){}

    bool read(const std::string &path);
    bool clkBias(MJD t, const std::string &prn, double &sck, double *drift=nullptr);
    void close();

private:
    bool update();

private:
    FILE *clkFile_;
    char buf_[256];
    MJD time_[2];
    double interval_;
    std::vector<std::string> prns_;
    //std::array<Coef_t, NMaxSat> coefs_[2];
    std::vector<Coef_t> coefs_[2];
};

#endif
