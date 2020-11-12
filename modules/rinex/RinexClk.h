#ifndef RINEXCLK_H
#define RINEXCLK_H

#include "../chrono/chrono.h"
// #include "../coord/coord.h"

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

    // const std::vector<CartCoor> &sta_poss()const { return sta_poss_; }
    const std::vector<std::string> &sta_names()const { return sta_names_; }

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

    // std::vector<CartCoor> sta_poss_; // postions for all stations
    std::vector<std::string> sta_names_;// list of staions, upper case
    // std::vector<std::vector<double>> sta_clks; // clk of all stations for all epoches
};

#endif
