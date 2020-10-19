#ifndef ANALYSECENTER_H
#define ANALYSECENTER_H

#include <string>
#include <vector>

#include <pppx/chrono.h>
#include <pppx/coord.h>
#include <pppx/rinex.h>
#include "Satellite.h"

const double None = 0.0;

class AnalyseCenter
{
public:
    explicit AnalyseCenter(const std::string &n) : name(n) {}
    AnalyseCenter(){}
    ~AnalyseCenter(){}

    bool read_orbit(const std::string &path);
    bool read_orbit(const std::vector<std::string> &paths);
    bool open_clock(const std::string &path);
    bool read_bias(const std::string &path, const std::vector<std::string> &prns,
                   const std::vector<Satellite> &sats);

    bool read_clock(MJD t, int length, int interval,
                    const std::vector<std::string> &prns, const RinexSp3 &refsp3);

    const RinexSp3 &rnxsp3()const { return rnxsp3_; }
    // void close() { rnxclk_.close(); }

private:
    bool read_grg_bias(const std::string &path, const std::vector<std::string> &prns);
    bool read_snx_bias(const std::string &path, const std::vector<std::string> &prns,
                       const std::vector<Satellite> &sats);

private:
    RinexSp3 rnxsp3_;
    RinexClk rnxclk_;

public:
    std::string name;
    std::string sp3_file;
    std::string clk_file;
    std::string bia_file;

    std::vector<std::vector<double>> sat_clks; //  clk of all satellites for all epoches
    std::vector<int>    have_bias;
    std::vector<double> wl_bias;
    std::vector<double> nl_bias;
};

#endif //ANALYSECENTER_H
