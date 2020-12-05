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
    bool read_sinex(const std::string &path);
    // bool read_orbit(const std::vector<std::string> &paths);
    bool open_atx(const std::string &path);
    bool open_clock(const std::string &path);
    bool read_bias(const std::string &path, const std::vector<std::string> &prns,
                   const std::vector<Satellite> &sats);

    bool read_clock(MJD t, int length, int interval,
                    const std::vector<std::string> &prns, const RinexSp3 &refsp3, RinexAtx &refatx);

    RinexAtx &rnxatx() { return rnxatx_; }
    const RinexSp3 &rnxsp3()const { return rnxsp3_; }
    const RinexSnx &rnxsnx()const { return rnxsnx_; }

    const std::vector<std::string> &sta_names() { return rnxclk_.sta_names(); }

    bool read_staclk(MJD t, int length, int interval,
                     const std::vector<std::string> &sta_list, const RinexSnx &refsnx);

private:
    bool read_grg_bias(const std::string &path, const std::vector<std::string> &prns);
    bool read_snx_bias(const std::string &path, const std::vector<std::string> &prns,
                       const std::vector<Satellite> &sats);

private:
    RinexAtx rnxatx_;
    RinexClk rnxclk_;
    RinexSnx rnxsnx_;
    RinexSp3 rnxsp3_;

public:
    std::string name;
    std::string atx_file;
    std::string clk_file;
    std::string bia_file;
    std::string snx_file;
    std::string sp3_file;

    std::vector<std::vector<double>> sta_clks; // clk of all stations for all epoches
    std::vector<std::vector<double>> sat_clks; // clk of all satellites for all epoches
    std::vector<int>    have_bias;
    std::vector<double> wl_bias;
    std::vector<double> nl_bias;
};

#endif //ANALYSECENTER_H
