// clkcomb - Clock and phase bias products Combination
// Copyright (C) 2021 Yuanxin Pan
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef ANALYSECENTER_H
#define ANALYSECENTER_H

#include <string>
#include <vector>

#include <chrono/chrono.h>
#include <coord/coord.h>
#include <rinex/rinex.h>
#include "Satellite.h"
#include "config.h"

const double None = 0.0;

class AnalyseCenter
{
public:
    explicit AnalyseCenter(const std::string &n) : name(n) {}
    AnalyseCenter(){}
    ~AnalyseCenter(){}

    bool read_att(const std::string &path)   { return rnxatt_.read(path); }
    bool open_atx(const std::string &path)   { return rnxatx_.open(path); }
    bool open_clock(const std::string &path) { return rnxclk_.read(path); }
    bool read_orbit(const std::string &path) { return rnxsp3_.read(path); }
    bool read_orbit(const std::vector<std::string> &paths) { return rnxsp3_.read(paths); }
    bool read_sinex(const std::string &path) { return rnxsnx_.read(path); }

    bool read_bias(const std::string &path, const std::vector<std::string> &prns,
                   const std::vector<Satellite> &sats);

    bool read_clock(const config_t &config, const RinexSp3 &refsp3,
                    const RinexAtx &refatx, const RinexAtt &refatt);

    bool read_staclk(MJD t, int length, int interval,
                     const std::vector<std::string> &sta_list, const RinexSnx &refsnx);

    bool att_corr(const config_t &config, const RinexSp3 &refsp3, const RinexAtt &refatt);

    const RinexAtt &rnxatt()const { return rnxatt_; }
    const RinexAtx &rnxatx()const { return rnxatx_; }
    const RinexSp3 &rnxsp3()const { return rnxsp3_; }
    const RinexSnx &rnxsnx()const { return rnxsnx_; }

    const std::vector<std::string> &sta_names() { return rnxclk_.sta_names(); }

#if 0
private:
    bool read_biassnx(const std::string &path, const std::vector<std::string> &prns,
                      const std::vector<Satellite> &sats);
    bool read_grg_bias(const std::string &path, const std::vector<std::string> &prns);
    bool read_snx_bias(const std::string &path, const std::vector<std::string> &prns,
                       const std::vector<Satellite> &sats);
    bool read_sgg_bias(const std::string &path, const std::vector<std::string> &prns,
                       const std::vector<Satellite> &sats);
    bool read_cnes_bias(const std::string &path,
                        const std::vector<std::string> &prns,
                        const std::vector<Satellite> &sats);
#endif

private:
    RinexAtt rnxatt_;
    RinexAtx rnxatx_;
    RinexClk rnxclk_;
    RinexSnx rnxsnx_;
    RinexSp3 rnxsp3_;

public:
    std::string name;
    double weight;
    std::string att_file;
    std::string atx_file;
    std::string clk_file;
    std::string bia_file;
    std::string snx_file;
    std::string sp3_file;

    std::vector<int>    have_bias;
    std::vector<double> wl_bias;
    std::vector<double> nl_bias;
    std::vector<std::vector<double>> sat_clks; // clk of all satellites for all epoches
    std::vector<std::vector<double>> sat_stds; // std of all satellites for all clocks
    std::vector<std::vector<double>> sta_clks; // clk of all stations for all epoches
    std::vector<std::vector<double>> sta_stds; // std of all stations for all clocks
    std::vector<std::vector<double>> init_sat_clks; // satellite clk before removing bias
    std::vector<std::vector<double>> init_sta_clks; // station clk before removing bias
};

#endif //ANALYSECENTER_H
