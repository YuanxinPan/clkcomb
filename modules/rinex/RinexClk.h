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

    // const std::vector<CartCoor> &sta_poss()const { return sta_poss_; }
    const std::vector<std::string> &sta_names()const { return sta_names_; }

private:
    bool update();

private:
    FILE *clkFile_;
    char buf_[256];
    MJD time_[2];
    double version_;
    double interval_;
    std::vector<std::string> prns_;
    //std::array<Coef_t, NMaxSat> coefs_[2];
    std::vector<Coef_t> coefs_[2];

    // std::vector<CartCoor> sta_poss_; // postions for all stations
    std::vector<std::string> sta_names_;// list of staions, upper case
    // std::vector<std::vector<double>> sta_clks; // clk of all stations for all epoches
};

#endif
