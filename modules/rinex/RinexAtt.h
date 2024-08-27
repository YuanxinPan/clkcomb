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

#ifndef RINEXATT_H
#define RINEXATT_H

#include "../chrono/chrono.h"

#include <map>
#include <string>
#include <vector>

class RinexAtt
{
private:
    struct Coef_t
    {
        MJD t;
        double q[4]; // ECEF => SV
        bool operator<(MJD mjd)const { return t<mjd; }
    };

public:
    RinexAtt(): attFile_(nullptr), interval_(0) {}
    ~RinexAtt() {}

    bool read(const std::string &path);
    bool sat_att(MJD t, const std::string &prn, double *q)const; // SV => ECEF
    bool empty()const { return atts_.empty(); }
    void close();
    double interval()const { return interval_; }
    double version()const { return version_; }

private:
    FILE *attFile_;
    double version_;
    double interval_;
    std::string creator_;
    std::map<std::string, std::string> svTypes_;
    std::vector<std::string> prns_;
    std::vector<std::vector<Coef_t>> atts_; // prn, epoch
};

// bool operator<(const RinexAtt::Coef_t &lhs, const MJD t) { return lhs.t < t; }

#endif // RINEXATT_H
