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

#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>

#include <const.h>
#include <chrono/chrono.h>

typedef struct TAGCONFIG
{
    // session
    MJD mjd;
    int interval = 0;        // sample rate: 1s, 30s, etc.
    int sta_interval = 300;  // sample rate for station clock
    int length = 0;          // time span of combination

    // constellation
    // int system = 0;       // sys |= 1<<GLS
    std::string strsys;      // GREC, etc.
    std::vector<enum GNSS_Tp> syss;
    std::vector<std::string>  prns;
    std::vector<std::string>  excluded_prns;

    std::vector<std::string> sta_list;

    // ac
    std::vector<std::string> ac_names;
    std::string orb_ac;      // combined orbit
    std::string weight_method;

    // product path
    bool combine_staclk = false;
    bool phase_clock = false;
    bool use_att = false;
    std::string product_path;
    std::string nav_pattern; // broadcast ephemeris
    std::string att_pattern; // satellite attitude
    std::string bia_pattern; // satellite bias
    std::string clk_pattern; // satellite clock
    std::string sp3_pattern; // precise orbit
    std::string snx_pattern; // SINEX
    std::string erp_pattern; // ERP file

    // table path
    std::string atx_pattern; // IGS atx
    // std::string eph_path; // JPL ephmeris
    // std::string lps_path; // leap second
    // std::string chn_path; // glonass channel

    // output
    bool align_brdc = false;
    std::string product_prefix;
    std::string cls_pattern; // summary file
    std::string dif_pattern; // clock diff
    std::string log_pattern; // log file

} config_t;

bool read_config(const std::string &path, config_t &config);

bool collect_valid_prns(config_t &config);

void print_config(FILE *fp, const config_t &config);

#endif  // CONFIG_H
