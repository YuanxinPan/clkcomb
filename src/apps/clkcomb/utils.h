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

#ifndef UTILS_H
#define UTILS_H

#include "config.h"
#include "AnalyseCenter.h"
#include "Satellite.h"

#include <vector>
#include <string>
#include <algorithm>

std::string replace_pattern(const std::string &pattern, MJD t,
                            const std::string &prefix=std::string(),
                            const std::string &suffix=std::string());

bool init_acs(config_t &config, const std::vector<Satellite> &sats,
              std::vector<AnalyseCenter> &acs, AnalyseCenter &combined_ac);

bool init_sats(const config_t &config, std::vector<Satellite> &sats);

bool construct_initclk(const config_t &config,
                       const std::vector<AnalyseCenter> &acs,
                       const std::vector<Satellite> &sats,
                       std::vector<std::vector<double>> &comb_clks,
                       int &refac);

bool construct_init_staclk(const config_t &config,
                           const std::vector<AnalyseCenter> &acs,
                           std::vector<std::vector<double>> &comb_clks);

bool align_widelane(const std::vector<Satellite> &sats,
                    std::vector<AnalyseCenter> &acs,
                    AnalyseCenter &combined_ac,
                    const std::vector<int> &prns);

void count_satclks(const config_t &config,
                   const std::vector<std::vector<double>> &sat_clks,
                   std::vector<size_t> &counts);

void remove_clock_datum(std::vector<std::vector<double>> &sat_clks,
                        const std::vector<int> &adjust_prns,
                        const std::vector<int> &common_prns);

void align_clock_datum(std::vector<std::vector<double>> &sat_clks,
                       std::vector<std::vector<double>> &ref_clks,
                       const std::vector<int> &adjust_prns,
                       const std::vector<int> &common_prns);

void align_clock_datum2(std::vector<std::vector<double>> &sat_clks,
                        std::vector<std::vector<double>> &ref_clks,
                        const std::vector<int> &adjust_prns,
                        const std::vector<int> &common_prns,
                        std::vector<std::vector<double>> &sta_clks, int sample_ratio);

double satclk_std(const std::vector<double> &clks, const std::vector<double> &refs, double T);

void remove_clock_bias(const std::string &ac_name,
                       const Satellite &sat, const int niter,
                       std::vector<double> &sat_clks,
                       const std::vector<double> &ref_clks,
                       bool phase_clock, double &std, bool edit);

void combine_one_epoch(const std::vector<double> &clks,
                       const std::vector<double> &wgts,
                       double &mean,
                       double &wrms,
                       std::vector<int> &deleted);

void write_satclks(FILE *fp, const config_t &config, const std::string &acn,
                   const std::vector<std::vector<double>> &sat_clks);

void write_clks_diff(FILE *fp, const std::vector<std::string> &name_list,
                     const std::string &prefix, const std::string &acn,
                     const std::vector<std::vector<double>> &src_clks,
                     const std::vector<std::vector<double>> &ref_clks);

bool write_fip(const std::string &path, MJD t,
               const std::vector<Satellite> &sats,
               const std::vector<double> &wl_bias,
               const std::vector<double> &nl_bias);

void write_clock_datum(const std::vector<std::string> &prns,
                       char sys, const std::string &acn,
                       const std::vector<std::vector<double>> &sat_clks,
                       const std::vector<int> &common_prns);

void compare_clks(const std::vector<std::string> &name_list,
                     const std::string &ac_name,
                     const std::vector<std::vector<double>> &src_clks,
                     const std::vector<std::vector<double>> &ref_clks,
                     std::vector<double> &biass,
                     std::vector<double> &stds,
                     bool epoch_output);

bool write_clkfile(const std::string &path, const config_t &config, const AnalyseCenter &combined);

bool write_bias(const std::string &path, MJD t,
                const std::vector<Satellite> &sats,
                const std::vector<int> &have_bias,
                const std::vector<double> &wl_bias);

bool write_att(const std::string &path, const config_t &config, const RinexAtt &rnxatt);

void clkfit(std::vector<double> clks, const int n, double &offset, double &drift);

bool write_summary(const std::string &path, config_t &config,
                   const std::vector<AnalyseCenter> &acs,
                   const std::vector<Satellite> &sats,
                   std::vector<std::vector<double>> &comb_clks,
                   std::vector<std::vector<double>> &comb_staclks);

double stable_mean(const std::vector<double> &v);

double weighted_mean(const std::vector<double> &vals, const std::vector<double> &wgts, double &wrms);

double median(const std::vector<double> &vals);

void unit_weight(std::vector<double> &wgts);

bool sat_freq(const std::string &prn, double f[]);

enum GNSS_Tp prn2sys(const std::string &prn);
char sys2char(enum GNSS_Tp tp);

#endif // UTILS_H
