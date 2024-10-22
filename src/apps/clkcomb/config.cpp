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

#include "config.h"
#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>

#include <io/io.h>
#include <const.h>

static void print_prns(FILE *fp, const std::vector<std::string> &prns)
{
    const size_t N = 12;
    for (size_t i=0; i<prns.size(); ++i) {
        if (i!=0 && i%N==0)
            fprintf(fp, "\n");
        printf(" %3s", prns[i].c_str());
    }
    if (prns.size() != 0u)
        fprintf(fp, "\n");
}

static void push_back(const char *value, std::vector<std::string> &vec)
{
    std::stringstream ss(value);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vstrings(begin, end);
    std::copy(vstrings.begin(), vstrings.end(), std::back_inserter(vec));
}

// static void intersection(const std::vector<std::string> &v1,
//                          const std::vector<std::string> &v2,
//                          std::vector<std::string> &out)
// {
//     std::vector<std::string> vec;
//     std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(vec));
//     out.swap(vec);
// }

static void difference(const std::vector<std::string> &v1,
                       const std::vector<std::string> &v2,
                       std::vector<std::string> &out)
{
    std::vector<std::string> vec;
    std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(vec));
    out.swap(vec);
}

static bool read_time(const char *str, MJD &t)
{
    int y, m, d;
    int h=0, min=0;
    double s=0.0;
    int n = sscanf(str, "%d%d%d%d%d%lf", &y, &m, &d, &h, &min, &s);
    if (n == 2) { // yr doy
        doy2date(y, m, &m, &d);
    }
    else if (n!=3 && n!=6) { // yr m d .or. yr m d h min sec
        return false;
    }
    t.d   = date2mjd(y, m, d);
    t.sod = hms2sod(h, min, s);
    return true;
}

static int handler(void* user, const char* section, const char* name, const char* value)
{
    config_t* pconfig = (config_t*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("session", "date")) {
        if (!read_time(value, pconfig->mjd)) return 0;
    } else if (MATCH("session", "interval")) {
        pconfig->interval = atoi(value);
    } else if (MATCH("session", "length")) {
        pconfig->length = atoi(value);
    } else if (MATCH("constellation", "system")) {
        pconfig->strsys.assign(value);
    } else if (MATCH("constellation", "prn")) {
        push_back(value, pconfig->prns);
    } else if (MATCH("constellation", "exclude")) {
        push_back(value, pconfig->excluded_prns);
    } else if (MATCH("ac", "name")) {
        push_back(value, pconfig->ac_names);
    } else if (MATCH("ac", "combined orbit")) {
        pconfig->orb_ac.assign(value);
    } else if (MATCH("ac", "weight method")) {
        pconfig->weight_method.assign(value);
    } else if (MATCH("product", "path")) {
        pconfig->product_path.assign(value);
    } else if (MATCH("product", "combine staclk")) {
        std::string v(value);
        if (v=="True" || v=="true")
            pconfig->combine_staclk = true;
    } else if (MATCH("product", "phase clock")) {
        std::string v(value);
        if (v=="True" || v=="true")
            pconfig->phase_clock = true;
    } else if (MATCH("product", "use att")) {
        std::string v(value);
        if (v=="True" || v=="true")
            pconfig->use_att = true;
    } else if (MATCH("product", "nav pattern")) {
        pconfig->nav_pattern.assign(value);
    } else if (MATCH("product", "att pattern")) {
        pconfig->att_pattern.assign(value);
    } else if (MATCH("product", "bia pattern")) {
        pconfig->bia_pattern.assign(value);
    } else if (MATCH("product", "clk pattern")) {
        pconfig->clk_pattern.assign(value);
    } else if (MATCH("product", "sp3 pattern")) {
        pconfig->sp3_pattern.assign(value);
    } else if (MATCH("product", "snx pattern")) {
        pconfig->snx_pattern.assign(value);
    // } else if (MATCH("product", "erp pattern")) {
    //     pconfig->erp_pattern.assign(value);
    } else if (MATCH("table", "atx")) {
        pconfig->atx_pattern.assign(value);
    // } else if (MATCH("table", "jpleph")) {
    //     pconfig->eph_path.assign(value);
    } else if (MATCH("output", "align brdc")) {
        std::string v(value);
        if (v=="True" || v=="true")
            pconfig->align_brdc = true;
    } else if (MATCH("output", "product prefix")) {
        pconfig->product_prefix.assign(value);
    } else if (MATCH("output", "cls pattern")) {
        pconfig->cls_pattern.assign(value);
    } else if (MATCH("output", "dif pattern")) {
        pconfig->dif_pattern.assign(value);
    } else if (MATCH("output", "log pattern")) {
        pconfig->log_pattern.assign(value);
    } else {
        //return 0;  /* unknown section/name, error */
    }

    return 1;
}

static bool check_config(const config_t &config)
{
    if (config.mjd.d == 0) {
        fprintf(stderr, MSG_ERR "check_config: [session] date not set\n");
        return false;
    } else if (config.interval == 0) {
        fprintf(stderr, MSG_ERR "check_config: [session] interval not set\n");
        return false;
    } else if (config.length == 0) {
        fprintf(stderr, MSG_ERR "check_config: [session] length not set\n");
        return false;
    } else if (config.strsys.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [constellation] system not set\n");
        return false;
    } else if (config.prns.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [constellation] prn not set\n");
        return false;
    } else if (config.ac_names.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [ac] name not set\n");
        return false;
    } else if (config.orb_ac.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [ac] combined orbit not set\n");
        return false;
    } else if (config.weight_method.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [ac] weight method not set\n");
        return false;
    } else if (config.product_path.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [product] path not set\n");
        return false;
    } else if (config.nav_pattern.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [product] nav pattern not set\n");
        return false;
    } else if (config.clk_pattern.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [product] clk pattern not set\n");
        return false;
    } else if (config.sp3_pattern.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [product] sp3 pattern not set\n");
        return false;
    } else if (config.combine_staclk && config.snx_pattern.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [product] snx pattern not set\n");
        return false;
    } else if (config.phase_clock && config.bia_pattern.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [product] bia pattern not set\n");
        return false;
    } else if (config.use_att && config.att_pattern.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [product] att pattern not set\n");
        return false;
    } else if (config.atx_pattern.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [table] atx pattern not set\n");
        return false;
    } else if (config.product_prefix.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [output] product prefix not set\n");
        return false;
    } else if (config.cls_pattern.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [output] cls pattern not set\n");
        return false;
    } else if (config.dif_pattern.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [output] dif pattern not set\n");
        return false;
    } else if (config.log_pattern.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [output] log pattern not set\n");
        return false;
    } else {
        return true;
    }
}

bool read_config(const std::string &path, config_t &config)
{
    if (ini_parse(path.c_str(), handler, &config) < 0) {
        fprintf(stderr, MSG_ERR "read_config: no such file: %s\n", path.c_str());
        return false;
    }

    if (!check_config(config))
        return false;

    // std::sort(config.strsys.begin(), config.strsys.end());
    std::sort(config.prns.begin(), config.prns.end());
    std::sort(config.excluded_prns.begin(), config.excluded_prns.end());

    if (!collect_valid_prns(config))
        return false;
    return true;
}

bool collect_valid_prns(config_t &config)
{
    // remove prns excluded
    difference(config.prns, config.excluded_prns, config.prns);
    config.excluded_prns.clear();

    // collect valid GNSS
    std::set<char> set;
    std::for_each(config.prns.begin(), config.prns.end(),
                  [&set] (std::string &prn) { set.insert(prn[0]); });
    std::string strtmp(config.strsys);
    config.strsys.clear();
    std::sort(strtmp.begin(), strtmp.end());
    std::set_intersection(strtmp.begin(), strtmp.end(),
                          set.begin(), set.end(), std::back_inserter(config.strsys));

    // remove prns without constellation request
    for (auto it=config.prns.begin(); it!=config.prns.end(); ++it) {
        if (config.strsys.find(it->front()) == std::string::npos)
            config.excluded_prns.push_back(*it);
    }
    difference(config.prns, config.excluded_prns, config.prns);

    // check prn number
    if (config.prns.empty() || config.strsys.empty()) {
        fprintf(stderr, MSG_ERR "collect_valid_prns: empty prns list\n");
        return false;
    }

    // fill vector:syss
    std::set<enum GNSS_Tp> tps;
    for_each(config.prns.begin(), config.prns.end(),
            [&tps] (std::string &prn) { tps.insert(prn2sys(prn)); });
    config.syss.assign(tps.begin(), tps.end());

    // auto igps = std::find(config.syss.begin(), config.syss.end(), _GPS_);
    // if (igps!=config.syss.end() && igps!=config.syss.begin()) {
    //     std::swap(*igps, config.syss.front());
    // }

    // map syss to strsys
    config.strsys.clear();
    for_each(config.syss.begin(), config.syss.end(),
            [&config](enum GNSS_Tp tp) { config.strsys.push_back(sys2char(tp)); });

    // std::map<char, int> sys{ {'G', _GPS_}, {'R', _GLS_}, {'E', _GAL_},
    //                          {'C', _BDS_}, {'J', _QZS_} };
    // for (auto it=config.strsys.begin(); it<config.strsys.end(); ++it)
    //     config.system |= 1<<sys[*it];

    return true;
}

// output configuration
void print_config(FILE *fp, const config_t &config)
{
    int week; double sow;
    mjd2wksow(config.mjd.d, config.mjd.sod, week, sow);
    int dow = static_cast<int>(sow/86400);

    fprintf(fp, "===============start of config===============\n");
    fprintf(fp, "[session]\n");
    fprintf(fp, "time    : %-d %-d (GPSWeek)\n", week, dow);
    fprintf(fp, "length  : %-d\n", config.length);
    fprintf(fp, "interval: %-d\n", config.interval);

    fprintf(fp, "\n[constellation]\n");
    fprintf(fp, "system   : %s\n", config.strsys.c_str());
    // fprintf(fp, "system: %3d\n", config.system);
    fprintf(fp, "satellite:\n");
    print_prns(fp, config.prns);

    fprintf(fp, "\n[ac]\n");
    fprintf(fp, "name: ");
    print_prns(fp, config.ac_names);
    fprintf(fp, "combined orbit: %3s\n", config.orb_ac.c_str());
    fprintf(fp, "weight method: %s\n", config.weight_method.c_str());


    fprintf(fp, "\n[product]\n");
    fprintf(fp, "path: %s\n", config.product_path.c_str());
    std::string str;
    str = config.phase_clock ? "true" : "false";
    fprintf(fp, "phase clock: %s\n", str.c_str());
    str = config.combine_staclk ? "true" : "false";
    fprintf(fp, "combine staclk: %s\n", str.c_str());
    str = config.use_att ? "true" : "false";
    fprintf(fp, "use att: %s\n", str.c_str());
    fprintf(fp, "nav pattern: %s\n", config.nav_pattern.c_str());
    fprintf(fp, "att pattern: %s\n", config.att_pattern.c_str());
    fprintf(fp, "bia pattern: %s\n", config.bia_pattern.c_str());
    fprintf(fp, "clk pattern: %s\n", config.clk_pattern.c_str());
    fprintf(fp, "snx pattern: %s\n", config.snx_pattern.c_str());
    fprintf(fp, "sp3 pattern: %s\n", config.sp3_pattern.c_str());

    fprintf(fp, "\n[table]\n");
    fprintf(fp, "atx: %s\n", config.atx_pattern.c_str());
    // fprintf(fp, "eph: %s\n", config.eph_path.c_str());
    // fprintf(fp, "lps: %s\n", config.lps_path.c_str());

    fprintf(fp, "\n[output]\n");
    str = config.align_brdc ? "true" : "false";
    fprintf(fp, "align brdc: %s\n", str.c_str());
    fprintf(fp, "product prefix: %s\n", config.product_prefix.c_str());
    fprintf(fp, "cls pattern: %s\n", config.cls_pattern.c_str());
    fprintf(fp, "dif pattern: %s\n", config.dif_pattern.c_str());
    fprintf(fp, "log pattern: %s\n", config.log_pattern.c_str());

    fprintf(fp, "================end of config================\n");
    fprintf(fp, "\n");
    fflush(fp);
}

