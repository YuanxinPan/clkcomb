#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>

#include <pppx/const.h>
#include <pppx/chrono.h>

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

} config_t;

bool read_config(const std::string &path, config_t &config);

bool collect_valid_prns(config_t &config);

void print_config(FILE *fp, const config_t &config);

#endif  // CONFIG_H
