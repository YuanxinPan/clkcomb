#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>

#include <pppx/chrono.h>

typedef struct TAGCONFIG
{
    // session
    MJD mjd;
    int interval = 0;        // sample rate: 1s, 30s, etc.
    int length = 0;          // time span of combination

    // constellation
    // int system = 0;       // sys |= 1<<GLS
    std::string strsys;      // GREC, etc.
    std::vector<std::string> prns;
    std::vector<std::string> excluded_prns;

    // ac
    std::vector<std::string> ac_names;
    std::string orb_ac;      // combined orbit
    std::string weight_method;

    // product path
    bool phase_clock = false;
    std::string product_path;
    std::string nav_pattern; // broadcast ephemeris
    std::string clk_pattern; // satellite clock
    std::string bia_pattern; // satellite bias
    std::string sp3_pattern; // precise orbit
    std::string erp_pattern; // ERP file

    // table path
    std::string atx_path; // IGS atx
    // std::string eph_path; // JPL ephmeris
    // std::string lps_path; // leap second
    // std::string chn_path; // glonass channel

} config_t;

bool read_config(const std::string &path, config_t &config);

bool collect_valid_prns(config_t &config);

void print_config(FILE *fp, const config_t &config);

#endif  // CONFIG_H
