#include "config.h"

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

#include <pppx/io.h>
#include <pppx/const.h>

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
    } else if (MATCH("product", "nav pattern")) {
        pconfig->nav_pattern.assign(value);
    } else if (MATCH("product", "clk pattern")) {
        pconfig->clk_pattern.assign(value);
    } else if (MATCH("product", "bia pattern")) {
        pconfig->bia_pattern.assign(value);
    } else if (MATCH("product", "sp3 pattern")) {
        pconfig->sp3_pattern.assign(value);
    } else if (MATCH("product", "snx pattern")) {
        pconfig->snx_pattern.assign(value);
    // } else if (MATCH("product", "erp pattern")) {
    //     pconfig->erp_pattern.assign(value);
    } else if (MATCH("table", "atx")) {
        pconfig->atx_path.assign(value);
    // } else if (MATCH("table", "jpleph")) {
    //     pconfig->eph_path.assign(value);
    // } else if (MATCH("table", "leapsec")) {
    //     pconfig->lps_path.assign(value);
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
    } else if (config.atx_path.empty()) {
        fprintf(stderr, MSG_ERR "check_config: [table] atx path not set\n");
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
    std::string strtmp, systmp=config.strsys;
    std::sort(systmp.begin(), systmp.end()); // to keep original order in config.strsys
    std::set_intersection(systmp.begin(), systmp.end(),
                          set.begin(), set.end(), std::back_inserter(strtmp));
    systmp.clear();
    for (auto it=config.strsys.begin(); it!= config.strsys.end(); ++it)
        if (std::string::npos != strtmp.find(*it))
            systmp.push_back(*it);
    config.strsys.swap(systmp);

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

    fprintf(fp, "[session]\n");
    fprintf(fp, "date    : %-d %-d (GPSWeek)\n", week, dow);
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
    fprintf(fp, "nav pattern: %s\n", config.nav_pattern.c_str());
    fprintf(fp, "sp3 pattern: %s\n", config.sp3_pattern.c_str());
    fprintf(fp, "clk pattern: %s\n", config.clk_pattern.c_str());
    fprintf(fp, "bia pattern: %s\n", config.bia_pattern.c_str());
    fprintf(fp, "snx pattern: %s\n", config.snx_pattern.c_str());

    fprintf(fp, "\n[table]\n");
    fprintf(fp, "atx: %s\n", config.atx_path.c_str());
    // fprintf(fp, "eph: %s\n", config.eph_path.c_str());
    // fprintf(fp, "lps: %s\n", config.lps_path.c_str());

    fflush(fp);
}

