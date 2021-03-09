#include "utils.h"

#include <math.h>
#include <stdlib.h>

#include <algorithm>
#include <map>
#include <numeric>

#include <const.h>
#include <io/io.h>
#include <coord/coord.h>
#include <rinex/rinex.h>

extern FILE *g_logfile;

const size_t MINAC = 3;

std::string replace_pattern(const std::string &pattern, MJD t,
                            const std::string &prefix,
                            const std::string &suffix)
{
    int y, m, d, h, min, week, dow, doy;
    double s;
    mjd2date(t.d, t.sod, &y, &m, &d, &h, &min, &s);
    mjd2doy(t.d, &y, &doy);
    date2gwk(y, m, d, &week, &dow);

    std::string file_name = pattern;
    size_t index;
    char str[64];
    if (std::string::npos != (index = file_name.find("_PREFIX_")))
        file_name.replace(index, 8, prefix);
    if (std::string::npos != (index = file_name.find("_YEAR_")))
        file_name.replace(index, 6, std::to_string(y));
    if (std::string::npos != (index = file_name.find("_YR_"))) {
        sprintf(str, "%02d", y%100);
        file_name.replace(index, 4, str);
    }
    if (std::string::npos != (index = file_name.find("_DOY_"))) {
        sprintf(str, "%03d", doy);
        file_name.replace(index, 5, str);
    }
    if (std::string::npos != (index = file_name.find("_WEEK_"))) {
        sprintf(str, "%04d", week);
        file_name.replace(index, 6, str);
    }
    if (std::string::npos != (index = file_name.find("_DOW_")))
        file_name.replace(index, 5, std::to_string(dow));
    if (std::string::npos != (index = file_name.find("_SUFFIX_")))
        file_name.replace(index, 8, suffix);

    return file_name;
}

bool init_acs(config_t &config, const std::vector<Satellite> &sats,
              std::vector<AnalyseCenter> &acs, AnalyseCenter &combined_ac)
{
    // AC of combined orbit
    fprintf(stdout, "(%s) read products of %3s...\n", run_time(), combined_ac.name.c_str()); fflush(stdout);
    combined_ac.atx_file = replace_pattern(config.atx_pattern, config.mjd, combined_ac.name);
    combined_ac.att_file = config.product_path + replace_pattern(config.att_pattern, config.mjd, combined_ac.name);
    combined_ac.bia_file = config.product_path + replace_pattern(config.bia_pattern, config.mjd, combined_ac.name);
    combined_ac.clk_file = config.product_path + replace_pattern(config.clk_pattern, config.mjd, combined_ac.name);
    combined_ac.snx_file = config.product_path + replace_pattern(config.snx_pattern, config.mjd, combined_ac.name);
    combined_ac.sp3_file = config.product_path + replace_pattern(config.sp3_pattern, config.mjd, combined_ac.name);
    std::string sp3_file;
    std::vector<std::string> sp3_files;
    for (int i=0; i!=2; ++i) {
        MJD t = config.mjd;
        t.d += i;
        sp3_files.push_back(config.product_path + replace_pattern(config.sp3_pattern, t, combined_ac.name));
    }

    if (config.use_att && !combined_ac.read_att(combined_ac.att_file))
        fprintf(stderr, MSG_WAR "use nominal attitude\n");
    // if (!combined_ac.read_orbit(sp3_files/*combined_ac.sp3_file*/) ||
    if (!combined_ac.read_orbit(sp3_files) ||
        !combined_ac.open_atx(combined_ac.atx_file) ||
        (config.combine_staclk && !combined_ac.read_sinex(combined_ac.snx_file))) {
        return false;
    }

    // each AC
    for (auto it=config.ac_names.begin(); it!=config.ac_names.end(); ++it) {
        fprintf(stdout, "(%s) read products of %3s...\n", run_time(), it->c_str()); fflush(stdout);

        acs.emplace_back(*it);
        AnalyseCenter &ac = acs.back();

        sp3_files.clear();
        for (int i=0; i!=2; ++i) {
            MJD t = config.mjd;
            t.d += i;
            sp3_file = config.product_path + replace_pattern(config.sp3_pattern, t, *it);
            sp3_files.push_back(sp3_file);
        }
        ac.atx_file = replace_pattern(config.atx_pattern, config.mjd, *it);
        ac.att_file = config.product_path + replace_pattern(config.att_pattern, config.mjd, *it);
        ac.bia_file = config.product_path + replace_pattern(config.bia_pattern, config.mjd, *it);
        ac.clk_file = config.product_path + replace_pattern(config.clk_pattern, config.mjd, *it);
        ac.snx_file = config.product_path + replace_pattern(config.snx_pattern, config.mjd, *it);
        ac.sp3_file = config.product_path + replace_pattern(config.sp3_pattern, config.mjd, *it);
        // if (!ac.read_orbit(sp3_files/*ac.sp3_file*/) || !ac.open_clock(ac.clk_file) || !ac.open_atx(ac.atx_file)
        if (!ac.read_orbit(sp3_files) || !ac.open_clock(ac.clk_file) || !ac.open_atx(ac.atx_file)
            || (config.use_att && !ac.read_att(ac.att_file))
            || (config.combine_staclk && !ac.read_sinex(ac.snx_file))
            || !ac.read_clock(config, combined_ac.rnxsp3(), combined_ac.rnxatx(), combined_ac.rnxatt())
            || (config.phase_clock && !ac.read_bias(ac.bia_file, config.prns, sats))) {
            acs.pop_back();
            fprintf(stderr, MSG_WAR "remove AC without products: %s\n", it->c_str());
        }
    }

    if (config.combine_staclk) {
        std::map<std::string, size_t> sta_counts;
        for (auto it=acs.begin(); it!=acs.end(); ++it)
        {
            auto sta_names = it->sta_names();
            for (auto name=sta_names.begin(); name!=sta_names.end(); ++name)
            {
                ++sta_counts[*name];
            }
        }
        for (auto it=sta_counts.begin(); it!=sta_counts.end(); ++it)
        {
            if (it->second < MINAC)
                continue;
            config.sta_list.push_back(it->first);
        }
        std::sort(config.sta_list.begin(), config.sta_list.end());

        fprintf(stdout, "\nCOMMON STATION: %3lu\n", config.sta_list.size());
        for (size_t i=0; i!=config.sta_list.size(); ++i)
        {
            if (i!=0 && i%10 == 0) printf("\n");
            fprintf(stdout, "%4s ", config.sta_list[i].c_str());
        }
        fprintf(stdout, "\n\n");
        fflush(stdout);
    }

    // if (!combined_ac.open_clock(combined_ac.clk_file))
    //     return false;
    // combined_ac.read_clock(config, combined_ac.rnxsp3(), combined_ac.rnxatx(), combined_ac.rnxatt());
    if (config.combine_staclk) {
        // combined_ac.read_staclk(config.mjd, config.length, config.sta_interval, config.sta_list, combined_ac.rnxsnx());
        for (auto it=acs.begin(); it!=acs.end(); ++it)
            it->read_staclk(config.mjd, config.length, config.sta_interval, config.sta_list, combined_ac.rnxsnx());
    }

    // Check whether we have enough ACs
    if (acs.size() < MINAC) {
        fprintf(stderr, MSG_ERR "too few ACs: %lu\n", acs.size());
        return false;
    } else {
        return true;
    }
}

bool sat_freq(const std::string &prn, double f[])
{
    switch (prn[0]) {
        case 'G':
        case 'J':
            f[0] = GPS_f1;
            f[1] = GPS_f2;
            break;
        case 'E':
            f[0] = GAL_f1;
            f[1] = GAL_f5;
            break;
        case 'C':
            f[0] = BDS_f1;
            f[1] = BDS_f3;
            break;
        case 'R':
            f[0] = GLS_f1;
            f[1] = GLS_f2;
            break;
        default:
            fprintf(stderr, MSG_ERR "sat_freq: unsupported system %c\n", prn[0]);
            return false;
    }
    return true;
}

bool sat_obstp(const std::string &prn, std::string obstp[])
{
    switch (prn[0]) {
        case 'G':
            obstp[0] = "L1C";
            obstp[1] = "L2W";
            obstp[2] = "C1W";
            obstp[3] = "C2W";
            break;
        case 'E':
            obstp[0] = "L1X"; // L1C L1X
            obstp[1] = "L5X"; // L5Q L5X
            obstp[2] = "C1X";
            obstp[3] = "C5X";
            break;
        case 'C':
            obstp[0] = "L2I";
            obstp[1] = "L6I";
            obstp[2] = "C2I";
            obstp[3] = "C6I";
            break;
        // case 'J':
        default:
            fprintf(stderr, MSG_ERR "sat_obstp: unsupported system %c\n", prn[0]);
            return false;
    }
    return true;
}

enum GNSS_Tp prn2sys(const std::string &prn)
{
    enum GNSS_Tp tp;
    int iprn;
    switch (prn[0]) {
    case 'G':
        tp = _GPS_; break;
    case 'R':
        tp = _GLS_; break;
    case 'E':
        tp = _GAL_; break;
    case 'C':
        iprn = atoi(prn.c_str()+1);
        tp = iprn > 16 ? _BD3_ : _BD2_;
        break;
    case 'J':
        tp = _QZS_; break;
    default:
        tp = _UNKNOWN_;
    }
    return tp;
}

char sys2char(enum GNSS_Tp tp)
{
    switch (tp) {
    case _GPS_:
        return 'G';
    case _GLS_:
        return 'R';
    case _GAL_:
        return 'E';
    case _BD2_:
        return 'C';
    case _BD3_:
        // return 'C';
        return 'B';
    case _QZS_:
        return 'J';
    default:
        return 'X';
    }
}

bool init_sats(const config_t &config, std::vector<Satellite> &sats)
{
    // RinexAtx rnxatx;
    // if (!rnxatx.open(config.atx_pattern))
    //     return false;

    MJD t = config.mjd;
    t.sod = 43200;
    const RinexAtx::atx_t *atx = nullptr;

    sats.clear();
    for (auto it=config.prns.begin(); it!=config.prns.end(); ++it) {
        sats.emplace_back(*it);
        Satellite &sat = sats.back();

        if (!sat_freq(*it, sat.f))
            return false;
        sat.waveLen[0] = LightSpeed/sat.f[0];
        sat.waveLen[1] = LightSpeed/sat.f[1];

        if (config.phase_clock && !sat_obstp(*it, sat.obstp))
            return false;

        // atx = rnxatx.atx(t, sat.prn);
        // if (nullptr != atx) {
        //     sat.svn = atx->svn();
        //     sat.blk = atx->blk();
        // } else {
        //     fprintf(stderr, MSG_WAR "init_sats: no SVN for %3s\n", sat.prn.c_str());
        // }
    }
    return true;
}

bool construct_initclk(const config_t &config,
                       const std::vector<AnalyseCenter> &acs,
                       const std::vector<Satellite> &sats,
                       std::vector<std::vector<double>> &comb_clks,
                       int &refac)
{
    size_t nac_total  = acs.size();
    size_t nsat_total = sats.size();
    size_t nepo_total = acs[0].sat_clks[0].size();

    comb_clks.clear();
    std::vector<double> stds(nac_total);
    std::vector<size_t> nulls;
    for (size_t iprn=0; iprn!=nsat_total; ++iprn)
    {
        // std::vector<double> vals(nac_total);
        bool null = false;
        for (size_t i=0; i<nac_total-1; ++i)
            for (size_t j=i+1; j<nac_total; ++j) {
                double T = config.phase_clock ? 1E9/(sats[iprn].f[0]+sats[iprn].f[1]) : 0.;
                double std = satclk_std(acs[i].sat_clks[iprn], acs[j].sat_clks[iprn], T);
                if (std == 0.0) {
                    null = true;
                    std = 0.099;
                }
                stds[i] += std*1E3;
                stds[j] += std*1E3;
                // fprintf(stderr, "%3s %3s %3s %8.3f\n", sats[iprn].prn.c_str(), acs[i].name.c_str(), acs[j].name.c_str(), std*1E3);
                // vals[i] += std*1E3;
                // vals[j] += std*1E3;
        }
        if (null) {
            nulls.push_back(iprn);
            // for (size_t i=0; i!=nac_total; ++i)
            //     stds[i] -= vals[i];
        }
    }

    int imin = std::min_element(stds.begin(), stds.end()) - stds.begin();
    printf("    reference AC: %3s\n", acs[imin].name.c_str());
    comb_clks = acs[imin].sat_clks;
    refac = imin;

    for (auto it=nulls.begin(); it!=nulls.end(); ++it)
    {
        int n = 0;
        for (size_t epo=0; epo!=nepo_total; ++epo)
            if (comb_clks[*it][epo] != None)
                ++n;

        if (n < nepo_total*0.7) {
            std::vector<double> vals(nac_total);
            for (size_t i=0; i<nac_total-1; ++i)
                for (size_t j=i+1; j<nac_total; ++j) {
                    double T = config.phase_clock ? 1E9/(sats[*it].f[0]+sats[*it].f[1]) : 0.;
                    double std = satclk_std(acs[i].sat_clks[*it], acs[j].sat_clks[*it], T);
                    if (std == 0.0)
                        std = 9.9;
                    vals[i] += std*1E3;
                    vals[j] += std*1E3;
            }
            auto itmin = std::min_element(vals.begin(), vals.end());
            int index = itmin - vals.begin();
            comb_clks[*it].assign(acs[index].sat_clks[*it].begin(), acs[index].sat_clks[*it].end());
            if (fabs(*itmin - 9.9*1E3*(nac_total-1)) > MaxWnd)
                printf("    null: %3s -> %3s\n", sats[*it].prn.c_str(), acs[index].name.c_str());
        }
    }
    return true;
}

bool construct_init_staclk(const config_t &config,
                           const std::vector<AnalyseCenter> &acs,
                           std::vector<std::vector<double>> &comb_clks)
{
    // size_t nac_total  = acs.size();
    size_t nsta_total = config.sta_list.size();
    size_t nepo_total = acs[0].sta_clks[0].size();

    comb_clks.clear();
    comb_clks.resize(nsta_total);
    for (size_t ista=0; ista!=nsta_total; ++ista)
    {
        for (size_t iepo=0; iepo!=nepo_total; ++iepo)
        {
            std::vector<double> clks;
            for (auto it=acs.begin(); it!=acs.end(); ++it)
            {
                if (it->sta_clks[ista][iepo] != None)
                    clks.push_back(it->sta_clks[ista][iepo]);
            }
            if (clks.size() < MINAC)
                comb_clks[ista].push_back(None);
            else
                comb_clks[ista].push_back(median(clks));
        }
    }

    return true;
}

void detect_outlier(const std::vector<double> &vals,
                    const std::vector<double> &wgts,
                    std::vector<int> &deleted)
{
    deleted.clear();
    if (vals.empty())
        return;

    // Distance to median
    double med = median(vals);
    std::vector<double> ress;
    for (auto it=vals.cbegin(), iw=wgts.cbegin(); it!=vals.cend(); ++it, ++iw) {
        ress.push_back(fabs(*it - med)*sqrt(*iw));
    }

    // Delete vals with large median-distance
    double res_med = median(ress);
    double sig = res_med/0.6745;
    for (auto it=ress.begin(); it!=ress.end(); ++it) {
        if (*it*1E3 > 20 && *it/sig > 6.0) { // ?
            deleted.push_back(it - ress.begin());
        }
    }
}

bool align_widelane(const std::vector<Satellite> &sats,
                    std::vector<AnalyseCenter> &acs,
                    AnalyseCenter &combined_ac,
                    const std::vector<int> &prns)
{
    size_t nsat = prns.size();
    size_t nac  = acs.size();

    int refsat = -1; // index of reference satellite used for alignment
    for (auto i=prns.begin(); i!=prns.end(); ++i) {
        size_t n = 0;
        for (auto it=acs.begin(); it!=acs.end(); ++it) {
            if (it->have_bias[*i] != 0)
                ++n;
        }
        if (n == nac) {
            refsat = *i;
            break;
        }
    }
    if (refsat == -1)
        return false;

    std::vector<std::vector<double>> diff(nac); // track change of WL
    // 1. coarse alignment
    for (auto it=acs.begin(); it!=acs.end(); ++it) {
        double d = acs[0].wl_bias[refsat] - it->wl_bias[refsat];
        // diff[it-acs.begin()].assign(nsat, d);
        diff[it-acs.begin()].assign(nsat, .0);
        for (auto i=prns.begin(); i!=prns.end(); ++i)
            if (it->have_bias[*i])
                it->wl_bias[*i] += d;
    }

    // 2. adjust WL to (-0.5, 0.5)
    for (auto i=prns.begin(); i!=prns.end(); ++i) {
        for (size_t j=0; j!=nac; ++j) {
            if (!acs[j].have_bias[*i])
                continue;
            while (acs[j].wl_bias[*i] > 0.5) {
                acs[j].wl_bias[*i] -= 1;
                diff[j][i-prns.begin()] -= 1;
            }
            while (acs[j].wl_bias[*i] < -0.5) {
                acs[j].wl_bias[*i] += 1;
                diff[j][i-prns.begin()] += 1;
            }
        }

        // check for WL cluster at two ends
        std::vector<double> bias;
        for (size_t j=0; j!=nac; ++j)
            if (acs[j].have_bias[*i])
                bias.push_back(acs[j].wl_bias[*i]);

        if (bias.size() < MINAC) {
            combined_ac.have_bias[*i] = 0;
            combined_ac.wl_bias[*i] = None;
            continue;
        }
        double max = *std::max_element(bias.begin(), bias.end());
        double min = *std::min_element(bias.begin(), bias.end());
        if (max - min > 0.75) {
            for (size_t j=0; j!=nac; ++j) {
                if (!acs[j].have_bias[*i])
                    continue;
                if (max - acs[j].wl_bias[*i] > 0.75) {
                    acs[j].wl_bias[*i] += 1;
                    diff[j][i-prns.begin()] += 1;
                }
            }
        }
    }

    // 3. adjust NL
    for (size_t i=0; i!=nac; ++i) {
        for (auto j=prns.begin(); j!=prns.end(); ++j) {
            if (!acs[i].have_bias[*j])
                continue;
            double coef = sats[*j].f[1]/(sats[*j].f[0] - sats[*j].f[1]);
            acs[i].nl_bias[*j] += coef*diff[i][j-prns.begin()];
        }
    }

    // align mean-WL
    int nhave_bias = 0;
    for (auto iprn=prns.begin(); iprn!=prns.end(); ++iprn)
    {
        size_t count = 0;
        for (auto iac=acs.begin(); iac!=acs.end(); ++iac) {
            if (iac->have_bias[*iprn])
                ++count;
        }
        if (count == nac) {
            combined_ac.have_bias[*iprn] = 1;
            ++nhave_bias;
        }
    }

    bool first = true;
    double mean = 0;
    for (auto it=acs.begin(); it!=acs.end(); ++it)
    {
        double sum = 0;
        for (auto j=prns.begin(); j!=prns.end(); ++j) {
            if (combined_ac.have_bias[*j])
                sum += it->wl_bias[*j];
        }
        sum /= nhave_bias;

        if (first) {
            mean = sum;
            first = false;
        } else {
            double diff = mean - sum;
            for (auto j=prns.begin(); j!=prns.end(); ++j) {
                if (it->have_bias[*j])
                    it->wl_bias[*j] += diff;
            }
        }
    }

    // 4. mean WL
    for (auto i=prns.begin(); i!=prns.end(); ++i)
    {
        printf("%3s WL ", sats[*i].prn.c_str());
        // double sum = 0;
        std::vector<double> vals;
        for (auto it=acs.begin(); it!=acs.end(); ++it) {
            if (!it->have_bias[*i]) {
                printf(" %7s", "");
                continue;
            }
            printf(" %7.3f", it->wl_bias[*i]);
            // sum += it->wl_bias[*i];
            if (it->name == "gbm") // gbm/esa only one
                continue;
            vals.push_back(it->wl_bias[*i]);
        }
        // double mean = sum/nac;
        double mean = stable_mean(vals);
        printf(" %7.3f\n", mean);

        if (mean == None) {
            combined_ac.have_bias[*i] = 0;
            combined_ac.wl_bias[*i] = None;
        } else {
            combined_ac.have_bias[*i] = 1;
            combined_ac.wl_bias[*i] = mean;
        }
    }

    return true;
}

void count_satclks(const config_t &config,
                   const std::vector<std::vector<double>> &sat_clks,
                   std::vector<size_t> &counts)
{
    const size_t nsat_total = config.prns.size();
    // const size_t nepo_total = sat_clks[0].size(); // valid epo

    counts.clear();
    counts.resize(nsat_total);
    for (size_t iprn=0; iprn!=nsat_total; ++iprn)
    {
        for (auto it=sat_clks[iprn].cbegin(), end=sat_clks[iprn].cend(); it!=end; ++it) {
            if (*it != None)
                ++counts[iprn];
        }
    }
}

void remove_clock_datum(std::vector<std::vector<double>> &sat_clks,
                        const std::vector<int> &adjust_prns,
                        const std::vector<int> &common_prns)
{
    // size_t nsat_total = sat_clks.size();
    size_t nepo_total = sat_clks[0].size();
    for (size_t epo=0; epo!=nepo_total; ++epo)
    {
        double sum = 0;
        for (auto it=common_prns.cbegin(); it!=common_prns.cend(); ++it)
            sum += sat_clks[*it][epo];
        double mean = sum/common_prns.size();

        for (auto it=adjust_prns.cbegin(); it!=adjust_prns.cend(); ++it)
            if (sat_clks[*it][epo] != None)
                sat_clks[*it][epo] -= mean;
    }
}

void align_clock_datum(std::vector<std::vector<double>> &sat_clks,
                       std::vector<std::vector<double>> &ref_clks,
                       const std::vector<int> &adjust_prns,
                       const std::vector<int> &common_prns)
{
    // size_t nsat_total = sat_clks.size();
    size_t nepo_total = sat_clks[0].size();
    size_t ncommon = common_prns.size();
    for (size_t epo=0; epo!=nepo_total; ++epo)
    {
        // std::vector<double> diff;
        double sum[2] = { 0, 0 };
        for (auto it=common_prns.cbegin(); it!=common_prns.cend(); ++it) {
            sum[0] += ref_clks[*it][epo];
            sum[1] += sat_clks[*it][epo];
            // diff.push_back(sat_clks[*it][epo] - ref_clks[*it][epo]);
        }
        double mean = (sum[1] - sum[0])/ncommon;
        // double mean = stable_mean(diff);

        for (auto it=adjust_prns.cbegin(); it!=adjust_prns.cend(); ++it)
            if (sat_clks[*it][epo] != None)
                sat_clks[*it][epo] -= mean;
    }
}

void align_clock_datum2(std::vector<std::vector<double>> &sat_clks,
                        std::vector<std::vector<double>> &ref_clks,
                        const std::vector<int> &adjust_prns,
                        const std::vector<int> &common_prns,
                        std::vector<std::vector<double>> &sta_clks, int sample_ratio)
{
    // size_t nsat_total = sat_clks.size();
    size_t nepo_total = sat_clks[0].size();
    size_t ncommon = common_prns.size();
    for (size_t epo=0; epo!=nepo_total; ++epo)
    {
        double sum[2] = { 0, 0 };
        // std::vector<double> diff;
        for (auto it=common_prns.cbegin(); it!=common_prns.cend(); ++it) {
            // diff.push_back(sat_clks[*it][epo] - ref_clks[*it][epo]);
            sum[0] += ref_clks[*it][epo];
            sum[1] += sat_clks[*it][epo];
        }
        double mean = (sum[1] - sum[0])/ncommon;
        // double mean = stable_mean(diff);

        for (auto it=adjust_prns.cbegin(); it!=adjust_prns.cend(); ++it)
            if (sat_clks[*it][epo] != None)
                sat_clks[*it][epo] -= mean;

        if (epo % sample_ratio)
            continue;
        for (auto it=sta_clks.begin(); it!=sta_clks.end(); ++it)
            if (it->at(epo/sample_ratio) != None)
                it->at(epo/sample_ratio) -= mean;
    }
}

double satclk_std(const std::vector<double> &clks, const std::vector<double> &refs, double T)
{
    std::vector<double> diffs;
    for (auto it=clks.begin(), ref=refs.begin(); it!=clks.end(); ++it, ++ref) {
        if (*it == None || *ref == None)
            continue;
        diffs.push_back(*it - *ref);
    }
    if (diffs.size() < 0.7*clks.size())
        return 0.0;

    // Bias & standard-deviation
    double bias = std::accumulate(diffs.begin(), diffs.end(), 0.0)/diffs.size();
    if (T != 0.0) {
        bias = floor(bias/T + 0.5)*T;
    }
    double sum = 0;
    for (auto it=diffs.begin(); it!=diffs.end(); ++it) {
        sum += pow(*it - bias, 2);
    }
    double std = sqrt(sum/diffs.size());
    return std;
}

static void remove_bias_seg(int beg, int end, std::vector<double> &biass,
                            std::vector<int> index, std::vector<double> &sat_clks,
                            const Satellite &sat, const std::string &ac_name,
                            const int niter, bool phase_clock)
{
    if (biass.empty())
        return;

    int n=-1, m=-1;
    // get beg(n) & end(m) index of biass
    for (int i=0, size=index.size(); i<size; ++i) {
        if (n==-1 && index[i]>=beg) n=i;
        if (m==-1 && index[i]>=end) m=i;
    }
    if (n==-1) return;
    if (m==-1) m=index.size();

    // Quality control
    std::vector<int> deleted;
    std::vector<int> id(index.begin()+n, index.begin()+m);
    std::vector<double> diffs(biass.begin()+n, biass.begin()+m);
    std::vector<double> wgts(diffs.size(), 1);
    detect_outlier(diffs, wgts, deleted);
    for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
        diffs[*it] = 0;
    }

    int epo_sum = diffs.size();
    int count = epo_sum - deleted.size();
    if (count==0)  return;

    // Bias & standard-deviation
    double bias = std::accumulate(diffs.begin(), diffs.end(), 0.0)/count;
    if (phase_clock) {
        double T = 1E9/(sat.f[0] + sat.f[1]);
        bias = floor(bias/T + 0.5)*T;
    }
    double sum = 0;
    for (auto it=diffs.begin(); it!=diffs.end(); ++it) {
        sum += pow(*it - bias, 2);
    }
    sum -= pow(bias, 2)*deleted.size();
    double std = sqrt(sum/count);

    for (auto i=beg; i!=end; ++i) {
        if (sat_clks[i] != None)
            sat_clks[i] -= bias;
    }

    fprintf(g_logfile, "bias %4d %4d %4d/%-4d %8.3f %8.3f\n", index[n], index[m-1], count, epo_sum, 1E3*bias, 1E3*std);

    // outliers
    //if (!deleted.empty()) {
    //    for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
    //        fprintf(g_logfile, "del %2d %3s %4s %4d\n", niter, ac_name.c_str(), sat.prn.c_str(), id[*it]);
    //    }
    //}
}

void remove_clkbias_seg(const std::string &ac_name, const Satellite &sat, const int niter,
                       std::vector<double> &sat_clks, const std::vector<double> &ref_clks, bool phase_clock)
{
    // Collect bias of each valid epoch
    std::vector<double> biass;
    std::vector<int> id;
    auto ref = ref_clks.cbegin();
    for (auto it=sat_clks.begin(); it!=sat_clks.end(); ++it, ++ref) {
        if (*it == None || *ref == None)
            continue;
        biass.push_back(*it - *ref);
        id.push_back(it - sat_clks.begin());
    }

    std::vector<int> segs; // end point of each bias segment
    for (int i=0, size=sat_clks.size(); i<size-1; ++i) {
        if (sat_clks[i]!=None && sat_clks[i+1]==None) {
            segs.push_back(i+1);
            fprintf(g_logfile, "gap %4d\n", i+1);
        }
    }

    std::vector<double> diffs;
    double rms = 0.;
    for (int i=0, size=biass.size(); i<size-1; ++i) {
        if (id[i+1]-id[i]==1) {
            diffs.push_back(biass[i]-biass[i+1]);
            rms += pow(diffs.back(), 2);
        } else {
            diffs.push_back(0.0);
        }
    }
    rms = sqrt(rms/diffs.size());
    for (int i=0, size=diffs.size(); i<size; ++i) {
        if (fabs(diffs[i]) > std::max(5*rms, 0.05)) {
            segs.push_back(id[i+1]);
            fprintf(g_logfile, "seg %4d\n", id[i+1]);
        }
    }
    segs.push_back(sat_clks.size());
    std::sort(segs.begin(), segs.end());
    std::unique(segs.begin(), segs.end());

    int beg = 0;
    for (auto it=segs.begin(); it!=segs.end(); ++it) {
        remove_bias_seg(beg, *it, biass, id, sat_clks, sat, ac_name, niter, phase_clock);
        beg = *it;
    }
}

void remove_clock_bias(const std::string &ac_name,
                       const Satellite &sat, const int niter,
                       std::vector<double> &sat_clks,
                       const std::vector<double> &ref_clks,
                       bool phase_clock, double &std, bool edit)
{
    // if (!edit && !phase_clock) {
    if (!edit) {
        remove_clkbias_seg(ac_name, sat, niter, sat_clks, ref_clks, phase_clock);
        return;
    }

    // Collect bias of each valid epoch
    std::vector<double> diffs;
    std::vector<int> id;
    auto ref = ref_clks.cbegin();
    for (auto it=sat_clks.begin(); it!=sat_clks.end(); ++it, ++ref) {
        if (*it == None || *ref == None)
            continue;
        diffs.push_back(*it - *ref);
        id.push_back(it - sat_clks.begin());
    }

    // Quality control
    std::vector<int> deleted;
    std::vector<double> wgts(diffs.size(), 1);
    detect_outlier(diffs, wgts, deleted);
    for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
        diffs[*it] = 0;
        if (edit) fprintf(g_logfile, "del %2d %3s %4s %4d\n", niter, ac_name.c_str(), sat.prn.c_str(), id[*it]);
        if (edit) sat_clks[id[*it]] = None; // delete
    }

    int epo_sum = diffs.size();
    int count = epo_sum - deleted.size();
    if (epo_sum == 0) {
        std = 99.9;
        return;
    }

    // Bias & standard-deviation
    double bias = std::accumulate(diffs.begin(), diffs.end(), 0.0)/count;
    double sum = 0;
    for (auto it=diffs.begin(); it!=diffs.end(); ++it) {
        sum += pow(*it - bias, 2);
    }
    sum -= pow(bias, 2)*deleted.size();
    std = sqrt(sum/count);

    double mean = bias;
    if (phase_clock) {
        double waveNL = 1E9/(sat.f[0] + sat.f[1]); // ns
        bias = floor(bias/waveNL + 0.5)*waveNL;
    }

    for (auto it=sat_clks.begin(); it!=sat_clks.end(); ++it) {
        if (*it != None)
            *it -= bias;
    }

    if (count < sat_clks.size()*0.5) {
        std += 20*1E-3;
    }

    double res_cycle = (mean-bias)/1E9*(sat.f[0]+sat.f[1]);
    fprintf(g_logfile, "bias  %4d/%-4d %8.3f %8.3f %8.3f\n", count, epo_sum, 1E3*bias, res_cycle, 1E3*std);
}

void combine_one_epoch(const std::vector<double> &clks,
                       const std::vector<double> &wgts,
                       double &mean,
                       double &wrms,
                       std::vector<int> &deleted)
{
    mean = weighted_mean(clks, wgts, wrms);

    std::vector<double> vals, w;
    if (wrms > 0.02) {
        detect_outlier(clks, wgts, deleted);
        for (int i=0, size=clks.size(); i!=size; ++i) {
            if (std::find(deleted.begin(), deleted.end(), i) != deleted.end())
                continue;
            vals.push_back(clks[i]);
            w.push_back(wgts[i]);
        }
        mean = weighted_mean(vals, w, wrms);

        // if (wrms > 0.02) {
        //     std::vector<int> dels;
        //     detect_outlier(vals, w, dels);
        //     std::vector<double> vs, ws;
        //     for (int i=0, size=vals.size(); i!=size; ++i) {
        //         if (std::find(dels.begin(), dels.end(), i) != dels.end())
        //             continue;
        //         vs.push_back(vals[i]);
        //         ws.push_back(w[i]);
        //     }
        //     mean = weighted_mean(vs, ws, wrms);
        // }

    }
}

double weighted_mean(const std::vector<double> &vals, const std::vector<double> &wgts, double &wrms)
{
    if (vals.size() < MINAC)
        return None;

    double sum = 0, sum_wgt = 0;
    for (auto it=vals.begin(), iw=wgts.begin(); it!=vals.end(); ++it, ++iw)
    {
        sum += *it*(*iw);
        sum_wgt += *iw;
    }
    double mean = sum/sum_wgt;

    // Residuals
    double res_sum = 0, wgt_sum = 0;
    // std::vector<double> ress;
    for (auto it=vals.begin(), iw=wgts.begin(); it!=vals.end(); ++it, ++iw)
    {
        // double res = fabs(*it - mean)*(*iw);
        // ress.push_back(res);

        wgt_sum += *iw;
        res_sum += pow(*it - mean, 2)*(*iw);
    }
    wrms = sqrt(res_sum/wgt_sum);

    return mean;
}

void compare_clks(const std::vector<std::string> &name_list, const std::string &ac_name,
                  const std::vector<std::vector<double>> &src_clks,
                  const std::vector<std::vector<double>> &ref_clks,
                  std::vector<double> &biass,
                  std::vector<double> &stds,
                  bool epoch_output)
{
    const size_t ntotal = name_list.size();

    for (size_t i=0; i!=ntotal; ++i) {
        int epo = 0;
        std::vector<double> diffs;
        for (auto it=src_clks[i].cbegin(), ref=ref_clks[i].cbegin(); it!=src_clks[i].cend(); ++ref, ++it, ++epo) {
            if (*it==None || *ref==None)
                continue;
            double diff = *it - *ref;
            diffs.push_back(diff);
            if (epoch_output)
                fprintf(g_logfile, "cmb-%3s %4d %4s %8.3f\n", ac_name.c_str(), epo, name_list[i].c_str(), 1E3*diff);
        }

        int count = diffs.size();
        double bias, std;
        if (count == 0){
            bias=None;
            std=None;
            biass.push_back(bias);
            stds.push_back(std);
            continue;
        }
        bias = std::accumulate(diffs.begin(), diffs.end(), 0.)/count;
        double sum = 0;
        for (auto it=diffs.cbegin(); it!=diffs.cend(); ++it)
            sum += pow(*it-bias, 2);
        std = sqrt(sum/count);
        //fprintf(g_logfile, "%3s %3s %4d %8.3f %8.3f\n", ac_name.c_str(), name_list[i].c_str(), count, 1E3*bias, 1E3*std);
        biass.push_back(bias);
        stds.push_back(std);
    }
}

void write_satclks(FILE *fp, const config_t &config, const std::string &acn,
                   const std::vector<std::vector<double>> &sat_clks)
{
    size_t nsat_total = sat_clks.size();
    size_t nepo_total = sat_clks[0].size();
    for (size_t epo=0; epo!=nepo_total; ++epo)
    {
        for (size_t prn=0; prn!=nsat_total; ++prn) {
            if (sat_clks[prn][epo] != None)
                fprintf(fp, "%4lu %3s %3s %16.3f\n", epo, acn.c_str(), config.prns[prn].c_str(), sat_clks[prn][epo]);
        }
    }
}

void write_clks_diff(FILE *fp, const std::vector<std::string> &name_list,
                     const std::string &prefix, const std::string &acn,
                     const std::vector<std::vector<double>> &src_clks,
                     const std::vector<std::vector<double>> &ref_clks)
{
    size_t name_total = name_list.size();
    size_t nepo_total = src_clks[0].size();
    for (size_t i=0; i!=name_total; ++i) {
        for (size_t epo=0; epo!=nepo_total; ++epo) {
            if (src_clks[i][epo]!=None && ref_clks[i][epo]!=None)
                fprintf(fp, "%s %3s %3s %4lu %8.3f\n", prefix.c_str(), acn.c_str(),
                        name_list[i].c_str(), epo, src_clks[i][epo] - ref_clks[i][epo]);
            else
                fprintf(fp, "%s %3s %3s %4lu %8s\n", prefix.c_str(), acn.c_str(), name_list[i].c_str(), epo, "NaN");
        }
    }
}

bool write_fip(const std::string &path, MJD t,
               const std::vector<Satellite> &sats,
               const std::vector<double> &wl_bias,
               const std::vector<double> &nl_bias)
{
    FILE *fp = fopen(path.c_str(), "w");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR"failed to create file: %s\n", path.c_str());
        return false;
    }

    int y, m, d, h, min;
    double s;
    mjd2date(t.d, t.sod, &y, &m, &d, &h, &min, &s);

    // header
    fprintf(fp, "Satellite fractional parts of hardware delays               COMMENT\n");
    fprintf(fp, "     Post-processing       undifferenced                    FIP TYPE\n");
    fprintf(fp, " %4d %2d %2d %2d %2d %9.6f                                 TIME OF START\n", y, m, d, 0, 0, 0.);
    fprintf(fp, " %4d %2d %2d %2d %2d %9.6f                                 TIME OF ENDING\n", y, m, d, 23, 59, 40.);
    fprintf(fp, "                                                            END OF HEADER\n");

    const size_t nsat = sats.size();
    for (size_t i=0; i!=nsat; ++i) {
        fprintf(fp, "*%3s %8.3f%8.3f  1\n", sats[i].prn.c_str(), wl_bias[i], 0.004);
        fprintf(fp, " %4d %2d %2d %2d %2d %9.6f  %4d %2d %2d %2d %2d %9.6f %7.3f %14.6f\n",
                y, m, d, 0, 0, 0., y, m, d, 23, 59, 40., nl_bias[i], 0.017);
    }

    fclose(fp);
    return true;
}

void write_clock_datum(const std::vector<std::string> &prns,
                       char sys, const std::string &acn,
                       const std::vector<std::vector<double>> &sat_clks,
                       const std::vector<int> &common_prns)
{
    size_t nepo_total = sat_clks[0].size();
    size_t ncommon = common_prns.size();
    for (size_t epo=0; epo!=nepo_total; ++epo)
    {
        double sum = 0;
        for (auto it=common_prns.cbegin(); it!=common_prns.cend(); ++it) {
            sum += sat_clks[*it][epo];
            fprintf(g_logfile, "raw %3s %3s %6lu %14.3f\n", acn.c_str(), prns[*it].c_str(), epo, sat_clks[*it][epo]);
        }
        double mean = sum/ncommon;
        fprintf(g_logfile, "datum %3s %c %6lu %14.3f\n", acn.c_str(), sys, epo, mean);
    }
}

//bool write_att(const std::string &path, const config_t &config, const RinexAtt &rnxatt)
//{
//    FILE *fp = fopen(path.c_str(), "w");
//    if (fp == nullptr) {
//        fprintf(stderr, "failed to create file: %s\n", path.c_str());
//        return false;
//    }
//
//    MJD t = config.mjd;
//    int interval = 300; // 5 min interval for obx output
//
//    fprintf(fp, "");
//    const std::vector<std::string> &prns = config.prns;
//    size_t nsat = prns.size();
//
//    for (size_t isat=0; isat!=nsat; ++isat)
//    {
//        for ();
//    }
//
//    fclose(fp);
//    return true;
//}

bool write_clkfile(const std::string &path, const config_t &config,
                   const std::vector<std::vector<double>> &satclks,
                   const std::vector<std::vector<double>> &staclks,
                   const RinexSnx &rnxsnx)
{
    FILE *fp = fopen(path.c_str(), "w");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR"failed to create file: %s\n", path.c_str());
        return false;
    }

    MJD t = config.mjd;
    int interval = config.interval;
    int interval_sta = config.sta_interval;
    const std::vector<std::string> &prns = config.prns;
    const std::vector<std::string> &site_list = config.sta_list;

    size_t nsat = prns.size();
    size_t nsta = site_list.size();
    size_t nepo = satclks[0].size();
    size_t nepo_sta = 0;
    if (config.combine_staclk)
        nepo_sta = staclks[0].size();
    size_t ratio = interval_sta/interval;

    // collect valid prns with clk
    std::vector<std::string> valid_prns;
    for (size_t i=0; i!=nsat; ++i) {
        for (size_t j=0; j!=nepo; ++j) {
            if (satclks[i][j] != None) {
                valid_prns.push_back(prns[i]);
                break;
            }
        }
    }

    // collect valid sites with clk
    std::vector<std::string> valid_sites;
    if (config.combine_staclk) {
        for (size_t i=0; i!=nsta; ++i) {
            for (size_t j=0; j!=nepo_sta; ++j) {
                if (staclks[i][j] != None) {
                    valid_sites.push_back(site_list[i]);
                    break;
                }
            }
        }
    }

    // leap second
    int init_leapsec = 21;
    std::vector<int> leap_mjds = { 45150, 45515, 46246, 47160, 47891, 48256, 48803,
                                   49168, 49533, 50082, 50629, 51178, 53735, 54831,
                                   56108, 57203, 57753, 59395 };
    auto it = std::lower_bound(leap_mjds.begin(), leap_mjds.end(), t.d);
    if (it == leap_mjds.end()) {
        fprintf(stderr, MSG_ERR "write_clkfile: %d out of leap_mjds range [%d, %d]\n", t.d, leap_mjds.front(), leap_mjds.back());
        return false;
    }
    int leap_sec = init_leapsec + (it-leap_mjds.begin()) - 1 - 19;

    // write header
    fprintf(fp, "     2.00           CLOCK DATA                              RINEX VERSION / TYPE\n");
    fprintf(fp, "    %2d                                                      LEAP SECONDS        \n", leap_sec);

    if (config.combine_staclk) {
        fprintf(fp, "   %3lu %73s\n", valid_sites.size(), "# OF SOLN STA / TRF ");
        for (size_t i=0; i!=valid_sites.size(); ++i) {
            std::string dome;
            CartCoor pos;
            rnxsnx.find_dome(valid_sites[i], dome);
            rnxsnx.find_pos(valid_sites[i], pos);
            fprintf(fp, "%4s %9s          %12.0f%12.0f%12.0f%s\n", valid_sites[i].c_str(),
                    dome.c_str(), pos.x*1E3, pos.y*1E3, pos.z*1E3, "SOLN STA NAME / NUM");
        }
    }

    fprintf(fp, "   %3lu %73s\n", valid_prns.size(), "# OF SOLN SATS      ");
    size_t n = (valid_prns.size()+14)/15*15;
    for (size_t i=0; i!=n; ++i) {
        if (i<valid_prns.size())
            fprintf(fp, "%3s ", valid_prns[i].c_str());
        else
            fprintf(fp, "    ");
        if ((i+1)%15 == 0)
            fprintf(fp, "PRN LIST            \n");
    }
    fprintf(fp, "%60sEND OF HEADER       \n", " ");

    // write clks
    int y, m, d, h, min;
    double s;
    for (size_t i=0; i!=nepo; ++i) {
        mjd2date(t.d, interval*i, &y, &m, &d, &h, &min, &s);

        if (config.combine_staclk && i%ratio == 0) {
            for (size_t j=0; j!=nsta; ++j) {
                if (staclks[j][i/ratio] == None)
                    continue;
                fprintf(fp, "AR %4s %4d %02d %02d %02d %02d %9.6f  1 %21.12E\n",
                        site_list[j].c_str(), y, m, d, h, min, s, staclks[j][i/ratio]);
            }
        }

        for (size_t j=0; j!=nsat; ++j) {
            if (satclks[j][i] == None)
                continue;
            fprintf(fp, "AS %3s  %4d %02d %02d %02d %02d %9.6f  1 %21.12E\n",
                    prns[j].c_str(), y, m, d, h, min, s, satclks[j][i]);
        }
    }

    fclose(fp);
    return true;
}

static void convert2raw_bias(double wl, const double f[], double raw[])
{
    double nl  = 0.0;
    double dcb = 0.0; // ns
    double g = f[0]/f[1];
    double dw = 1/(f[0] - f[1]);
    double dn = 1/(f[0] + f[1]);
    double alpha = g*g/(g*g - 1);
    double beta = 1 - alpha;

    raw[0] = (-wl/g*dw + (1+1/g)*nl*dn)*1E9 - beta*dcb; // L1W
    raw[1] = (-wl*g*dw + (1 + g)*nl*dn)*1E9 + alpha*dcb ; // L2W
    raw[2] = beta*dcb; // C1W
    raw[3] = -alpha*dcb; // C2W
}

bool write_bias(const std::string &path, MJD t,
                const std::vector<Satellite> &sats,
                const std::vector<int> &have_bias,
                const std::vector<double> &wl_bias)
{
    FILE *fp = fopen(path.c_str(), "w");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR"failed to create file: %s\n", path.c_str());
        return false;
    }

    // header
    fprintf(fp, "+BIAS/SOLUTION\n");
    fprintf(fp, "*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___ __ESTIMATED_SLOPE____ _STD_DEV___\n");

    const size_t nsat = sats.size();
    int y, doy;
    double raw[4];
    // const char *types[4] = { "L1C", "L2W", "C1W", "C2W" };
    mjd2doy(t.d, &y, &doy);
    for (size_t i=0; i!=nsat; ++i) {
        if (!have_bias[i])
            continue;
        convert2raw_bias(wl_bias[i], sats[i].f, raw);
        for (int j=0; j!=4; ++j)
            fprintf(fp, " OSB  %4s %3s %9s %-4s %-4s %4d:%03d:%05d %4d:%03d:%05d %-4s %21.3f %11.3f\n",
                    sats[i].svn.c_str(), sats[i].prn.c_str(), "", sats[i].obstp[j].c_str(), "", y, doy, 0, y, doy+1, 0, "ns", raw[j], 0.0);
    }

    fprintf(fp, "-BIAS/SOLUTION\n");
    fprintf(fp, "%%=ENDBIA\n");

    fclose(fp);
    return true;
}

void clkfit(std::vector<double> clks, const int n, double &a, double &b)
{
    int deleted = 1;
    while (deleted != 0) {
        deleted = 0;
        double xsum = 0, ysum = 0, xysum = 0, x2sum = 0;
        int count = 0;
        for (int iepo = 0; iepo != n; ++iepo) {
            if (clks[iepo] != None) {
                xsum += iepo;             ysum += clks[iepo];
                xysum += iepo*clks[iepo]; x2sum += iepo*iepo;
                ++count;
            }
        }
        a = (x2sum*ysum - xsum * xysum) / (x2sum*count - xsum * xsum);
        b = (count*xysum - xsum * ysum) / (count*x2sum - xsum * xsum);

        std::vector<double> res;
        for (int iepo = 0; iepo != n; ++iepo) {
            if (clks[iepo] == None)
                continue;
            res.push_back(clks[iepo] - (a + b*iepo));
            //fprintf(stderr, "res %12.2f\n", clks[iepo] - (a + b * iepo));
        }
        double med = median(res);

        std::vector<double> ress;
        for (int iepo = 0; iepo != n; ++iepo) {
            if (clks[iepo] == None)
                continue;
            ress.push_back(fabs(clks[iepo] - (a + b*iepo) - med));
            //fprintf(stderr, "ress %20.12f %12.5f\n", clks[iepo],fabs(clks[iepo] - (a + b * iepo) - med));
        }
        double mad = median(ress)/0.6745;
        //fprintf(stderr, "%12.5f %12.5f\n", med, mad);

        for (int iepo = 0; iepo != n; ++iepo) {
            if (clks[iepo] != None && fabs(clks[iepo] - (a + b*iepo) - med) > 6*mad) {
                fprintf(g_logfile, "outliers in brdc alignment %4d %20.4f\n", iepo, clks[iepo] - (a + b*iepo));
                clks[iepo] = None;
                ++deleted;
            }
        }
    }
}

bool write_summary(const std::string &path, config_t &config,
                   const std::vector<AnalyseCenter> &acs,
                   const std::vector<Satellite> &sats,
                   std::vector<std::vector<double>> &comb_clks,
                   std::vector<std::vector<double>> &comb_staclks)
{
    FILE *fp = fopen(path.c_str(), "w");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR"failed to creat file: %s\n", path.c_str());
        return false;
    }

    size_t nac_total = acs.size();
    size_t nsat_total = sats.size();
    size_t nsta_total = config.sta_list.size();
    size_t nepo_total = config.length/config.interval + 1;

    // Compare with combined clock
    std::vector<std::vector<double>> sta_bias(nac_total); // ac station RMS
    std::vector<std::vector<double>> sta_std(nac_total);  // ac station RMS
    if (config.combine_staclk) {
        for (size_t i=0; i!=nac_total; ++i) {
            std::vector<double> tmp_sta_bias;
            std::vector<double> tmp_sta_std;
            compare_clks(config.sta_list, acs[i].name, acs[i].init_sta_clks, comb_staclks, sta_bias[i], tmp_sta_std, false);
            compare_clks(config.sta_list, acs[i].name, acs[i].sta_clks, comb_staclks, tmp_sta_bias, sta_std[i], false);

        }
    }
    std::vector<std::vector<double>> bias(nac_total);  // ac sattelite RMS
    std::vector<std::vector<double>> std(nac_total);   // ac sattelite RMS
    for (size_t i=0; i!=nac_total; ++i) {
        std::vector<double> tmp_bias;
        std::vector<double> tmp_std;
        compare_clks(config.prns, acs[i].name, acs[i].init_sat_clks, comb_clks, bias[i], tmp_std, false);
        compare_clks(config.prns, acs[i].name, acs[i].sat_clks, comb_clks, tmp_bias, std[i], false);
    }

    std::vector<double> rms_total(nac_total);
    std::vector<int> count_total(nac_total);

    std::vector<double> rms_sat(nsat_total);
    for (size_t isat=0; isat!=nsat_total; ++isat) {
        rms_sat[isat] = 0.0;
        int count_sat = 0;
        for (size_t i=0; i!=nac_total; ++i){
            if (std[i][isat] == None) continue;
            rms_sat[isat] += pow(std[i][isat]*1E3,2);
            rms_total[i] += pow(std[i][isat]*1E3,2);
            count_total[i]++;
            count_sat++;
        }
        if (count_sat == 0){
            rms_sat[isat] = None;
            continue;
        }
        rms_sat[isat] = sqrt(rms_sat[isat]/count_sat);
    }

    std::vector<double> rms_sta(nsta_total);
    if(config.combine_staclk){
        for (size_t ista=0; ista!=nsta_total; ++ista)
        {
            rms_sta[ista] = 0.0;
            int count_sta = 0;
            for (size_t i=0; i!=nac_total; ++i){
                if (sta_std[i][ista] == None) continue;
                rms_sta[ista] += pow(sta_std[i][ista]*1E3,2);
                rms_total[i] += pow(sta_std[i][ista]*1E3,2);
                count_total[i]++;
                count_sta++;
            }
            if (count_sta == 0){
                rms_sta[ista] = None;
                continue;
            }
            rms_sta[ista] = sqrt(rms_sta[ista]/count_sta);
        }
    }

    for (size_t i=0; i!=nac_total; ++i){
        if (count_total[i] == 0){
            rms_total[i] = None;
        }
        rms_total[i] = sqrt(rms_total[i]/count_total[i]);
    }

    std::vector<int> nepo_sat(nsat_total);
    for (size_t isat=0; isat!=nsat_total; ++isat){
        nepo_sat[isat]=0;
        for (size_t epo=0; epo!=nepo_total; ++epo){
            if (comb_clks[isat][epo]!=None){
                nepo_sat[isat] += 1;
            }
        }
    }

    std::vector<int> nepo_sta(nsta_total);
    for (size_t ista=0; ista!=nsta_total; ++ista){
        nepo_sta[ista]=0;
        for (size_t epo=0, total=config.length/config.sta_interval+1; epo!=total; ++epo){
            if ( comb_staclks[ista][epo]!=None){
                nepo_sta[ista] += 1;
            }
        }
    }

    int week; double sow;
    mjd2wksow(config.mjd.d, config.mjd.sod, week, sow);
    int dow = int(sow/86400);

    // Summary
    fprintf(fp, "\n");
    fprintf(fp, "%s\n", " cls-file Version 2.0");
    fprintf(fp, "%10s %4d %6s %d %6s %5d\n", " GPS week:", week, "Day:", dow, "MJD:", config.mjd.d);
    fprintf(fp, "\n");
    fprintf(fp, "%s\n", " RESULTS OF FINAL WEIGHTED COMBINATION");
    fprintf(fp, "%s\n", " -------------------------------------");
    fprintf(fp, "\n");

    // Summary--WEIGHT
    fprintf(fp, "%13s\n", " CEN | WEIGHT");
    fprintf(fp, "%13s\n", "     |    (%)");
    fprintf(fp, "%13s\n", " ----+-------");
    for (size_t i=0; i!=nac_total; ++i){
        fprintf(fp, " %3s | %6.2f\n", acs[i].name.c_str(), 100.0*acs[i].weight);
    }
    fprintf(fp, "\n");

    // Summary--MEAN bias
    //fprintf(stderr,"%s\n", "summary-bias");
    fprintf(fp, "      | MEAN (ps) OF AC CLOCK COMPARED TO COMBINATION\n");
    fprintf(fp, "      |");
    for (size_t i=0; i!=nac_total; ++i){
        fprintf(fp, "%6s", acs[i].name.c_str());
    }
    fprintf(fp," |\n");
    fprintf(fp, " -----+");
    for (size_t i=0; i!=nac_total; ++i){
        fprintf(fp, "------");
    }
    fprintf(fp, "-+\n");
    for (size_t isat=0; isat!=nsat_total; ++isat) {
        fprintf(fp, "%5s |", config.prns[isat].c_str());
        for (size_t i=0; i!=nac_total; ++i){
            if (bias[i][isat] == None) {
                fprintf(fp, "%6s", "X");
                continue;
            }
            fprintf(fp, "%6d", int(bias[i][isat]*1E3));
        }
        fprintf(fp," |\n");
    }
    fprintf(fp, " -----+");
    for (size_t i=0; i!=nac_total; ++i){
        fprintf(fp, "------");
    }
    fprintf(fp, "-+\n");

    if(config.combine_staclk){
        for (size_t ista=0; ista!=nsta_total; ++ista)
        {
            fprintf(fp, "%5s |", config.sta_list[ista].c_str());
            for (size_t i=0; i!=nac_total; ++i){
                if (sta_bias[i][ista] == None) {
                    fprintf(fp, "%6s", "X");
                    continue;
                }
                fprintf(fp, "%6d", int(sta_bias[i][ista]*1E3));
            }
            fprintf(fp, " |\n");
        }
        fprintf(fp, " -----+");
        for (size_t i=0; i!=nac_total; ++i){
            fprintf(fp, "------");
        }
        fprintf(fp, "-+\n");
    }

    fprintf(fp, "\n");

    // Summary--RMS
    //fprintf(stderr,"%s\n", "summary-rms");
    fprintf(fp, "%s", "      | RMS (ps) OF AC CLOCK COMPARED TO COMBINATION\n");
    fprintf(fp, "%s", "      |");
    for (size_t i=0; i!=nac_total; ++i){
        fprintf(fp, "%6s", acs[i].name.c_str());
    }
    fprintf(fp, " | NEPO    RMS\n");
    fprintf(fp, " -----+");
    for (size_t i=0; i!=nac_total; ++i){
        fprintf(fp, "------");
    }
    fprintf(fp, "-+-------------\n");

    for (size_t isat=0; isat!=nsat_total; ++isat) {
        fprintf(fp, "%5s |", config.prns[isat].c_str());
        for (size_t i=0; i!=nac_total; ++i){
            if (std[i][isat] == None) {
                fprintf(fp, "%6s", "X");
                continue;
            }
            fprintf(fp, "%6d", int(std[i][isat]*1E3));
        }
        if (rms_sat[isat] == None) {
            fprintf(fp," | %4d%7s\n", nepo_sat[isat],"X");
            continue;
        }
        fprintf(fp," | %4d%7d\n", nepo_sat[isat],int(rms_sat[isat]));
    }

    fprintf(fp, " -----+");
    for (size_t i=0; i!=nac_total; ++i){
        fprintf(fp, "------");
    }
    fprintf(fp, "-+-------------\n");

    if(config.combine_staclk){
        for (size_t ista=0; ista!=nsta_total; ++ista)
        {
            fprintf(fp, "%5s |", config.sta_list[ista].c_str());
            for (size_t i=0; i!=nac_total; ++i){
                if (sta_std[i][ista] == None) {
                    fprintf(fp, "%6s", "X");
                    continue;
                }
                fprintf(fp, "%6d", int(sta_std[i][ista]*1E3));
            }
            if (rms_sta[ista] == None) {
                fprintf(fp," | %4d%7s\n", nepo_sta[ista],"X");
                continue;
            }
            fprintf(fp," | %4d%7d\n", nepo_sta[ista],int(rms_sta[ista]));
        }
        fprintf(fp, " -----+");
        for (size_t i=0; i!=nac_total; ++i){
            fprintf(fp, "------");
        }
        fprintf(fp, "-+-------------\n");
    }
    fprintf(fp, " TOT  |");
    for (size_t i=0; i!=nac_total; ++i){
        if (rms_total[i] == None){
            fprintf(fp, "%6s", "X");
            continue;
        }
        fprintf(fp, "%6d", int(rms_total[i]));
    }
    fprintf(fp, " |\n");

    fclose(fp);
    return true;
}

double median(const std::vector<double> &vals)
{
    size_t s = vals.size();
    std::vector<double> v(vals);
    std::sort(v.begin(), v.end());

    if (s == 0) {
        return 0;
    } else if (s%2 == 0) {
        return (v[s/2-1] + v[s/2])/2;
    } else {
        return v[s/2];
    }
}

double stable_mean(const std::vector<double> &vals)
{
    // std::vector<double> vals(v);
    std::vector<double> wgts(vals.size(), 1.0);
    std::vector<int> deleted;
    detect_outlier(vals, wgts, deleted);
    if (!deleted.empty()) {
        for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
            // vals[*it] = 0;
            wgts[*it] = 0;
        }
    }
    double wrms;
    double mean = weighted_mean(vals, wgts, wrms);
    return mean;
}

void unit_weight(std::vector<double> &wgts)
{
    double sum = 0;
    for (auto it=wgts.begin(); it!=wgts.end(); ++it)
        sum += *it;
    for (auto it=wgts.begin(); it!=wgts.end(); ++it)
        *it /= sum;
}
