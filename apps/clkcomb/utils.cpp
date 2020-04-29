#include "utils.h"

#include <math.h>
#include <numeric>
#include <pppx/io.h>
#include <pppx/const.h>

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

bool init_acs(const config_t &config, std::vector<AnalyseCenter> &acs, AnalyseCenter &combined_ac)
{
    for (auto it=config.ac_names.begin(); it!=config.ac_names.end(); ++it) {
        acs.emplace_back(*it);
        AnalyseCenter &ac = acs.back();

        std::vector<std::string> sp3_files;
        for (int i=0; i!=3; ++i) {
            MJD t = config.mjd;
            t.d += (i-1);
            sp3_files.push_back(config.product_path + replace_pattern(config.sp3_pattern, t, *it));
        }
        ac.sp3_file = config.product_path + replace_pattern(config.sp3_pattern, config.mjd, *it);
        ac.clk_file = config.product_path + replace_pattern(config.clk_pattern, config.mjd, *it);
        ac.bia_file = config.product_path + replace_pattern(config.bia_pattern, config.mjd, *it);
        // if (!ac.read_orbit(sp3_files) || !ac.open_clock(ac.clk_file)
        if (!ac.read_orbit(ac.sp3_file) || !ac.open_clock(ac.clk_file)
            || (config.phase_clock && !ac.read_bias(ac.bia_file, config.prns))) {
            acs.pop_back();
            fprintf(stderr, MSG_WAR "remove AC without products: %s\n", it->c_str());
        }
    }

    std::vector<std::string> sp3_files;
    for (int i=0; i!=3; ++i) {
        MJD t = config.mjd;
        t.d += (i-1);
        sp3_files.push_back(config.product_path + replace_pattern(config.sp3_pattern, t, combined_ac.name));
    }
    combined_ac.sp3_file = config.product_path + replace_pattern(config.sp3_pattern, config.mjd, combined_ac.name);
    combined_ac.clk_file = config.product_path + replace_pattern(config.clk_pattern, config.mjd, combined_ac.name);
    combined_ac.bia_file = config.product_path + replace_pattern(config.bia_pattern, config.mjd, combined_ac.name);
    // if (!combined_ac.read_orbit(sp3_files) || !combined_ac.open_clock(combined_ac.clk_file)) {
    if (!combined_ac.read_orbit(combined_ac.sp3_file) || !combined_ac.open_clock(combined_ac.clk_file)) {
        fprintf(stderr, MSG_ERR "no combined sp3 file: %s\n", combined_ac.sp3_file.c_str());
        return false;
    }
    combined_ac.read_clock(config.mjd, config.length, config.interval, config.prns, combined_ac.rnxsp3());

    // Check whether we have enough ACs
    if (acs.size() < 3u) {
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
            f[1] = GAL_f8; // ?
            break;
        case 'C':
            f[0] = BDS_f1;
            f[1] = BDS_f2;
            break;
        default:
            fprintf(stderr, MSG_ERR "sat_freq: unsupported system %c\n", prn[0]);
            return false;
    }
    return true;
}

bool init_sats(const std::vector<std::string> &prns, std::vector<Satellite> &sats)
{
    sats.clear();
    for (auto it=prns.begin(); it!=prns.end(); ++it) {
        sats.emplace_back(*it);
        Satellite &sat = sats.back();

        if (!sat_freq(*it, sat.f))
            return false;
        sat.waveLen[0] = LightSpeed/sat.f[0];
        sat.waveLen[1] = LightSpeed/sat.f[1];
    }
    return true;
}

bool construct_initclk(const config_t &config,
                       const std::vector<AnalyseCenter> &acs,
                       const std::vector<Satellite> &sats,
                       std::vector<std::vector<double>> &comb_clks)
{
    size_t nac_total  = acs.size();
    size_t nsat_total = sats.size();
    size_t nepo_total = acs[0].sat_clks[0].size();

    comb_clks.clear();
    std::vector<double> stds(nac_total);
    std::vector<size_t> nulls;
    for (size_t iprn=0; iprn!=nsat_total; ++iprn)
    {
        printf("%3s std", sats[iprn].prn.c_str());
        // std::vector<double> stds(nac_total);
        std::vector<double> vals(nac_total);
        bool null = false;
        for (size_t i=0; i<nac_total-1; ++i)
            for (size_t j=i+1; j<nac_total; ++j) {
                double T = config.phase_clock ? 1E9/(sats[iprn].f[0]+sats[iprn].f[1]) : 0.;
                double std = satclk_std(acs[i].sat_clks[iprn], acs[j].sat_clks[iprn], T);
                if (std == 0.0)
                    null = true;
                stds[i] += std*1E3;
                stds[j] += std*1E3;
                vals[i] += std*1E3;
                vals[j] += std*1E3;
                printf(" %12.3f", std*1E3);
        }
        if (null) {
            nulls.push_back(iprn);
            for (size_t i=0; i!=nac_total; ++i)
                stds[i] -= vals[i];
        }

        printf("\n>>>     %12.3f %12.3f %12.3f\n", vals[0], vals[1], vals[2]);
    }
    auto imin = std::min_element(stds.begin(), stds.end());
    printf("reference: %12.3f %12.3f %12.3f => %3s\n", stds[0], stds[1], stds[2], acs[imin-stds.begin()].name.c_str());
    comb_clks = acs[imin-stds.begin()].sat_clks;

    for (auto it=nulls.begin(); it!=nulls.end(); ++it)
    {
        int n = 0;
        for (size_t epo=0; epo!=nepo_total; ++epo)
            if (comb_clks[*it][epo] != None)
                ++n;

        if (n < nepo_total*0.8) {
            std::vector<double> vals(nac_total);
            for (size_t i=0; i<nac_total-1; ++i)
                for (size_t j=i+1; j<nac_total; ++j) {
                    double T = config.phase_clock ? 1E9/(sats[*it].f[0]+sats[*it].f[1]) : 0.;
                    double std = satclk_std(acs[i].sat_clks[*it], acs[j].sat_clks[*it], T);
                    if (std == 0.0)
                        std = 99.9;
                    vals[i] += std*1E3;
                    vals[j] += std*1E3;
            }
            auto itmin = std::min_element(vals.begin(), vals.end());
            int index = itmin - vals.begin();
            comb_clks[*it].assign(acs[index].sat_clks[*it].begin(), acs[index].sat_clks[*it].end());
            printf("null %3s %3s\n", sats[*it].prn.c_str(), acs[index].name.c_str());
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

    std::vector<double> diff[nac]; // track change of WL
    // 1. coarse alignment
    for (auto it=acs.begin(); it!=acs.end(); ++it) {
        double d = acs[0].wl_bias[refsat] - it->wl_bias[refsat];
        diff[it-acs.begin()].assign(nsat, d);
        for (auto i=prns.begin(); i!=prns.end(); ++i)
            it->wl_bias[*i] += d;
    }

    // 2. adjust WL to (-0.5, 0.5)
    for (auto i=prns.begin(); i!=prns.end(); ++i) {
        for (size_t j=0; j!=nac; ++j) {
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
        bool pos = false, neg = false;
        for (size_t j=0; j!=nac; ++j) {
            if (acs[j].wl_bias[*i] >  0.45) pos=true;
            if (acs[j].wl_bias[*i] < -0.45) neg=true;
        }
        if (pos && neg) {
            for (size_t j=0; j!=nac; ++j) {
                if (acs[j].wl_bias[*i] > 0.45)
                    continue;
                acs[j].wl_bias[*i] += 1;
                diff[j][i-prns.begin()] += 1;
            }
        }
    }

    // 3. mean WL
    for (auto i=prns.begin(); i!=prns.end(); ++i)
    {
        printf("%3s WL ", sats[*i].prn.c_str());
        // double sum = 0;
        std::vector<double> vals;
        std::vector<double> wgts;
        for (auto it=acs.begin(); it!=acs.end(); ++it) {
            printf(" %7.3f", it->wl_bias[*i]);
            // sum += it->wl_bias[*i];
            vals.push_back(it->wl_bias[*i]);
            wgts.push_back(1.0);
        }
        // double mean = sum/nac;
        std::vector<int> deleted;
        detect_outlier(vals, wgts, deleted);
        if (!deleted.empty()) {
            for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
                vals[*it] = 0;
                wgts[*it] = 0;
            }
        }
        double wrms;
        double mean = weighted_mean(vals, wgts, wrms);
        printf(" %7.3f\n", mean);

        for (size_t j=0; j!=nac; ++j) {
            double d = mean - acs[j].wl_bias[*i];
            acs[j].wl_bias[*i] += d;
            diff[j][i-prns.begin()] += d;
        }
    }

    // 4. adjust NL
    for (size_t i=0; i!=nac; ++i) {
        for (auto j=prns.begin(); j!=prns.end(); ++j) {
            double coef = sats[*j].f[1]/(sats[*j].f[0] - sats[*j].f[1]);
            acs[i].nl_bias[*j] += coef*diff[i][j-prns.begin()];
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
        fprintf(stdout, "%3s %4lu\n", config.prns[iprn].c_str(), counts[iprn]);
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

double satclk_std(const std::vector<double> &clks, const std::vector<double> &refs, double T)
{
    std::vector<double> diffs;
    for (auto it=clks.begin(), ref=refs.begin(); it!=clks.end(); ++it, ++ref) {
        if (*it == None || *ref == None)
            continue;
        diffs.push_back(*it - *ref);
    }
    if (diffs.size() < 0.9*clks.size())
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
                            std::vector<int> index, std::vector<double> &sat_clks, bool phase_clock)
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
    fprintf(stdout, "%4d %4d => %4d %4d\n", n, m, index[n], index[m-1]);

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
        double T = 1E9/(GPS_f1 + GPS_f2);
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

    fprintf(stdout, "bias:  %4d/%4d %8.3f %8.3f\n", count, epo_sum, 1E3*bias, 1E3*std);

    // Debug
    if (!deleted.empty()) {
        for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
            fprintf(stdout, " %4d", id[*it]);
        }
        fprintf(stdout, "\n");
    }
}

void remove_clkbias_seg(const std::string &name, const std::string &prn,
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
            printf("fault %4d\n", i+1);
        }
    }

    std::vector<double> diffs;
    double rms = 0.;
    for (int i=0, size=biass.size(); i<size-1; ++i) {
        // printf("diff %3s %3s %7d %12.3f\n", name.c_str(), prn.c_str(), id[i], biass[i]-biass[i+1]);
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
            printf("seg %4d\n", id[i+1]);
        }
    }
    segs.push_back(sat_clks.size());
    std::sort(segs.begin(), segs.end());
    std::unique(segs.begin(), segs.end());

    int beg = 0;
    for (auto it=segs.begin(); it!=segs.end(); ++it) {
        remove_bias_seg(beg, *it, biass, id, sat_clks, phase_clock);
        beg = *it;
    }
}

// void remove_clock_bias(const std::string &prn,
void remove_clock_bias(const std::string &name, const std::string &prn,
                       std::vector<double> &sat_clks,
                       const std::vector<double> &ref_clks,
                       bool phase_clock, double &std, bool edit)
{
    // if (!edit && !phase_clock) {
    if (!edit) {
        remove_clkbias_seg(name, prn, sat_clks, ref_clks, phase_clock);
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

    if (phase_clock) {
        double f[2];
        sat_freq(prn, f);
        double waveNL = 1E9/(f[0] + f[1]); // ns
        bias = floor(bias/waveNL + 0.5)*waveNL;
    }
    for (auto it=sat_clks.begin(); it!=sat_clks.end(); ++it) {
        if (*it != None)
            *it -= bias;
    }

    if (count < sat_clks.size()*0.5) {
        std += 20*1E-3;
    }

    fprintf(stdout, "bias:  %4d/%4d %8.3f %8.3f\n", count, epo_sum, 1E3*bias, 1E3*std);

    // Debug
    if (!deleted.empty()) {
        for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
            fprintf(stdout, " %4d", id[*it]);
        }
        fprintf(stdout, "\n");
    }
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
    if (vals.size() < 2u)
        return None;

    double sum = 0, sum_wgt = 0;
    for (auto it=vals.begin(), iw=wgts.begin(); it!=vals.end(); ++it, ++iw)
    {
        sum += *it*(*iw);
        sum_wgt += *iw;
    }
    double mean = sum/sum_wgt;

    // Residuals
    // double res;
    double res_sum = 0, wgt_sum = 0;
    // std::vector<double> ress;
    for (auto it=vals.begin(), iw=wgts.begin(); it!=vals.end(); ++it, ++iw)
    {
        // res = fabs(*it - mean)*(*iw);
        // ress.push_back(res);

        wgt_sum += *iw;
        res_sum += pow(*it - mean, 2)*(*iw);
    }
    wrms = sqrt(res_sum/wgt_sum);

    return mean;
}

void compare_satclks(const config_t &config,
                     const std::vector<std::vector<double>> &sat_clks,
                     const std::vector<std::vector<double>> &ref_clks)
{
    const size_t nsat_total = config.prns.size();
    const size_t nepo_total = sat_clks[0].size();

    double diff;
    for (size_t iprn=0; iprn!=nsat_total; ++iprn)
    {
        for (size_t epo=0; epo!=nepo_total; ++epo)
        {
            if (sat_clks[iprn][epo]==None || ref_clks[iprn][epo]==None)
                continue;
            diff = sat_clks[iprn][epo] - ref_clks[iprn][epo];
            fprintf(stdout, "%4lu %3s com-igs %8.3f\n", epo, config.prns[iprn].c_str(), 1E3*diff);
        }
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

bool write_clkfile(const std::string &path, MJD t, int length, int interval,
                   const std::vector<std::string> &prns,
                   const std::vector<std::vector<double>> &satclks)
{
    FILE *fp = fopen(path.c_str(), "w");
    if (fp == nullptr) {
        fprintf(stderr, "failed to create file: %s\n", path.c_str());
        return false;
    }

    size_t nsat = prns.size();
    size_t nepo = length/interval + 1;

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

    // write header
    fprintf(fp, "     2.00           CLOCK DATA                              RINEX VERSION / TYPE\n");
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

    raw[0] = beta*dcb; // C1W
    raw[1] = -alpha*dcb; // C2W
    raw[2] = (-wl/g*dw + (1+1/g)*nl*dn)*1E9 - beta*dcb; // L1W
    raw[3] = (-wl*g*dw + (1 + g)*nl*dn)*1E9 + alpha*dcb ; // L2W
}

bool write_bias(const std::string &path, MJD t,
                const std::vector<std::string> &prns,
                const std::vector<double> &wl_bias)
{
    FILE *fp = fopen(path.c_str(), "w");
    if (fp == nullptr) {
        fprintf(stderr, "failed to create file: %s\n", path.c_str());
        return false;
    }

    // header
    fprintf(fp, "+BIAS/SOLUTION\n");
    fprintf(fp, "*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___ __ESTIMATED_SLOPE____ _STD_DEV___\n");

    const size_t nsat = prns.size();
    int y, doy;
    double f[2], raw[4];
    const char *types[4] = { "C1W", "C2W", "L1C", "L2W" };
    mjd2doy(t.d, &y, &doy);
    for (size_t i=0; i!=nsat; ++i) {
        sat_freq(prns[i], f);
        convert2raw_bias(wl_bias[i], f, raw);
        for (int j=0; j!=4; ++j)
            fprintf(fp, " OSB  %4s %3s %9s %-4s %-4s %4d:%03d:%05d %4d:%03d:%05d %-4s %21.3f %11.3f\n",
                    "SVN_", prns[i].c_str(), "", types[j], "", y, doy, 0, y, doy+1, 0, "ns", raw[j], 0.0);
    }

    fprintf(fp, "-BIAS/SOLUTION\n");
    fprintf(fp, "%%=ENDBIA\n");

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

void unit_weight(std::vector<double> &wgts)
{
    double sum = 0;
    for (auto it=wgts.begin(); it!=wgts.end(); ++it)
        sum += *it;
    for (auto it=wgts.begin(); it!=wgts.end(); ++it)
        *it /= sum;
}
