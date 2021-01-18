#include <stdio.h>
#include <math.h>

#include <numeric>
#include <pppx/io.h>

#include "config.h"
#include "utils.h"
#include "AnalyseCenter.h"

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

bool write_fip(const std::string &path, MJD t,
               const std::vector<Satellite> &sats,
               const std::vector<double> &wl_bias,
               const std::vector<double> &nl_bias)
{
    FILE *fp = fopen(path.c_str(), "w");
    if (fp == nullptr) {
        fprintf(stderr, "failed to create file: %s\n", path.c_str());
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
            fprintf(stdout, "raw %3s %3s %6lu %14.3f\n", acn.c_str(), prns[*it].c_str(), epo, sat_clks[*it][epo]);
        }
        double mean = sum/ncommon;
        fprintf(stdout, "datum %3s %c %6lu %14.3f\n", acn.c_str(), sys, epo, mean);
    }
}

int main(int argc, char *argv[])
{
    if (argc > 2) {
        printf("usage clkcomb [comb.ini]\n");
        return 1;
    }

    std::string config_path = "comb.ini";
    if (argc == 2)
        config_path.assign(argv[1]);
    config_t config;
    if (!read_config(config_path, config)) {
        return 1;
    } else {
        print_config(stdout, config);
    }

    std::vector<Satellite> sats;
    if (!init_sats(config, sats))
        return 1;

    AnalyseCenter combined_ac(config.orb_ac);
    std::vector<AnalyseCenter> acs;
    acs.reserve(config.ac_names.size()); // avoid copy
    if (!init_acs(config, sats, acs, combined_ac))
        return 1;

    const size_t nac_total  = acs.size();
    const size_t nsat_total = config.prns.size();
    const size_t nsys_total = config.syss.size();
    const size_t nepo_total = config.length/config.interval + 1;
    const size_t nsta_total = config.sta_list.size();
    const size_t nepo_total_sta = config.length/config.sta_interval + 1;

    fprintf(stdout, "\nconstellation: %lu\n", nsys_total);
    fprintf(stdout, "satellite: %lu\n", nsat_total);
    fprintf(stdout, "epo_total: %-lu\n\n", nepo_total);

    // if (config.phase_clock) {
    //     for (auto it=acs.begin(); it!=acs.end(); ++it)
    //     {
    //         std::string bia_file = replace_pattern(config.bia_pattern, config.mjd, it->name);
    //         write_fip(bia_file, config.mjd, sats, it->wl_bias, it->nl_bias);
    //     }
    // }

    std::vector<std::vector<int>> prns_per_system(nsys_total);
    for (size_t i=0; i<nsat_total; ++i) {
        enum GNSS_Tp tp = prn2sys(config.prns[i]);
        int isys = std::find(config.syss.begin(), config.syss.end(), tp) - config.syss.begin();
        prns_per_system[isys].push_back(i);
    }

    if (config.phase_clock) {
        combined_ac.have_bias.resize(sats.size());
        combined_ac.wl_bias.resize(sats.size());
        for (size_t i=0; i!=nsys_total; ++i)
            if (!align_widelane(sats, acs, combined_ac, prns_per_system[i]))
                return 1;
    }

    // Adust sat_clks with NL
    if (config.phase_clock) {
        for (size_t i=0; i!=nsat_total; ++i) {
            double freq = sats[i].f[0] + sats[i].f[1]; // NL
            for (auto it=acs.begin(); it!=acs.end(); ++it) {
                if (!it->have_bias[i] || !combined_ac.have_bias[i]) {
                    it->sat_clks[i].assign(nepo_total, None);
                    continue;
                }
                double val = it->nl_bias[i]*(1E9/freq); // ns
                val = it->name=="grg" ? val : -val;
                // fprintf(stdout, "NL %3s %4s %8.3f %8.3f %8.3f\n", sats[i].prn.c_str(), it->name.c_str(), it->wl_bias[i], it->nl_bias[i], val);
                for (size_t epo=0; epo!=nepo_total; ++epo) {
                    if (it->sat_clks[i][epo] == None)
                        continue;
                    it->sat_clks[i][epo] += val;
                }
            }
        }
    }

    // Find common PRNs
    std::vector<size_t> prn_counts(nsat_total);
    for (auto it=acs.begin(); it!=acs.end(); ++it) {
        fprintf(stdout, "\n## %3s\n", it->name.c_str());
        std::vector<size_t> clk_count;
        count_satclks(config, it->sat_clks, clk_count);
        for (size_t iprn=0; iprn!=nsat_total; ++iprn)
        {
            fprintf(stdout, "%3s %4lu\n", sats[iprn].prn.c_str(), clk_count[iprn]);
            if (clk_count[iprn] == nepo_total)
                ++prn_counts[iprn];
        }
    }

    // Align clock datum with common PRNs
    std::vector<std::vector<int>> common_prns(nsys_total);
    for (size_t isys=0; isys<nsys_total; ++isys)
    {
        int begin = prns_per_system[isys].front();
        int end   = prns_per_system[isys].back() + 1; // one past last
        size_t count_max = *std::max_element(&prn_counts[begin], &prn_counts[end]);
        if (count_max==0) {
            fprintf(stdout, MSG_ERR "bad [session]:length or [constellation]:system setting\n");
            return 1;
        }

        for (auto iprn=begin; iprn!=end; ++iprn) {
            if (prn_counts[iprn] == count_max)
                common_prns[isys].push_back(iprn);

            fprintf(stdout, "%3s %2lu\n", config.prns[iprn].c_str(), prn_counts[iprn]);
        }
        fprintf(stdout, "common_prns: %c %lu\n", config.strsys[isys], common_prns[isys].size());
    }

    // for (size_t i=0; i!=nac_total; ++i) {
    //     for (size_t isys=0; isys!=nsys_total; ++isys) {
    //         write_clock_datum(config.prns, config.strsys[isys], acs[i].name, acs[i].sat_clks, common_prns[isys]);
    //     }
    // }

    // for (size_t i=0; i!=nac_total; ++i) {
    for (size_t i=1; i!=nac_total; ++i) {
        for (size_t isys=0; isys!=nsys_total; ++isys) {
            // remove_clock_datum(acs[i].sat_clks, prns_per_system[isys], common_prns[isys]);
            if (isys == 0 && config.combine_staclk)
                align_clock_datum2(acs[i].sat_clks, acs[0].sat_clks, prns_per_system[isys],
                        common_prns[isys], acs[i].sta_clks, config.sta_interval/config.interval);
            else
                align_clock_datum(acs[i].sat_clks, acs[0].sat_clks, prns_per_system[isys], common_prns[isys]);
        }
        // write_satclks(stdout, config, acs[i].name, acs[i].sat_clks);
    }

    std::vector<double> datum_clks;
    for (size_t iepo=0; iepo!=nepo_total; ++iepo) {
        double sum = 0;
        for (auto it=common_prns[0].cbegin(); it!=common_prns[0].cend(); ++it) {
            sum += acs[0].sat_clks[*it][iepo];
        }
        double mean = sum/common_prns[0].size();
        datum_clks.push_back(mean);
    }

    // Output aligned station clks
    // for (size_t iac=0; iac!=acs.size(); ++iac)
    //     for (size_t ista=0; ista!=config.sta_list.size(); ++ista)
    //         for (size_t epo=0, total=config.length/config.sta_interval+1; epo!=total; ++epo)
    //         {
    //             if (acs[iac].sta_clks[ista][epo] == None)
    //                 continue;
    //             fprintf(stdout, "%3s %4s %6lu %12.3f\n", acs[iac].name.c_str(),
    //                     config.sta_list[ista].c_str(), epo, acs[iac].sta_clks[ista][epo]);
    //         }

    // for (size_t isys=0; isys!=nsys_total; ++isys) {
    //     if (isys==0 && config.combine_staclk)
    //         align_clock_datum2(combined_ac.sat_clks, acs[0].sat_clks, prns_per_system[isys],
    //                 common_prns[isys], combined_ac.sta_clks, config.sta_interval/config.interval);
    //     else
    //         align_clock_datum(combined_ac.sat_clks, acs[0].sat_clks, prns_per_system[isys], common_prns[isys]);
    // }

    // Construct initial combined clock
    std::vector<std::vector<double>> comb_clks;
    std::vector<std::vector<double>> comb_staclks;
    construct_initclk(config, acs, sats, comb_clks);
    if (config.combine_staclk)
        construct_init_staclk(config, acs, comb_staclks);
    std::vector<std::vector<double>> comb_clks_backup = comb_clks;;

    // Remove fractional clock difference for each constellation
    if (config.phase_clock) {
        for (size_t iac=0; iac!=nac_total; ++iac) {
            for (size_t isys=0; isys!=nsys_total; ++isys)
            {
                std::vector<double> bias;
                int syssat = prns_per_system[isys].front();
                double freq = sats[syssat].f[0] + sats[syssat].f[1]; // NL
                for (auto iprn=common_prns[isys].begin(); iprn!=common_prns[isys].end(); ++iprn)
                {
                    std::vector<double> diff;
                    for (size_t iepo=0; iepo!=nepo_total; ++iepo)
                    {
                        double src = acs[iac].sat_clks[*iprn][iepo];
                        double ref = comb_clks[*iprn][iepo];
                        if (src==None || ref==None)
                            continue;
                        diff.push_back(src - ref);
                    }
                    if (diff.size() == 0u)
                        continue;

                    double mean = stable_mean(diff);
                    // double mean = std::accumulate(diff.begin(), diff.end(), 0.0)/diff.size();
                    double cycle = mean/1E9*freq;
                    cycle = cycle - int(cycle);
                    if (cycle < 0.0) cycle += 1;
                    fprintf(stdout, "%3s %3s %8.3f\n", acs[iac].name.c_str(), sats[*iprn].prn.c_str(), cycle);
                    bias.push_back(cycle);
                }

                double max = *std::max_element(bias.begin(), bias.end());
                for (auto it=bias.begin(); it!=bias.end(); ++it)
                    if (max - *it > 0.6)
                        *it += 1;
                double cycle = stable_mean(bias);
                // double cycle = std::accumulate(bias.begin(), bias.end(), 0.0)/bias.size();
                double mean = cycle*1E9/freq;
                fprintf(stdout, "%3s %c =>%8.3f\n", acs[iac].name.c_str(), config.strsys[isys], cycle);
                for (auto iprn=prns_per_system[isys].begin(); iprn!=prns_per_system[isys].end(); ++iprn)
                {
                    for (size_t iepo=0; iepo!=nepo_total; ++iepo)
                        if (acs[iac].sat_clks[*iprn][iepo] != None)
                            acs[iac].sat_clks[*iprn][iepo] -= mean;
                }
            }
        }
    }

    fprintf(stderr, "Write clock diff...\n");
    for (size_t i=0; i!=nac_total; ++i)
        write_satclks_diff(stdout, config, "1st", acs[i].name, acs[i].sat_clks, comb_clks);

    bool weight_ac = config.weight_method == "ac" || config.combine_staclk;
    // Iteration
    std::vector<std::vector<double>> clk_wgts(nac_total, std::vector<double>(nsat_total));
    int niter = 0;
    const int max_niter = 3;
    while (niter < max_niter)
    {
        ++niter;
        fprintf(stdout, "** iter #%2d **\n", niter);

        // Remove clock bias between ACs & Determine weight
        std::vector<std::vector<double>> clk_stds(nac_total);
        for (size_t i=0; i!=nac_total; ++i) {
            double clk_std; // std_sum = 0;
            for (size_t iprn=0; iprn!=nsat_total; ++iprn) {
                fprintf(stdout, "## %3s %3s  \n", acs[i].name.c_str(), config.prns[iprn].c_str());
                remove_clock_bias(acs[i].name, sats[iprn], acs[i].sat_clks[iprn], comb_clks[iprn],
                                  config.phase_clock, clk_std, niter!=1);
                // if (niter == 1 && clk_std != 0.0) {
                if (niter == 1) {
                    clk_std = 1;
                }
                // std_sum += clk_std;

                clk_stds[i].push_back(clk_std);
                clk_wgts[i][iprn] = 1/clk_stds[i][iprn];
            }
            if (niter == 1) write_satclks_diff(stdout, config, "2nd", acs[i].name, acs[i].sat_clks, comb_clks);

            if (config.phase_clock) {
                for (size_t isys=0; isys!=nsys_total; ++isys)
                    for (size_t epo=0; epo!=nepo_total; ++epo) {
                        int n=0; double diff=0;
                        for (auto iprn=prns_per_system[isys].begin(); iprn!=prns_per_system[isys].end(); ++iprn) {
                            if (acs[i].sat_clks[*iprn][epo]!=None && comb_clks[*iprn][epo]!=None) {
                                ++n;
                                diff += comb_clks[*iprn][epo] - acs[i].sat_clks[*iprn][epo];
                            }
                        }
                        if (n!=0) diff/=n;
                        for (auto iprn=prns_per_system[isys].begin(); iprn!=prns_per_system[isys].end(); ++iprn)
                            if (acs[i].sat_clks[*iprn][epo]!=None)
                                acs[i].sat_clks[*iprn][epo] += diff;
                    }
                if (niter == 1) write_satclks_diff(stdout, config, "3rd", acs[i].name, acs[i].sat_clks, comb_clks);
            }

            // if (niter==1) write_satclks(stdout, config, acs[i].name, acs[i].sat_clks);

            // Station clock
            for (size_t isit=0; isit!=nsta_total; ++isit) {
                fprintf(stdout, "## %3s %4s  \n", acs[i].name.c_str(), config.sta_list[isit].c_str());
                remove_clock_bias(acs[i].name, sats[0], acs[i].sta_clks[isit], comb_staclks[isit],
                        false, clk_std, niter!=1);
            }
        }

        // AC weight
        if (weight_ac) {
            std::vector<double> std_sums(nac_total);
            for (size_t isat=0; isat!=nsat_total; ++isat)
            {
                bool usable = true;
                for (size_t iac=0; iac!=nac_total; ++iac) {
                    if (clk_stds[iac][isat] == 99.9) {
                        usable = false;
                        break;
                    }
                }
                if (!usable)
                    continue;
                for (size_t iac=0; iac!=nac_total; ++iac)
                    std_sums[iac] += clk_stds[iac][isat];
            }

            double std_sum = std::accumulate(std_sums.begin(), std_sums.end(), 0.0);
            for (size_t iac=0; iac!=nac_total; ++iac) {
                double wgt = std_sum/std_sums[iac];
                clk_wgts[iac].assign(nsat_total, wgt);
                // fprintf(stdout, "AC weight: %3s %8.3f\n", acs[iac].name.c_str(), wgt);
            }
        }

        // ouput weight
        for (size_t iprn=0; iprn!=nsat_total; ++iprn) {
            fprintf(stdout, "weight  %3s ", config.prns[iprn].c_str());
            double wgt_sum = 0;
            for (size_t i=0; i!=nac_total; ++i)
                wgt_sum += clk_wgts[i][iprn];

            for (size_t i=0; i!=nac_total; ++i) {
                fprintf(stdout, " %8.3f", clk_wgts[i][iprn]/wgt_sum);
            }
            fprintf(stdout, "\n");
        }

        // Combination
        for (size_t iprn=0; iprn!=nsat_total; ++iprn)
        {
            for (size_t epo=0; epo!=nepo_total; ++epo)
            {
                double wrms;
                std::vector<size_t> ac_indexs;
                std::vector<int> deleted;
                std::vector<double> clks, wgts;
                for (size_t i=0; i!=nac_total; ++i)
                {
                    if (acs[i].sat_clks[iprn][epo] == None)
                        continue;
                    ac_indexs.push_back(i);
                    clks.push_back(acs[i].sat_clks[iprn][epo]);
                    wgts.push_back(clk_wgts[i][iprn]);
                }
                combine_one_epoch(clks, wgts, comb_clks[iprn][epo], wrms, deleted);

                fprintf(stdout, "iter #%2d %4lu %3s %2lu %8.3f\n", niter, epo, config.prns[iprn].c_str(), clks.size(), 1E3*wrms);

                // remove outliers
                if (!deleted.empty())
                    fprintf(stdout, "del %3s %4lu ", config.prns[iprn].c_str(), epo);
                for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
                    acs[ac_indexs[*it]].sat_clks[iprn][epo] = None;
                    fprintf(stdout, " %3s", acs[ac_indexs[*it]].name.c_str());
                }
                if (!deleted.empty()) fprintf(stdout, "\n");
            }
        }
        // if (niter == 1)
            // write_satclks(stdout, config, "com", comb_clks);
        // Debug
        if (niter == 1) {
            write_satclks_diff(stdout, config, "3rd", "xxx", comb_clks, comb_clks_backup);
        }

        if (!config.combine_staclk)
            continue;
        // Station clock
        for (size_t isit=0; isit!=nsta_total; ++isit)
        {
            for (size_t epo=0; epo!=nepo_total_sta; ++epo)
            {
                double wrms;
                std::vector<size_t> ac_indexs;
                std::vector<int> deleted;
                std::vector<double> clks, wgts;
                for (size_t i=0; i!=nac_total; ++i)
                {
                    if (acs[i].sta_clks[isit][epo] == None)
                        continue;
                    ac_indexs.push_back(i);
                    clks.push_back(acs[i].sta_clks[isit][epo]);
                    wgts.push_back(clk_wgts[i][0]);
                }
                combine_one_epoch(clks, wgts, comb_staclks[isit][epo], wrms, deleted);

                fprintf(stdout, "iter #%2d %4lu %4s %2lu %8.3f\n", niter, epo, config.sta_list[isit].c_str(), clks.size(), 1E3*wrms);

                // remove outliers
                if (!deleted.empty())
                    fprintf(stdout, "del %3s %4lu ", config.sta_list[isit].c_str(), epo);
                for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
                    acs[ac_indexs[*it]].sta_clks[isit][epo] = None;
                    fprintf(stdout, " %3s", acs[ac_indexs[*it]].name.c_str());
                }
                if (!deleted.empty()) fprintf(stdout, "\n");
            }
        }

    } // while

    // write_satclks(stdout, config, "com", comb_clks);

    // Compare with IGS clocks
    if (config.combine_staclk)
        for (size_t i=0; i!=nac_total; ++i)
            compare_staclks(config, acs[i].name, acs[i].sta_clks, comb_staclks, false);

    for (size_t i=0; i!=nac_total; ++i)
        compare_satclks(config, acs[i].name, acs[i].sat_clks, comb_clks, false);

    // compare_satclks(config, "xxx", comb_clks, combined_ac.sat_clks, true);

    // align to broadcast ephemeris
    if (config.align_brdc) {
        std::string nav_file = config.product_path + replace_pattern(config.nav_pattern, config.mjd);
        RinexNav rnxnav;
        if (!rnxnav.read(nav_file))
            return 1;

        MJD t = config.mjd;
        double xsum=0, ysum=0, xysum=0, x2sum=0;
        for (size_t iepo=0; iepo!=nepo_total; ++iepo) {
            t.sod = iepo*config.interval;
            double sum = 0, clk;
            for (auto it=common_prns[0].cbegin(); it!=common_prns[0].cend(); ++it) {
                const Ephemeris *eph = rnxnav.eph(t, config.prns[*it]);
                if (eph == nullptr) {
                    fprintf(stderr, "no nav: %3s %5d %7.0f\n", config.prns[*it].c_str(), t.d, t.sod);
                    return 1;
                }
                eph->satClk(t, clk);
                sum += clk*1E9;
            }
            double mean = sum/common_prns[0].size();
            double diff = mean - datum_clks[iepo];
            xsum += iepo;       ysum += diff;
            xysum += iepo*diff; x2sum += iepo*iepo;
        }
        double a = (x2sum*ysum - xsum*xysum)/(x2sum*nepo_total - xsum*xsum);
        double b = (nepo_total*xysum - xsum*ysum)/(nepo_total*x2sum - xsum*xsum);
        fprintf(stdout, "brdc fit: %12.3f %12.3f\n", a*1E3, b*1E3*86400/config.interval);

        size_t sample_ratio = config.sta_interval/config.interval;
        for (size_t iepo=0; iepo!=nepo_total; ++iepo)
        {
            double corr = a + b*iepo;
            for (size_t isat=0; isat!=nsat_total; ++isat)
                if (comb_clks[isat][iepo] != None)
                    comb_clks[isat][iepo] += corr;

            if (!config.combine_staclk || iepo%sample_ratio != 0)
                continue;
            for (auto it=comb_staclks.begin(); it!=comb_staclks.end(); ++it)
                if (it->at(iepo/sample_ratio) != None)
                    it->at(iepo/sample_ratio) += corr;
        }
    }

    // ns => s
    for (size_t sat=0; sat!=nsat_total; ++sat)
        for (size_t epo=0; epo!=nepo_total; ++epo)
            comb_clks[sat][epo] /= 1E9;

    if (config.combine_staclk)
        for (size_t sta=0; sta!=nsta_total; ++sta)
            for (size_t epo=0; epo!=nepo_total_sta; ++epo)
                comb_staclks[sta][epo] /= 1E9;

    std::string clk_file = replace_pattern(config.clk_pattern, config.mjd, "xxx");
    write_clkfile(clk_file, config, comb_clks, comb_staclks, combined_ac.rnxsnx());

    if (config.phase_clock) {
        std::string bia_file = replace_pattern(config.bia_pattern, config.mjd, "xxx");
        write_bias(bia_file, config.mjd, sats, combined_ac.have_bias, combined_ac.wl_bias);

        // bia_file = replace_pattern(config.bia_pattern, config.mjd, "fip");
        // combined_ac.nl_bias.resize(nsat_total);
        // write_fip(bia_file, config.mjd, sats, acs[0].wl_bias, combined_ac.nl_bias);
    }

    return 0;
}
