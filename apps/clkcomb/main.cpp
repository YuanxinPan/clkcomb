#include <stdio.h>
#include <math.h>

#include <numeric>
#include <pppx/io.h>

#include "config.h"
#include "utils.h"
#include "AnalyseCenter.h"

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
    const size_t nsys_total = config.strsys.size();
    const size_t nepo_total = config.length/config.interval + 1;

    fprintf(stdout, "\nconstellation: %lu\n", nsys_total);
    fprintf(stdout, "satellite: %lu\n", nsat_total);
    fprintf(stdout, "epo_total: %-lu\n\n", nepo_total);

    if (config.phase_clock) {
        // std::string bia_file = replace_pattern(config.bia_pattern, config.mjd, "cod");
        // write_fip(bia_file, config.mjd, sats, acs[0].wl_bias, acs[0].nl_bias);

        // bia_file = replace_pattern(config.bia_pattern, config.mjd, "grg");
        // write_fip(bia_file, config.mjd, sats, acs[1].wl_bias, acs[1].nl_bias);

        // bia_file = replace_pattern(config.bia_pattern, config.mjd, "whu");
        // write_fip(bia_file, config.mjd, sats, acs[2].wl_bias, acs[2].nl_bias);
    }

    std::vector<std::vector<int>> prns_per_system(nsys_total);
    for (size_t i=0; i<nsat_total; ++i) {
        int isys = config.strsys.find(config.prns[i][0]);
        prns_per_system[isys].push_back(i);
    }

    std::vector<size_t> prn_counts(nsat_total);
    for (auto it=acs.begin(); it!=acs.end(); ++it) {
        it->read_clock(config.mjd, config.length, config.interval, config.prns, combined_ac.rnxsp3());

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

    if (config.phase_clock) {
        for (size_t i=0; i!=nsys_total; ++i)
            if (!align_widelane(sats, acs, prns_per_system[i]))
                return false;
    }

    // Adust sat_clks with NL
    if (config.phase_clock) {
        for (size_t i=0; i!=nsat_total; ++i) {
            double freq = sats[i].f[0] + sats[i].f[1]; // NL
            for (auto it=acs.begin(); it!=acs.end(); ++it) {
                if (!it->have_bias[i]) {
                    it->sat_clks[i].assign(nepo_total, None);
                    continue;
                }
                double val = it->nl_bias[i]*(1E9/freq); // ns
                val = it->name=="grg" ? val : -val;
                for (size_t epo=0; epo!=nepo_total; ++epo) {
                    if (it->sat_clks[i][epo] == None)
                        continue;
                    it->sat_clks[i][epo] += val;
                }
            }
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
    for (size_t i=1; i!=nac_total; ++i) {
        for (size_t isys=0; isys!=nsys_total; ++isys) {
            // remove_clock_datum(acs[i].sat_clks, prns_per_system[isys], common_prns[isys]);
            align_clock_datum(acs[i].sat_clks, acs[0].sat_clks, prns_per_system[isys], common_prns[isys]);
        }
        // write_satclks(stdout, config, acs[i].name, acs[i].sat_clks);
    }

    // for (size_t isys=0; isys!=nsys_total; ++isys)
    //     align_clock_datum(combined_ac.sat_clks, acs[0].sat_clks, prns_per_system[isys], common_prns[isys]);

    // Construct initial combined clock
    std::vector<std::vector<double>> comb_clks;
    construct_initclk(config, acs, sats, comb_clks);

    // for (size_t iac=0; iac!=nac_total; ++iac)
    //     for (size_t isys=0; isys!=nsys_total; ++isys)
    //         for (auto it=prns_per_system[isys].begin()+1; it!=prns_per_system[isys].end(); ++it)
    //             for (size_t iepo=0; iepo!=nepo_total; ++iepo)
    //             {
    //                 double diff = acs[iac].sat_clks[*it][iepo] - acs[iac].sat_clks[0][iepo] - (comb_clks[*it][iepo] - comb_clks[0][iepo]);
    //                 if (comb_clks[*it][iepo]!=None && acs[iac].sat_clks[*it][iepo]!=None)
    //                     fprintf(stdout, "diff %4lu %3s %3s %16.3f\n", iepo, acs[iac].name.c_str(), sats[*it].prn.c_str(), diff);
    //             }

    for (size_t i=0; i!=nac_total; ++i)
        write_satclks_diff(stdout, config, acs[i].name, acs[i].sat_clks, comb_clks);

    bool weight_ac = config.weight_method == "ac";
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
                // remove_clock_bias(acs[i].sat_clks[iprn], comb_clks[iprn], clk_std, true);
                // if (niter == 1 && clk_std != 0.0) {
                if (niter == 1) {
                    clk_std = 1;
                }
                // std_sum += clk_std;
                clk_stds[i].push_back(clk_std);
                clk_wgts[i][iprn] = 1/clk_stds[i][iprn];
            }

            // if (config.phase_clock) {
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

                    // fprintf(stdout, "%3s %7lu %8.1f\n", acs[i].name.c_str(), epo, 1E3*diff);
                }
            // }

            // if (niter==1) write_satclks(stdout, config, acs[i].name, acs[i].sat_clks);
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
    } // while

    // write_satclks(stdout, config, "com", comb_clks);

    // Compare with IGS clocks
    // compare_satclks(config, comb_clks, combined_ac.sat_clks, false);

    for (size_t i=0; i!=nac_total; ++i)
        compare_satclks(config, acs[i].name, acs[i].sat_clks, comb_clks, false);

    // ns => s
    for (size_t sat=0; sat!=nsat_total; ++sat)
        for (size_t epo=0; epo!=nepo_total; ++epo)
            comb_clks[sat][epo] /= 1E9;

    // align to broadcast ephemeris
    std::string nav_file = config.product_path + replace_pattern(config.nav_pattern, config.mjd);
    RinexNav rnxnav;

    if (!rnxnav.read(nav_file))
        return 1;

    // 1. calculate ISB
    // std::vector<double> isb_bias;
    // std::vector<double> clk_means(nsys_total);
    // for (size_t sys=0; sys!=nsys_total; ++sys)
    // {
    //     double clk;
    //     const Ephemeris *eph;
    //     int n = 0;
    //     for (auto it=common_prns[sys].begin(); it!=common_prns[sys].end(); ++it)
    //     {
    //         if (comb_clks[*it][0] == None)
    //             continue;
    //         ++n;
    //         eph = rnxnav.eph(config.mjd, config.prns[*it]);
    //         eph->satClk(config.mjd, clk);
    //         clk_means[sys] += clk;
    //         // fprintf(stdout, "%3s %12.5E %12.5E %12.5E\n", config.prns[*it].c_str(),
    //         //         clk, comb_clks[*it][0]*1E-9, clk-comb_clks[*it][0]*1E-9);
    //     }
    //     clk_means[sys] /= n;
    //     isb_bias.push_back(clk_means[sys] - clk_means[0]);
    // }

    // 2. align
    MJD t = config.mjd;
    for (size_t epo=0; epo!=nepo_total; ++epo, t+=config.interval)
    {
        const Ephemeris *eph;
        double clk, datum=0, clk_sum[2] = { 0, 0 }; // nav, combined
        for (size_t sys=0; sys!=nsys_total; ++sys)
        {
            if (sys == 0) {
                int n = 0;
                for (auto it=common_prns[sys].begin(); it!=common_prns[sys].end(); ++it)
                {
                    if (comb_clks[*it][epo] == None)
                        continue;
                    eph = rnxnav.eph(t, config.prns[*it]);
                    if (eph == nullptr) {
                        fprintf(stderr, "no nav: %3s %4lu\n", config.prns[*it].c_str(), epo);
                        continue;
                    }
                    ++n;
                    eph->satClk(t, clk);
                    clk_sum[0] += clk;
                    clk_sum[1] += comb_clks[*it][epo];
                }
                datum = (clk_sum[0] - clk_sum[1])/n;
            }

            for (auto it=prns_per_system[sys].begin(); it!=prns_per_system[sys].end(); ++it)
            {
                if (comb_clks[*it][epo] != None)
                    comb_clks[*it][epo] += datum; // + isb_bias[sys];
            }
        }
    }

    std::string clk_file = replace_pattern(config.clk_pattern, config.mjd, "xxx");
    write_clkfile(clk_file, config.mjd, config.length, config.interval, config.prns, comb_clks);

    if (config.phase_clock) {
        std::string bia_file = replace_pattern(config.bia_pattern, config.mjd, "xxx");
        write_bias(bia_file, config.mjd, sats, acs[0].wl_bias);

        // bia_file = replace_pattern(config.bia_pattern, config.mjd, "fip");
        // combined_ac.nl_bias.resize(nsat_total);
        // write_fip(bia_file, config.mjd, sats, acs[0].wl_bias, combined_ac.nl_bias);
    }

    return 0;
}
