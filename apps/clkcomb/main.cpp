#include <stdio.h>
#include <math.h>

#include <io/io.h>
#include <numeric>

#include "config.h"
#include "utils.h"
#include "AnalyseCenter.h"

FILE *g_logfile = nullptr;

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

    AnalyseCenter combined(config.orb_ac);
    std::vector<AnalyseCenter> acs;
    acs.reserve(config.ac_names.size()); // avoid copy
    if (!init_acs(config, sats, acs, combined))
        return 1;

    const size_t nac_total  = acs.size();
    const size_t nsat_total = config.prns.size();
    const size_t nsys_total = config.syss.size();
    const size_t nepo_total = config.length/config.interval + 1;
    const size_t nsta_total = config.sta_list.size();
    const size_t nepo_total_sta = config.length/config.sta_interval + 1;

    fprintf(stdout, "    constellation: %4lu\n", nsys_total);
    fprintf(stdout, "    satellite    : %4lu\n", nsat_total);
    fprintf(stdout, "    epoch        : %4lu\n", nepo_total);
    fflush(stdout);

    std::string dif_file = replace_pattern(config.dif_pattern, config.mjd);
    FILE *diffile = fopen(dif_file.c_str(), "w");
    if (diffile == nullptr) {
        fprintf(stderr, MSG_ERR"unable to create: %s\n", dif_file.c_str());
        return 1;
    }

    std::string log_file = replace_pattern(config.log_pattern, config.mjd);
    g_logfile = fopen(log_file.c_str(), "w");
    if (g_logfile == nullptr) {
        fprintf(stderr, MSG_ERR"unable to create: %s\n", log_file.c_str());
        return 1;
    }

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
        combined.have_bias.resize(sats.size());
        combined.wl_bias.resize(sats.size());
        for (size_t i=0; i!=nsys_total; ++i)
            if (!align_widelane(sats, acs, combined, prns_per_system[i]))
                return 1;
    }

    // Adust sat_clks with NL
    if (config.phase_clock) {
        for (size_t i=0; i!=nsat_total; ++i) {
            double freq = sats[i].f[0] + sats[i].f[1]; // NL
            for (auto it=acs.begin(); it!=acs.end(); ++it) {
                if (!it->have_bias[i] || !combined.have_bias[i]) {
                    it->sat_clks[i].assign(nepo_total, None);
                    continue;
                }
                double val = it->nl_bias[i]*(1E9/freq); // ns
                // fprintf(g_logfile, "NL %3s %4s %8.3f %8.3f %8.3f\n", sats[i].prn.c_str(), it->name.c_str(), it->wl_bias[i], it->nl_bias[i], val);
                for (size_t epo=0; epo!=nepo_total; ++epo) {
                    if (it->sat_clks[i][epo] == None)
                        continue;
                    it->sat_clks[i][epo] -= val;
                }
            }
        }
    }

    fprintf(stdout, "(%s) %s\n", run_time(), "align clock datum..."); fflush(stdout);

    // Count epochs for each PRN
    std::vector<std::vector<size_t>> clk_counts(nac_total);
    std::vector<size_t> prn_counts(nsat_total);
    for (size_t iac=0; iac!=nac_total; ++iac) {
        count_satclks(config, acs[iac].sat_clks, clk_counts[iac]);
        for (size_t iprn=0; iprn!=nsat_total; ++iprn) {
            if (clk_counts[iac][iprn] == nepo_total)
                ++prn_counts[iprn];
        }
    }

    // output count table
    fprintf(g_logfile, "# Epoch count\nPRN");
    for (size_t iac=0; iac!=nac_total; ++iac)
        fprintf(g_logfile, " %5s", acs[iac].name.c_str());
    fprintf(g_logfile, " valid\n");
    for (size_t iprn=0; iprn!=nsat_total; ++iprn) {
        fprintf(g_logfile, "%3s", sats[iprn].prn.c_str());
        for (size_t iac=0; iac!=nac_total; ++iac) {
            fprintf(g_logfile, " %5lu", clk_counts[iac][iprn]);
        }
        fprintf(g_logfile, " %3lu\n", prn_counts[iprn]);
    }

    // Find common PRNs
    std::vector<std::vector<int>> common_prns(nsys_total);
    for (size_t isys=0; isys<nsys_total; ++isys)
    {
        int begin = prns_per_system[isys].front();
        int end   = prns_per_system[isys].back() + 1; // one past last
        size_t count_max = *std::max_element(&prn_counts[begin], &prn_counts[end]);
        if (count_max==0) {
            fprintf(stderr, MSG_ERR "bad [session]:length or [constellation]:system setting\n");
            return 1;
        }

        for (auto iprn=begin; iprn!=end; ++iprn) {
            if (prn_counts[iprn] == count_max)
                common_prns[isys].push_back(iprn);
        }
        fprintf(g_logfile, "\ncommon PRNs: %c %lu\n", config.strsys[isys], common_prns[isys].size());
    }

    // for (size_t i=0; i!=nac_total; ++i) {
    //     for (size_t isys=0; isys!=nsys_total; ++isys) {
    //         write_clock_datum(config.prns, config.strsys[isys], acs[i].name, acs[i].sat_clks, common_prns[isys]);
    //     }
    // }

    // Find datum ref ac
    size_t idatum = 0; // AC used for align_datum
    for (size_t iac=1; iac!=nac_total; ++iac) {
        size_t iprn = common_prns[0][0];
        if (clk_counts[iac][iprn] > clk_counts[idatum][iprn]
            && acs[iac].name!="gfz" && acs[iac].name!="esa") // nonlinear clk datum
            idatum = iac;
    }
    fprintf(g_logfile, "datum reference: %3s %4lu\n", acs[idatum].name.c_str(), clk_counts[idatum][common_prns[0][0]]);

    // Align clock datum with common PRNs
    for (size_t i=0; i!=nac_total; ++i) {
        if (i == idatum)
            continue;
        for (size_t isys=0; isys!=nsys_total; ++isys) {
            // remove_clock_datum(acs[i].sat_clks, prns_per_system[isys], common_prns[isys]);
            if (isys == 0 && config.combine_staclk)
                align_clock_datum2(acs[i].sat_clks, acs[idatum].sat_clks, prns_per_system[isys],
                        common_prns[isys], acs[i].sta_clks, config.sta_interval/config.interval);
            else
                align_clock_datum(acs[i].sat_clks, acs[idatum].sat_clks, prns_per_system[isys], common_prns[isys]);
        }
        // write_satclks(stdout, config, acs[i].name, acs[i].sat_clks);
    }

    //std::vector<double> datum_clks;
    //for (size_t iepo=0; iepo!=nepo_total; ++iepo) {
    //    double sum = 0;
    //    for (auto it=common_prns[0].cbegin(); it!=common_prns[0].cend(); ++it) {
    //        sum += acs[idatum].sat_clks[*it][iepo];
    //    }
    //    double mean = sum/common_prns[0].size();
    //    datum_clks.push_back(mean);
    //}

    // For compare IGS clock
    // for (size_t isys=0; isys!=nsys_total; ++isys) {
    //     if (isys==0 && config.combine_staclk)
    //         align_clock_datum2(combined.sat_clks, acs[0].sat_clks, prns_per_system[isys],
    //                 common_prns[isys], combined.sta_clks, config.sta_interval/config.interval);
    //     else
    //         align_clock_datum(combined.sat_clks, acs[0].sat_clks, prns_per_system[isys], common_prns[isys]);
    // }

    fprintf(stdout, "(%s) %s\n", run_time(), "construct reference clock..."); fflush(stdout);

    // Construct initial combined clock
    int refac;
    construct_initclk(config, acs, sats, combined.sat_clks, refac);
    if (config.combine_staclk)
        construct_init_staclk(config, acs, combined.sta_clks);

    // Allocate stds vec
    combined.sat_stds.assign(nsat_total, std::vector<double>(nepo_total));
    combined.sta_stds.assign(nsta_total, std::vector<double>(nepo_total_sta));

    // Backup clocks before removing bias
    for (size_t iac=0; iac!=nac_total; ++iac) {
        acs[iac].init_sat_clks = acs[iac].sat_clks;
        acs[iac].init_sta_clks = acs[iac].sta_clks;
    }

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
                        double ref = combined.sat_clks[*iprn][iepo];
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
                    fprintf(g_logfile, "%3s %3s %8.3f\n", acs[iac].name.c_str(), sats[*iprn].prn.c_str(), cycle);
                    bias.push_back(cycle);
                }

                double max = *std::max_element(bias.begin(), bias.end());
                for (auto it=bias.begin(); it!=bias.end(); ++it)
                    if (max - *it > 0.75)
                        *it += 1;
                double cycle = stable_mean(bias);
                // double cycle = std::accumulate(bias.begin(), bias.end(), 0.0)/bias.size();
                double mean = cycle*1E9/freq;
                fprintf(g_logfile, "%3s %c =>%8.3f\n", acs[iac].name.c_str(), config.strsys[isys], cycle);
                for (auto iprn=prns_per_system[isys].begin(); iprn!=prns_per_system[isys].end(); ++iprn)
                {
                    for (size_t iepo=0; iepo!=nepo_total; ++iepo)
                        if (acs[iac].sat_clks[*iprn][iepo] != None)
                            acs[iac].sat_clks[*iprn][iepo] -= mean;
                }
            }
        }
    }

    fprintf(stdout, "(%s) %s\n", run_time(), "iterative combination..."); fflush(stdout);

    bool weight_ac = config.weight_method == "ac" || config.combine_staclk;
    // Iteration
    std::vector<std::vector<double>> clk_wgts(nac_total, std::vector<double>(nsat_total));
    int niter = 0;
    const int max_niter = 3;
    while (niter < max_niter)
    {
        ++niter;
        fprintf(stdout, "(%s) iter #%1d...\n", run_time(), niter); fflush(stdout);
        fprintf(g_logfile, "\n** iter #%1d **\n", niter);

        // Remove clock bias between ACs & Determine weight
        std::vector<std::vector<double>> clk_stds(nac_total);
        for (size_t i=0; i!=nac_total; ++i) {
            double clk_std;
            for (size_t iprn=0; iprn!=nsat_total; ++iprn) {
                fprintf(g_logfile, "## %3s %3s  \n", acs[i].name.c_str(), config.prns[iprn].c_str());
                remove_clock_bias(acs[i].name, sats[iprn], niter, acs[i].sat_clks[iprn],
                                  combined.sat_clks[iprn], config.phase_clock, clk_std, niter!=1);
                // if (niter == 1 && clk_std != 0.0) {
                if (niter == 1) {
                    clk_std = 1;
                }
                clk_stds[i].push_back(clk_std);
                clk_wgts[i][iprn] = 1/clk_stds[i][iprn];
            }

            // Station clock
            if (config.combine_staclk) {
                for (size_t isit=0; isit!=nsta_total; ++isit) {
                    fprintf(g_logfile, "## %3s %4s  \n", acs[i].name.c_str(), config.sta_list[isit].c_str());
                    remove_clock_bias(acs[i].name, Satellite(config.sta_list[isit]), niter,
                                      acs[i].sta_clks[isit], combined.sta_clks[isit], false, clk_std, niter!=1);
                }
            }

            if (config.phase_clock) {
                for (size_t isys=0; isys!=nsys_total; ++isys)
                    for (size_t epo=0; epo!=nepo_total; ++epo) {
                        int n=0; double diff=0;
                        for (auto iprn=prns_per_system[isys].begin(); iprn!=prns_per_system[isys].end(); ++iprn) {
                            if (acs[i].sat_clks[*iprn][epo]!=None && combined.sat_clks[*iprn][epo]!=None) {
                                ++n;
                                diff += combined.sat_clks[*iprn][epo] - acs[i].sat_clks[*iprn][epo];
                            }
                        }
                        if (n!=0) diff/=n;
                        for (auto iprn=prns_per_system[isys].begin(); iprn!=prns_per_system[isys].end(); ++iprn)
                            if (acs[i].sat_clks[*iprn][epo]!=None)
                                acs[i].sat_clks[*iprn][epo] += diff;
                    }
            }

            if (niter == 1) {
                //write_satclks(stdout, config, acs[i].name, acs[i].sat_clks);
                write_clks_diff(diffile, config.prns, "satclk", acs[i].name, acs[i].sat_clks, combined.sat_clks);
                if(config.combine_staclk) {
                    write_clks_diff(diffile, config.sta_list, "staclk", acs[i].name, acs[i].sta_clks, combined.sta_clks);
                }
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
                if (!usable || sats[isat].prn[0] == 'R') // don't use GLONASS for AC weight
                    continue;
                for (size_t iac=0; iac!=nac_total; ++iac)
                    std_sums[iac] += clk_stds[iac][isat];
            }

            double wgt_sum = 0;
            double std_sum = std::accumulate(std_sums.begin(), std_sums.end(), 0.0);
            for (size_t iac=0; iac!=nac_total; ++iac) {
                double wgt = std_sum/std_sums[iac];
                wgt_sum += wgt;
                clk_wgts[iac].assign(nsat_total, wgt);
            }
            for (size_t iac=0; iac!=nac_total; ++iac) {
                fprintf(stdout, "    weight: %3s %6.2f%%\n", acs[iac].name.c_str(), 100*clk_wgts[iac][0]/wgt_sum);
                acs[iac].weight = clk_wgts[iac][0]/wgt_sum;
            }
        } else { // ouput satellite weight
            for (size_t iprn=0; iprn!=nsat_total; ++iprn) {
                fprintf(stdout, "    weight: %3s", config.prns[iprn].c_str());
                double wgt_sum = 0;
                for (size_t i=0; i!=nac_total; ++i)
                    wgt_sum += clk_wgts[i][iprn];

                for (size_t i=0; i!=nac_total; ++i) {
                    fprintf(stdout, " %6.2f", 100*clk_wgts[i][iprn]/wgt_sum);
                }
                fprintf(stdout, "\n");
            }
        }
        fflush(stdout);

        // Combination
        for (size_t iprn=0; iprn!=nsat_total; ++iprn)
        {
            for (size_t epo=0; epo!=nepo_total; ++epo)
            {
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
                combine_one_epoch(clks, wgts, combined.sat_clks[iprn][epo], combined.sat_stds[iprn][epo], deleted);
                //fprintf(g_logfile, "iter #%2d %4lu %3s %2lu %8.3f\n", niter, epo, config.prns[iprn].c_str(), clks.size(), 1E3*wrms);

                // remove outliers
                for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
                    fprintf(g_logfile, "del %2d %3s %4s %4lu\n", niter, acs[ac_indexs[*it]].name.c_str(), config.prns[iprn].c_str(), epo);
                    acs[ac_indexs[*it]].sat_clks[iprn][epo] = None;
                }
            }
        }

        // if (niter == 1)
            // write_satclks(stdout, config, "com", combined.sat_clks);

        if (!config.combine_staclk)
            continue;
        // Station clock
        for (size_t isit=0; isit!=nsta_total; ++isit)
        {
            for (size_t epo=0; epo!=nepo_total_sta; ++epo)
            {
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
                combine_one_epoch(clks, wgts, combined.sta_clks[isit][epo], combined.sta_stds[isit][epo], deleted);

                //fprintf(g_logfile, "iter #%2d %4lu %4s %2lu %8.3f\n", niter, epo, config.sta_list[isit].c_str(), clks.size(), 1E3*wrms);

                // remove outliers
                for (auto it=deleted.begin(); it!=deleted.end(); ++it) {
                    acs[ac_indexs[*it]].sta_clks[isit][epo] = None;
                    fprintf(g_logfile, "del %2d %3s %4s %4lu\n", niter, acs[ac_indexs[*it]].name.c_str(), config.sta_list[isit].c_str(), epo);
                }
            }
        }
    } // while

    // write_satclks(stdout, config, "com", combined.sat_clks);

    // Compare with comb_clks
    // for (size_t i=0; i!=nac_total; ++i)
    //     compare_satclks(config, acs[i].name, acs[i].sat_clks, combined.sat_clks, false);
    // if (config.combine_staclk)
    //     for (size_t i=0; i!=nac_total; ++i)
    //         compare_staclks(config, acs[i].name, acs[i].sta_clks, combined.sta_clks, false);

    // Compare with IGS clocks
    // compare_satclks(config, "xxx", combined.sat_clks, combined.sat_clks, true);

    // Output CLS file
    std::string cls_file = replace_pattern(config.cls_pattern, config.mjd, config.product_prefix);
    fprintf(stdout, "(%s) %s\n", run_time(), "write cls file..."); fflush(stdout);
    write_summary(cls_file, config, acs, sats, combined.sat_clks, combined.sta_clks);
    fprintf(stdout, "    -> %s\n", cls_file.c_str()); fflush(stdout);

    // align to broadcast ephemeris
    if (config.align_brdc) {
        fprintf(stdout, "(%s) %s\n", run_time(), "align to broadcast ephemeris..."); fflush(stdout);

        std::string nav_file = config.product_path + replace_pattern(config.nav_pattern, config.mjd);
        RinexNav rnxnav;
        if (!rnxnav.read(nav_file))
            return 1;

        fprintf(g_logfile, "\n# Align broadcast ephemeris\n");
        MJD t = config.mjd;
        std::vector<double> diffs;
        for (size_t iepo = 0; iepo != nepo_total; ++iepo) {
            size_t n = 0;
            double sum_diff = 0, clk;
            t.sod = iepo*config.interval;

            for (auto it = common_prns[0].cbegin(); it != common_prns[0].cend(); ++it) {
                const Ephemeris *eph = rnxnav.eph(t, config.prns[*it]);
                if (eph == nullptr) {
                    fprintf(stderr, MSG_WAR"no nav: %3s %5d %7.0f\n", config.prns[*it].c_str(), t.d, t.sod);
                    break;
                }
                eph->satClk(t, clk);

                if (combined.sat_clks[*it][iepo] != None) {
                    ++n;
                    sum_diff += clk*1E9 - combined.sat_clks[*it][iepo];
                }
            }
            double diff = sum_diff/n;
            diffs.push_back(n==common_prns[0].size() ? diff : None);
            if (n == common_prns[0].size())
                fprintf(g_logfile, "cmb-brdc %4lu %8.3f\n", iepo, diff);
        }

        double offset, drift;
        clkfit(diffs, nepo_total, offset, drift);
        fprintf(stdout, "    brdc fit: %12.3f(ps) %12.3f(ps/d)\n", offset*1E3, drift*1E3*86400/config.interval); fflush(stdout);
        fprintf(g_logfile,  "brdc fit: %12.3f(ps) %12.3f(ps/d)\n", offset*1E3, drift*1E3*86400/config.interval);

        // double xsum=0, ysum=0, xysum=0, x2sum=0;
        // for (size_t iepo=0; iepo!=nepo_total; ++iepo) {
        //     t.sod = iepo*config.interval;
        //     double sum = 0, clk;
        //     for (auto it=common_prns[0].cbegin(); it!=common_prns[0].cend(); ++it) {
        //         const Ephemeris *eph = rnxnav.eph(t, config.prns[*it]);
        //         if (eph == nullptr) {
        //             fprintf(stderr, "no nav: %3s %5d %7.0f\n", config.prns[*it].c_str(), t.d, t.sod);
        //             return 1;
        //         }
        //         eph->satClk(t, clk);
        //         sum += clk*1E9;
        //     }
        //     double mean = sum/common_prns[0].size();
        //     double diff = mean - datum_clks[iepo];
        //     xsum += iepo;       ysum += diff;
        //     xysum += iepo*diff; x2sum += iepo*iepo;
        // }
        // double a = (x2sum*ysum - xsum*xysum)/(x2sum*nepo_total - xsum*xsum);
        // double b = (nepo_total*xysum - xsum*ysum)/(nepo_total*x2sum - xsum*xsum);
        // fprintf(g_logfile, "brdc fit: %12.3f %12.3f\n", a*1E3, b*1E3*86400/config.interval);

        size_t sample_ratio = config.sta_interval/config.interval;
        for (size_t iepo=0; iepo!=nepo_total; ++iepo)
        {
            double corr = offset + drift*iepo;
            for (size_t isat=0; isat!=nsat_total; ++isat)
                if (combined.sat_clks[isat][iepo] != None)
                    combined.sat_clks[isat][iepo] += corr;

            if (!config.combine_staclk || iepo%sample_ratio != 0)
                continue;
            for (auto it=combined.sta_clks.begin(); it!=combined.sta_clks.end(); ++it)
                if (it->at(iepo/sample_ratio) != None)
                    it->at(iepo/sample_ratio) += corr;
        }
    }

    fprintf(stdout, "(%s) %s\n", run_time(), "write clk file..."); fflush(stdout);

    // ns => s
    for (size_t sat=0; sat!=nsat_total; ++sat)
        for (size_t epo=0; epo!=nepo_total; ++epo) {
            combined.sat_clks[sat][epo] /= 1E9;
            combined.sat_stds[sat][epo] /= 1E9;
        }

    if (config.combine_staclk)
        for (size_t sta=0; sta!=nsta_total; ++sta)
            for (size_t epo=0; epo!=nepo_total_sta; ++epo) {
                combined.sta_clks[sta][epo] /= 1E9;
                combined.sta_stds[sta][epo] /= 1E9;
            }

    std::string clk_file = replace_pattern(config.clk_pattern, config.mjd, config.product_prefix);
    write_clkfile(clk_file, config, combined);
    fprintf(stdout, "    -> %s\n", clk_file.c_str()); fflush(stdout);

    if (config.phase_clock) {
        std::string bia_file = replace_pattern(config.bia_pattern, config.mjd, config.product_prefix);
        write_bias(bia_file, config.mjd, sats, combined.have_bias, combined.wl_bias);

        // bia_file = replace_pattern(config.bia_pattern, config.mjd, "fip");
        // combined.nl_bias.resize(nsat_total);
        // write_fip(bia_file, config.mjd, sats, acs[0].wl_bias, combined.nl_bias);
    }

    fclose(diffile);
    fclose(g_logfile);
    return 0;
}
