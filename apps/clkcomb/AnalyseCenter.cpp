#include "AnalyseCenter.h"
#include "utils.h"

#include <pppx/const.h>
#include <pppx/io.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <algorithm>
#include <array>
#include <cctype>

bool AnalyseCenter::read_clock(const config_t &config, const RinexSp3 &refsp3,
                               const RinexAtx &refatx, const RinexAtt &refatt)
{
    MJD t = config.mjd;
    const int length = config.length;
    const int interval = config.interval;
    const std::vector<std::string> &prns = config.prns;
    const size_t nsat_total = prns.size();
    const size_t nepo_total = length/interval + 1;

    sat_clks.clear();
    sat_clks.resize(nsat_total);

    // clock & orbit
    double clk, corr;
    CartCoor pos, ref;

    // Atx
    const RinexAtx::atx_t *atx_src = nullptr;
    const RinexAtx::atx_t *atx_dst = nullptr;
    double psrc[3], pdst[3];
    std::string f1, f2;
    std::vector<const RinexAtx::atx_t *> src_atxs;
    std::vector<const RinexAtx::atx_t *> dst_atxs;
    for (size_t isat=0; isat!=nsat_total; ++isat) {
        src_atxs.push_back(rnxatx_.atx(t, prns[isat]));
        dst_atxs.push_back(refatx.atx(t, prns[isat]));
    }

    // Att
    double qr[4], qs[4];
    double xsun[3], f[2];
    std::vector<double> diffs[2];
    diffs[0].resize(prns.size());
    diffs[1].resize(prns.size());

    for (size_t epo=0; epo!=nepo_total; ++epo, t+=interval) {
        for (size_t isat=0; isat!=nsat_total; ++isat) {
            if (!rnxclk_.clkBias(t, prns[isat], clk)) {
                // fprintf(stderr, "%3s %5lu %3s: no clk\n", name.c_str(), epo, prns[isat].c_str());
                clk = None;
            }

            // Antenna correction
            atx_src = src_atxs[isat];
            atx_dst = dst_atxs[isat];
            if (atx_src == nullptr || atx_dst == nullptr)
                clk = None;
            else if (clk != None) {
                if (prns[isat][0] == 'G') {
                    f1 = "G01"; f2 = "G02";
                } else if (prns[isat][0] == 'E') {
                    f1 = "E01"; f2 = "E05";
                } else if (prns[isat][0] == 'R') {
                    f1 = "R01"; f2 = "R02";
                } else if (prns[isat][0] == 'C') {
                    f1 = "C02"; f2 = "C06";
                }
                atx_src->pco(psrc, f1, f2);
                atx_dst->pco(pdst, f1, f2);
                clk += (psrc[2] - pdst[2])/LightSpeed;
            }

            // Orbit correction
            if (!rnxsp3_.satPos(t, prns[isat], pos) ||
                !refsp3.satPos( t, prns[isat], ref)) {
                // fprintf(stderr, "%3s %5lu %3s: no orb\n", name.c_str(), epo, prns[isat].c_str());
                clk = None;
            } else if (clk != None) {
                corr = (pos - ref).dot(pos)/pos.norm();
                clk -= corr/LightSpeed;
            }

            // Attitude correction
            if (config.use_att) {
                if (!rnxatt_.sat_att(t, prns[isat], qs) || (!refatt.empty() && !refatt.sat_att(t, prns[isat], qr))) {
                    clk = None;  // try remove this, issue when read igs clock without att
                } else {
                    if (refatt.empty()) { // use nominal attitude
                        SunPosition(t.d, t.sod, xsun);
                        nominal_att(ref.data(), xsun, qr);
                    }
                    diffs[1][isat] = qAngularDist(qs, qr);
                    while (diffs[1][isat] - diffs[0][isat] > 180)
                        diffs[1][isat] -= 360;
                    while (diffs[1][isat] - diffs[0][isat] < -180)
                        diffs[1][isat] += 360;

                    sat_freq(prns[isat], f);
                    double pwu = diffs[1][isat]/360/(f[0] + f[1]); // phase windup correction
                    if (clk != None)
                        clk += pwu;
                    diffs[0][isat] = diffs[1][isat];
                }
            }

            sat_clks[isat].push_back(clk*1E9); // sec => ns
        }
    }

    rnxclk_.close();
    rnxatx_.close();
    if (config.use_att) rnxatt_.close();
    return true;
}

bool AnalyseCenter::read_staclk(MJD t, int length, int interval,
                                const std::vector<std::string> &sta_list,
                                const RinexSnx &refsnx)
{
    const size_t nsta_total = sta_list.size();
    const size_t nepo_total = length/interval + 1;

    FILE *fp = fopen(clk_file.c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR "no such file: %s\n", clk_file.c_str());
        return false;
    }

    sta_clks.clear();
    sta_clks.resize(nsta_total);
    for (auto it=sta_clks.begin(); it!=sta_clks.end(); ++it)
        it->resize(nepo_total);

    skip_header(fp);
    char buf[BUFSIZ];
    int y, m, d, h, min, n, epo;
    double s, sod, bias;
    std::string site;
    CartCoor src, dst;
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strncmp(buf, "AR ", 3) != 0)
            continue;
        site.assign(buf+3, 4);
        std::transform(site.begin(), site.end(), site.begin(),
                        [](unsigned char c) { return std::toupper(c); });
        if (!std::binary_search(sta_list.begin(), sta_list.end(), site))
            continue;
        sscanf(buf+8, "%d %d %d %d %d %lf %d %lf", &y, &m, &d, &h, &min, &s, &n, &bias);
        if (n > 2) {
            fprintf(stderr, MSG_ERR "RinexClk::open: n>2\n");
            return false;
        }
        sod = h*3600 + min*60 + s;
        if (fmod(sod, interval)>MaxWnd || length-sod<0)
            continue;

        auto it = std::lower_bound(sta_list.begin(), sta_list.end(), site);
        int index = it - sta_list.begin();
        if (!rnxsnx_.find_pos(site, src) || !refsnx.find_pos(site, dst))
            bias = None; // 0
        else
            bias -= (src - dst).dot(src)/src.norm()/LightSpeed;

        epo = int(sod/interval);
        sta_clks[index][epo] = bias*1E9; // sec => ns
    }

    fclose(fp);
    return true;
}

static void convert_bias(const double f[], const std::array<double, 4> &bias, double &wl, double &nl)
{
    double g = f[0]/f[1];
    double fw = f[0] - f[1];
    double fn = f[0] + f[1];
    double w1 =  g/(g - 1);
    double w2 = -1/(g - 1);
    double n1 = -g/(g + 1);
    double n2 = -1/(g + 1);
    double alpha = g*g/(g*g - 1);
    double beta = 1 - alpha;

    // bias: L1 L2 P1 P2 (ns)
    wl = (bias[0]*w1 + bias[1]*w2 + bias[2]*n1 + bias[3]*n2)*fw/1E9;
    nl = (bias[0]*alpha + bias[1]*beta)*fn/1E9;
}

bool AnalyseCenter::read_bias(const std::string &path, const std::vector<std::string> &prns,
                              const std::vector<Satellite> &sats)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR "AnalyseCenter::read_bias: no such file: %s\n", path.c_str());
        return false;
    }

    have_bias.assign(prns.size(), 0);
    wl_bias.resize(prns.size());
    nl_bias.resize(prns.size());
    const int interval = 30; // sec
    const size_t nepo_total = sat_clks[0].size();
    const size_t nsat_total = prns.size();
    std::vector<std::vector<std::array<double, 4>>> bias(nsat_total); // L1 L2 P1 P2
    for (size_t i=0; i!=nsat_total; ++i)
        bias[i].resize(nepo_total);

    char buf[BUFSIZ];
    std::string prn;
    int mjd=0, mjd0, mjd1, sod0, sod1, year, doy;
    int index;
    double val;
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strncmp(buf, " OSB ", 5) || buf[15]!=' ')
            continue;

        prn.assign(buf+11, 3);
        val = atof(buf+70); // ns
        if (!std::binary_search(prns.begin(), prns.end(), prn))
            continue;
        auto it = std::lower_bound(prns.begin(), prns.end(), prn);
        index =  it - prns.begin();

        sscanf(buf+35, "%d:%d:%d", &year, &doy, &sod0);
        mjd0 = ydoy2mjd(year, doy);
        sscanf(buf+50, "%d:%d:%d", &year, &doy, &sod1);
        mjd1 = ydoy2mjd(year, doy);

        if (mjd == 0) mjd = mjd0;
        sod0 = (mjd0 - mjd)*86400 + sod0;
        sod1 = (mjd1 - mjd)*86400 + sod1;

        // have_bias[index] = 1;
        for (int i=0; i!=4; ++i) { // L1 L2 P1 P2
            if (strncmp(buf+25, sats[index].obstp[i].c_str(), 3))
                continue;

            size_t beg = sod0/interval;
            size_t end = std::min(static_cast<size_t>(sod1/interval), nepo_total);
            for (size_t epo=beg; epo<end; ++epo)
                bias[index][epo][i] = val;
            break;
        }
    }

    std::vector<std::vector<double>> wls(nsat_total), nls(nsat_total);
    for (size_t i=0; i!=nsat_total; ++i) {
        wls[i].resize(nepo_total);
        nls[i].resize(nepo_total);
    }
    for (size_t isat=0; isat!=nsat_total; ++isat) {
        for (size_t iepo=0; iepo!=nepo_total; ++iepo) {
            if (bias[isat][iepo][0]!=0.0 || bias[isat][iepo][1]!=0.0) {
                convert_bias(sats[isat].f, bias[isat][iepo], wls[isat][iepo], nls[isat][iepo]);
            }
        }
    }
    // Average WLs
    for (size_t isat=0; isat!=nsat_total; ++isat)
    {
        int count=0;
        double sum=0;
        for (size_t iepo=0; iepo!=nepo_total; ++iepo) {
            if (wls[isat][iepo]==0.0 && nls[isat][iepo]==0.0)
                continue;
            sum += wls[isat][iepo];
            ++count;
        }
        if (count != 0) {
            wl_bias[isat] = sum/count;
            have_bias[isat] = 1;
        }
    }

    for (size_t isat=0; isat!=nsat_total; ++isat) {
        double freq = sats[isat].f[0] + sats[isat].f[1]; // NL
        for (size_t iepo=0; iepo!=nepo_total; ++iepo) {
            if (sat_clks[isat][iepo] != None)
                sat_clks[isat][iepo] -= nls[isat][iepo]*(1E9/freq); // ns
        }
    }

    fclose(fp);
    return true;
}

#if 0
bool AnalyseCenter::read_bias(const std::string &path, const std::vector<std::string> &prns,
                              const std::vector<Satellite> &sats)
{
    // if (name == "grg") {
    //     return read_grg_bias(clk_file, prns);
    // } else if (name == "cos" || name == "esa" || name == "gfs" || name == "grs" ) {
    //     return read_sgg_bias(bia_file, prns, sats);
    // } else if (name == "gfz") {
    //     return read_cnes_bias(bia_file, prns, sats);
    // } else if (name == "cod" || name == "whu") {
    //     return read_snx_bias(bia_file, prns, sats);
    // } else {
    //     fprintf(stdout, MSG_ERR "AnalyseCenter::read_bias: unsupported ac: %3s\n", name.c_str());
    //     return false;
    // }
    return read_biassnx(bia_file, prns, sats);
}

bool AnalyseCenter::read_snx_bias(const std::string &path, const std::vector<std::string> &prns,
                                  const std::vector<Satellite> &sats)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR "AnalyseCenter::read_bias: no such file: %s\n", path.c_str());
        return false;
    }

    have_bias.assign(prns.size(), 0);
    wl_bias.resize(prns.size());
    nl_bias.resize(prns.size());
    std::vector<std::array<double, 4>> bias(prns.size()); // L1 L2 P1 P2

    char buf[BUFSIZ];
    std::string prn;
    int index;
    double val;
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strncmp(buf, " OSB ", 5) || buf[15]!=' ')
            continue;

        prn.assign(buf+11, 3);
        val = atof(buf+70); // ns
        if (!std::binary_search(prns.begin(), prns.end(), prn))
            continue;
        auto it = std::lower_bound(prns.begin(), prns.end(), prn);
        index =  it - prns.begin();

        have_bias[index] = 1;
        for (int i=0; i!=4; ++i) // L1 L2 P1 P2
            if (!strncmp(buf+25, sats[index].obstp[i].c_str(), 3)) {
                bias[index][i] = val;
                break;
            }
    }

    for (size_t i=0; i!=prns.size(); ++i) {
        if (!have_bias[i])
            continue;
        convert_bias(sats[i].f, bias[i], wl_bias[i], nl_bias[i]);
        fprintf(stdout, "%3s %3s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", name.c_str(), prns[i].c_str(),
                bias[i][0], bias[i][1], bias[i][2], bias[i][3], wl_bias[i], nl_bias[i]);
    }

    fclose(fp);
    return true;
}

bool AnalyseCenter::read_grg_bias(const std::string &path, const std::vector<std::string> &prns)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR "AnalyseCenter::read_bias: no such file: %s\n", path.c_str());
        return false;
    }

    have_bias.assign(prns.size(), 0);
    wl_bias.resize(prns.size());
    nl_bias.resize(prns.size());

    char buf[BUFSIZ];
    std::string prn;
    int index;
    double val;
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strncmp(buf+60, "END OF HEADER", 13) == 0)
            break;
        else if (strncmp(buf, "WL ", 3) != 0)
            continue;

        prn.assign(buf+3, 3);
        val = atof(buf+40);
        if (!std::binary_search(prns.begin(), prns.end(), prn))
            continue;
        auto it = std::lower_bound(prns.begin(), prns.end(), prn);
        index =  it - prns.begin();
        have_bias[index] = 1;
        wl_bias[index] = -val; // cycle

        fprintf(stdout, "%3s %3s %8.3f\n", name.c_str(), prns[index].c_str(), wl_bias[index]);
    }

    fclose(fp);
    return true;
}

bool AnalyseCenter::read_sgg_bias(const std::string &path,
                                  const std::vector<std::string> &prns,
                                  const std::vector<Satellite> &sats)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR "AnalyseCenter::read_bias: no such file: %s\n", path.c_str());
        return false;
    }

    have_bias.assign(prns.size(), 0);
    wl_bias.resize(prns.size());
    nl_bias.resize(prns.size());

    char buf[BUFSIZ];
    std::string prn;
    int index;
    double val;
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strncmp(buf, "WL ", 3)==0) {
            prn.assign(buf+4, 3);
            val = atof(buf+11); // WL cycle
            if (!std::binary_search(prns.begin(), prns.end(), prn))
                continue;
            auto it = std::lower_bound(prns.begin(), prns.end(), prn);
            index =  it - prns.begin();
            have_bias[index] = 1;
            wl_bias[index] = val;
            fprintf(stdout, "%3s %3s %8.3f\n", name.c_str(), prns[index].c_str(), wl_bias[index]);
        } else if (strncmp(buf+60, "END OF HEADER", 13)==0) {
            break;
        }
    }

    int y, m, d, h, min;
    int nepo=0, nepo_total=sat_clks[0].size();
    double sec, sod, nl_time;
    while (fgets(buf, sizeof(buf), fp))
    {
        if (buf[0] == '*') {
            sscanf(buf+1, "%d%d%d%d%d%lf", &y, &m, &d, &h, &min, &sec);
            sod = h*3600 + min*60 + sec;
            nepo = sod/30.0; // assumed interval: 30s
        } else if (buf[0] == 'P') {
            prn.assign(buf+1, 3);
            if (!std::binary_search(prns.begin(), prns.end(), prn))
                continue;
            auto it = std::lower_bound(prns.begin(), prns.end(), prn);
            index =  it - prns.begin();
            val = atof(buf+20);
            // fprintf(stdout, "%2d %2d %8.3f %3s %8.3f\n", h, min, sec, sats[index].prn.c_str(), val);
            nl_time = val*1E9/(sats[index].f[0] + sats[index].f[1]); // ns
            for (int i=0; i<900/30; ++i)
            {
                if (i+nepo < nepo_total && sat_clks[index][i+nepo] != None)
                    sat_clks[index][i+nepo] -= nl_time;
            }
        }
    }

    fclose(fp);
    return true;
}

bool AnalyseCenter::read_cnes_bias(const std::string &path,
                                   const std::vector<std::string> &prns,
                                   const std::vector<Satellite> &sats)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR "AnalyseCenter::read_bias: no such file: %s\n", path.c_str());
        return false;
    }

    have_bias.assign(prns.size(), 0);
    wl_bias.resize(prns.size());
    nl_bias.resize(prns.size());
    std::vector<std::array<double, 4>> bias(prns.size()); // L1 L2 P1 P2
    // std::array<double, 4> bias;

    char buf[BUFSIZ];
    std::string prn; //, last_prn("XXX");
    int index;
    int year=0, doy, sod, sod0=0, mjd0=0, mjd;
    double val;
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strncmp(buf, " OSB ", 5) || buf[15]!=' ')
            continue;

        prn.assign(buf+11, 3);

        val = atof(buf+70); // ns
        if (year == 0) {
            year = atoi(buf+35);
            doy = atoi(buf+40);
            mjd0 = ydoy2mjd(year, doy);
        }
        year = atoi(buf+35);
        doy = atoi(buf+40);
        mjd = ydoy2mjd(year, doy);
        sod = atoi(buf+44) + (mjd - mjd0)*86400;
        if (!std::binary_search(prns.begin(), prns.end(), prn))
            continue;
        auto it = std::lower_bound(prns.begin(), prns.end(), prn);
        index =  it - prns.begin();

        if (sod > 86370)
            continue;

        if (sod!=sod0 && sod!=0) {
            int ref = 0;
            for (size_t i=0; i!=prns.size(); ++i) {
                if (!have_bias[i])
                    continue;
                if (prns[i][0] != 'C')
                    continue;
                // ref = std::lower_bound(prns.begin(), prns.end(), "C07") - prns.begin();
                if (sats[i].blk[7] == '2')
                    ref = std::lower_bound(prns.begin(), prns.end(), "C06") - prns.begin();
                else
                    ref = std::lower_bound(prns.begin(), prns.end(), "C19") - prns.begin();
                fprintf(stdout, "diff %3s %8d %8.3f %8.3f\n",  prns[i].c_str(), sod,
                        wl_bias[i]-wl_bias[ref], nl_bias[i]-nl_bias[ref]);
            }
            have_bias.assign(prns.size(), 0);
        }

        int i=0;
        for (i=0; i!=4; ++i) // L1 L2 P1 P2
            if (!strncmp(buf+25, sats[index].obstp[i].c_str(), 3)) {
                // bias[i] = val;
                bias[index][i] = val;
                break;
            }


        if (i == 1) {
            have_bias[index] = 1;
            convert_bias(sats[index].f, bias[index], wl_bias[index], nl_bias[index]);
            fprintf(stdout, "%3s %3s %8d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", name.c_str(), prns[index].c_str(),
                    sod, bias[index][0], bias[index][1], bias[index][2], bias[index][3], wl_bias[index], nl_bias[index]);
                    // sod, bias[0], bias[1], bias[2], bias[3], wl_bias[index], nl_bias[index]);
        }
        sod0 = sod;
    }

    fclose(fp);
    return true;
}
#endif
