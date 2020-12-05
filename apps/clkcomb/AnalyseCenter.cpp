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

bool AnalyseCenter::read_orbit(const std::string &path)
{
    return rnxsp3_.read(path);
}

// bool AnalyseCenter::read_orbit(const std::vector<std::string> &paths)
// {
//     if (rnxsp3_.read(paths)) {
//         return true;
//     } else {
//         return false;
//     }
// }

bool AnalyseCenter::read_sinex(const std::string &path)
{
    return rnxsnx_.read(path);
}

bool AnalyseCenter::open_atx(const std::string &path)
{
    return rnxatx_.open(path);
}

bool AnalyseCenter::open_clock(const std::string &path)
{
    return rnxclk_.read(path);
}

bool AnalyseCenter::read_clock(MJD t, int length, int interval,
                               const std::vector<std::string> &prns, const RinexSp3 &refsp3, RinexAtx &refatx)
{
    const size_t nsat_total = prns.size();
    const size_t nepo_total = length/interval + 1;

    sat_clks.clear();
    sat_clks.resize(nsat_total);

    double clk, corr;
    const RinexAtx::atx_t *atx_src = nullptr;
    const RinexAtx::atx_t *atx_dst = nullptr;
    double psrc[3], pdst[3];
    std::string f1, f2;
    CartCoor pos, ref;
    for (size_t epo=0; epo!=nepo_total; ++epo, t+=interval) {
        for (size_t isat=0; isat!=nsat_total; ++isat) {
            if (!rnxclk_.clkBias(t, prns[isat], clk))
                clk = None;
            atx_src = rnxatx_.atx(t, prns[isat]);
            atx_dst = refatx.atx(t, prns[isat]);
            if (atx_src == nullptr || atx_dst == nullptr)
                clk = None;
            else {
                if (prns[isat][0] == 'G') {
                    f1 = "G01"; f2 = "G02";
                } else if (prns[isat][0] == 'E') {
                    f1 = "E01"; f2 = "E05";
                } else if (prns[isat][0] == 'R') {
                    f1 = "R01"; f2 = "R02";
                } else if (prns[isat][0] == 'C') {
                    f1 = "C01"; f2 = "C07";
                }
                atx_src->pco(psrc, f1, f2);
                atx_dst->pco(pdst, f1, f2);
            }
            if (!rnxsp3_.satPos(t, prns[isat], pos) ||
                !refsp3.satPos( t, prns[isat], ref)) {
                clk = None;
            } else if (clk != None) {
                corr = (pos - ref).dot(pos)/pos.norm();
                clk -= corr/LightSpeed;
                clk += (psrc[2] - pdst[2])/LightSpeed;
            }

            sat_clks[isat].push_back(clk*1E9); // sec => ns
        }
    }

    rnxclk_.close();
    rnxatx_.close();
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

bool AnalyseCenter::read_bias(const std::string &path, const std::vector<std::string> &prns,
                              const std::vector<Satellite> &sats)
{
    if (name == "grg") {
        return read_grg_bias(clk_file, prns);
    } else if (name == "cod" || name == "whu") {
        return read_snx_bias(bia_file, prns, sats);
    } else {
        fprintf(stdout, MSG_ERR "AnalyseCenter::read_bias: unsupported ac: %3s\n", name.c_str());
        return false;
    }
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
