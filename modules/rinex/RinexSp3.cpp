#include "RinexSp3.h"
#include "../io/io.h"
#include "../const.h"
#include "../coord/coord.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>

bool operator<(const RinexSp3::Sp3_t &lhs, const RinexSp3::Sp3_t &rhs)
{
    return lhs.t < rhs.t;
}

bool RinexSp3::read(const std::string &path)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexSp3::read: no such file: %s\n", path.c_str());
        return false;
    }

    size_t ns = 0u;
    int nline = 0;
    char buf[128];
    std::string prn;

    while (fgets(buf, sizeof(buf), fp) && buf[0]!='*') {
        ++nline;
        if (nline == 3) {
            ns = static_cast<size_t>(atoi(buf+3));
            prns_.reserve(ns);
            for (int j=0; j<17&&prns_.size()<ns; ++j) {
                prn.assign(buf+9+3*j, 3);
                if (prn[0] == ' ') prn[0] = 'G';
                if (prn[1] == ' ') prn[1] = '0';
                prns_.push_back(prn);
            }
        } else if (nline>3 && prns_.size()<ns) {
            for (int j=0; j<17&&prns_.size()<ns; ++j) {
                prn.assign(buf+9+3*j, 3);
                if (prn[0] == ' ') prn[0] = 'G';
                if (prn[1] == ' ') prn[1] = '0';
                prns_.push_back(prn);
            }
        }
    }

    std::sort(prns_.begin(), prns_.end());
    satposs_.assign(ns, std::vector<Sp3_t>());

    int y, m, d, h, min;
    double s;
    Sp3_t sp3;
    do {
        if (buf[0] == '*') {
            int n = sscanf(buf+1, "%d%d%d%d%d%lf", &y, &m, &d, &h, &min, &s);
            if (n != 6) {
                fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                        "RinexSp3::read: %s", buf);
                return false;
            }
            sp3.t.d = date2mjd(y, m, d);
            sp3.t.sod = hms2sod(h, min, s);
        } else {
            int n = sscanf(buf+4, "%lf%lf%lf", &sp3.pos.x, &sp3.pos.y, &sp3.pos.z);
            if (n != 3) {
                fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                        "RinexSp3::read: %s", buf);
                return false;
            }
            if (sp3.pos.norm() < MaxWnd || isnan(sp3.pos.x+sp3.pos.y+sp3.pos.z))
                continue;

            prn.assign(buf+1, 3);
            if (prn[0] == ' ') prn[0] = 'G';
            if (prn[1] == ' ') prn[1] = '0';
            if (std::binary_search(prns_.begin(), prns_.end(), prn)) {
                auto it = std::lower_bound(prns_.begin(), prns_.end(), prn);
                int iprn = it - prns_.begin();
                satposs_[iprn].push_back(sp3);
            }
        }
    } while (fgets(buf, sizeof(buf), fp) && strncmp(buf, "EOF", 3) != 0);

    fclose(fp);
    return true;
}

bool RinexSp3::read(const std::vector<std::string> &paths)
{
    if (paths.size() != 2u) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexSp3::read: 2-day SP3 files needed\n");
        return false;
    }

    FILE *fp = fopen(paths[0].c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexSp3::read: no such file: %s\n", paths[0].c_str());
        return false;
    }

    size_t ns = 0u;
    int nline = 0;
    char buf[128];
    std::string prn;

    // Read PRNs
    while (fgets(buf, sizeof(buf), fp) && buf[0]!='*') {
        ++nline;
        if (nline == 3) {
            ns = static_cast<size_t>(atoi(buf+3));
            prns_.reserve(ns);
            for (int j=0; j<17&&prns_.size()<ns; ++j) {
                prn.assign(buf+9+3*j, 3);
                if (prn[0] == ' ') prn[0] = 'G';
                if (prn[1] == ' ') prn[1] = '0';
                prns_.push_back(prn);
            }
        } else if (nline>3 && prns_.size()<ns) {
            for (int j=0; j<17&&prns_.size()<ns; ++j) {
                prn.assign(buf+9+3*j, 3);
                if (prn[0] == ' ') prn[0] = 'G';
                if (prn[1] == ' ') prn[1] = '0';
                prns_.push_back(prn);
            }
        }
    }
    std::sort(prns_.begin(), prns_.end());
    satposs_.assign(ns, std::vector<Sp3_t>());
    fclose(fp);

    for (auto it=paths.begin(); it!=paths.end(); ++it)
    {
        fp = fopen(it->c_str(), "r");
        if (fp == nullptr) {
            fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                    "RinexSp3::read: no such file: %s\n", it->c_str());
            return false;
        }

        // skip header
        while (fgets(buf, sizeof(buf), fp) && buf[0]!='*') { }

        bool end_of_day_record = false; // have SP3 record for 86400 s
        int nrec = 0;
        int y, m, d, h, min;
        double s;
        MJD t;
        Sp3_t sp3;
        do {
            if (buf[0] == '*') {
                if (it!=paths.begin() && ++nrec>1)
                    break;
                int n = sscanf(buf+1, "%d%d%d%d%d%lf", &y, &m, &d, &h, &min, &s);
                if (n != 6) {
                    fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                            "RinexSp3::read: %s", buf);
                    return false;
                }
                sp3.t.d = date2mjd(y, m, d);
                sp3.t.sod = hms2sod(h, min, s);

                if (t.d == 0) {
                    t = sp3.t;
                } else if (sp3.t-t > 86400.0-1E-3) { // no need for SP3 of the next day
                    end_of_day_record = true;
                }
            } else {
                int n = sscanf(buf+4, "%lf%lf%lf", &sp3.pos.x, &sp3.pos.y, &sp3.pos.z);
                if (n != 3) {
                    fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                            "RinexSp3::read: %s", buf);
                    return false;
                }
                if (sp3.pos.norm() < MaxWnd)
                    continue;

                prn.assign(buf+1, 3);
                if (prn[0] == ' ') prn[0] = 'G';
                if (prn[1] == ' ') prn[1] = '0';
                if (std::binary_search(prns_.begin(), prns_.end(), prn)) {
                    auto it = std::lower_bound(prns_.begin(), prns_.end(), prn);
                    int iprn = it - prns_.begin();
                    satposs_[iprn].push_back(sp3);
                }
            }
        } while (fgets(buf, sizeof(buf), fp) && strncmp(buf, "EOF", 3) != 0);

        fclose(fp);
        if (end_of_day_record)
            break;
    } // for in paths

    for (auto it=satposs_.begin(); it!=satposs_.end(); ++it)
        std::sort(it->begin(), it->end());

    return true;
}

static double interp_orbit(const double *x, double *y, int n)
{
    int i, j;
    for (j=1; j<n; ++j) {
        for (i=0; i<n-j; ++i) {
            y[i] = (x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
        }
    }
    return y[0];
}

bool RinexSp3::satPos(MJD t, const std::string &prn, CartCoor &pos)const
{
    return satPos(t, prn, pos.data());
}

bool RinexSp3::satPos(MJD t, const std::string &prn, double *pos)const
{
    if (!std::binary_search(prns_.begin(), prns_.end(), prn))
        return false;
    auto it = std::lower_bound(prns_.begin(), prns_.end(), prn);
    int iprn = it - prns_.begin();
    if (satposs_[iprn].size() == 0u)
        return false;

    // no extrapolation
    if (satposs_[iprn].front().t-t>0.1 || t-satposs_[iprn].back().t>0.1)
        return false;

    const int NMAX = 10; // polynomial degree
    int i, j, k, n=satposs_[iprn].size();
    double x[NMAX+1], y[3][NMAX+1];
    // double theta, sinl, cosl;

    /* binary search */
    for (i=0, j=n-1; i<j;) {
        k = (i+j)/2;
        if (satposs_[iprn][k].t < t) i=k+1; else j=k;
    }
    int index = i<=0?0:i-1;
    i = index - (NMAX+1)/2;
    if (i<0) i=0; else if (i+NMAX >= n) i=n-NMAX-1;

    CartCoor p;
    for (j=0; j<=NMAX; ++j) {
        x[j] = satposs_[iprn][i+j].t - t;
        p = satposs_[iprn][i+j].pos;
        y[0][j] = p.x;
        y[1][j] = p.y;
        y[2][j] = p.z;
    }

    // unit: m
    pos[0] = 1E3*interp_orbit(x, y[0], NMAX+1);
    pos[1] = 1E3*interp_orbit(x, y[1], NMAX+1);
    pos[2] = 1E3*interp_orbit(x, y[2], NMAX+1);
    return true;
}

bool RinexSp3::satPosVel(MJD t, const std::string &prn, CartCoor &pos, CartCoor &vel)const
{
    return satPosVel(t, prn, pos.data(), vel.data());
}

bool RinexSp3::satPosVel(MJD t, const std::string &prn, double *pos, double *vel)const
{
    const double dt = 1E-3;
    double poss[2][3];
    if (satPos(t-dt, prn, poss[0]) && satPos(t+dt, prn, poss[1])) {
        pos[0] = (poss[0][0] + poss[1][0])/2.0;
        pos[1] = (poss[0][1] + poss[1][1])/2.0;
        pos[2] = (poss[0][2] + poss[1][2])/2.0;

        vel[0] = (poss[1][0] - poss[0][0])/(2.0*dt);
        vel[1] = (poss[1][1] - poss[0][1])/(2.0*dt);
        vel[2] = (poss[1][2] - poss[0][2])/(2.0*dt);
        return true;
    }
    return false;
}

