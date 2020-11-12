#include "RinexClk.h"
#include "../io/io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cctype>
#include <algorithm>
#include <pppx/const.h>

void RinexClk::close()
{
    if (clkFile_ != nullptr)
        fclose(clkFile_);
    clkFile_ = nullptr;
}

bool RinexClk::read(const std::string &path)
{
    clkFile_ = fopen(path.c_str(), "r");
    if (clkFile_ == nullptr) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexClk::open: no such file: %s\n", path.c_str());
        return false;
    }

    size_t ns = 0u;
    std::string prn;
    while (fgets(buf_, 256, clkFile_))
    {
        if (strncmp(buf_+60, "END OF HEADER", 13) == 0)
            break;
        else if (strncmp(buf_+60, "# OF SOLN STA / TRF", 19) == 0) {
            int n = atoi(buf_);
            sta_names_.reserve(n);
            // sta_poss_.reserve(n);
        } else if (strncmp(buf_+60, "SOLN STA NAME / NUM", 19) == 0) {
            std::string name(buf_, 4);
            std::transform(name.begin(), name.end(), name.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            sta_names_.push_back(name);
            // double x, y, z;
            // sscanf(buf_+14, "%lf %lf %lf", &x, &y, &z);
            // sta_poss_.emplace_back(x/1E3, y/1E3, z/1E3);
        } else if (strncmp(buf_+60, "# OF SOLN SATS", 14) == 0) {
            int n = atoi(buf_+3);
            ns = static_cast<size_t>(n);
            prns_.reserve(ns);
        } else if (strncmp(buf_+60, "PRN LIST", 8) == 0) {
            for (int i=0; i<15 && prns_.size()<ns; ++i) {
                //if (buf[] == ' ') buf[] = '0';
                //if (buf[] == ' ') buf[] = 'G';
                prn.assign(buf_+4*i, 3);
                prns_.push_back(prn);
            }
        }
    }

    // printf("%s: %3lu\n", path.c_str(), sta_names_.size());
    // for (size_t i=0, size=sta_names_.size(); i!=size; ++i) {
    //     printf(" %4s %14.3f %14.3f %14.3f\n", sta_names_[i].c_str(), sta_poss_[i].x, sta_poss_[i].y, sta_poss_[i].z);
    // }
    // printf("\n");
    // fflush(stdout);

    if (prns_.size() == 0u) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexClk::open: no PRN LIST\n");
        return false;
    }
    std::sort(prns_.begin(), prns_.end());
    coefs_[0].assign(prns_.size(), Coef_t());
    coefs_[1].assign(prns_.size(), Coef_t());

    MJD cur;
    int y, m, d, h, min, i=0, n;
    double s, bias;
    auto beg = prns_.begin();
    while (fgets(buf_, 256, clkFile_))
    {
        if (strncmp(buf_, "AS ", 3) != 0)
            continue;
        sscanf(buf_+7, "%d %d %d %d %d %lf %d %lf", &y, &m, &d, &h, &min, &s, &n, &bias);
        if (n > 2) {
            fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                    "RinexClk::open: n>2\n");
            exit(1);
        }
        cur.d = date2mjd(y, m, d);
        cur.sod = hms2sod(h, min, s);
        if (time_[i].d == 0)
            time_[i] = cur;
        else if (cur - time_[i] > 0) {
            ++i;
            if (i==2) break;
        }
        prn.assign(buf_+3, 3);
        if (!std::binary_search(beg, prns_.end(), prn))
            continue;
        auto it = std::lower_bound(beg, prns_.end(), prn);
        coefs_[i][it-beg].a = bias;
        coefs_[i][it-beg].valid = true;
    }

    interval_ = time_[1] - time_[0];
    return true;
}

bool RinexClk::clkBias(MJD t, const std::string &prn, double &sck, double *_drift)
{
    if (time_[0]-t > MaxWnd)
        return false;
    while (t-time_[1] > MaxWnd) {
        if (!update())
            return false;
    }

    if (!std::binary_search(prns_.begin(), prns_.end(), prn))
        return false;
    auto it = std::lower_bound(prns_.begin(), prns_.end(), prn);
    int iprn = it - prns_.begin();
    if (!coefs_[0][iprn].valid || !coefs_[1][iprn].valid)
        return false;

    double a[2] = { coefs_[0][iprn].a, coefs_[1][iprn].a };
    double drift = (a[1] - a[0])/interval_;
    sck = a[0] + drift*(t - time_[0]);
    //fprintf(logFile, "%6.0f %3s %6.0f %6.0f %16.12f\n",
    //        t.sod, prns[iprn].c_str(), time_[0].sod, time_[1].sod, sck);
    if (_drift) *_drift = drift;

    return true;
}

bool RinexClk::update()
{
    if (feof(clkFile_))
        return false;

    coefs_[0].swap(coefs_[1]);
    std::for_each(coefs_[1].begin(), coefs_[1].end(),
                  [](Coef_t &coef) { coef.valid=false; });
    time_[0] = time_[1];
    time_[1].sod += interval_;

    bool first = true;
    MJD cur;
    std::string prn;
    int y, m, d, h, min, n;
    double s, bias;
    auto beg = prns_.begin();
    do {
        if (strncmp(buf_, "AS ", 3) != 0)
            continue;
        sscanf(buf_+7, "%d %d %d %d %d %lf %d %lf", &y, &m, &d, &h, &min, &s, &n, &bias);
        if (n > 2) {
            fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                    "RinexClk::open: n>2\n");
            exit(1);
        }
        cur.d = date2mjd(y, m, d);
        cur.sod = hms2sod(h, min, s);
        if (first) {
            first = false;
            time_[1] = cur;
        }
        else if (cur - time_[1] > MaxWnd) {
            break;
        }
        prn.assign(buf_+3, 3);
        // remove if use IGS clk
        if (!std::binary_search(beg, prns_.end(), prn))
            continue;
        auto it = std::lower_bound(beg, prns_.end(), prn);
        coefs_[1][it-beg].a = bias;
        coefs_[1][it-beg].valid = true;
    } while (fgets(buf_, 256, clkFile_));

    return true;
}
