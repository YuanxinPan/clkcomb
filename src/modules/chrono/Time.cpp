// clkcomb - Clock and phase bias products Combination
// Copyright (C) 2021 Yuanxin Pan
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "chrono.h"
#include "../io/io.h"

#include <stdio.h>
#include <math.h>
#include <algorithm>

bool CalendTime::read(const char buf[])
{
    if (6 == sscanf(buf, "%d %d %d %d %d %lf", &y, &m, &d, &h, &min, &sec))
        return true;
    return false;
}

void CalendTime::write(FILE *fp)
{
    fprintf(fp, "%4d-%02d-%02d-%02d:%02d:%4.1f", yr2year(y), m, d, h, min, sec);
}

int CalendTime::mjd()const
{
    return date2mjd(y, m, d);
}

double CalendTime::sod()const
{
    return hms2sod(h, min, sec);
}

//std::ostream &operator<<(std::ostream &out, const GPST &t)
//{
//    out << std::setw(8) << t.week << std::setw(10) << t.sow;
//    return out;
//}

void GPST::write(FILE *fp)
{
    fprintf(fp, "%4d%10.1f", week, sow);
}

void MJD::set(int week, double sow)
{
    d =  week*7 + 44244 + static_cast<int>(sow/86400);
    sod = fmod(sow, 86400.0);
}

void MJD::adjust()
{
    if (sod < 0.0) {
        d -= 1;
        sod += 86400.0;
    }
    else if (sod >= 86400.0) {
        d += 1;
        sod -= 86400.0;
    }
}

std::vector<MJD::leapsec_t> MJD::leapsecs_;

bool MJD::read_leapsec(const std::string &path)
{
    if (leapsecs_.size() != 0u)
        return true;
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR "MJD::read_leapsec: no such file: %s\n", path.c_str());
        return false;
    }

    leapsec_t leapsec;
    char buf[128];
    while (fgets(buf, sizeof(buf), fp)) {
        int n = sscanf(buf, "%lf %d", &leapsec.t, &leapsec.leap);
        if (n == 2) leapsecs_.push_back(leapsec);
    }
    std::sort(leapsecs_.begin(), leapsecs_.end(),
              [](const leapsec_t &lhs, const leapsec_t &rhs) { return lhs.t < rhs.t; });

    fclose(fp);
    return true;
}

int MJD::leapsec(double t)
{
    t = static_cast<int>(t);  // only doy-part matters
	if (t < leapsecs_.front().t || t > leapsecs_.back().t) {
        fprintf(stderr, MSG_ERR "CoordConverter::leapsec: request time %7.1f "
                "is out of range [%7.1f, %7.1f]\n", t,
                leapsecs_.front().t, leapsecs_.back().t);
        return 0;
	}
    // return the leap seconds at mjd
    auto it = std::lower_bound(leapsecs_.begin(), leapsecs_.end(), t,
              [](const leapsec_t &leap, double t) { return leap.t < t; });
    //printf("%10.1f %3d\n", t, it->leap - 1);
    return it->leap - 1;
}
