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

#ifndef _CHRONO_H_
#define _CHRONO_H_

#include "Time.h"

inline int yr2year(int y)
{
    if (y > 99) // 2017
        return y;
    else if (y < 80) // 17
        return y + 2000;
    else // 98
        return y + 1900;
}

inline double hms2sod(int h, int m, double s)
{
    return h*3600 + m*60 + s;
}

void sod2hms(double sod, int *h, int *m, double *s);

int date2mjd(int y, int m, int d);

void mjd2date(int jd, double sod, int *iy, int *imon, int *id, int *ih, int *im, double *is);

void wksow2mjd(int week, double sow, int *mjd, double *sod);

void mjd2wksow(int mjd, double sod, int &week, double &sow);

void mjd2doy(int mjd, int *yr, int *doy);

int date2doy(int y, int m, int d);

void doy2date(int iyear, int idoy, int *imonth, int *iday);

int ydoy2mjd(int year, int doy);

void date2gwk(int year, int month, int day, int *week, int *wd);

const char *run_time();

#endif
