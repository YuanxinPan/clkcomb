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

#ifndef JPLEPH_H
#define JPLEPH_H

#include <string>
#include <stdio.h>

class JplEph
{
public:
    enum {
        Earth = 3,
        Moon  = 10,
        Sun   = 11
    };

private:
    JplEph(const JplEph &);
    JplEph &operator=(const JplEph &);

public:
    JplEph();
    ~JplEph();
    // open jpleph file before interpolation
    bool open(const std::string &path);
    // interpolate position & velocity of ntarg in ICRS coordinate system
    bool pleph(double et, int ntarg, int ncent, double rrd[]);

private:
    void constan(char nam[][6], double val[], double sss[], int *n);
    void state(double et2[], int list[], double pv[][6], double nut[]);
    void split(double tt, double fr[]);
    void interp(double buf[], double t[], int ncf, int ncm, int na, int ifl, double pv[]);

private:
    FILE *ephFile_;
    int KM, BARY;
    double PVSUN[6];
};

#endif
