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

#ifndef CLKCOMB_CONST_HPP
#define CLKCOMB_CONST_HPP

const double PI = 3.141592653589793238462643383;

const double R2D = 180.0/PI;

const double D2R = PI/180.0;

const double LightSpeed = 299792458;

// earth gravitational constant
const double GM = 3.986005E+14;

// earth rotate angle-velocity
const double We = 7.2921151467E-05;

// relative effect
const double F = -4.442807633E-10;

enum GNSS_Tp {
    _GPS_ = 0,
    _GLS_,
    _GAL_,
    _BD2_,
    _BD3_,
    _QZS_,
    _UNKNOWN_ = 99
};

// GPS frequency
const double GPS_f1 = 1.57542E+09;
const double GPS_f2 = 1.22760E+09;
const double GPS_f5 = 1.17645E+09;
// GLS frequency
const double GLS_f1  = 1.6020E+09;
const double GLS_f2  = 1.2460E+09;
const double GLS_df1 = 0.5625E+06;
const double GLS_df2 = 0.4375E+06;
// BDS frequency
const double BDS_f1 = 1.561098E+09;
const double BDS_f2 = 1.207140E+09;
const double BDS_f3 = 1.268520E+09;
// GAL frequency
const double GAL_f1 = 1.575420E+09;
const double GAL_f5 = 1.176450E+09;
const double GAL_f6 = 1.278750E+09;
const double GAL_f7 = 1.207140E+09;
const double GAL_f8 = 1.191795E+09;

// iono-free coef
const double IF_coef1 = (GPS_f1*GPS_f1)/(GPS_f1*GPS_f1 - GPS_f2*GPS_f2);
const double IF_coef2 = -(GPS_f2*GPS_f2)/(GPS_f1*GPS_f1 - GPS_f2*GPS_f2);

const double MaxWnd = 1E-8;

#endif
