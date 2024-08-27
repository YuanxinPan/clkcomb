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

#ifndef _COORD_H_
#define _COORD_H_

#include "CartCoor.h"
#include "Ellipsoid.h"
#include "JplEph.h"

// 3-d vector only
double dot(const double *p, const double *q);

// out can share the same address of p/q
void cross(const double *p, const double *q, double *out);

// out can share the same address of p
void unit_vector(const double *p, double *out);

double qAngularDist(const double *qs, const double *qr);

void xyz2blh(const double xyz[], double blh[]);

void blh2xyz(const double blh[], double xyz[]);

// convert vector in ECEF to ENU in local-coordinate-system
// vec and enu can share the same address
void ecef2enu(const double vec[], const double blh_ref[], double enu[]);

// convert vector in ENU to ECEF
// vec and ecef can share the same address
void enu2ecef(const double vec[], const double blh_ref[], double ecef[]);

void elevazim(const double *site_pos, const double *sat_pos, double &elev, double &azim);

// los: line-of-sight in original coordinate system
// R:   rotation matrix from original to destination coordinate system
void elevazim(const double *los, const double *R, double *elev, double *azim);

// convert displacements in ECEF to ENU
// disp and enu can share the same address
void disp2enu(const double *disp, const double *R, double *enu);

// R: rotation matrix, column major
void rot_enu2ecef(const double *xyz, double *R);

// R: rotation matrix, column major
void rot_ecef2enu(const double *xyz, double *R);

void SunPosition(int mjd, double sod, double rs[]);
void MoonPosition(int mjd, double sod, double rm[]);

void nominal_att(const double *xsat, const double *xsun, double *q);

#endif
