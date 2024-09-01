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

#include "coord.h"
#include "../const.h"

#include <math.h>

double dot(const double *p, const double *q)
{
    //double d = 0.0;
    //while (n > 0) {
    //    d += *(p)**(q);
    //    ++p; ++q; --n;
    //}
    return p[0]*q[0] + p[1]*q[1] + p[2]*q[2];
}

void cross(const double *p, const double *q, double *out)
{
    double tmp[3];
    tmp[0] = p[1]*q[2] - p[2]*q[1];
    tmp[1] = p[2]*q[0] - p[0]*q[2];
    tmp[2] = p[0]*q[1] - p[1]*q[0];
    out[0] = tmp[0];
    out[1] = tmp[1];
    out[2] = tmp[2];
}

void unit_vector(const double *p, double *out)
{
    double norm = sqrt(dot(p, p));
    out[0] = p[0]/norm;
    out[1] = p[1]/norm;
    out[2] = p[2]/norm;
}

static void rotmat2quater(const double *R, double *q)
{
    q[0] = sqrt(1 + R[0] + R[4] + R[8])/2;
    q[1] = (R[5] - R[7])/4/q[0];
    q[2] = (R[6] - R[2])/4/q[0];
    q[3] = (R[1] - R[3])/4/q[0];
}

static void mat3_transpose(const double *m, double *n)
{
    n[0] = m[0]; n[4] = m[4]; n[8] = m[8];
    n[1] = m[3]; n[2] = m[6]; n[3] = m[1];
    n[5] = m[7]; n[6] = m[2]; n[7] = m[5];
}

static void qcross(const double *q1, const double *q2, double *q)
{
    q[0] = q1[0]*q2[0] - dot(q1+1, q2+1);
    double c[3];
    cross(q2+1, q1+1, c);
    for (int i=0; i<3; ++i)
        q[i+1] = q2[0]*q1[i+1] + q1[0]*q2[i+1] + c[i];
}

static void qinv(const double *q, double *qi)
{
    qi[0] =  q[0];
    qi[1] = -q[1];
    qi[2] = -q[2];
    qi[3] = -q[3];
}

double qAngularDist(const double *qs, const double *qr)
{
    // double cosa = 0;
    // for (int i=0; i<4; i++) {
    //     cosa += qs[i]*qr[i];
    // }
    // if (cosa < 0.0) {
    //     for (int i=0; i<4; i++) qs[i] = -qs[i]
    // }

    double qi[4], qc[4];
    qinv(qs, qi);
    qcross(qi, qr, qc);
    double diff = 2*atan2(qc[3], qc[0])*R2D;
    return diff;
}

static void matmpy(const double *mat, const double *vec, double *out, bool transpose=false)
{
    double v[3] = { 0 };
    int istep = transpose ?  3 : 1;
    int jstep = transpose ?  1 : 3;
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            v[i] += mat[i*istep + j*jstep]*vec[j];
    out[0] = v[0];
    out[1] = v[1];
    out[2] = v[2];
}

void xyz2blh(const double xyz[], double blh[])
{
    double x=xyz[0], y=xyz[1], z=xyz[2];
    double r = sqrt(x*x + y*y);
    blh[1] = atan2(y, x);

    double B=0.7, B0, H, N;
    double e2 = WGS84.e2;
    do {
        B0 = B;
        double sinB = sin(B);
        N = WGS84.a/sqrt(1 - e2*sinB*sinB);
        H = r/cos(B) - N;
        B = atan(z/r/(1 - e2*N/(N+H)));
    } while(fabs(B-B0) > 1.E-12);   // 1E-10 -> 0.64 mm

    blh[0]=B;
    blh[2]=H;
}

void blh2xyz(const double blh[], double xyz[])
{
    double e2 = WGS84.e2;
    double sinB = sin(blh[0]);
    double cosB = cos(blh[0]);
    double N = WGS84.a/sqrt(1 - e2*sinB*sinB);

    xyz[0] = (N + blh[2])*cosB*cos(blh[1]);
    xyz[1] = (N + blh[2])*cosB*sin(blh[1]);
    xyz[2] = ((1-e2)*N + blh[2])*sinB;
}

static void _ecef2enu(const double xyz[], const double xyz_ref[], double enu[])
{
    double blh[3];
    xyz2blh(xyz_ref, blh);
    double sinB = sin(blh[0]);
    double cosB = cos(blh[0]);
    double sinL = sin(blh[1]);
    double cosL = cos(blh[1]);

    double dx = xyz[0]-xyz_ref[0];
    double dy = xyz[1]-xyz_ref[1];
    double dz = xyz[2]-xyz_ref[2];

    enu[0] = -sinL*dx + cosL*dy;
    enu[1] = -cosL*sinB*dx - sinL*sinB*dy + cosB*dz;
    enu[2] =  cosL*cosB*dx + sinL*cosB*dy + sinB*dz;
}

void ecef2enu(const double vec[], const double blh_ref[], double enu[])
{
    double sinB = sin(blh_ref[0]);
    double cosB = cos(blh_ref[0]);
    double sinL = sin(blh_ref[1]);
    double cosL = cos(blh_ref[1]);

    double e = -sinL*vec[0] + cosL*vec[1];
    double n = -cosL*sinB*vec[0] - sinL*sinB*vec[1] + cosB*vec[2];
    double u =  cosL*cosB*vec[0] + sinL*cosB*vec[1] + sinB*vec[2];
    enu[0] = e; enu[1] = n; enu[2] = u;
}

void enu2ecef(const double vec[], const double blh_ref[], double ecef[])
{
    double sinB = sin(blh_ref[0]);
    double cosB = cos(blh_ref[0]);
    double sinL = sin(blh_ref[1]);
    double cosL = cos(blh_ref[1]);

    double x = -sinL*vec[0] - cosL*sinB*vec[1] + cosL*cosB*vec[2];
    double y =  cosL*vec[0] - sinL*sinB*vec[1] + sinL*cosB*vec[2];
    double z =  cosB*vec[1] + sinB*vec[2];
    ecef[0] = x; ecef[1] = y; ecef[2] = z;
}

void elevazim(const double *site_pos, const double *sat_pos, double &elev, double &azim)
{
    double enu[3];
    _ecef2enu(sat_pos, site_pos, enu);
    double hori = sqrt(enu[0]*enu[0] + enu[1]*enu[1]);
    elev = atan2(enu[2], hori);
    azim = atan2(enu[0], enu[1]);
}

void elevazim(const double *los, const double *R, double *elev, double *azim)
{
    double enu[3];
    matmpy(R, los, enu);
    double hori = sqrt(enu[0]*enu[0] + enu[1]*enu[1]);
    *elev = atan2(enu[2], hori);
    *azim = atan2(enu[0], enu[1]);
}

void disp2enu(const double *disp, const double *R, double *enu)
{
    matmpy(R, disp, enu);
}

void rot_enu2ecef(const double *xyz, double *R)
{
    double blh[3];
    xyz2blh(xyz, blh);
    double sinB = sin(blh[0]);
    double cosB = cos(blh[0]);
    double sinL = sin(blh[1]);
    double cosL = cos(blh[1]);

    // column major
    R[0 + 0*3] = -sinL;
    R[0 + 1*3] = -cosL*sinB;
    R[0 + 2*3] = cosL*cosB;

    R[1 + 0*3] = cosL;
    R[1 + 1*3] = -sinL*sinB;
    R[1 + 2*3] = sinL*cosB;

    R[2 + 0*3] = 0;
    R[2 + 1*3] = cosB;
    R[2 + 2*3] = sinB;
}

void rot_ecef2enu(const double *xyz, double *R)
{
    double blh[3];
    xyz2blh(xyz, blh);
    double sinB = sin(blh[0]);
    double cosB = cos(blh[0]);
    double sinL = sin(blh[1]);
    double cosL = cos(blh[1]);

    // column major
    R[0*3 + 0] = -sinL;
    R[0*3 + 1] = -cosL*sinB;
    R[0*3 + 2] = cosL*cosB;

    R[1*3 + 0] = cosL;
    R[1*3 + 1] = -sinL*sinB;
    R[1*3 + 2] = sinL*cosB;

    R[2*3 + 0] = 0;
    R[2*3 + 1] = cosB;
    R[2*3 + 2] = sinB;
}

// Functions to compute Sun & Moon position

static void getghar(int jd, double sod, double *ghar)
{
    double D2R=PI/180.0;
    /* need UT to get sidereal time */
    // double tsecgps = sod;                                  /* GPS time (sec of day) */
    double tsecutc = sod-18;                               /* UTC time (sec of day) */
    double fmjdutc = tsecutc / 86400.0;                    /* UTC time (fract. day) */
    double d = (jd - 51544) + (fmjdutc - 0.50);            /* days since J2000      */
    double ghad = 280.46061837504 + 360.9856473662862 * d; /* corrn.   (+digits)    */
    int i = (int)(ghad / 360.0);
    *ghar = (ghad - i * 360.0) * D2R;
    while (*ghar >= 2 * PI) {
        *ghar = *ghar - 2 * PI;
    }
    while (*ghar < 0.0) {
        *ghar = *ghar + 2 * PI;
    }
}

static void rot1(double theta, double x, double y, double z, double *u, double *v, double *w)
{
    double s = sin(theta);
    double c = cos(theta);

    *u = x;
    *v = c*y+s*z;
    *w = c*z-s*y;
}

static void rot3(double theta, double x, double y, double z, double *u, double *v, double *w)
{
    double s = sin(theta);
    double c = cos(theta);
    *u = c * x + s * y;
    *v = c * y - s * x;
    *w = z;
}

extern void SunPosition(int jd, double sod, double rs[])
{
    /* mean elements for year 2000, sun ecliptic orbit wrt. Earth */
    double obe     = 23.43929111 * D2R;         /* obliquity of the J2000 ecliptic */
    double sobe    = sin(obe);
    double cobe    = cos(obe);
    double opod    = 282.94;                    /* RAAN + arg.peri. (deg.) */

    /* use TT for solar ephemerides */
    // double tsecgps = sod;                       /* GPS time (sec of day)   */
    double tsectt  = sod+19.0+32.184;
                                                /* TT  time (sec of day) */
    double fmjdtt  = tsectt / 86400.0;          /* TT  time (fract. day) */

    /* julian centuries since 1.5 january 2000 (J2000) */
    /*    (note also low precision use of mjd --> tjd) */
    double tjdtt   = jd + fmjdtt + 2400000.5;           /* Julian Date, TT      */
    double t       = (tjdtt - 2451545.0) / 36525.0;     /* Julian centuries, TT */
    double emdeg   = 357.5256 + 35999.049 * t;          /* degrees              */
    double em      = emdeg * D2R;                       /* radians              */
    double em2     = em + em;                           /* radians              */

    /* series expansions in mean anomaly, em */
    double r       = (149.619 - 2.499 * cos(em) - 0.021 * cos(em2)) * 1.0e9;        /* m */
    double slond   = opod + emdeg + (6892.0 * sin(em) + 72.0 * sin(em2)) / 3600.0;

    /* precession of equinox */
    slond         += 1.3972 * t; /* degrees */

    /* position vector of sun (mean equinox & ecliptic of J2000) (EME2000, ICRF) */
    /*                                    (plus long. advance due to precession) */
    double slon    = slond * D2R;      /* radians */
    double sslon   = sin(slon);
    double cslon   = cos(slon);

    double rs1     = r * cslon;        /* meters */
    double rs2     = r * sslon * cobe; /* meters */
    double rs3     = r * sslon * sobe; /* meters */

    /* convert position vector of sun to ECEF (ignore polar motion / LOD) */
    double ghar = 0.0;
    getghar(jd, sod,&ghar);
    rot3(ghar, rs1, rs2, rs3, rs, rs+1, rs+2);
}

extern void MoonPosition(int jd, double sod, double rm[])
{
    /* use TT for lunar ephemerides */
    // double tsecgps = sod;                                        /* GPS time (sec of day) */
    double tsectt  = sod+19.0+32.184;                            /* TT  time (sec of day) */
    double fmjdtt  = tsectt / 86400.0;                           /* TT  time (fract. day) */
    /* julian centuries since 1.5 january 2000 (J2000) */
    /*    (note also low precision use of mjd --> tjd) */
    double tjdtt   = jd + fmjdtt + 2400000.5;                    /* Julian Date, TT      */
    double t       = (tjdtt - 2451545.0) / 36525.0;              /* Julian centuries, TT */

    /* el0 -- mean longitude of Moon (deg)                             */
    /* el  -- mean anomaly of Moon (deg)                               */
    /* elp -- mean anomaly of Sun (sun)                                */
    /* f   -- mean angular distance of Moon from ascending node (deg)  */
    /* d   -- difference between mean longitudes of Sun and Moon (deg) */
    double el0     = 218.31617 + 481267.88088 * t - 1.3972 * t;
    double el      = 134.96292 + 477198.86753 * t;
    double elp     = 357.52543 +  35999.04944 * t;
    double f       =  93.27283 + 483202.01873 * t;
    double d       = 297.85027 + 445267.11135 * t;

    /* longitude w.r.t. equinox and ecliptic of year 2000 */
    double selond  = el0
                   + ((22640.0 * sin((el              ) * D2R))
                   +  (  769.0 * sin((el + el         ) * D2R))
                   -  ( 4586.0 * sin((el - d - d      ) * D2R))
                   +  ( 2370.0 * sin((d + d           ) * D2R))
                   -  (  668.0 * sin((elp             ) * D2R))
                   -  (  412.0 * sin((f + f           ) * D2R))
                   -  (  212.0 * sin((el + el - d - d ) * D2R))
                   -  (  206.0 * sin((el + elp - d - d) * D2R))
                   +  (  192.0 * sin((el + d + d      ) * D2R))
                   -  (  165.0 * sin((elp - d - d     ) * D2R))
                   +  (  148.0 * sin((el - elp        ) * D2R))
                   -  (  125.0 * sin((d               ) * D2R))
                   -  (  110.0 * sin((el + elp        ) * D2R))
                   -  (   55.0 * sin((f + f - d - d   ) * D2R))) / 3600.0;

    /* latitude w.r.t. equinox and ecliptic of year 2000 */
    double q      = (412.0 * sin((f + f) * D2R) + 541.0 * sin((elp) * D2R)) / 3600.0;
    double selatd =
                  + ((18520.0 * sin((f + selond - el0 + q) * D2R))
                  - (   526.0 * sin((f - d - d           ) * D2R))
                  + (    44.0 * sin((el + f - d - d      ) * D2R))
                  - (    31.0 * sin((-el + f - d - d     ) * D2R))
                  - (    25.0 * sin((-el - el + f        ) * D2R))
                  - (    23.0 * sin((elp + f - d - d     ) * D2R))
                  + (    21.0 * sin((-el + f             ) * D2R))
                  + (    11.0 * sin((-elp + f - d - d    ) * D2R))) / 3600.0;

    /* distance from Earth center to Moon (m) */
    double rse    = (385000.0
                  - ( 20905.0 * cos((el              ) * D2R))
                  - (  3699.0 * cos((d + d - el      ) * D2R))
                  - (  2956.0 * cos((d + d           ) * D2R))
                  - (   570.0 * cos((el + el         ) * D2R))
                  + (   246.0 * cos((el + el - d - d ) * D2R))
                  - (   205.0 * cos((elp - d - d     ) * D2R))
                  - (   171.0 * cos((el + d + d      ) * D2R))
                  - (   152.0 * cos((el + elp - d - d) * D2R))) * 1000.0;

    /* convert spherical ecliptic coordinates to equatorial cartesian */
    /* precession of equinox wrt. J2000 */
    selond += 1.3972 * t;    /* degrees */

    /* position vector of moon (mean equinox & ecliptic of J2000) (EME2000, ICRF) */
    /*                                     (plus long. advance due to precession) */
    double oblir  = 23.43929111 * D2R;         /* obliquity of the J2000 ecliptic */

    double sselat = sin(selatd * D2R);
    double cselat = cos(selatd * D2R);
    double sselon = sin(selond * D2R);
    double cselon = cos(selond * D2R);

    double t1 = rse * cselon * cselat;      /* meters */
    double t2 = rse * sselon * cselat;      /* meters */
    double t3 = rse *          sselat;      /* meters */

    double rm1 = 0.0;
    double rm2 = 0.0;
    double rm3 = 0.0;
    rot1(-oblir, t1, t2, t3, &rm1, &rm2, &rm3);

    /* convert position vector of moon to ECEF (ignore polar motion / LOD) */
    double ghar = 0.0;
    getghar(jd, sod,&ghar);
    rot3(ghar, rm1, rm2, rm3, rm, rm+1, rm+2);
}

// ECEF => SV
void nominal_att(const double *xsat, const double *xsun, double *q)
{
    double R[9]; // sv => ecef
    double *xscf = R;
    double *yscf = R+3;
    double *zscf = R+6;

    //sc-fixed z-axis, from sc to earth
    unit_vector(xsat, zscf);
    zscf[0] = -zscf[0];
    zscf[1] = -zscf[1];
    zscf[2] = -zscf[2];

    //unint vector from sc to sun
    double u_sc2sun[3];
    u_sc2sun[0] = xsun[0] - xsat[0];
    u_sc2sun[1] = xsun[1] - xsat[1];
    u_sc2sun[2] = xsun[2] - xsat[2];
    unit_vector(u_sc2sun, u_sc2sun);

    //the sc-fixed y-axis
    cross(zscf, u_sc2sun, yscf);
    unit_vector(yscf, yscf);

    //the sc-fixed x-axis
    cross(yscf, zscf, xscf);

    double Ri[9]; // ecef => sv
    mat3_transpose(R, Ri);
    rotmat2quater(Ri, q);
}
