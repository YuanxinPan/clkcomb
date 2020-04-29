#include "CoordConverter.h"
#include "coord.h"
#include "../io/io.h"
#include "../chrono/chrono.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <Eigen/Dense>
#include <pppx/const.h>

using Eigen::Matrix3d;

double ut1ut1r(double mjd);
double ERA2000(double DJ1, double DJ2);
double EE2000(double DATE1, double DATE2, double DPSI);
double EECT2000(double DATE1, double DATE2);
double GMST2000(double UTA, double UTB, double TTA, double TTB);
double GST2000(double UTA, double UTB, double TTA, double TTB, double DPSI);
double SP2000(double DATE1, double DATE2);

void T2C2000(const Matrix3d &RPOM, double THETA, const Matrix3d &RBPN, Matrix3d &RMAT, Matrix3d &RT2C);
Matrix3d CBPN2000(double DATE1, double DATE2, double DPSI, double DEPS);
void NU2000B(double DATE1, double DATE2, double *DPSI, double *DEPS);
void ortho_eop(double time, double *dx, double *dy, double *dut1);
Matrix3d POM2000(double XP, double YP, double SP);

enum Axis {
    XAXIS = 0,
    YAXIS = 1,
    ZAXIS = 2
};

std::vector<CoordConverter::erp_t> CoordConverter::erps_;

void *CoordConverter::mat_enu2ecef(const double *xyz)
{
    double blh[3];
    xyz2blh(xyz, blh);
    double sinB = sin(blh[0]);
    double cosB = cos(blh[0]);
    double sinL = sin(blh[1]);
    double cosL = cos(blh[1]);

    Matrix3d *p = new Matrix3d(Matrix3d::Identity());
    (*p)(0, 0) = -sinL;
    (*p)(0, 1) = -cosL*sinB;
    (*p)(0, 2) = cosL*cosB;

    (*p)(1, 0) = cosL;
    (*p)(1, 1) = -sinL*sinB;
    (*p)(1, 2) = sinL*cosB;

    (*p)(2, 1) = cosB;
    (*p)(2, 2) = sinB;
    return p;
}

CartCoor CoordConverter::convert(void *pmat, CartCoor pos, enum Type tp)
{
    Eigen::Vector3d v(pos.x, pos.y, pos.z);
    Matrix3d *p = reinterpret_cast<Matrix3d*>(pmat);
    switch (tp)
    {
    case TOPO2ECEF:
        v = (*p)*v;
        break;
    case ECEF2TOPO:
        v = p->transpose()*v;
        break;
    default:
        break;
    }
    return CartCoor(v[0], v[1], v[2]);
}

bool CoordConverter::read_igserp(const std::string &path)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "CoordConverter::read_igserp: no such file: %s\n", path.c_str());
        return false;
    }

    char buf[256];
    fgets(buf, sizeof(buf), fp);
    if (strncmp(buf, "VERSION 2", 9)!=0 && strncmp(buf, "version 2", 9)!=0) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "CoordConverter::read_igserp: unknown version: %s", buf);
        return false;
    }

    bool but1utc = true; // true:UT1UTC; false:UT1TAI
    while (fgets(buf, sizeof(buf), fp)) {
        if (strncmp(buf+2, "MJD", 3) == 0) {
            std::string str(buf);
            if (str.find("UTC") != std::string::npos) {
                but1utc = true;
                if (!MJD::read_leapsec("leap.sec")) return false;
            }
            else if (str.find("TAI") != std::string::npos) {
                but1utc = false;
            }
            else {
                fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                        "CoordConverter::read_igserp: UTC/TAI not found: %s", buf);
                return false;
            }
            fgets(buf, sizeof(buf), fp);
            break;
        }
    }
    if (feof(fp)) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "CoordConverter::read_igserp: MJD not found\n");
        return false;
    }

    erp_t erp;
    while (fgets(buf, sizeof(buf), fp)) {
        int n = sscanf(buf, "%lf%lf%lf%lf", &erp.t, &erp.dx, &erp.dy, &erp.ut1tai);
        if (n != 4) continue;
        erp.dx *= 1E-6;
        erp.dy *= 1E-6;
        erp.ut1tai *= 1E-7;
        erp.ut1tai -= ut1ut1r(erp.t+GPST2TT/86400.0);
        if (but1utc) {
            int sec = MJD::leapsec(erp.t);
            if (sec == 0) return false;
            //printf("leapsec: %8.1f%3d\n", erp.t, sec);
            erp.ut1tai -= sec;
        }
        erps_.push_back(erp);
    }
    std::sort(erps_.begin(), erps_.end(),
              [](const erp_t &lhs, const erp_t &rhs) { return lhs.t < rhs.t; });

    fclose(fp);
    return true;
}

static Eigen::Matrix3d ecef2iner;
static Eigen::Matrix3d mrate;
static Eigen::Matrix3d iner2ecef;

void ef2int(double t, double ut1tai, double *xhelp, Matrix3d &mate2j,
            Matrix3d &rmte2j, double *xpole, double *ypole);

bool CoordConverter::update(double t, double *xpole, double *ypole)
{
    if (t < erps_.front().t || t > erps_.back().t) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "CoordConverter::update: %7.1f "
                "is out of range [%7.1f, %7.1f]\n", t,
                erps_.front().t, erps_.back().t);
        return 0;
    }

    double pole[2];
    double ut1tai;
    auto ie = std::lower_bound(erps_.begin(), erps_.end(), t,
              [](const erp_t &erp, double t) { return erp.t < t; });
    if (ie->t == t) {
        pole[0] = ie->dx;
        pole[1] = ie->dy;
        ut1tai = ie->ut1tai;
    }
    else {
        auto ib = ie - 1;
        //printf("update: %7.1f use %7.1f %7.1f\n", t, ib->t, ie->t);
        double dx = ie->dx - ib->dx;
        double dy = ie->dy - ib->dy;
        double dut1tai = ie->ut1tai - ib->ut1tai;
        double ratio = (t - ib->t)/(ie->t - ib->t);
        pole[0] = ib->dx + dx*ratio;
        pole[1] = ib->dy + dy*ratio;
        ut1tai = ib->ut1tai + dut1tai*ratio;
        //printf("ratio: %8.6f\n", ratio);
    }
    //printf("%10.5E %10.5E %10.6E\n", pole[0], pole[1], ut1tai+37.0);
    double xp, yp;
    ef2int(t, ut1tai, pole, ecef2iner, mrate, &xp, &yp);
    iner2ecef = ecef2iner.transpose();
    if (xpole) {
        *xpole = xp;
        *ypole = yp;
    }
    return true;
}

CartCoor CoordConverter::convert(CartCoor pos, enum Type tp)
{
    Eigen::Vector3d v(pos.x, pos.y, pos.z);
    switch (tp) {
    case ITRS2ICRS:
        v = ecef2iner*v;
        break;
    case ICRS2ITRS:
        v = iner2ecef*v;
        break;
    default:
        break;
    }
    return CartCoor(v[0], v[1], v[2]);
}

Matrix3d rot_matrix(double theta, enum Axis axis)
{
    Matrix3d R(Matrix3d::Identity());
    switch (axis) {
    case XAXIS:
        R(1,1) = R(2,2) = cos(theta);
        R(1,2) =  sin(theta);
        R(2,1) = -sin(theta);
        break;
    case YAXIS:
        R(0,0) = R(2,2) = cos(theta);
        R(0,2) = -sin(theta);
        R(2,0) =  sin(theta);
        break;
    case ZAXIS:
        R(0,0) = R(1,1) = cos(theta);
        R(0,1) =  sin(theta);
        R(1,0) = -sin(theta);
        break;
    default:
        break;
    }
    return R;
}

// purpose  : rotation from earth - fixed system to inertial system(IERS 2003)
// input :
// jd, sod -- GPS time
// output : mate2j -- rotation matrix
// rmte2j -- first derivative of rotation matrix
// gast   -- GAST angle(radian)
// xpole, ypole -- pole movement(arcsec)
void ef2int(double t, double ut1tai, double *xhelp, Matrix3d &mate2j,
            Matrix3d &rmte2j, double *xpole, double *ypole)
{
    double dphi, depsi;
    double dx, dy, dut1, dels;

    //printf("XPOLE: %20.15f%20.15f%20.15f\n", xhelp[0], xhelp[1], ut1tai);

    double tt = t + GPST2TT/86400.0;
    // Tidal variations in Earth rotation
    ortho_eop(tt, &dx, &dy, &dut1);
    *xpole = xhelp[0] + dx;
    *ypole = xhelp[1] + dy;
    //printf("DXPOLE: %20.15f%20.15f%20.15f\n", dx, dy, dut1);

    // Quantity s'
    dels = SP2000(tt, MJD2JD);
    //printf("Quantity: %20.15f\n", dels);

    // polar motion
    Matrix3d xmat = POM2000(*xpole/3600.0/180.0*PI, *ypole/3600.0/180.0*PI, dels);
    //printf("Polar Motion:\n");
    //for (int i=0; i<3; ++i)
    //    printf("%20.15f%20.15f%20.15f\n", xmat(i,0), xmat(i,1), xmat(i,2));

    // nutation
    NU2000B(tt, MJD2JD, &dphi, &depsi);
    //printf("Nutation: %20.15f%20.15f\n", dphi, depsi);

    double dut1_long = ut1ut1r(tt);
    double ut1 = t + (GPST2TAI + ut1tai + dut1_long + dut1)/86400.0;
    //printf("%24.15f%24.15f%24.15f\n", dut1_long, ut1tai+dut1_long, dut1);

    //int d = int(ut1);
    //double sod = (ut1 - d)*86400.0;
    //printf("UT1: %8d %20.15f\n", d, sod);
    // Greenwich
    double gast = GST2000(ut1, MJD2JD, tt, MJD2JD, dphi);
    //printf("GAST: %20.15f\n", gast);

    // Q matrix
    Matrix3d pmat = CBPN2000(tt, MJD2JD, dphi, depsi);
    //printf("Q matrix:\n");
    //for (int i=0; i<3; ++i)
    //    printf("%20.15f%20.15f%20.15f\n", pmat(i,0), pmat(i,1), pmat(i,2));

    // transformation
    T2C2000(xmat, gast, pmat, rmte2j, mate2j);
    //printf("Finally:\n");
    //for (int i=0; i<3; ++i)
    //    printf("%20.15f%20.15f%20.15f\n", mate2j(i,0), mate2j(i,1), mate2j(i,2));
}

// Form the TRS - to - CRS matrix, IAU 2000, from components.
// Annex to IERS Conventions 2000, Chapter 5
//
// Given:
// RPOM     d(3, 3)   polar motion matrix(W)
//
// followed by either(for the CEO - based transformation) :
// THETA      d      Earth Rotation Angle(radians, giving matrix R)
// RBPN     d(3, 3)   intermediate - to - celestial matrix(Q)
//
// or alternatively(for the classical, equinox - based, transformation) :
// THETA      d      Greenwich Sidereal Time(radians)
// RBPN     d(3, 3)   true - to - celestial matrix
// RMAT     d(3, 3)   rate matrix
//
// Returned:
// RT2C     d(3, 3)   terrestrial - to - celestial matrix
//
// Calls the SOFA routines iau_CR, iau_RZ, iau_RXR.
// revision : 2002 November 25
void T2C2000(const Matrix3d &RPOM, double THETA, const Matrix3d &RBPN, Matrix3d &RMAT, Matrix3d &RT2C)
{
    double  SG = sin(THETA), CG = cos(THETA);
    // Rate matrix
    RMAT(2,0) = 0; RMAT(0,2) = 0;
    RMAT(1,2) = 0; RMAT(2,1) = 0;
    RMAT(2,2) = 0;
    RMAT(0,0) = -SG*We; RMAT(0,1) = -CG*We;
    RMAT(1,0) =  CG*We; RMAT(1,1) = -SG*We;

    // Polar motion.
    //R = RPOM
    //MATMPY(RMAT[0], R[0], RMAT[0], 3, 3, 3);
    RMAT = RMAT*RPOM;

    // Earth rotation.
    //ROT_Z(-THETA, R);
    Matrix3d R = rot_matrix(-THETA, ZAXIS)*RPOM;

    // CIP motion.
    //MATMPY(RBPN[0],    R[0], RT2C[0], 3, 3, 3);
    //MATMPY(RBPN[0], RMAT[0], RMAT[0], 3, 3, 3);
    RT2C = RBPN*R;
    RMAT = RBPN*RMAT;
}

// Classical bias - precession - nutation matrix.
// Annexe to IERS Conventions 2000, Chapter 5
//
// Given:
// DATE1, DATE2   d      TT date(JD = DATE1 + DATE2)
// DPSI, DEPS     d      nutation(luni - solar + planetary, radians)
//
// Returned :
// RBPNC       d(3, 3)   true - to - celestial matrix
// revision : 2002 November 26
Matrix3d CBPN2000(double DATE1, double DATE2, double DPSI, double DEPS)
{
    // Arcseconds to radians
    const double DAS2R = 4.848136811095359935899141E-6;

    // Reference epoch(J2000), JD
    const double DJ0 = 2451545.0;

    // Days per Julian century
    const double DJC = 36525.0;

    // J2000 obliquity(Lieske et al. 1977)
    const double EPS0 = 84381.448 * DAS2R;

    // The ICRS RA of the J2000 equinox(Chapront et al., 2002)
    const double DRA0 = -0.0146 * DAS2R;

    // The precession and obliquity corrections(radians per century)
    const double PRECOR = -0.29965*DAS2R;
    const double OBLCOR = -0.02524*DAS2R;

    // The frame bias corrections in longitude and obliquity
    const double DPSIBI = -0.041775*DAS2R;
    const double DEPSBI = -0.0068192*DAS2R;

    double DPSIPR, DEPSPR, EPSA80, PSIA77, OMA77, CHIA, PSIA, OMA, EPSA;

    // Interval between fundamental epoch J2000.0 and given date(JC).
    double T = ((DATE1 - DJ0) + DATE2) / DJC;

    // Precession rate contributions with respect to IAU 1976 / 80.
    DPSIPR = PRECOR*T;
    DEPSPR = OBLCOR*T;

    // IAU 1980 mean obliquity of date.
    EPSA80 = EPS0 + (-46.8150 + (-0.00059 + (0.001813)*T)*T)*T*DAS2R;

    // Precession angles(Lieske et al. 1977)
    PSIA77 = (5038.7784 + (-1.07259 + (-0.001147)*T)*T)*T*DAS2R;
    OMA77 = EPS0 + ((0.05127 + (-0.007726)*T)*T)*T*DAS2R;
    CHIA = (10.5526 + (-2.38064 + (-0.001125)*T)*T)*T*DAS2R;

    // Apply IAU 2000A precession corrections.
    PSIA = PSIA77 + DPSIPR;
    OMA = OMA77 + DEPSPR;
    EPSA = EPSA80 + DEPSPR;

    Matrix3d RBPNC = rot_matrix(-EPSA, XAXIS) * rot_matrix(DPSI, ZAXIS)
                     * rot_matrix(EPSA+DEPS, XAXIS);

    RBPNC = rot_matrix(-EPS0, XAXIS) * rot_matrix(PSIA, ZAXIS)
            * rot_matrix(OMA, XAXIS) * rot_matrix(-CHIA, ZAXIS) * RBPNC;

    RBPNC = rot_matrix(-DRA0, ZAXIS) * rot_matrix(-DPSIBI*sin(EPS0), YAXIS)
            * rot_matrix(DEPSBI, XAXIS) * RBPNC;

    return RBPNC;
}

// Greenwich Mean Sidereal Time(model consistent with IAU 2000
// resolutions).
// Annexe to IERS Conventions 2000, Chapter 5
//
// Given:
// UTA, UTB     d      UT1 date(JD = UTA + UTB)
// TTA, TTB     d      TT date(JD = TTA + TTB)
//
// The result is the Greenwich Mean Sidereal Time(radians), in the
// range 0 to 2pi.
//
// revision : 2002 December 2
double GMST2000(double UTA, double UTB, double TTA, double TTB)
{
    // Arcseconds to radians
    const double DAS2R = 4.848136811095359935899141E-6;
    // Reference epoch(J2000), JD
    const double DJ0 = 2451545.0;
    // Days per Julian century
    const double DJC = 36525.0;
    const double D2PI = 6.283185307179586476925287;

    // TT Julian centuries since J2000.0.
    double T = ((TTA - DJ0) + TTB)/DJC;

    // Greenwich Mean Sidereal Time, IAU 2000.
    double ret = ERA2000(UTA, UTB) + (0.0145060 + (4612.157399660 + (1.396677210 +
                 (-0.000093440 + 0.000018820*T)*T)*T)*T)*DAS2R;
    ret = fmod(ret, D2PI);
    if (ret < 0.0)
        ret = ret + D2PI;
    return ret;
}

// Nutation, IAU 2000B (truncated)model.
// Annexe to IERS Conventions 2000, Chapter 5
//
// Given:
// DATE1, DATE2    d   TT date(JD = DATE1 + DATE2)
//
// Returned :
// DPSI, DEPS      d   nutation(luni - solar + planetary, radians)
//
//  revision : 2002 November 25
#define NLS  77
void NU2000B(double DATE1, double DATE2, double *DPSI, double *DEPS)
{
    // Arcseconds to radians
    const double DAS2R = 4.848136811095359935899141E-6;

    // Milliarcseconds to radians
    const double DMAS2R = DAS2R/1.E3;

    // Arc seconds in a full circle
    const double TURNAS = 1296000;

    // 2Pi
    const double D2PI = 6.283185307179586476925287;

    // Units of 0.1 microarcsecond to radians
    const double U2R = DAS2R/1.E7;

    // Reference epoch(J2000), JD
    const double DJ0 = 2451545;

    // Days per Julian century
    const double DJC = 36525;

    // Miscellaneous
    double T, EL, ELP, F, D, OM, ARG, DP, DE, SARG, CARG, DPSILS, DEPSLS, DPSIPL, DEPSPL;
    // ------------------------ -
    // Luni - Solar nutation model
    // ------------------------ -

    // Number of terms in the luni - solar nutation model
    // const int NLS = 77;

    // Coefficients for fundamental arguments
    static int NALS[5][NLS] = {
        0,  0,  0, 0,  0,  0, 1,  0, 1, 0, 0, -1, -1, 1, -1, -1, 1, -2, 0, 0,
        0, -2,  2, 1, -1,  2, 0,  0,-1, 0, 0,  1,  0,-1,  0,  1,-2,  0, 0, 0,
        0,  1,  2,-2,  2,  0, 0, -1, 2, 1, 0,  1, -2, 3,  0,  1, 0, -1,-1, 0,
        -2, 1,  2,-1,  1,  1,-1,  1,-1, 0, -1,-1,  0, 1, -2, -1, 1,  0, 0, 0,
        0,  1,  1, 0,  0,  0,-1,  0, 0, 0, 0,  0,  0, 0,  0,  0, 0, -2, 0, 0,
        0,  0,  0, 0,  1,  0, 2,  0, 0,-1, 0,  2,  0, 0,  1,  0,-1,  0, 0, 0,
        0,  0, -1, 0, -1,  0, 0,  1,-1, 0, 0, -1, -1, 0, -1,  0,-1,  0, 1, 0,
        1,  1,  0, 0,  0,  0, 0,  0, 1,-2, 0,  0,  0, 1,  0,  2, 2,  0, 0, 2,
        0,  2,  2, 2,  2,  2, 0,  0, 0, 2, 2,  2,  0, 2,  2,  0, 2,  2, 2, 0,
        2,  0,  0, 2, -2,  0, 0,  2, 0, 2, 2,  2,  2, 2,  0,  2, 2,  0, 2, 2,
        0,  0,  0, 0,  2,  0, 2,  2, 0, 2, 0,  2,  2, 2,  0,  2, 0,  0, 0, 2,
        2,  0,  0, 2,  2,  0, 2,  2, 2, 0, 2,  0, -2, 0,  0,  0,-2,  0, 0, 0,
        -2,-2,  0, 2,  0,  0, 2,  0, 0, 2, 2, -2,  2, 0, -2,  0, 0,  0, 0, 2,
        -2, 2, -2, 0,  2,  0, 2,  0, 0, 2, 0,  2, -2,-2,  2,  0,-2, -2, 2,-2,
        2, -2,  0, 0,  0,  2, 0,  1, 2, 0, 2,  0,  0, 0,  1,  0, 0, -2, 0, 1,
        1,  4,  1,-2,  2,  2, 0, -2, 1, 2, 2,  2,  0, 2,  0,  1, 2,  2, 1, 2,
        0,  1,  1, 2,  1,  1, 0,  2, 2, 0, 2,  2,  1, 0,  0,  1, 1,  2, 0, 1,
        1,  1,  0, 2,  0,  2, 1,  2, 1, 1, 2,  1,  1, 1,  1,  0, 1,  0, 1, 0,
        2,  2,  0, 2,  0,  2, 0,  2, 1, 2, 1,  0,  0, 0,  1,  2, 0,  2, 2, 1,
        1,  1,  2, 2,  2 };

    // Longitude and obliquity coefficients
    // double CLS[6][NLS];

    // ------------------
    // Planetary nutation(radians)
    // ------------------

    const double DPPLAN = -0.135*DMAS2R, DEPLAN = +0.388*DMAS2R;

    // n.b.The above fixed terms account for the omission of the
    // long - period planetary terms in the truncated model.

    // ----------------------------------------
    // Tables of argument and term coefficients
    // ----------------------------------------

    //
    // Luni - Solar argument multipliers :
    // L     L'    F     D     Om

    // Luni - Solar nutation coefficients, unit 1e-7 arcsec:
    // longitude(sin, t*sin, cos), obliquity(cos, t*cos, sin)
    //
    // Each row of coefficients in CLS belongs with the corresponding row of
    // fundamental - argument multipliers in NALS.
    static const double CLS[6][NLS] = {
    -172064161, -13170906, -2276413, 2074554, 1475877, -516821, 711159, -387298, -301461, 215829,
        128227,    123457,   156994,   63110,  -57976,  -59641, -51613,   45893,   63384, -38571,
         32481,    -47722,   -31046,   28593,   20441,   29243,  25887,  -14053,   15164, -15794,
         21783,    -12873,   -12654,  -10204,   16707,   -7691, -11024,    7566,   -6637,  -7141,
         -6302,      5800,     6443,   -5774,   -5350,   -4752,  -4940,    7350,    4065,   6579,
          3579,      4725,    -3075,   -2904,    4348,   -2878,  -4230,   -2819,   -4056,  -2647,
         -2294,      2481,     2179,    3276,   -3389,    3339,  -1987,   -1981,    4026,   1660,
         -1521,      1314,    -1283,   -1331,    1383,    1405,   1290, -174666,   -1675,   -234,
           207,     -3633,     1226,      73,    -367,     -36,   -494,     137,      11,     10,
            63,       -63,      -11,     -42,      50,      11,     -1,       0,       0,     -1,
             0,        21,        0,       0,     -25,      10,     72,       0,     -10,     11,
             0,       -85,        0,       0,     -21,     -11,     21,     -11,      10,      0,
           -11,         0,      -11,     -11,       0,       0,      0,       0,       0,      0,
             0,         0,        0,       0,       0,       0,      0,       0,       0,      0,
             0,         0,        0,       0,       0,       0,      0,       0,       0,      0,
             0,         0,        0,       0,   33386,  -13696,   2796,    -698,   11817,   -524,
          -872,       380,      816,     111,     181,      19,   -168,      27,    -189,    149,
           129,        31,     -150,     158,       0,     -18,    131,      -1,      10,    -74,
           -66,        79,       11,     -16,      13,     -37,     63,      25,     -10,     44,
           -14,       -11,       25,       8,       2,       2,     -7,     -15,      21,     -3,
           -21,        -8,        6,     -24,       5,      -6,     -2,      15,     -10,      8,
             5,         7,        5,      11,     -10,      -7,     -2,       1,       5,    -13,
            -6,         0,     -353,      -5,       9,       0,      0,       8,      -2,      4,
             0,  92052331,  5730336,  978459, -897492,   73871, 224386,   -6750,  200728, 129025,
        -95929,    -68982,   -53311,   -1235,  -33228,   31429,  25543,   26366,  -24236,  -1220,
         16452,    -13870,      477,   13238,  -12338,  -10758,   -609,    -550,    8551,  -8001,
          6850,      -167,     6953,    6415,    5222,     168,   3268,     104,   -3250,   3353,
          3070,      3272,    -3045,   -2768,    3041,    2695,   2719,    2720,     -51,  -2206,
          -199,     -1900,      -41,    1313,    1233,     -81,   1232,     -20,    1207,     40,
          1129,      1266,    -1062,   -1129,      -9,      35,   -107,    1073,     854,   -553,
          -710,       647,     -700,     672,     663,    -594,   -610,    -556,    9086,  -3015,
          -485,       470,     -184,    -677,       0,      18,    -63,     299,      -9,     32,
             0,         0,        0,     -11,       0,     -10,      0,     -11,       0,      0,
           -11,        10,        0,       0,       0,      -2,      0,     -42,       0,      0,
             0,         0,       -1,       0,       0,       0,      0,       0,       0,      0,
             0,         0,        0,       0,       0,       0,      0,       0,       0,      0,
             0,         0,        0,       0,       0,       0,      0,       0,       0,      0,
             0,         0,        0,       0,       0,       0,      0,       0,       0,      0,
             0,         0,        0,       0,       0,   15377,  -4587,    1374,    -291,  -1924,
          -174,       358,      318,     367,     132,      39,     -4,      82,      -9,    -75,
            66,        78,       20,      29,      68,       0,    -25,      59,      -3,     -3,
            13,        11,      -45,      -1,      -5,      13,    -14,      26,      15,     10,
            19,         2,       -5,      14,       4,       4,     -1,      -4,      -5,     12,
            -3,        -9,        4,       1,       2,       1,      3,      -1,       7,      2,
             4,        -2,        3,      -2,       5,      -4,     -3,      -2,       0,     -2,
             1,        -2,        0,    -139,      -2,       4,      0,       0,       4,     -2,
             2,         0 };

    // Interval between fundamental epoch J2000.0 and given date(JC).
    T = ((DATE1 - DJ0) + DATE2)/DJC;

    // ------------------ -
    // LUNI - SOLAR NUTATION
    // ------------------ -

    // Fundamental(Delaunay) arguments from Simon et al. (1994)

    // Mean anomaly of the Moon.
    EL = fmod(485868.249036 + T*1717915923.2178, TURNAS)*DAS2R;

    // Mean anomaly of the Sun.
    ELP = fmod(1287104.79305 + T*129596581.0481, TURNAS)*DAS2R;

    // Mean argument of the latitude of the Moon.
    F = fmod(335779.526232 + T*1739527262.8478, TURNAS)*DAS2R;

    // Mean elongation of the Moon from the Sun.
    D = fmod(1072260.70369 + T*1602961601.2090, TURNAS)*DAS2R;

    // Mean longitude of the ascending node of the Moon.
    OM = fmod(450160.398036 - T*6962890.5431, TURNAS)*DAS2R;

    // Initialize the nutation values.
    DP = 0;
    DE = 0;

    // Summation of luni - solar nutation series(in reverse order).
    //DO I = NLS, 1, -1
    for (int I = NLS; I >= 1; --I) {
        // Argument and functions.
        ARG = fmod((double)(NALS[0][I-1])*EL + (double)(NALS[1][I-1])*ELP +
              (double)(NALS[2][I-1])*F + (double)(NALS[3][I-1])*D +
              (double)(NALS[4][I-1])*OM, D2PI);
        SARG = sin(ARG);
        CARG = cos(ARG);

        // Term.
        DP = DP + (CLS[0][I-1] + CLS[1][I-1]*T)*SARG + CLS[2][I-1]*CARG;
        DE = DE + (CLS[3][I-1] + CLS[4][I-1]*T)*CARG + CLS[5][I-1]*SARG;
    }

    // Convert from 0.1 microarcsec units to radians.
    DPSILS = DP*U2R;
    DEPSLS = DE*U2R;

    // ------------------
    // PLANETARY NUTATION
    // ------------------

    // Fixed terms to allow for long - period nutation.
    DPSIPL = DPPLAN;
    DEPSPL = DEPLAN;

    // Add planetary and luni - solar components.
    *DPSI = DPSIPL + DPSILS;
    *DEPS = DEPSPL + DEPSLS;
}

void cnmtx(double dmjd, double h[])
{
    const int  nlines = 71;
    const double dt = 2.0;
    const double nmax = 2;
    const double d1960 = 37076.5;
    //the orthotide weight factors
    static const double sp[6][2] = {
        0.0298, 0.0200,  0.1408,  0.0905, +0.0805, +0.0638,
        0.6002, 0.3476, +0.3025, +0.1645,  0.1517,  0.0923 };
    //for (i= 0; i<41; ++i) { nj[i]=2; mj[i]=1; }
    //for (i=41; i<71; ++i) { nj[i]=2; mj[i]=2; }
    //static const int nj[nlines] = {
    //  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    //  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    //  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    //  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };
    static const int mj[nlines] = {
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };
    //tidal potential model for 71 diurnal and semidiurnal lines
    static const double freq[nlines]= {
         5.18688050,  5.38346657,  5.38439079,  5.41398343,
         5.41490765,  5.61149372,  5.61241794,  5.64201057,
         5.64293479,  5.83859664,  5.83952086,  5.84044508,
         5.84433381,  5.87485066,  6.03795537,  6.06754801,
         6.06847223,  6.07236095,  6.07328517,  6.10287781,
         6.24878055,  6.26505830,  6.26598252,  6.28318449,
         6.28318613,  6.29946388,  6.30038810,  6.30131232,
         6.30223654,  6.31759007,  6.33479368,  6.49789839,
         6.52841524,  6.52933946,  6.72592553,  6.75644239,
         6.76033111,  6.76125533,  6.76217955,  6.98835826,
         6.98928248, 11.45675174, 11.48726860, 11.68477889,
        11.71529575, 11.73249771, 11.89560406, 11.91188181,
        11.91280603, 11.93000800, 11.94332289, 11.96052486,
        12.11031632, 12.12363121, 12.13990896, 12.14083318,
        12.15803515, 12.33834347, 12.36886033, 12.37274905,
        12.37367327, 12.54916865, 12.56637061, 12.58357258,
        12.59985198, 12.60077620, 12.60170041, 12.60262463,
        12.82880334, 12.82972756, 13.06071921 };
    static const double phase[nlines] = {
         9.0899831,  8.8234208, 12.1189598,  1.4425700,  4.7381090,
         4.4715466,  7.7670857, -2.9093042,  0.3862349, -3.1758666,
         0.1196725,  3.4152116, 12.8946194,  5.5137686,  6.4441883,
        -4.2322016, -0.9366625,  8.5427453, 11.8382843,  1.1618945,
         5.9693878, -1.2032249,  2.0923141, -1.7847596,  8.0679449,
         0.8953321,  4.1908712,  7.4864102, 10.7819493,  0.3137975,
         6.2894282,  7.2198478, -0.1610030,  3.1345361,  2.8679737,
        -4.5128771,  4.9665307,  8.2620698, 11.5576089,  0.6146566,
         3.9101957, 20.6617051, 13.2808543, 16.3098310,  8.9289802,
         5.0519065, 15.8350306,  8.6624178, 11.9579569,  8.0808832,
         4.5771061,  0.7000324, 14.9869335, 11.4831564,  4.3105437,
         7.6060827,  3.7290090, 10.6350594,  3.2542086, 12.7336164,
        16.0291555, 10.1602590,  6.2831853,  2.4061116,  5.0862033,
         8.3817423, 11.6772814, 14.9728205,  4.0298682,  7.3254073,
         9.1574019 };
    static const double hs[nlines] = {
         -1.94,  -1.25,  -6.64,  -1.51,   -8.02,   -9.47,
        -50.20,  -1.80,  -9.54,   1.52,  -49.45, -262.21,
          1.70,   3.43,   1.94,   1.37,    7.41,   20.62,
          4.14,   3.94,  -7.14,   1.37, -122.03,    1.02,
          2.89,  -7.30, 368.78,  50.01,   -1.08,    2.93,
          5.25,   3.95,  20.62,   4.09,    3.42,    1.69,
         11.29,   7.23,   1.51,   2.16,    1.38,    1.80,
          4.67,  16.01,  19.32,   1.30,   -1.02,   -4.51,
        120.99,   1.13,  22.98,   1.06,   -1.90,   -2.18,
        -23.58, 631.92,   1.92,  -4.66,  -17.86,    4.47,
          1.97,  17.20, 294.00,  -2.46,   -1.02,   79.96,
         23.83,   2.59,   4.47,   1.95,    1.17 };

    int j, k, m, n;
    double  anm[2][4][3],  bnm[2][4][3];
    double pinm, alpha, dt60, ap, bp, am, bm, p[3][2], q[3][2];

    //compute the time dependent potential matrix
    for (k = -1; k <=1; ++k)
    {
        dt60 = dmjd - k*dt - d1960;
        anm[0][1][k+1] = 0.0; anm[0][2][k+1] = 0.0;
        bnm[0][1][k+1] = 0.0; bnm[0][2][k+1] = 0.0;
        for (j = 0; j < nlines; ++j)
        {
            //n = 2;        //n = nj[j];
            m = mj[j];      //m = j<41?1:2;
            //pinm = fmod((n+m),2)*2*PI / 4.0;
            pinm = ((2+m)%2)*2*PI / 4.0;
            alpha = phase[j] + freq[j]*dt60 - pinm;
            alpha = fmod(alpha, 2*PI);
            //anm[n - 2][m][k + 1] = anm[n - 2][m][k + 1] + hs[j] * cos(alpha);
            //bnm[n - 2][m][k + 1] = bnm[n - 2][m][k + 1] - hs[j] * sin(alpha);
            anm[0][m][k + 1] = anm[0][m][k + 1] + hs[j] * cos(alpha);
            bnm[0][m][k + 1] = bnm[0][m][k + 1] - hs[j] * sin(alpha);
        }
    }
    //orthogonalize the response terms
    for (m = 1; m <= 2; ++m)
    {
        ap = anm[0][m][2] + anm[0][m][0];
        am = anm[0][m][2] - anm[0][m][0];
        bp = bnm[0][m][2] + bnm[0][m][0];
        bm = bnm[0][m][2] - bnm[0][m][0];
        p[0][m-1] = sp[0][m-1]*anm[0][m][1];
        p[1][m-1] = sp[1][m-1]*anm[0][m][1] - sp[2][m-1]*ap;
        p[2][m-1] = sp[3][m-1]*anm[0][m][1] - sp[4][m-1]*ap + sp[5][m-1]*bm;
        q[0][m-1] = sp[0][m-1]*bnm[0][m][1];
        q[1][m-1] = sp[1][m-1]*bnm[0][m][1] - sp[2][m-1]*bp;
        q[2][m-1] = sp[3][m-1]*bnm[0][m][1] - sp[4][m-1]*bp - sp[5][m-1]*am;
        anm[0][m][0] = p[0][m-1];
        anm[0][m][1] = p[1][m-1];
        anm[0][m][2] = p[2][m-1];
        bnm[0][m][0] = q[0][m-1];
        bnm[0][m][1] = q[1][m-1];
        bnm[0][m][2] = q[2][m-1];
    }

    //fill partials vector
    j = 0;
    for (n=2; n<=nmax; ++n)
        for (m=1; m<= n; ++m)
            for (k=-1; k<=1; k++) {
                h[j] = anm[n-2][m][k+1];
                h[j+1] = bnm[n-2][m][k+1];
                j = j + 2;
            }
    return;
}

//subroutines to compute the diurnal and semidiurnal variations
//in the earth orientation (x, y, UT1) from the version of Richard
//Ray's ocean tide model that was listed in IERS Technical Note
// 21, July 1996. This code includes the variations from 71 diurnal
//and semidiurnal terms instead of the 8 that are listed in the report.
//
//...Coded by : Richard Eanes, UT / CSR, Feb 1997
//
//...Input : time = Modified Julian Date
//
//...Output : dx, dy, dut1 = (delta_x, delta_y, delta_UT1)
//microarcsec for x and y, microsec for UT1
//(transferred to arcsec and time sec)
void ortho_eop(double time, double *dx, double *dy, double *dut1)
{
    static const double orthow[12][3]= {
       -6.77832, 14.86283, -1.76335, -14.86323, -6.77846,  1.03364,
        0.47884,  1.45234, -0.27553,  -1.45303,  0.47888,  0.34569,
        0.16406, -0.42056, -0.12343,   0.42030,  0.16469, -0.10146,
        0.09398, 15.30276, -0.47119,  25.73054, -4.30615,  1.28997,
       -4.77974,  0.07564, -0.19336,   0.28080,  2.28321,  0.02724,
        1.94539, -0.45717,  0.08955,  -0.73089, -1.62010,  0.04726 };

    int k, j;
    double eop[3], h[12];

    cnmtx(time, h);
    for (k = 0; k < 3; ++k) {
        eop[k] = 0.0;
        for (j = 0; j < 12; ++j) {
            eop[k] = eop[k] + h[j] * orthow[j][k];
        }
    }
    *dx = eop[0] * 1.0E-6;
    *dy = eop[1] * 1.0E-6;
    *dut1 = eop[2] * 1.0E-6;
}

// Form the matrix of polar motion, IAU 2000.
//
// Annex to IERS Conventions 2000, Chapter 5
//
// Given:
// XP, YP      d      coordinates of the pole(radians)
// SP         d      the quantity s' (radians)
//
// Returned :
// RPOM     d(3, 3)   polar - motion matrix
//
// The returned rotation matrix is the first to be applied when
// transforming a TRS vector into a CRS vector.
//
// revision : 2002 November 25
Matrix3d POM2000(double XP, double YP, double SP)
{
    //ROT_X(YP, RPOM);
    //ROT_Y(XP, RPOM);
    //ROT_Z(-SP, RPOM);
    Matrix3d R = rot_matrix(-SP, ZAXIS)*rot_matrix(XP, YAXIS)*rot_matrix(YP, XAXIS);
    return R;
}

// The quantity s', positioning the Terrestrial Ephemeris Origin on the
// equator of the Celestial Intermediate Pole.
//
// Annex to IERS Conventions 2000, Chapter 5
//
// Given:
// DATE1, DATE2    d      TT date(JD = DATE1 + DATE2)
//
// Returned :
// SP2000         d      the quantity s' in radians
double SP2000(double DATE1, double DATE2)
{
    // Arcseconds to radians
    const double DAS2R = 4.848136811095359935899141E-6;

    // Reference epoch(J2000), JD
    const double DJ0 = 2451545.0;

    // Days per Julian century
    const double DJC = 36525.0;

    // Time since J2000, in Julian centuries
    double T;

    // Interval between fundamental epoch J2000.0 and current date(JC).
    T = ((DATE1 - DJ0) + DATE2)/DJC;

    // Approximate S'.
    return -47E-6*T*DAS2R;
}


// Greenwich Sidereal Time(model consistent with IAU 2000 resolutions).
//
// Annexe to IERS Conventions 2000, Chapter 5
//
// Given:
// UTA, UTB     d      UT1 date(JD = UTA + UTB)
// TTA, TTB     d      TT date(JD = TTA + TTB)
// DPSI         d      nutation in longitude(radians)
//
// The result is the Greenwich Mean(apparent) Sidereal Time(radians),
// in the range 0 to 2pi.
//
// revision : 2002 December 9
double  GST2000(double UTA, double UTB, double TTA, double TTB, double DPSI)
{
    const double D2PI = 6.283185307179586476925287;
    // Greenwich Sidereal Time, IAU 2000.
    double ret = GMST2000(UTA, UTB, TTA, TTB) + EE2000(TTA, TTB, DPSI);
    ret = fmod(ret, D2PI);
    if (ret < 0.0)
        ret += D2PI;

    return ret;
}

// The equation of the equinoxes, compatible with IAU 2000
// resolutions, given the nutation in longitude.
//
// Annex to IERS Conventions 2000, Chapter 5
//
// Given:
// DATE1, DATE2    d      TT date(JD = DATE1 + DATE2)
// DPSI           d      nutation in longitude(radians)
//
// Returned :
// EE2000         d      equation of the equinoxes
//
// Calls the IERS routine EECT2000.
// revision : 2002 November 26
double EE2000(double DATE1, double DATE2, double DPSI)
{
    // Arcseconds to radians
    const double DAS2R = 4.848136811095359935899141E-6;

    // Reference epoch(J2000), JD
    const double DJ0 = 2451545.0;

    // Days per Julian century
    const double DJC = 36525.0;

    // J2000 obliquity(Lieske et al. 1977)
    const double EPS0 = 84381.448*DAS2R;


    // Interval between fundamental epoch J2000.0 and given date(JC).
    double T = ((DATE1 - DJ0) + DATE2)/DJC;

    // Mean obliquity from Chapter 5, expression(32).
    double EPSA = EPS0 + (-46.8402 + (-0.00059 + (0.001813)*T)*T)*T*DAS2R;

    // Equation of the equinoxes.
    return DPSI*cos(EPSA) + EECT2000(DATE1, DATE2);
}

double ANPM(double A)
{
    const double D2PI = 6.283185307179586476925287;
    const double DPI  = 3.141592653589793238462643;
    double W = fmod(A, D2PI);
    double sign = 1.0;
    if (A < 0.0)
        sign = -1.0;
    if (fabs(W) >= DPI)
        W -= D2PI*sign;
    return W;
}

// Equation of the equinoxes complementary terms,
// consistent with IAU 2000 resolutions.
//
// Annexe to IERS Conventions 2000, Chapter 5
//
// Given:
// DATE1, DATE2   d    TT date(JD = DATE1 + DATE2)
//
// Returned :
// EECT2000      d    complementary terms(radians)
// revision : 2003 May 4
double EECT2000(double DATE1, double DATE2)
{
    // Number of terms in the series
    const int NE0 = 33;
    const int NE1 = 1;

    const double D2PI = 6.283185307179586476925287;
    // Arcseconds to radians
    const double DAS2R = 4.848136811095359935899141E-6;

    // Reference epoch(J2000), JD
    const double DJ0 = 2451545;

    // Days per Julian century
    const double DJC = 36525;

    // Time since J2000, in Julian centuries
    double  T;

    // Miscellaneous
    double  A, S0, S1;
    int i, j;

    // Fundamental arguments
    double  FA[14];

    // ---------------------------------------- -
    // The series for the EE complementary terms
    // ---------------------------------------- -

    // Coefficients of l, l',F,D,Om,LMe,LVe,LE,LMa,LJu,LSa,LU,LN,pA
    int  KE1[14][NE1] = { 0 };

    // Sine and cosine coefficients for t ^ 1
    double SE1[2][NE1] = { -0.87E-6, +0.00E-6 };

    // Argument coefficients for t ^ 0
    int KE0[14][NE0]= {
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,
      0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2,  1,  0,  1,  0,
      0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,
      1,  1,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,
      1,  0,  0,  0,  0,  0,  0,  0,  2,  2,  2,  2,  2,  0,  0,
      0,  0,  0,  2,  2,  4,  1,  2,  2,  2,  2,  2, -2, -2,  0,
      0, -2,  0,  2,  0,  4,  2, -2, -2,  0,  0, -2, -2, -2,  0,
      0,  0,  0,  0,  0,  0, -2, -2, -4, -1,  0,  0,  0,  0, -2,
      2,  2,  0,  2,  0, -2, -2, -2, -2, -2,  0,  0,  1,  2,  3,
      1,  2,  3,  1,  3,  1, -1, -1,  1,  3,  1,  4,  1,  0,  2,
      3,  1,  0, -3, -1,  0,  0, -1,  1,  2, -1,  4,  4, -3, -1,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0, -8,  0,  0,  0,  0,  0,  0,  0,  8,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0, 12,  0,  0,  0,  0,  0,  0,  0,-13,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 };

    // Argument coefficients for t ^ 1
    KE1[4][0] = 1;

    // Sine and cosine coefficients for t ^ 0
    double SE0[2][NE0] = {
    2640.96E-6, +63.52E-6, +11.75E-6, +11.21E-6,  -4.55E-6,  +2.02E-6,  +1.98E-6,  -1.72E-6,
      -1.41E-6,  -1.26E-6,  -0.63E-6,  -0.63E-6,  +0.46E-6,  +0.45E-6,  +0.36E-6,  -0.24E-6,
      +0.32E-6,  +0.28E-6,  +0.27E-6,  +0.26E-6,  -0.21E-6,  +0.19E-6,  +0.18E-6,  -0.10E-6,
      +0.15E-6,  -0.14E-6,  +0.14E-6,  -0.14E-6,  +0.14E-6,  +0.13E-6,  -0.11E-6,  +0.11E-6,
      +0.11E-6,  -0.39E-6,  -0.02E-6,  +0.01E-6,  +0.01E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,
      +0.00E-6,  -0.01E-6,  -0.01E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,
      -0.12E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,
      +0.05E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,  +0.00E-6,
      +0.00E-6,  +0.00E-6 };

    // Interval between fundamental epoch J2000.0 and current date(JC).
    T = ((DATE1 - DJ0) + DATE2) / DJC;

    // Fundamental Arguments(from IERS Conventions 2000)

    // Mean Anomaly of the Moon.
    FA[0] = ANPM((485868.249036 +
            (715923.2178 + (31.8792 + (0.051635 + (-0.00024470)
            *T)*T)*T)*T)*DAS2R + fmod(1325*T, 1)*D2PI);

    // Mean Anomaly of the Sun.
    FA[1] = ANPM((1287104.793048 +
            (1292581.0481 + (-0.5532 + (+0.000136 + (-0.00001149)
            *T)*T)*T)*T)*DAS2R + fmod(99*T, 1)*D2PI);

    // Mean Longitude of the Moon minus Mean Longitude of the Ascending
    // Node of the Moon.
    FA[2] = ANPM((335779.526232 +
            (295262.8478 + (-12.7512 + (-0.001037 + (0.00000417)
            *T)*T)*T)*T)*DAS2R + fmod(1342*T, 1)*D2PI);

    // Mean Elongation of the Moon from the Sun.
    FA[3] = ANPM((1072260.703692 +
                (1105601.2090 + (-6.3706 + (0.006593 + (-0.00003169)
                *T)*T)*T)*T)*DAS2R + fmod(1236*T, 1)*D2PI);

    // Mean Longitude of the Ascending Node of the Moon.
    FA[4] = ANPM((450160.398036 +
            (-482890.5431 + (7.4722 + (0.007702 + (-0.00005939)
            * T) * T) * T) * T) * DAS2R + fmod(-5 * T, 1) * D2PI);

    FA[5] = ANPM(4.402608842 + 2608.7903141574 * T);
    FA[6] = ANPM(3.176146697 + 1021.3285546211 * T);
    FA[7] = ANPM(1.753470314 + 628.3075849991 * T);
    FA[8] = ANPM(6.203480913 + 334.0612426700 * T);
    FA[9] = ANPM(0.599546497 + 52.9690962641 * T);
    FA[10] = ANPM(0.874016757 + 21.3299104960 * T);
    FA[11] = ANPM(5.481293872 + 7.4781598567 * T);
    FA[12] = ANPM(5.311886287 + 3.8133035638 * T);
    FA[13] = (0.024381750 + 0.00000538691 * T) * T;

    // Evaluate the EE complementary terms.
    S0 = 0;
    S1 = 0;

    for (i = NE0; i >= 1; --i) {
        A = 0;
        for (j = 0; j < 14; ++j) {
            A += (double)(KE0[j][i-1])*FA[j];
        }
        S0 += SE0[0][i-1]*sin(A) + SE0[1][i-1]*cos(A);
    }

    for (i = NE1; i >= 1; --i) {
        A = 0;
        for (j = 0; j < 14; ++j) {
            A += (double)(KE1[j][i-1])*FA[j];
        }
        S1 += SE1[0][i-1]*sin(A) + SE1[1][i-1]*cos(A);
    }

    return (S0 + S1*T)*DAS2R;
}

// Earth rotation angle(IAU 2000 model).
// Annexe to IERS Conventions 2000, Chapter 5
//
// Given:
// DJ1, DJ2     d      UT1 date(JD = DJ1 + DJ2)
//
// The result is the Earth Rotation Angle(radians), in the range [0, 2pi]
// revision : 2003 May 4
double ERA2000(double DJ1, double DJ2)
{
    const double D2PI = 6.283185307179586476925287;

    // Reference epoch(J2000), JD
    double DJ0 = 2451545.0;

    double  D1, D2, T, F;
    double ret;

    // Days since fundamental epoch.
    if (DJ1 < DJ2) {
        D1 = DJ1;
        D2 = DJ2;
    }
    else {
        D1 = DJ2;
        D2 = DJ1;
    }

    T = D1 + (D2 - DJ0);

    // Fractional part of T(days).
    F = fmod(D1, 1.0) + fmod(D2, 1.0);

    // Earth rotation angle at this UT1.
    ret = D2PI * (F + 0.7790572732640 + 0.00273781191135448 * T);
    ret = fmod(ret, D2PI);
    if (ret < 0.0)
        ret += D2PI;
    return ret;
}

// purpose   : fundamental arguments of nutation
// parameters :
// rmjd_tdt -- mjd in TDT system
// arg      -- five nutation arguments
void fund_arg_nutation(double rmjd_tdt, double *arg)
{
    // IERS resolution 2003
    const double deg2rad = atan(1.0)/45.0;
    const double sec2rad = atan(1.0)/162000.0;
    static double  coef[5][5] = {
       134.96340251*deg2rad,   357.52910918*deg2rad,     93.27209062*deg2rad,    297.85019547*deg2rad,  125.04455501*deg2rad,
    1717915923.2178*sec2rad, 129596581.0481*sec2rad, 1739527262.8478*sec2rad, 1602961601.2090*sec2rad, -6962890.5431*sec2rad,
            31.8792*sec2rad,        -0.5532*sec2rad,        -12.7512*sec2rad,         -6.3706*sec2rad,        7.4722*sec2rad,
           0.051635*sec2rad,       0.000136*sec2rad,       -0.001037*sec2rad,        0.006593*sec2rad,      0.007702*sec2rad,
        -0.00024470*sec2rad,    -0.00001149*sec2rad,      0.00000417*sec2rad,     -0.00003169*sec2rad,   -0.00005939*sec2rad };

    // calculate the angles
    double ttc = (rmjd_tdt - 51544.50) / 36525.0;
    for (int i = 0; i < 5; ++i)
        arg[i] = coef[0][i] + (coef[1][i] + (coef[2][i] +
                 (coef[3][i] + coef[4][i] * ttc)*ttc)*ttc)*ttc;
}

// purpose  : compute tide variation on UT1 due to Earth rotation(IERS2003)
// parameter :
// mjd -- mod.julian date
// dut1 -- ut1 correction
// UT1 = UT1R + dut1
//
// notice : if ut1.table is UT1, dut1 must be subtracted for interpolation.
// dut1 must be added to the interpolated UT1R for coordinate trans.
// Zontal tide terms (IERS 2003)
double ut1ut1r(double mjd)
{
    static char mi[5][62] = {
     1, 2, 2, 0, 0, 1, 1, 1, 3,-1,-1, 1, 2, 0, 0, 0, 0, 2, 2, 2,
     0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1, 0, 1,-1,-1,-1, 1,
    -1, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 1, 2,-2,-1,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
     0, 0,-1, 0, 0, 0,-1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0,
     0, 0,-1, 2, 1, 1, 0, 0, 0, 2, 0, 0, 0,-1, 1,-1, 1, 1, 0, 0,
     0, 1, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 0, 2, 2, 2, 2,
     2, 0, 0, 0, 2, 0, 0, 0, 0, 2, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0,
     0, 0, 0,-2, 0, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 2, 0, 2, 0, 0,
     0,-2, 2, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 0, 0, 2, 2, 2,-2, 0,
     0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2,-2,-2, 0, 0, 0, 0, 0, 0, 0,
     1, 0, 2, 2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 0,-2,-2,-2,-2, 0,-2,
     0, 0,-1, 0, 0, 1, 0, 0, 2, 1, 2, 1, 2, 0, 1, 2, 0, 1, 2, 0,
     2, 2, 0, 1, 2,-1, 0, 1, 2,-1, 0, 1, 0, 1, 2, 0, 0, 1, 2,-1,
     0, 1, 0, 0,-1, 0, 1,-1, 0, 2, 1, 2, 0, 1, 2, 0,-1, 0, 1, 1,
    -1, 2, 0, 1, 0, 0, 1, 0, 2, 1 };

    static double sc[2][62] = {
     -0.02, -0.04, -0.10, -0.05, -0.12,   -0.04,  -0.41, -1.00,
     -0.02, -0.08, -0.20, -0.08,  0.02,    0.03,  -0.30, -3.22,
     -7.79,  0.02, -0.34,  0.02, -0.02,    0.05,  -0.74, -0.05,
     -0.05,  0.05,  0.10,  0.04,  0.05,    0.18,   0.44,  0.54,
     -8.33,  0.55,  0.05, -0.06,  0.12,   -1.84,   0.13,  0.02,
     -0.09, -0.06,  0.03, -1.91,  0.26,    1.18, -49.06, -0.20,
      0.05, -0.56,  0.04, -0.05,  0.09,    0.82, -15.65, -0.14,
      0.03, -0.14,  0.43, -0.04,  8.20,-1689.54,   0.00,  0.00,
      0.00,  0.00,  0.00,  0.00,  0.00,    0.01,   0.00,  0.00,
      0.00,  0.00,  0.00,  0.00,  0.00,    0.02,   0.05,  0.00,
      0.00,  0.00,  0.00,  0.00,  0.00,    0.00,   0.00,  0.00,
      0.00,  0.00,  0.00,  0.00,  0.00,    0.00,   0.06,  0.00,
      0.00,  0.00,  0.00,  0.01,  0.00,    0.00,   0.00,  0.00,
      0.00,  0.02,  0.00, -0.01,  0.43,    0.00,   0.00,  0.01,
      0.00,  0.00,  0.00, -0.01,  0.15,    0.00,   0.00,  0.00,
     -0.01,  0.00,  0.11,-25.04 };

    int i, j;
    double argu[5], arg;
    // get arguement
    fund_arg_nutation(mjd, argu);
    double dut1 = 0.0;
    for (i=0; i<62; ++i) {
        arg = 0.0;
        for (j=0; j<5; ++j)
            arg += mi[j][i]*argu[j];
        dut1 += sc[0][i]*sin(arg) + sc[1][i]*cos(arg);
    }
    // 0.1 mimiseconds to rad.
    dut1 *= 1E-4;
    return dut1;
}


#if 0
void timinc(int jd, double sec, double delt, int *jd1, double *sec1)
{
    int inc, nss = 86400;

    *sec1 = sec + delt;
    inc = (int)(*sec1 / nss);
    *sec1 -= inc*nss;
    *jd1 = jd + inc;
    if (*sec1 >= 0)
        return;
    --*jd1;
    *sec1 = nss + *sec1;
}

void ROT_X(double PHI, double R[][3])
{
    // Matrix representing new rotation.
    double S = sin(PHI);
    double C = cos(PHI);

    double A[3][3] = {0};
    for (int i = 0; i < 3; ++i)
            A[i][i] = 1.0;

    A[1][1] = C;
    A[2][1] = -S;
    A[1][2] = S;
    A[2][2] = C;

    // Rotate.
    MATMPY(A[0], R[0], R[0], 3, 3, 3);
}

void ROT_Y(double THETA, double R[][3])
{
    // Matrix representing new rotation.
    double S = sin(THETA);
    double C = cos(THETA);

    double A[3][3] = { 0 };
    for (int i = 0; i < 3; ++i)
            A[i][i] = 1.0;

    A[0][0] = C;
    A[2][0] = S;
    A[0][2] = -S;
    A[2][2] = C;

    // Rotate.
    MATMPY(A[0], R[0], R[0], 3, 3, 3);
}

void ROT_Z(double PSI, double R[][3])
{
    // Matrix representing new rotation.
    double S = sin(PSI);
    double C = cos(PSI);

    double A[3][3] = {0};
    for (int i = 0; i < 3; ++i)
            A[i][i] = 1.0;

    A[0][0] = C;
    A[1][0] = -S;
    A[0][1] = S;
    A[1][1] = C;

    // Rotate.
    MATMPY(A[0], R[0], R[0], 3, 3, 3);
}

void unit_vector(int n, double *v, double *u, double *length)
{
    *length = 0;
    for (int i = 0; i < n; ++i)
        *length += v[i] * v[i];

    *length = sqrt(*length);
    for (int i = 0; i < n; ++i)
        u[i] = v[i] / *length;
}

void MATMPY(double *A, double *B, double *C, int NROW, int NCOLA, int NCOLB)
{
    // REAL * 8  A(NROW, NCOLA), B(NCOLA, NCOLB), C(NROW, NCOLB);

    // ALLOCATE(C1(1:NROW, 1 : NCOLB));
    // enable B, C to be the same array
    double *C1 = (double*)malloc(NROW*NCOLB*sizeof(double));
    int i,j,k;
    for (i = 0; i < NROW; ++i) {
        for (j = 0; j < NCOLB; ++j) {
            double sum = 0.0;
            for (k = 0; k < NCOLA; ++k)
                sum += *(A + i*NCOLA + k) * *(B + k*NCOLB + j);
            *(C1 + i*NCOLB + j) = sum;
        }
    }

    for (i = 0; i < NROW*NCOLB; ++i)
        *(C + i) = *(C1 + i);

    free(C1);
}

void cross(double *v1, double *v2, double *vout)
{
    vout[0] = v1[1]*v2[2] - v1[2]*v2[1];
    vout[1] = v1[2]*v2[0] - v1[0]*v2[2];
    vout[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

double dot(int n, double *v1, double *v2)
{
    double ret = 0.0;
    for (int i = 0; i < n; ++i)
        ret += v1[i] * v2[i];

    return ret;
}
#endif
