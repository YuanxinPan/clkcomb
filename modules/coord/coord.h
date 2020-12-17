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

void SunPosition(int jd, double sod, double rs[]);
void MoonPosition(int jd, double sod, double rm[]);

void nominal_att(const double *xsat, const double *xsun, double *q);

#endif
