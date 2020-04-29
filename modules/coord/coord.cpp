#include "coord.h"

#include <pppx/const.h>
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
