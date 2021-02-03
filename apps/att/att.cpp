#include <pppx/const.h>
#include <pppx/rinex.h>
#include <pppx/coord.h>
#include <math.h>

double q_diff(const double *q1, const double*q2)
{
    double dot=0;
    for (int i=0; i<4; ++i)
        dot += q1[i]*q2[i];
    double theta = acos(2*dot*dot-1);
    return theta;
}

void qcross(const double *q1, const double *q2, double *q)
{
    q[0] = q1[0]*q2[0] - dot(q1+1, q2+1);
    double c[3];
    cross(q2+1, q1+1, c);
    for (int i=0; i<3; ++i)
        q[i+1] = q2[0]*q1[i+1] + q1[0]*q2[i+1] + c[i];
}

void qinv(const double *q, double *qi)
{
    qi[0] =  q[0];
    qi[1] = -q[1];
    qi[2] = -q[2];
    qi[3] = -q[3];
}

int main(int argc, char *argv[])
{
    if (argc == 1) {
        fprintf(stdout, "att srcatt refatt sp3 mjd prn\n");
        return 0;
    }

    RinexAtt srcatt;
    RinexAtt refatt;
    RinexSp3 rnxsp3;

    if (!srcatt.read(argv[1]) )
        return 1;

    bool have_ref = true;
    if (!refatt.read(argv[2]))
        have_ref = false;
    if (!have_ref && !rnxsp3.read(argv[3]))
        return 1;

    MJD t(atoi(argv[4]), 0);
    std::string prn(argv[5]);
    // double qs[4], qr[4], qi[4], qc[4];
    double diff[2] = { 0 };
    //
    double interval = 30;
    // MJD t(56658, 0);
    // std::string prn("G03");
    double qs[4], qr[4]={0}, qn[4]={0};

    for (int iepo=0; iepo!=2850; ++iepo) 
    { 
        t.sod = interval*iepo;
        srcatt.sat_att(t, prn, qs);

        if (have_ref)
            refatt.sat_att(t, prn, qr);
        else {
            double xsat[3], xsun[3];
            rnxsp3.satPos(t, prn, xsat);
            SunPosition(t.d, t.sod, xsun);
            nominal_att(xsat, xsun, qn);
        }

        //printf("nominal %4d %20.16f %20.16f %20.16f %20.16f\n", iepo, qn[0], qn[1], qn[2], qn[3]);
        //printf("src     %4d %20.16f %20.16f %20.16f %20.16f\n", iepo, qs[0], qs[1], qs[2], qs[3]);
        //printf("ref     %4d %20.16f %20.16f %20.16f %20.16f\n", iepo, qr[0], qr[1], qr[2], qr[3]);

        if((dot(qn+1,qs+1)+qn[0]*qs[0])<0) {
            qn[0] = -qn[0];
            qn[1] = -qn[1];
            qn[2] = -qn[2];
            qn[3] = -qn[3];
        }

        double qi[4], qc[4];
        qinv(qs, qi);
        if (have_ref)
            qcross(qi, qr, qc);
        else
            qcross(qi, qn, qc);
        diff[1] = 2*atan2(qc[3], qc[0])*R2D;

    //     qinv(qs, qi);
    //     qcross(qi, qr, qc);
    //     // norm = sqrt(dot(qc+1, qc+1));
    //     // diff = 2*atan2(norm, qc[0])*R2D;
    //     diff[1] = 2*atan2(qc[3], qc[0])*R2D;
        while (diff[1] - diff[0] > 180) {
            diff[1] -= 360;
        }
        while (diff[1] - diff[0] < -180) {
            diff[1] += 360;
        }

    //     double pwu = diff[1]/360*1E12/(GPS_f1 + GPS_f2);
    //     fprintf(stdout, "%8.0f %8.3f %8.3f\n", t.sod, diff[1], pwu);
    //
        diff[0] = diff[1];
        printf("%4d %8.3f\n", iepo, diff[1]);
    }

    return 0;
}
