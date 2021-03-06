#include "Ephemeris.h"
#include "../const.h"
#include "../io/io.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <algorithm>

Ephemeris::Ephemeris() {}

Ephemeris::~Ephemeris() {}

double Ephemeris::orbit_period()const
{
    double n0 = sqrt(GM) / pow(sqrtA, 3);
    double n = n0 + dn;
    return PI*2/n;
}

double Ephemeris::calcEt(double Mt)const
{
    int iter = 0;
    double Et = Mt, tmp;
    do{
        tmp = Et;
        Et = Mt + e*sin(Et);
    } while (fabs(Et - tmp) > 1E-13 && ++iter<9);
    return Et;
}

bool Ephemeris::satPosClk(MJD t, CartCoor &pos, double &sck)const
{
    double n0 = sqrt(GM) / pow(sqrtA, 3);
    double n = n0 + dn;
    double dt = t-toe;    // time interval
    //if (dt > 302400.0)
    //    dt -= 86400.0*7;
    //else if (dt < -302400.0)
    //    dt += 86400.0*7;

    double Mt = M0 + n*dt;
    double Et = calcEt(Mt);

    //double ft = atan(sqrt(1 - e*e)*sin(Et) / (cos(Et) - e));
    double cosf = cos(Et) - e;
    double sinf = sqrt(1 - e*e) * sin(Et);
    double ft = atan2(sinf, cosf);

    double ut_ = omega + ft;
    double rt_ = pow(sqrtA, 2)*(1 - e*cos(Et));

    double ut = ut_ + Cuc*cos(2 * ut_) + Cus*sin(2 * ut_);
    double rt = rt_ + Crc*cos(2 * ut_) + Crs*sin(2 * ut_);
    double it = i0 + iDot*dt + Cic*cos(2 * ut_) + Cis*sin(2 * ut_);

    double xt = rt*cos(ut);
    double yt = rt*sin(ut);

    double Lt = OMEGA + (omegaDot - We)*(dt /*+ TOE*/) - We*TOE;
    // L = omig0 + omig1*(t - toe) - omigae*t;

    pos.x = xt*cos(Lt) - yt*cos(it)*sin(Lt);
    pos.y = xt*sin(Lt) + yt*cos(it)*cos(Lt);
    pos.z = yt*sin(it);

    // double deltaTr = F*e*sqrtA*sin(Et);
    sck = clkBias + dt*clkDrift + dt*dt*clkDriftRate; // + deltaTr;// - TGD; dual freq obs
    return true;
}

bool Ephemeris::satPos(MJD t, CartCoor &pos)const
{
    double n0 = sqrt(GM) / pow(sqrtA, 3);
    double n = n0 + dn;
    double dt = t-toe;    // time interval

    double Mt = M0 + n*dt;
    double Et = calcEt(Mt);

    double cosf = cos(Et) - e;
    double sinf = sqrt(1 - e*e) * sin(Et);
    double ft = atan2(sinf, cosf);

    double ut_ = omega + ft;
    double rt_ = pow(sqrtA, 2)*(1 - e*cos(Et));

    double ut = ut_ + Cuc*cos(2 * ut_) + Cus*sin(2 * ut_);
    double rt = rt_ + Crc*cos(2 * ut_) + Crs*sin(2 * ut_);
    double it = i0 + iDot*dt + Cic*cos(2 * ut_) + Cis*sin(2 * ut_);

    double xt = rt*cos(ut);
    double yt = rt*sin(ut);

    double Lt = OMEGA + (omegaDot - We)*(dt /*+ TOE*/) - We*TOE;

    pos.x = xt*cos(Lt) - yt*cos(it)*sin(Lt);
    pos.y = xt*sin(Lt) + yt*cos(it)*cos(Lt);
    pos.z = yt*sin(it);
    return true;
}

bool Ephemeris::satClk(MJD t, double &sck)const
{
    // double n0 = sqrt(GM) / pow(sqrtA, 3);
    // double n = n0 + dn;
    double dt = t-toe;    // time interval

    // double Mt = M0 + n*dt;
    // double Et = calcEt(Mt);
    // double deltaTr = F*e*sqrtA*sin(Et);
    sck = clkBias + dt*clkDrift + dt*dt*clkDriftRate; // + deltaTr // - TGD; dual freq obs
    return true;
}

bool Ephemeris::write(FILE *fp)const
{
    fprintf(fp, "%3s %20.12E%20.12E%20.12E\n", prn.c_str(), /*toc, */clkBias, clkDrift, clkDriftRate);
    fprintf(fp, "    %20.12E%20.12E%20.12E%20.12E\n", IODE, Crs, dn, M0);
    fprintf(fp, "    %20.12E%20.12E%20.12E%20.12E\n", Cuc, e, Cus, sqrtA);
    fprintf(fp, "    %20.12E%20.12E%20.12E%20.12E\n", /* toe*/0.0, Cic, OMEGA, Cis);
    fprintf(fp, "    %20.12E%20.12E%20.12E%20.12E\n", i0, Crc, omega, omegaDot);
    fprintf(fp, "    %20.12E%20.12E%20.12E%20.12E\n", iDot, codesOnL2, GPSWeek, L2PDataFlag);
    fprintf(fp, "    %20.12E%20.12E%20.12E%20.12E\n", SVAccuracy, SVHealth, TGD, IODC);
    fprintf(fp, "    %20.12E%20.12E\n", transTimeOfMsg, fitInterval);
    return true;
}

static bool parse(char *buf, double *val)
{
    buf[18]='E'; buf[37]='E';
    buf[56]='E'; buf[75]='E';
    val[0] = atof(buf + 3);
    val[1] = atof(buf + 22);
    val[2] = atof(buf + 41);
    val[3] = atof(buf + 60);
    return true;
}

bool Ephemeris::read(FILE *fp, char tp)
{
    char buf[128];
    if (!fgets(buf, sizeof(buf), fp))
        return false;

    int shift = 0;  // field shift diff between n/p RINEX-NAV (n:0, p:1)
    CalendTime t;
    if (tp == 'n') {
        prn = 'G';
        prn.append(buf, 2);
        t.read(buf+3);
    }
    else if (tp == 'p') {
        shift = 1;
        prn.assign(buf, 3);
        t.read(buf+6);
        if (prn[0] == 'I') {
            skip_nline(fp, 8-1);
            return true;
        }
        else if (prn[0]=='R' || prn[0]=='S') {
            skip_nline(fp, 4-1);
            return true;
        }
    }
    if (prn[1] == ' ') prn[1] = '0';
    toc = t;

    // cvt D to E, otherwise atof() failed on Linux
    buf[37+shift]='E'; buf[56+shift]='E'; buf[75+shift]='E';
    sscanf(buf+22+shift, "%lf%lf%lf", &clkBias, &clkDrift, &clkDriftRate);

    double *p[6] = { &IODE, &Cuc, &TOE, &i0, &iDot, &SVAccuracy };
    for (int i=0; i<6; ++i) {
        fgets(buf, sizeof(buf), fp);
        parse(buf+shift, p[i]);
    }

    fgets(buf, sizeof(buf), fp);
    buf[18+shift]='E'; buf[37+shift]='E';
    //buf[56+shift]='E'; buf[75+shift]='E';
    transTimeOfMsg = atof(buf + shift + 3);
    fitInterval    = atof(buf + shift + 22);

    if(!quality_check())
        toc.d = 0;

    int dWeek = prn[0]=='C' ? 1356 : 0;
    double dSow = prn[0]=='C' ? 14 : 0;
    toe.set(static_cast<int>(GPSWeek)+dWeek, TOE+dSow);  // for convinence to calculate dt
    return true;
}

bool Ephemeris::quality_check()const
{
    int empty = 0;
    const double *q = &IODE;
    for (int i=0; i<20; ++i) {
        if (fabs(*q++) < DBL_EPSILON)
            ++empty;
    }
    return empty>5 ? false:true;
}
