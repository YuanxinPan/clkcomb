#include "RinexAtt.h"
#include "../io/io.h"
#include "../const.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

bool RinexAtt::read(const std::string &path)
{
    attFile_ = fopen(path.c_str(), "r");
    if (attFile_ == nullptr) {
        fprintf(stderr, MSG_ERR "RinexAtt::open: no such file: %s\n", path.c_str());
        return false;
    }

    Coef_t coef;
    int y, m, d, h, min, n;
    double sec;
    std::string prn;
    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), attFile_))
    {
        if (strncmp(buf, " ATT ", 5) == 0) {
            prn.assign(buf+5, 3);
            n = sscanf(buf+24, "%lf %lf %lf %lf", coef.q, coef.q+1, coef.q+2, coef.q+3);
            if (n!=4) {
                fprintf(stderr, MSG_ERR "RinexAtt::read:: %s", buf);
                return false;
            } else if (coef.q[0] == 0.0) {
                buf[9] = '\0';
                fprintf(stderr, MSG_WAR "## %4d %2d %2d %2d %2d %8.2f %s empty\n", y, m, d, h, min, sec, buf);
                continue;
            }
            if (!std::binary_search(prns_.begin(), prns_.end(), prn)) {
                fprintf(stderr, MSG_WAR "## %4d %2d %2d %2d %2d %8.2f %3s not in SATELLITE LIST\n", y, m, d, h, min, sec, prn.c_str());
                continue;
            }

            // Special case for JPL BLOCK IIR satellites
            if (creator_ == "JPL" && svTypes_[prn] == "IIR") {
                std::swap(coef.q[0], coef.q[3]);
                std::swap(coef.q[1], coef.q[2]);
                coef.q[0] = -coef.q[0];
                coef.q[2] = -coef.q[2];
            }

            int iprn = std::lower_bound(prns_.begin(), prns_.end(), prn) - prns_.begin();
            atts_[iprn].push_back(coef);
            // printf(" %3s %19.16f %19.16f %19.16f %19.16f\n", prn.c_str(), coef.q[0], coef.q[1], coef.q[2], coef.q[3]);
        } else if (strncmp(buf, "## ", 3) == 0) {
            sscanf(buf+3, "%d %d %d %d %d %lf", &y, &m, &d, &h, &min, &sec);
            coef.t.d = date2mjd(y, m, d);
            coef.t.sod = hms2sod(h, min, sec);
            // printf("## %4d %2d %2d %2d %2d %10.2f\n", y, m, d, h, min, sec);
        } else if (strncmp(buf, "%=ORBEX", 7) == 0) {
            version_ = atof(buf+7);
            // printf("VERSION %6.3f\n", version_);
        } else if (strncmp(buf, " CREATED_BY", 11) == 0) {
            creator_.assign(buf+21, 3);
            // printf("CREATOR %s\n", creator_.c_str());
        } else if (strncmp(buf, " EPOCH_INTERVAL", 15) == 0) {
            interval_ = atof(buf+21);
            // printf("INTERVAL %8.3f\n", interval_);
        } else if (strncmp(buf, "+SATELLITE/ID_AND_DESCRIPTION", 29) == 0) {
            while (fgets(buf, sizeof(buf), attFile_)) {
                if (strncmp(buf, "-SATELLITE/ID_AND_DESCRIPTION", 29) == 0)
                    break;
                else if (buf[0] != ' ')
                    continue;
                prn.assign(buf+1, 3);
                prns_.push_back(prn);
                if (creator_ == "JPL") {
                    svTypes_[prn] = std::string(buf+13, 3);
                    // printf("SATELLITE %3s %s\n", prns_.back().c_str(), svTypes_[prn].c_str());
                }
            }
            std::sort(prns_.begin(), prns_.end());
            atts_.resize(prns_.size());
            // fprintf(stderr, "PRN SIZE %lu\n", prns_.size());
            // for (auto it=prns_.begin(); it!=prns_.end(); ++it)
            //     fprintf(stderr, "%3s\n", it->c_str());
        }
    }

    // for (size_t i=0; i!=prns_.size(); ++i) {
    //     fprintf(stderr, "%3s %8lu\n", prns_[i].c_str(), atts_[i].size());
    // }

    return true;
}

void RinexAtt::close()
{
    if (attFile_ != nullptr) {
        fclose(attFile_);
        attFile_ = nullptr;
    }
}

static void quatern_interp(const double *qs, const double *qe, double dist, double *q)
{
    double cosa = 0, sign = 1;
    for (int i=0; i<4; i++) {
        cosa += qs[i]*qe[i];
    }
    if (cosa < 0.0) {
        cosa = -cosa;
        sign = -1;
    }

    double k0, k1;
    if (cosa > 0.9995) { // linear interpolation
        k0 = 1 - dist;
        k1 = dist;
    } else { // spherical linear interpolation
        double sina = sqrt(1 - cosa*cosa);
        double a = atan(sina / cosa);
        k0 = sin((1-dist)*a)/sina;
        k1 = sin(dist*a)/sina;
    }

    double norm = 0;
    for (int i=0; i<4; ++i) {
        q[i] = qs[i]*k0 + sign*qe[i]*k1;
        norm += q[i]*q[i];
    }
    norm = sqrt(norm);
    for (int i=0; i<4; ++i)
        q[i] /= norm;
}

void quatern2rotmat(const double *q, double *R)
{
    double q00 = q[0]*q[0], q01 = q[0]*q[1], q02 = q[0]*q[2], q03 = q[0]*q[3];
    double q11 = q[1]*q[1], q12 = q[1]*q[2], q13 = q[1]*q[3];
    double q22 = q[2]*q[2], q23 = q[2]*q[3];
    double q33 = q[3]*q[3];
    
    // column major
    // row 1
    R[0+3*0] = q00 + q11 - q22 - q33;
    R[0+3*1] = 2*(q12 - q03);
    R[0+3*2] = 2*(q13 + q02);

    // row 2
    R[1+3*0] = 2*(q12 + q03);
    R[1+3*1] = q00 - q11 + q22 - q33;
    R[1+3*2] = 2*(q23 - q01);

    // row 3
    R[2+3*0] = 2*(q13 - q02);
    R[2+3*1] = 2*(q23 + q01);
    R[2+3*2] = q00 - q11 - q22 + q33;
}

bool RinexAtt::sat_att(MJD t, const std::string &prn, double *_q)const
{
    if (!std::binary_search(prns_.begin(), prns_.end(), prn))
        return false;

    int iprn = std::lower_bound(prns_.begin(), prns_.end(), prn) - prns_.begin();
    auto it = std::lower_bound(atts_[iprn].begin(), atts_[iprn].end(), t);
    if (it == atts_[iprn].end())
        return false;

    double q[4];
    if (fabs(it->t - t) < 1E-3) { // 1 ms
        memcpy(q, it->q, sizeof(q));
    } else {
        if (it+1 == atts_[iprn].end())
            return false;
        double dist = (t - it->t)/((it+1)->t - it->t);
        quatern_interp(it->q, (it+1)->q, dist, q);
    }

    // quatern2rotmat(q, R); // ECEF => SV
    memcpy(_q, q, sizeof(q));
    return true;
}
