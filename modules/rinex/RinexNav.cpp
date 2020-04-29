#include "RinexNav.h"
#include "../io/io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

bool RinexNav::Header::read(FILE *fp)
{
    char buf[256];
    fgets(buf, sizeof(buf), fp);
    if (strlen(buf) == 0 || strncmp(buf+60, "RINEX VERSION / TYPE", 20) || buf[20] != 'N')
    {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexNav::Header::read: Not a RINEX-NAV file\n");
        return false;
    }
    version = atof(buf+5);
    if (version < 2.00 || version > 3.03) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexNav::Header::read: unknown version: %4.2f\n", version);
        return false;
    }

    while (fgets(buf, sizeof(buf), fp))
    {
        if (strncmp(buf+60, "END OF HEADER", 13) == 0) {
            break;
        }
        else if (strncmp(buf+60, "ION ALPHA", 9) == 0) {
            buf[10] = 'E'; buf[22] = 'E';
            buf[34] = 'E'; buf[46] = 'E';
            sscanf(buf, "%lf%lf%lf%lf", ionoCoefs, ionoCoefs+1, ionoCoefs+2, ionoCoefs+3);
        }
        else if (strncmp(buf+60, "ION BETA", 8) == 0) {
            buf[10] = 'E'; buf[22] = 'E';
            buf[34] = 'E'; buf[46] = 'E';
            sscanf(buf, "%lf%lf%lf%lf", ionoCoefs+4, ionoCoefs+5, ionoCoefs+6, ionoCoefs+7);
        }
        else if (strncmp(buf+60, "DELTA-UTC: A0,A1,T,W", 20) == 0) {
            buf[18] = 'E'; buf[37] = 'E';
            sscanf(buf, "%lf%lf%d%d", &A0, &A1, &T, &W);
        }
        else if (strncmp(buf+60, "LEAP SECONDS", 12) == 0) {
            leapSec = atoi(buf);
        }
    }
    return true;
}

bool RinexNav::Header::write(FILE *fp)
{
    for (int i=0; i<8; ++i) {
        fprintf(fp, "%13.4E", ionoCoefs[i]);
        if (i%4==0) fprintf(fp, "\n");
    }
    fprintf(fp, "%20.12E%20.12E%9d%9d\n", A0, A1, T, W);

    return true;
}

const Ephemeris *RinexNav::eph(MJD t, const std::string &prn)const
{
    auto it = ephs_.find(prn);
    if (it == ephs_.end())
        return nullptr;

    //auto targ = std::lower_bound(it->second.begin(), it->second.end(), t,
    //            [](const Ephemeris &eph, MJD t) { return eph.toe < t; });
    //if (targ == it->second.begin()) {
    //    return &*(it->second.begin());
    //}
    //else if (targ == it->second.end()) {
    //    return &*(--it->second.end());
    //}
    //else {
    //    auto candi = --targ;
    //    double diff[2] = { fabs(t-targ->toe), fabs(t-candi->toe) };
    //    return diff[0]<diff[1] ? &*targ : &*candi;
    //}

    double diff;
    double minVal = 86400.0*7.0;
    auto i =  it->second.begin();
    auto target = i;
    for (; i != it->second.end(); ++i)
    {
        diff = fabs(t-i->toe);
        if (diff < minVal) {
            minVal = diff;
            target = i;
        }
    }
    return &*target;
}

bool RinexNav::read(const std::string &path)//, const std::vector<std::string> &prns)
{
    FILE *navFile = fopen(path.c_str(), "r");
    if (navFile == nullptr) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexNav::read: no such file: %s\n", path.c_str());
        return false;
    }
    if (!header_.read(navFile))
        return false;

    char tp = path.back();
    //FILE *fp = fopen("nav", "w");
    //if (fp == nullptr) {
    //    fprintf(stderr, "error: fail to create nav\n");
    //}
    if (header_.version > 2.99) tp = 'p';
    else if (header_.version > 1.99) tp = 'n';
    Ephemeris eph;
    while (eph.read(navFile, tp))
    {
        //if (std::binary_search(prns.begin(), prns.end(), eph.prn)) {
        if (eph.toc.d != 0) {
            ephs_[eph.prn].push_back(eph);
            //eph.write(fp);
        }
        //}
    }
    //fclose(fp);

    //for (auto it=ephs_.begin(); it !=ephs_.end(); ++it) {
    //    std::sort(it->second.begin(), it->second.end(),
    //              [](const Ephemeris &lhs, const Ephemeris &rhs)
    //              { return lhs.toe < rhs.toe; });
    //}
    return true;
}

void RinexNav::repeat_time(FILE *fp)const
{
    for (auto &pair : ephs_)
    {
        int count = 0;
        double avg = 0;
        std::string prn = pair.first;
        for (auto &eph : pair.second)
        {
            ++count;
            avg += eph.orbit_period();
            fprintf(fp, "%3s %8.2f\n", prn.c_str(), eph.orbit_period());
        }
        avg/=count;

        double sum = 0;
        for (auto &eph : pair.second)
        {
            double t = eph.orbit_period();
            if (fabs(avg-t) > 3) {
                --count;
                continue;
            }
            sum += t;
        }
        avg = sum/count;
        fprintf(fp, "%3s avg => %8.2f\n", prn.c_str(), avg);

        char sys = prn.front();
        if (sys == 'G')
            fprintf(fp, "%3s timeshift => %4.0f\n", prn.c_str(), avg*2-86400);
        else if (sys == 'C') {
            int No = atoi(prn.c_str()+1);
            if (No<11||No==13)
                fprintf(fp, "%3s timeshift => %4.0f\n", prn.c_str(), avg-86400);
            else
                fprintf(fp, "%3s timeshift => %4.0f\n", prn.c_str(), avg*13-86400*7);
        }
        else if (sys == 'J')
            fprintf(fp, "%3s timeshift => %4.0f\n", prn.c_str(), avg-86400);
        else if (sys == 'E') {
            int No = atoi(prn.c_str()+1);
            if (No==14 || No==18) {
                //fprintf(fp, "%3s timeshift => %4.0f\n", avg*13-86400*7);
            }
            else
                fprintf(fp, "%3s timeshift => %4.0f\n", prn.c_str(), avg*17-86400*10);
        }
    }
}

#if 0
bool RinexNav::Header::read(std::istream &in)
{
    std::string buf, substr;
    std::getline(in, buf);
    if (buf.empty() || buf.substr(60, 20) != "RINEX VERSION / TYPE" || buf[20] != 'N')
    {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexNav::Header::read: Not a RINEX-NAV file\n");
        return false;
    }

    while (std::getline(in, buf))
    {
        substr = buf.substr(60, buf.find_last_not_of(' ')-59);
        if (substr == "END OF HEADER") {
            break;
        }
        else if (substr == "ION ALPHA") {
            buf[10] = 'E'; buf[22] = 'E';
            buf[34] = 'E'; buf[46] = 'E';
            for (int i=0; i<4; ++i)
                ionoCoefs[i] = atof(buf.c_str() + 2 + 12*i);
        }
        else if (substr == "ION BETA") {
            buf[10] = 'E'; buf[22] = 'E';
            buf[34] = 'E'; buf[46] = 'E';
            for (int i=4; i<8; ++i)
                ionoCoefs[i] = atof(buf.c_str() + 2 + 12*i);
        }
        else if (substr == "DELTA-UTC: A0,A1,T,W") {
            buf[18] = 'E'; buf[37] = 'E';
            A0 = atof(buf.c_str() + 4);
            A1 = atof(buf.c_str() + 23);
            T  = atoi(buf.c_str() + 41);
            W  = atoi(buf.c_str() + 50);
        }
        else if (substr == "LEAP SECONDS") {
            leapSec = atoi(buf.c_str());
        }
    }

    return true;
}
#endif
