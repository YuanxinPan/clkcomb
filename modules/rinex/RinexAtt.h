#ifndef RINEXATT_H
#define RINEXATT_H

#include "../chrono/chrono.h"

#include <map>
#include <string>
#include <vector>

class RinexAtt
{
private:
    struct Coef_t
    {
        MJD t;
        double q[4]; // ECEF => SV
        bool operator<(MJD mjd) { return t<mjd; }
    };

public:
    RinexAtt(): attFile_(nullptr), interval_(0) {}
    ~RinexAtt() {}

    bool read(const std::string &path);
    bool satAtt(MJD t, const std::string &prn, double *R); // SV => ECEF
    void close();
    double interval()const { return interval_; }
    double version()const { return version_; }

private:
    FILE *attFile_;
    double version_;
    double interval_;
    std::string creator_;
    std::map<std::string, std::string> svTypes_;
    std::vector<std::string> prns_;
    std::vector<std::vector<Coef_t>> atts_; // prn, epoch
};

// bool operator<(const RinexAtt::Coef_t &lhs, const MJD t) { return lhs.t < t; }

#endif // RINEXATT_H
