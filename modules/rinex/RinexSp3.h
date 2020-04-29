#ifndef RINEXSP3_H
#define RINEXSP3_H

#include "../coord/coord.h"
#include "../chrono/chrono.h"

#include <vector>
#include <string>

class RinexSp3
{
public:
    struct Sp3_t {
        MJD t;
        CartCoor pos;
    };

public:
    RinexSp3(){}
    ~RinexSp3(){}

    //static std::vector<std::string> satPRN(const std::string &path);
    bool read(const std::string &path);
    bool read(const std::vector<std::string> &paths);
    bool satPos(MJD t, const std::string &prn, double *pos)const;
    bool satPos(MJD t, const std::string &prn, CartCoor &pos)const;
    bool satPosVel(MJD t, const std::string &prn, double *pos, double *vel)const;
    bool satPosVel(MJD t, const std::string &prn, CartCoor &pos, CartCoor &vel)const;

    const std::vector<std::string> &sat_prns() { return prns_; }

private:
    std::vector<std::string> prns_;
    std::vector<std::vector<Sp3_t>> satposs_;
};

bool operator<(const RinexSp3::Sp3_t &lhs, const RinexSp3::Sp3_t &rhs);

#endif
