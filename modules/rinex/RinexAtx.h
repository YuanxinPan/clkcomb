#ifndef RINEXATX_H
#define RINEXATX_H

#include "../chrono/chrono.h"

#include <array>
#include <string>
#include <vector>
#include <stdio.h>

class RinexAtx
{
public:
    struct atx_t
    {
    private:
        friend class RinexAtx;
        std::string name;
        std::string svn_;
        std::string svType_;
        int nfreq;
        std::vector<std::string> freq;
        std::vector<std::array<double,3>> pcos;  // nfreq
        double dazi;
        double dzen, zen1, zen2;
        std::vector<std::vector<std::vector<double>>> pcvs;   // nfreq, nazi, nzen
        mutable double R_[9];    // rotation matrix: spacecraft to ecef

    public:
        void pco(double *pco, const std::string &f1, const std::string &f2)const;
        void pco(const double *xsat, const double *vsat, const double *xsun, double *pco, double *R=nullptr)const;
        double pcv(double zen, double azi)const;
        inline const std::string &svn()const { return svn_; }
        inline const std::string &svType()const { return svType_; }
        inline const double *R()const { return R_; }
        inline friend bool operator==(const atx_t &atx, const std::string &ant) {
            return atx.name == ant;
        }
    };

#if 0
    static const int MaxSatFreq = 4;
    static const int MaxSatZen = 18;
    struct SatAtx {
        std::string name;
        int nfreq;
        double pco[MaxSatFreq][3];  // NEU
        double dzen, zen1, zen2;
        double pcv[MaxSatFreq][MaxSatZen]; // 

        inline friend bool operator<(const SatAtx &lhs, const SatAtx &rhs) {
            return lhs.name < rhs.name;
        }
    };
    static const int MaxRcvFreq = 5;
    static const int MaxRcvZen = 19;
    struct RcvAtx {
        std::string name;
        int nfreq;
        double pco[MaxRcvFreq][3]; // NEU
        double dazi;
        double dzen, zen1, zen2;
        double pcv[MaxRcvFreq][73][MaxRcvZen];

        inline friend bool operator<(const RcvAtx &lhs, const RcvAtx &rhs) {
            return lhs.name < rhs.name;
        }
    };
#endif

// private:
//     RinexAtx(const RinexAtx &);
//     RinexAtx &operator=(const RinexAtx &);

public:
    RinexAtx(): atxFile_(nullptr) {}
    ~RinexAtx(){}

    bool open(const std::string &path);
    void close();
    const atx_t *atx(MJD t, const std::string &ant)const;

private:
    // find require ant in atxFile_
    bool find_atx(MJD t, const std::string &ant, atx_t &atx)const;

    // read the found atx in atxFile_
    void read_atx(FILE *fp, atx_t &atx)const;

    // skip current atx
    void skip_atx(FILE *fp)const;

private:
    FILE *atxFile_;
    mutable std::vector<atx_t> atxs_;
};

#endif
