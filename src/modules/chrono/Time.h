// clkcomb - Clock and phase bias products Combination
// Copyright (C) 2021 Yuanxin Pan
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef TIME_H
#define TIME_H

#include <stdio.h>
#include <string>
#include <vector>

class GPST;
class MJD;

const double GPST2TAI = 19.0;
const double TAI2TT   = 32.184;
const double GPST2TT  = GPST2TAI+TAI2TT;
const double MJD2JD   = 2400000.5;

class CalendTime
{
public:
    CalendTime() {}
    CalendTime(const char *str) { read(str); }
    ~CalendTime() {}

    //friend std::istream& operator>>(std::istream &in, CalendTime &t);
    //friend std::ostream& operator<<(std::ostream &out, const CalendTime &t);
    bool read(const char *str); // const char *buf: avoid to use substr()
    void write(FILE *fp);

    int mjd()const;
    double sod()const;

private:
    int y, m, d;
    int h, min;
    double sec;
};


class GPST
{
public:
    GPST(int w=0, double s=0) :week(w), sow(s) {}
    ~GPST() {}

    //friend std::ostream &operator<<(std::ostream &out, const GPST &t);
    void write(FILE *fp);

private:
    int week;
    double sow;
};


class MJD
{
private:
    struct leapsec_t {
        double t;  // MJD
        int leap;
        //inline bool operator<(int t) { return leap.t < t; }
    };
    static std::vector<leapsec_t> leapsecs_;

public:
    static bool read_leapsec(const std::string &path);

    // function: return leap seconds at MJD::t
    // return 0: no leapsec was found at MJD::t
    //           otherwise leapsec was returned.
    static int  leapsec(double t);

public:
    MJD(const CalendTime &t) {
        d = t.mjd();
        sod = t.sod();
    }

    explicit MJD(double t) {
        d = static_cast<int>(t);
        sod = (t - d)*86400.0;
    }

    MJD(int _d=0, double _sod=0): d(_d), sod(_sod) {}

    ~MJD() {}

    void set(int week, double sow);
    void adjust();
    inline double mjd() { return d + sod/86400.0; }
    inline double cvt2tt()  { return d+(sod+GPST2TT)/86400; }
    inline double cvt2tai() { return d+(sod+GPST2TAI)/86400; }

    //inline operator double() { return sod/86400.0 + d; }

    inline const MJD &operator+=(double t) { sod+=t; return *this; }
    inline const MJD &operator-=(double t) { sod-=t; return *this; }

    //friend MJD operator+(MJD lhs, double s);
    //friend MJD operator-(MJD mjd, double s);
    //friend double operator-(MJD lhs, MJD rhs);
    //friend bool operator<(MJD lhs, MJD rhs);
    //friend bool operator>(MJD lhs, MJD rhs);

public:
    int d;
    double sod;
};


inline MJD operator+(MJD mjd, double s)
{
    mjd.sod += s;
    return mjd;
}

inline MJD operator-(MJD mjd, double s)
{
    mjd.sod -= s;
    return mjd;
}

inline double operator-(const MJD lhs, const MJD rhs)
{
    return 86400*(lhs.d - rhs.d) + lhs.sod - rhs.sod;
}

inline bool operator<(MJD lhs, MJD rhs)
{
    return lhs-rhs < 0;
}

inline bool operator>(MJD lhs, MJD rhs)
{
    return lhs-rhs > 0;
}

inline bool operator<=(MJD lhs, MJD rhs)
{
    return lhs-rhs <= 0;
}

inline bool operator>=(MJD lhs, MJD rhs)
{
    return lhs-rhs >= 0;
}

#endif

