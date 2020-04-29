#ifndef JPLEPH_H
#define JPLEPH_H

#include <string>
#include <stdio.h>

class JplEph
{
public:
    enum {
        Earth = 3,
        Moon  = 10,
        Sun   = 11
    };

private:
    JplEph(const JplEph &);
    JplEph &operator=(const JplEph &);

public:
    JplEph();
    ~JplEph();
    // open jpleph file before interpolation
    bool open(const std::string &path);
    // interpolate position & velocity of ntarg in ICRS coordinate system
    bool pleph(double et, int ntarg, int ncent, double rrd[]);

private:
    void constan(char nam[][6], double val[], double sss[], int *n);
    void state(double et2[], int list[], double pv[][6], double nut[]);
    void split(double tt, double fr[]);
    void interp(double buf[], double t[], int ncf, int ncm, int na, int ifl, double pv[]);

private:
    FILE *ephFile_;
    int KM, BARY;
    double PVSUN[6];
};

#endif
