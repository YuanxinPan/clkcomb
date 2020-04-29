#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include "../chrono/chrono.h"
#include "../coord/CartCoor.h"

#include <string>
#include <stdio.h>

// GPS broadcast-ephmeris
class Ephemeris
{
public:
    Ephemeris();
    ~Ephemeris();
    // read n/p RINEX-NAV from stream
    bool read(FILE *fp, char tp);
    bool write(FILE *fp)const;
    bool satClk(MJD t, double &sck)const;
    bool satPos(MJD t, CartCoor &pos)const;
    bool satPosClk(MJD t, CartCoor &pos, double &sck)const;
    double orbit_period()const;
    
private:
    double calcEt(double Mt)const;
    bool quality_check()const;

public:
    std::string prn;
    MJD toc;            // time from clk reference epoach
    MJD toe;

    double clkBias;      // f0 
    double clkDrift;     // f1
    double clkDriftRate; // f2 = 0.5*acc

    // channle-1
    double IODE;         //
    double Crs;          //
    double dn;           //
    double M0;           //

    // channle-2
    double Cuc;          //
    double e;            //
    double Cus;          //
    double sqrtA;        //

    // channle-3
    double TOE;          //
    double Cic;          //
    double OMEGA;        //
    double Cis;          //

    // channle-4
    double i0;           //
    double Crc;          //
    double omega;        //
    double omegaDot;     //

    // channle-5
    double iDot;         //
    double codesOnL2;    //
    double GPSWeek;      //
    double L2PDataFlag;  //

    // channle-6
    double SVAccuracy;   //
    double SVHealth;     //
    double TGD;          //
    double IODC;         //

    // channle-7
    double transTimeOfMsg;//
    double fitInterval;   //
    double spare1;
    double spare2;
};

#endif
