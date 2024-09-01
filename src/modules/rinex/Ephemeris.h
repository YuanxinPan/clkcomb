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
