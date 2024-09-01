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

#ifndef RINEXNAV_H
#define RINEXNAV_H

#include "Ephemeris.h"
#include "../chrono/chrono.h"

#include <map>
#include <string>
#include <vector>
#include <stdio.h>

class RinexNav
{
private:
	class Header
	{
	public:
        Header() {}
        ~Header() {}
	    
	    bool read(FILE *fp);
	    bool write(FILE *fp);

	public:
        // inonsphere correction coef: 0-3:alpha, 4-7:beta
        double version;
	    double ionoCoefs[8];
	    double A0, A1;
	    int T, W, leapSec;
	};

private:
    RinexNav(const RinexNav &);
    RinexNav &operator=(const RinexNav &);

public:
    RinexNav() {}
    ~RinexNav() {}

    bool read(const std::string &path);//, const std::vector<std::string> &prns);
    const Ephemeris *eph(MJD t, const std::string &prn)const;
    void  repeat_time(FILE *fp)const;

private:
    Header header_;
    std::map<std::string, std::vector<Ephemeris>> ephs_;
};

#endif
