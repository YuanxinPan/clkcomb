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
