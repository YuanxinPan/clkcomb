#ifndef _CHRONO_H_
#define _CHRONO_H_

#include "Time.h"

inline int yr2year(int y)
{
    if (y > 99) // 2017
        return y;
    else if (y < 80) // 17
        return y + 2000;
    else // 98
        return y + 1900;
}

inline double hms2sod(int h, int m, double s)
{
    return h*3600 + m*60 + s;
}

void sod2hms(double sod, int *h, int *m, double *s);

int date2mjd(int y, int m, int d);

void mjd2date(int jd, double sod, int *iy, int *imon, int *id, int *ih, int *im, double *is);

void wksow2mjd(int week, double sow, int *mjd, double *sod);

void mjd2wksow(int mjd, double sod, int &week, double &sow);

void mjd2doy(int mjd, int *yr, int *doy);

int date2doy(int y, int m, int d);

void doy2date(int iyear, int idoy, int *imonth, int *iday);

int ydoy2mjd(int year, int doy);

void date2gwk(int year, int month, int day, int *week, int *wd);

const char *run_time();

#endif
