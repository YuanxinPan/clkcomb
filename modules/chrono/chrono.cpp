#include "chrono.h"

#include <time.h>
#include <math.h>
#include <string.h>

void sod2hms(double sod, int *h, int *m, double *sec)
{
    *h = (int)(sod / 3600.0);
    *m = (int)((sod - *h*3600.0) / 60.0);
    *sec = sod - *h*3600.0 - *m*60.0;
}

int date2mjd(int y, int m, int d)
{
    y = yr2year(y);
    int mjd = -678987 + 367*y + d + static_cast<int>(275*m/9.0)
              - static_cast<int>(1.75*(y + static_cast<int>((m+9.0)/12)));
    return mjd;
}

void mjd2date(int jd, double sod, int *iy, int *imon, int *id, int *ih, int *im, double *is)
{
    int mjd, doy;
    double  msod;

    // check type of modified julday
    if (jd != 0) {
        mjd = jd;
        msod = sod;
    }
    else {
        mjd = (int)sod;
        msod = (sod - mjd)*86400.0;
    }

    // transformation
    mjd2doy(mjd, iy, &doy);
    doy2date(*iy, doy, imon, id);
    sod2hms(msod, ih, im, is);
}

int ydoy2mjd(int year, int doy)
{
    int m, d=0;
    year = yr2year(year);
    doy2date(year, doy, &m, &d);
    return date2mjd(year, m, d);
}

void wksow2mjd(int week, double sow, int *mjd, double *sod)
{
    *mjd = week * 7 + 44244 + (int)(sow / 86400.0);
    *sod = fmod(sow, 86400.0);
}

void mjd2wksow(int mjd, double sod, int &week, double &sow)
{
    week = (mjd - 44244)/7;
    sow  = (mjd - 44244 - week*7)*86400 + sod;
}

// parameters : day month year-- if month == 0, day = doy
void date2gwk(int year, int month, int day, int *week, int *wd)
{
    int mjd;

    if (year != 0){
        if (month != 0)
            mjd = date2mjd(year, month, day);
        else
            mjd = date2mjd(year, 1, 1) + day - 1;
    }
    else
        mjd = day;

    *week = (mjd - 44244) / 7;
    *wd = mjd - 44244 - *week * 7;
}

void  mjd2doy(int jd, int *iyear, int *idoy)
{
    *iyear = (jd + 678940) / 365;
    *idoy = jd - date2mjd(*iyear, 1, 1);
    while (*idoy <= 0)
    {
        --*iyear;
        *idoy = jd - date2mjd(*iyear, 1, 1) + 1;
    }
}

int date2doy(int y, int m, int d)
{
    return date2mjd(y, m, d) - date2mjd(y, 1, 1) + 1;
}

void doy2date(int iyear, int idoy, int *imonth, int *iday)
{
    int days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

    iyear = yr2year(iyear);
    if ((iyear%4 == 0 && iyear%100 != 0) || iyear%400 == 0)
        days_in_month[1] = 29;

    for (*imonth = 0; *imonth < 12; ++*imonth)
    {
        idoy -= days_in_month[*imonth];
        if (idoy > 0)
            continue;
        *iday = idoy + days_in_month[*imonth];
        break;
    }
    ++*imonth;
}

const char *run_time()
{
    static char ret[12];
    time_t t;
    char   ct[32];

    time(&t);
    strcpy(ct, ctime(&t));
    //ct[19] = '\0';
    strncpy(ret, ct+11, 8);
    ret[9] = '\0';

    return ret;
}
