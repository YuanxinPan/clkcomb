/***************************************************************************
*******                  JPLBIN.H   v.1.2                          *********
****************************************************************************
**  This header file is used both by ASC2EPH and TESTEPH programs.        **
**  It is NECESSARY TO ADJUST IT MANUALLY for different ephemerides.      **
****************************************************************************
**  Written: May 28, 1997 by PAD   **  Last modified: June 23,1997 by PAD **
****************************************************************************
**  PAD: dr. Piotr A. Dybczynski,          e-mail: dybol@amu.edu.pl       **
**   Astronomical Observatory of the A.Mickiewicz Univ., Poznan, Poland   **
***************************************************************************/
#include "JplEph.h"
#include "../io/io.h"

#include <math.h>

/* UNCOMMENT ONE AND ONLY ONE OF THE FOLLOWING DENUM DEFINITIONS: */
/*#define DENUM 200*/
/*#define DENUM 403*/
/*#define DENUM 404*/
#define DENUM 405
/*#define DENUM 406*/

#if   DENUM==200
#define KSIZE 1652
#elif DENUM==403 || DENUM==405
#define KSIZE 2036
#elif DENUM==404 || DENUM==406
#define KSIZE 1456
#endif

#define NRECL 4
#define RECSIZE (NRECL*KSIZE)
#define NCOEFF (KSIZE/2)
#define TRUE 1
#define FALSE 0

/* we use long int instead of int for DOS-Linux compatibility */

struct rec1{
  char ttl[3][84];
  char cnam[400][6];
  double ss[3];
  long int ncon;
  double au;
  double emrat;
  long int ipt[12][3];
  long int numde;
  long int lpt[3];
};
struct{
  struct rec1 r1;
  char spare[RECSIZE-sizeof(struct rec1)];
} R1;

struct rec2{
  double cval[400];
};
struct{
  struct rec2 r2;
  char spare[RECSIZE-sizeof(struct rec2)];
} R2;


JplEph::JplEph()
{
    KM   = TRUE;
    BARY = FALSE;
    ephFile_ = nullptr;
}

JplEph::~JplEph()
{
    if (ephFile_ != nullptr) {
        fclose(ephFile_);
    }
}

bool JplEph::open(const std::string &path)
{
    ephFile_ = fopen(path.c_str(), "rb");
    if (ephFile_ == nullptr) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "JplEph::open: no such file: %s\n", path.c_str());
        return false;
    }

    int nvs;
    char nams[400][6];
    double vals[400];
    double ss[3];  // beg/end/step
    constan(nams,vals,ss,&nvs);

    //printf("%14.2f  %14.2f  %14.2f %d\n",ss_[0],ss_[1],ss_[2], nvs);
    //for(int i=0;i<nvs;++i)
    //    printf("%.6s  %24.16E\n",nams[i],vals[i]);
    return true;
}

/****************************************************************************/
/*****************************************************************************
**                         pleph(et,ntar,ncent,rrd)                         **
******************************************************************************
**                                                                          **
**    This subroutine reads the jpl planetary ephemeris                     **
**    and gives the position and velocity of the point 'ntarg'              **
**    with respect to 'ncent'.                                              **
**                                                                          **
**    Calling sequence parameters:                                          **
**                                                                          **
**      et = (double) julian ephemeris date at which interpolation          **
**           is wanted.                                                     **
**                                                                          **
**    ntarg = integer number of 'target' point.                             **
**                                                                          **
**    ncent = integer number of center point.                               **
**                                                                          **
**    The numbering convention for 'ntarg' and 'ncent' is:                  **
**                                                                          **
**            1 = mercury           8 = neptune                             **
**            2 = venus             9 = pluto                               **
**            3 = earth            10 = moon                                **
**            4 = mars             11 = sun                                 **
**            5 = jupiter          12 = solar-system barycenter             **
**            6 = saturn           13 = earth-moon barycenter               **
**            7 = uranus           14 = nutations (longitude and obliq)     **
**                                 15 = librations, if on eph. file         **
**                                                                          **
**            (If nutations are wanted, set ntarg = 14.                     **
**             For librations, set ntarg = 15. set ncent= 0)                **
**                                                                          **
**     rrd = output 6-element, double array of position and velocity        **
**           of point 'ntarg' relative to 'ncent'. The units are au and     **
**           au/day. For librations the units are radians and radians       **
**           per day. In the case of nutations the first four words of      **
**           rrd will be set to nutations and rates, having units of        **
**           radians and radians/day.                                       **
**                                                                          **
**           The option is available to have the units in km and km/sec.    **
**           for this, set km=TRUE at the beginning of the program.         **
*****************************************************************************/
bool JplEph::pleph(double et,int ntarg,int ncent,double rrd[] )
{
  if (et<R1.r1.ss[0] || et >R1.r1.ss[1]) {
      fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
              "JplEph::pleph: epoch%10.1f out of range [%8.1f, %8.1f]\n",
              et, R1.r1.ss[0], R1.r1.ss[0]);
      return false;
  }
  double et2[2],pv[13][6];/* pv is the position/velocity array
                             NUMBERED FROM ZERO: 0=Mercury,1=Venus,...
                             8=Pluto,9=Moon,10=Sun,11=SSBary,12=EMBary
                             First 10 elements (0-9) are affected by state(),
                             all are adjusted here.                         */

  int bsave,i,k;
  int list[12];          /* list is a vector denoting, for which "body"
                            ephemeris values should be calculated by state():
                            0=Mercury,1=Venus,2=EMBary,...,8=Pluto,
                            9=geocentric Moon, 10=nutations in long. & obliq.
                            11= lunar librations  */

  /*  initialize et2 for 'state' and set up component count   */

  et2[0]=et;
  et2[1]=0.0;

  for(i=0;i<6;++i) rrd[i]=0.0;

  if(ntarg == ncent) {
      fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
              "JplEph::pleph: same tar & cent\n");
      return false;
  }

  for(i=0;i<12;++i) list[i]=0;

  /*   check for nutation call    */
  if(ntarg == 14)
  {
    if(R1.r1.ipt[11][1] > 0) /* there is nutation on ephemeris */
    {
      list[10]=2;
      state(et2,list,pv,rrd);
    }
    else {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "JplEph::pleph: no nutations on the ephemeris file\n");
    }
    return false;
  }

  /*  check for librations   */
  if(ntarg == 15)
  {
    if(R1.r1.lpt[1] > 0) /* there are librations on ephemeris file */
    {
      list[11]=2;
      state(et2,list,pv,rrd);
      for(i=0;i<6;++i)  rrd[i]=pv[10][i]; /* librations */
    }
    else {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "JplEph::pleph: no librations on the ephemeris file\n");
    }
    return false;
  }

  /*  force barycentric output by 'state'     */
  bsave=BARY;
  BARY= TRUE;

  /*  set up proper entries in 'list' array for state call     */
  for(i=0;i<2;++i) /* list[] IS NUMBERED FROM ZERO ! */
  {
    k=ntarg-1;
    if(i == 1) k=ncent-1;   /* same for ntarg & ncent */
    if(k <= 9) list[k]=2;   /* Major planets */
    if(k == 9) list[2]=2;   /* for moon state earth state is necessary*/
    if(k == 2) list[9]=2;   /* for earth state moon state is necessary*/
    if(k == 12) list[2]=2;  /* EMBary state additionaly */
  }

  /*   make call to state   */

  state(et2,list,pv,rrd);
  /* Solar System barycentric Sun state goes to pv[10][] */
  if(ntarg == 11 || ncent == 11) for(i=0;i<6;++i) pv[10][i]=PVSUN[i];

  /* Solar System Barycenter coordinates & velocities equal to zero */
  if(ntarg == 12 || ncent == 12) for(i=0;i<6;++i) pv[11][i]=0.0;

  /* Solar System barycentric EMBary state:  */
  if(ntarg == 13 || ncent == 13) for(i=0;i<6;++i) pv[12][i]=pv[2][i];

  /* if moon from earth or earth from moon ..... */
  if( (ntarg*ncent) == 30 && (ntarg+ncent) == 13)
  for(i=0;i<6;++i) pv[2][i]=0.0;
  else
  {
    if(list[2] == 2) /* calculate earth state from EMBary */
    for(i=0;i<6;++i) pv[2][i] -= pv[9][i]/(1.0+R1.r1.emrat);

    if(list[9] == 2) /* calculate Solar System barycentric moon state */
    for(i=0;i<6;++i) pv[9][i] += pv[2][i];
  }

  for(i=0;i<6;++i)
  {
    rrd[i] = (pv[ntarg-1][i]-pv[ncent-1][i])*1E3; // unit: m
    if(i>2) rrd[i]*=86400.0;
  }
  BARY=bsave;

  return true;
}

/*****************************************************************************
**                     interp(buf,t,ncf,ncm,na,ifl,pv)                      **
******************************************************************************
**                                                                          **
**    this subroutine differentiates and interpolates a                     **
**    set of chebyshev coefficients to give position and velocity           **
**                                                                          **
**    calling sequence parameters:                                          **
**                                                                          **
**      input:                                                              **
**                                                                          **
**        buf   1st location of array of d.p. chebyshev coefficients        **
**              of position                                                 **
**                                                                          **
**          t   t[0] is double fractional time in interval covered by       **
**              coefficients at which interpolation is wanted               **
**              (0 <= t[0] <= 1).  t[1] is dp length of whole               **
**              interval in input time units.                               **
**                                                                          **
**        ncf   # of coefficients per component                             **
**                                                                          **
**        ncm   # of components per set of coefficients                     **
**                                                                          **
**         na   # of sets of coefficients in full array                     **
**              (i.e., # of sub-intervals in full interval)                 **
**                                                                          **
**         ifl  integer flag: =1 for positions only                         **
**                            =2 for pos and vel                            **
**                                                                          **
**                                                                          **
**      output:                                                             **
**                                                                          **
**        pv   interpolated quantities requested.  dimension                **
**              expected is pv(ncm,ifl), dp.                                **
**                                                                          **
*****************************************************************************/
void JplEph::interp(double coef[],double t[2],int ncf,int ncm,int na,int ifl,
                    double posvel[6])
{
  static double pc[18] = { 1.0 };
  static double vc[18] = { 0.0, 1.0 };
  static int np=2, nv=3;// first=1;
  static double twot=0.0;
  double dna,dt1,temp,tc,vfac,temp1;
  int l,i,j;

  //if(first){           /* initialize static vectors when called first time */
  //  pc[0]=1.0;
  //  pc[1]=0.0;
  //  vc[1]=1.0;
  //  first=0;
  //}

  /*  entry point. get correct sub-interval number for this set
  of coefficients and then get normalized chebyshev time
  within that subinterval.                                             */

  dna=(double)na;
  modf(t[0],&dt1);
  temp=dna*t[0];
  l=(int)(temp-dt1);

  /*  tc is the normalized chebyshev time (-1 <= tc <= 1)    */

  tc=2.0*(modf(temp,&temp1)+dt1)-1.0;

  /*  check to see whether chebyshev time has changed,
  and compute new polynomial values if it has.
  (the element pc[1] is the value of t1[tc] and hence
  contains the value of tc on the previous call.)     */

  if(tc != pc[1])
  {
    np=2;
    nv=3;
    pc[1]=tc;
    twot=tc+tc;
  }

  /*  be sure that at least 'ncf' polynomials have been evaluated
  and are stored in the array 'pc'.    */

  if(np < ncf)
  {
    for(i=np;i<ncf;++i)  pc[i]=twot*pc[i-1]-pc[i-2];
    np=ncf;
  }

  /*  interpolate to get position for each component  */

  for(i=0;i<ncm;++i) /* ncm is a number of coordinates */
  {
    posvel[i]=0.0;
    for(j=ncf-1;j>=0;--j)
    posvel[i]=posvel[i]+pc[j]*coef[j+i*ncf+l*ncf*ncm];
  }

  if(ifl <= 1) return;


  /*  if velocity interpolation is wanted, be sure enough
  derivative polynomials have been generated and stored.    */

  vfac=(dna+dna)/t[1];
  vc[2]=twot+twot;
  if(nv < ncf)
  {
    for(i=nv;i<ncf;++i) vc[i]=twot*vc[i-1]+pc[i-1]+pc[i-1]-vc[i-2];
    nv=ncf;
  }

  /*  interpolate to get velocity for each component    */

  for(i=0;i<ncm;++i)
  {
    posvel[i+ncm]=0.0;
    for(j=ncf-1;j>0;--j)
    posvel[i+ncm]=posvel[i+ncm]+vc[j]*coef[j+i*ncf+l*ncf*ncm];
    posvel[i+ncm]=posvel[i+ncm]*vfac;
  }
  return;
}

/****************************************************************************
****                       split(tt,fr)                                  ****
*****************************************************************************
****  this subroutine breaks a d.p. number into a d.p. integer           ****
****  and a d.p. fractional part.                                        ****
****                                                                     ****
****  calling sequence parameters:                                       ****
****                                                                     ****
****    tt = d.p. input number                                           ****
****                                                                     ****
****    fr = d.p. 2-word output array.                                   ****
****         fr(1) contains integer part                                 ****
****         fr(2) contains fractional part                              ****
****                                                                     ****
****         for negative input numbers, fr(1) contains the next         ****
****         more negative integer; fr(2) contains a positive fraction.  ****
****************************************************************************/
void JplEph::split(double tt, double fr[2])
{
  /*  main entry -- get integer and fractional parts  */

  fr[1]=modf(tt,&fr[0]);

  if(tt >= 0.0 || fr[1] == 0.0) return;

  /*  make adjustments for negative input number   */

  fr[0]=fr[0]-1.0;
  fr[1]=fr[1]+1.0;

  return;
}

/*****************************************************************************
**                        state(et2,list,pv,nut)                            **
******************************************************************************
** This subroutine reads and interpolates the jpl planetary ephemeris file  **
**                                                                          **
**    Calling sequence parameters:                                          **
**                                                                          **
**    Input:                                                                **
**                                                                          **
**        et2[] double, 2-element JED epoch at which interpolation          **
**              is wanted.  Any combination of et2[0]+et2[1] which falls    **
**              within the time span on the file is a permissible epoch.    **
**                                                                          **
**               a. for ease in programming, the user may put the           **
**                  entire epoch in et2[0] and set et2[1]=0.0               **
**                                                                          **
**               b. for maximum interpolation accuracy, set et2[0] =        **
**                  the most recent midnight at or before interpolation     **
**                  epoch and set et2[1] = fractional part of a day         **
**                  elapsed between et2[0] and epoch.                       **
**                                                                          **
**               c. as an alternative, it may prove convenient to set       **
**                  et2[0] = some fixed epoch, such as start of integration,**
**                  and et2[1] = elapsed interval between then and epoch.   **
**                                                                          **
**       list   12-element integer array specifying what interpolation      **
**              is wanted for each of the "bodies" on the file.             **
**                                                                          **
**                        list[i]=0, no interpolation for body i            **
**                               =1, position only                          **
**                               =2, position and velocity                  **
**                                                                          **
**              the designation of the astronomical bodies by i is:         **
**                                                                          **
**                        i = 0: mercury                                    **
**                          = 1: venus                                      **
**                          = 2: earth-moon barycenter                      **
**                          = 3: mars                                       **
**                          = 4: jupiter                                    **
**                          = 5: saturn                                     **
**                          = 6: uranus                                     **
**                          = 7: neptune                                    **
**                          = 8: pluto                                      **
**                          = 9: geocentric moon                            **
**                          =10: nutations in longitude and obliquity       **
**                          =11: lunar librations (if on file)              **
**                                                                          **
**                                                                          **
**    output:                                                               **
**                                                                          **
**    pv[][6]   double array that will contain requested interpolated       **
**              quantities.  The body specified by list[i] will have its    **
**              state in the array starting at pv[i][0]  (on any given      **
**              call, only those words in 'pv' which are affected by the    **
**              first 10 'list' entries (and by list(11) if librations are  **
**              on the file) are set.  The rest of the 'pv' array           **
**              is untouched.)  The order of components in pv[][] is:       **
**              pv[][0]=x,....pv[][5]=dz.                                   **
**                                                                          **
**              All output vectors are referenced to the earth mean         **
**              equator and equinox of epoch. The moon state is always      **
**              geocentric; the other nine states are either heliocentric   **
**              or solar-system barycentric, depending on the setting of    **
**              global variables (see below).                               **
**                                                                          **
**              Lunar librations, if on file, are put into pv[10][k] if     **
**              list[11] is 1 or 2.                                         **
**                                                                          **
**        nut   dp 4-word array that will contain nutations and rates,      **
**              depending on the setting of list[10].  the order of         **
**              quantities in nut is:                                       **
**                                                                          **
**                       d psi  (nutation in longitude)                     **
**                       d epsilon (nutation in obliquity)                  **
**                       d psi dot                                          **
**                       d epsilon dot                                      **
**                                                                          **
**    Global variables:                                                     **
**                                                                          **
**         KM   logical flag defining physical units of the output          **
**              states. KM = TRUE, km and km/sec                            **
**                         = FALSE, au and au/day                           **
**              default value = FALSE.  (Note that KM determines time unit  **
**              for nutations and librations. Angle unit is always radians.)**
**                                                                          **
**       BARY   logical flag defining output center.                        **
**              only the 9 planets (entries 0-8) are affected.              **
**                       BARY = TRUE =\ center is solar-system barycenter   **
**                            = FALSE =\ center is sun                      **
**              default value = FALSE                                       **
**                                                                          **
**      PVSUN   double, 6-element array containing                          **
**              the barycentric position and velocity of the sun.           **
**                                                                          **
*****************************************************************************/
void JplEph::state(double et2[2],int list[12],double pv[][6],double nut[4])
{
  int i,j;
  static int ipt[13][3], first=TRUE;
  /* local copy of R1.r1.ipt[] is extended for R1.r1.lpt[[] */
  long int nr;
  static long int nrl=0;
  double pjd[4];
  static double buf[NCOEFF];
  double s,t[2],aufac;
  double pefau[6];

  if(first)
  {
    first=FALSE;
    for(i=0;i<3;++i)
    {
      for(j=0;j<12;++j) ipt[j][i]=(int)R1.r1.ipt[j][i];
      ipt[12][i] = (int)R1.r1.lpt[i];
    }
  }

  /*  ********** main entry point **********  */

  s=et2[0] - 0.5;
  split(s,&pjd[0]);
  split(et2[1],&pjd[2]);
  pjd[0]=pjd[0]+pjd[2]+0.5;
  pjd[1]=pjd[1]+pjd[3];
  split(pjd[1],&pjd[2]);
  pjd[0]=pjd[0]+pjd[2];
  /* here pjd[0] contains last midnight before epoch desired (in JED: *.5)
  and pjd[3] contains the remaining, fractional part of the epoch         */

  /*  error return for epoch out of range  */

  //if( (pjd[0]+pjd[3]) < R1.r1.ss[0] || (pjd[0]+pjd[3]) > R1.r1.ss[1] )
  //{
  //  fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
  //          "JplEph::state: Requested JED not within ephemeris limits.\n");
  //  return;
  //}

  /*   calculate record # and relative time in interval   */

  nr=(long)((pjd[0]-R1.r1.ss[0])/R1.r1.ss[2])+2;
  /* add 2 to adjus for the first two records containing header data */
  if(pjd[0] == R1.r1.ss[1]) nr=nr-1;
  t[0]=( pjd[0]-( (1.0*nr-2.0)*R1.r1.ss[2]+R1.r1.ss[0] ) +
  pjd[3] )/R1.r1.ss[2];

  /*   read correct record if not in core (static vector buf[])   */

  if(nr != nrl)
  {
    nrl=nr;
    fseek(ephFile_,nr*RECSIZE,SEEK_SET);
    fread(buf,sizeof(buf),1,ephFile_);
  }

  if(KM)
  {
    t[1]=R1.r1.ss[2]*86400.0;
    aufac=1.0;
  }
  else
  {
    t[1]=R1.r1.ss[2];
    aufac=1.0/R1.r1.au;
  }

  /*  every time interpolate Solar System barycentric sun state   */

  interp(&buf[ipt[10][0]-1],t,ipt[10][1],3,ipt[10][2],2,pefau);

  for(i=0;i<6;++i)  PVSUN[i]=pefau[i]*aufac;

  /*  check and interpolate whichever bodies are requested   */

  for(i=0;i<10;++i)
  {
    if(list[i] == 0) continue;

    interp(&buf[ipt[i][0]-1],t,ipt[i][1],3,ipt[i][2],list[i],pefau);

    for(j=0;j<6;++j)
    {
      if(i < 9 && !BARY)   pv[i][j]=pefau[j]*aufac-PVSUN[j];
      else                 pv[i][j]=pefau[j]*aufac;
    }
  }

  /*  do nutations if requested (and if on file)    */

  if(list[10] > 0 && ipt[11][1] > 0)
  interp(&buf[ipt[11][0]-1],t,ipt[11][1],2,ipt[11][2],list[10],nut);

  /*  get librations if requested (and if on file)    */

  if(list[11] > 0 && ipt[12][1] > 0)
  {
    interp(&buf[ipt[12][0]-1],t,ipt[12][1],3,ipt[12][2],list[11],pefau);
    for(j=0;j<6;++j) pv[10][j]=pefau[j];
  }
  return;
}

/****************************************************************************
**                        constan(nam,val,sss,n)                           **
*****************************************************************************
**                                                                         **
**    this function obtains the constants from the ephemeris file          **
**    calling sequence parameters (all output):                            **
**      char nam[][6] = array of constant names (max 6 characters each)    **
**      double val[] = array of values of constants                        **
**      double sss[3] = JED start, JED stop, step of ephemeris             **
**      int n = number of entries in 'nam' and 'val' arrays                **
** Note: we changed name of this routine because const is a reserved word  **
**       in the C language.                                                **
*****************************************************************************
**    external variables:                                                  **
**         struct R1 and R2 (first two records of ephemeris)               **
**         defined in file:     jplbin.h                                   **
**         ephFile_ = ephemeris binary file pointer (obtained from fopen())**
****************************************************************************/
void JplEph::constan(char nam[][6], double val[], double sss[], int *n)
{
  int i,j;

  fread(&R1,sizeof(R1),1,ephFile_);
  fread(&R2,sizeof(R2),1,ephFile_);

  *n =(int)R1.r1.ncon;

  for(i=0;i<3;++i) sss[i]=R1.r1.ss[i];

  for(i=0;i<*n;++i) {
    val[i] = R2.r2.cval[i];
    for(j=0;j<6;++j) nam[i][j] = R1.r1.cnam[i][j];
  }
  return;
}
/*************************** THE END ***************************************/
