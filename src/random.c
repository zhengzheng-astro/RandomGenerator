#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"

/*
  The Mersenne Twister pseudorandom number generator (by Makoto Matsumoto 
  & Takuji Nishimura) is adopted here to produce random deviates uniformly 
  distributed in [0,1): 
    seed_ran_mt() - set the initial state based on the input seed
    ran_mt()      - generate the random deviates (double) in [0,1)

  The code follows the three functions in GNU Scientific Library's rng/mt.c,
    static void mt_set (void *vstate, unsigned long int s)
    static inline unsigned long mt_get (void *vstate)
    static double mt_get_double (void * vstate)
*/

#define N 624   /* Period parameters */
#define M 397

/* most significant w-r bits */
const unsigned long UPPER_MASK = 0x80000000UL;   

/* least significant r bits */
const unsigned long LOWER_MASK = 0x7fffffffUL;   

/** RNG with the MT algorithm **/
double ran_mt(void *vstate)
{
  mt_state_t *state = (mt_state_t *) vstate;

  unsigned long k ;
  unsigned long int *const mt = state->mt;

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

  if (state->mti >= N)
    {   /* generate N words at one time */
      int kk;

      for (kk = 0; kk < N - M; kk++)
        {
          unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
          mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
        }
      for (; kk < N - 1; kk++)
        {
          unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
          mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
        }

      {
        unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
      }

      state->mti = 0;
    }

  /* Tempering */
  
  k = mt[state->mti];
  k ^= (k >> 11);
  k ^= (k << 7) & 0x9d2c5680UL;
  k ^= (k << 15) & 0xefc60000UL;
  k ^= (k >> 18);

  state->mti++;

  return k / 4294967296.0;
}


/** set the intial state of the MT RNG **/
void seed_ran_mt(void *vstate, unsigned long int s)
{
  mt_state_t *state = (mt_state_t *) vstate;
  int i;

  if (s == 0)
    s = 4357;   /* the default seed is 4357 */

  state->mt[0]= s & 0xffffffffUL;

  for (i = 1; i < N; i++)
    {
      /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
         Ed. p.106 for multiplier. */

      state->mt[i] =
        (1812433253UL * (state->mt[i-1] ^ (state->mt[i-1] >> 30)) + i);
      
      state->mt[i] &= 0xffffffffUL;
    }

  state->mti = i;
}


/** uz sampling method in the Zheng&Miralda-Escude(2002; ZM2002) paper **/

void get_uz0(double *uz, double a, double x, void *vstate)
{
  double u,r,d,t,t1,emd2,P,ratio,sign;

  sign=1.0;
  if(x<0) {x=-x; sign=-1.0;}
  d=x/(1.01+x/210.0+x*x/105.0); /*experience*/
  emd2=exp(-d*d);
  t1=atan((d-x)/a);
  P=(t1+0.5*M_PI)/(0.5*M_PI*(1+emd2)+t1*(1-emd2));

  do {
      if(ran_mt(vstate)<P)
        {
         t=-0.5*M_PI+ran_mt(vstate)*(t1+0.5*M_PI);
         u=a*tan(t)+x;
        }
       else
        {
         t=t1+ran_mt(vstate)*(0.5*M_PI-t1);
         u=a*tan(t)+x;
        }

       r=ran_mt(vstate);
       ratio=exp(-u*u);
       if(u>=d) ratio=ratio/emd2;
     }
  while(r > ratio);

  *uz=sign*u;
}


/** 
  do the v^2 <-> exp(-u^2) mapping 
  sampling > thres points for exp(-u^2) distribution 
**/

double thres_expu2(void *vstate, double thres) 
{
 //double thres;                  /** pick up > thres points for exp(-u^2) distribution **/
 double crit;                   /** probability boundary for > thres, exp(-thres^2) **/
                                /** we only need random[0,crit] **/
 double sqrtcrit;               /** sqrt(crit) **/
 double y1,y2,y,fac,r,v1,v2;

 y1=0.0;
 y2=0.0;
 crit=exp(-thres*thres);
 sqrtcrit=sqrt(crit);

 do {
     do {
         v1=sqrtcrit*ran_mt(vstate);
         v2=sqrtcrit*ran_mt(vstate);
         r=v1*v1+v2*v2;
     } while (r >= crit);
     fac=sqrt(-log(r)/r);
     y1=v1*fac;
     y2=v2*fac;
 } while (y1<thres && y2<thres);

 if(y1>=thres) y=y1;
 else y=y2;

 return y;
}


/** 
  do the v^2 <-> exp(-u^2) mapping 
  thres1 < thres2 
  sampling (thres1, thres2) points for exp(-u^2) distribution 
  adopt bin2_expu2 in the end 
**/
double bin_expu2(void *vstate, double thres1, double thres2) 
{
 //double thres;                  /** pick up > thres points for exp(-u^2) distribution **/
 double crit;              /** probability boundary for > thres, exp(-thres^2) **/
                                  /** we only need random[crit1,crit2] **/
 double sqrtcrit;      /** sqrt(crit) **/
 double y1,y2,y,fac,r,v1,v2;

 y1=0.0;
 y2=0.0;
 crit=exp(-thres1*thres1);
 sqrtcrit=sqrt(crit);

 do {
     do {
         v1=sqrtcrit*ran_mt(vstate);
         v2=sqrtcrit*ran_mt(vstate);
         r=v1*v1+v2*v2;
     } while (r >= crit);
     fac=sqrt(-log(r)/r);
     y1=v1*fac;
     y2=v2*fac;
   } while (!( (y1>thres1 && y1<=thres2) || (y2>thres1 && y2<=thres2) ) );

 if(y1>thres1 && y1<=thres2) y=y1;
 else y=y2;

 return y;
}

double bin2_expu2(void *vstate, double thres1, double thres2) 
{
 double r1,r2,r;
 double t1,t2,t,tmax,R2,ratio;

 r1=thres1;
 r2=thres2;
 tmax=acos(r1/r2);
 do {
    r=sqrt(r1*r1-log(ran_mt(vstate)));
    if(r<r2) t1=0.0;
    else     t1=acos(r2/r);
    t2=acos(r1/r);
    ratio=(t2-t1)/tmax;
    R2=ran_mt(vstate);
 } while (R2>ratio);

 t=t1+ran_mt(vstate)*(t2-t1);
  
 return r*cos(t);
}


/** 
  making the uz sampling far more efficient by choosing better comparison 
  functions for rejection method, in large |x| regime; 
  (ZZ20190912)
  [This one is not optimal at x~0.] 
  We end up using the optimized one (ZM02+LZ26): get_uz_opt()
  See Appendix A of Li&Zheng(2026)
**/
void get_uz1(double *uz, double a, double xt, void *vstate)
{
  double halfsqrtpi=0.886226925452757940959713778284;
  double x,ratio,sign;
  double I1t,I1,I2,I3,I4,I5,Itot,C1,C2,C3,C4;
  double invx,xx,aa,p,q,u,u1,u2,u3,umx,u1mx,u2mx,u3mx;
  double xx_p_aa,xx_p_aa_p_1;
  double up,up2,up3,dup,upmx;
  double t0,t1,t2,t3;
  double exp_upsq,A,Acoef;
  double Rc,R;
  int    i;

  x=xt;
  sign=1.0;
  if(x<0) {x=-x; sign=-1.0;}

  invx=1.0/x;
  aa = a*a;
  xx = x*x;
  xx_p_aa = xx+aa;
  t0 = atan(x/a);

  I1t=halfsqrtpi/xx_p_aa;

  if(x<=1.0) {
    q =1.1;
    u1=q*x;
    u1mx=u1-x;
    t1=atan(u1mx/a);
    I1=I1t;
    I2=(t1+t0)/a;
    I3=halfsqrtpi*(1.0-erf(u1))/(u1mx*u1mx+aa);
    Itot=I1+I2+I3;
    C1=I1/Itot;
    C2=(I1+I2)/Itot;
    //printf("I1 I2 I3 = %e %e %e C1 C2 C3 = %e %e %e\n",I1,I2,I3,C1,C2,1.0-(C1+C2)); exit(0);
    do {
      Rc=ran_mt(vstate);
      if(Rc<=C1) {
        u=-thres_expu2(vstate,0.0);
        umx=u-x;
        ratio = xx_p_aa/(umx*umx+aa);
      } 
      else if(Rc<=C2) {
        R=ran_mt(vstate);
        u = x + a*tan(R*t1+(R-1.0)*t0);
        ratio = exp(-u*u);
      } 
      else {
        u=thres_expu2(vstate,u1);
        umx=u-x;
        ratio = (u1mx*u1mx+aa)/(umx*umx+aa);
      }
      R=ran_mt(vstate);
    } while(R>ratio);
  }

  if(1.0<x && x<=2.0) {
    p=0.7;
    q=1.1;
    u1=p*x; 
    u2=q*x;
    u1mx=u1-x;
    u2mx=u2-x;
    t1=atan(u1mx/a);
    t2=atan(u2mx/a);
    I1=I1t;
    I2=(t1+t0)/a;
    I3=(t2-t1)/(a*exp(u1*u1));
    I4=halfsqrtpi*(1.0-erf(u2))/(u2mx*u2mx+aa);
    //printf("I1 I2 I3 I4 = %e %e %e %e\n",I1,I2,I3,I4); exit(0);
    Itot=I1+I2+I3+I4;
    C1=I1/Itot;
    C2=(I1+I2)/Itot;
    C3=(I1+I2+I3)/Itot;
 
    do {
      Rc=ran_mt(vstate);
      if(Rc<=C1) {
        u=-thres_expu2(vstate,0.0);
        umx=u-x;
        ratio = xx_p_aa/(umx*umx+aa);
      }
      else if(Rc<=C2) {
        R=ran_mt(vstate);
        u = x + a*tan(R*t1+(R-1.0)*t0);
        ratio = exp(-u*u);
      }
      else if(Rc<=C3) {
        R=ran_mt(vstate);
        u = x + a*tan(R*t2+(1.0-R)*t1);
        ratio = exp(-(u*u-u1*u1));
      }
      else {
         u=thres_expu2(vstate,u2);
         umx=u-x;
         ratio = (u2mx*u2mx+aa)/(umx*umx+aa);
      }
      R=ran_mt(vstate);
    } while(R>ratio);
  }

  if(x>2.0) {
    /** find u_p where f(u) peaks, based on Newton's method of finding roots **/
    xx_p_aa_p_1=xx_p_aa + 1.0;
    up=x/xx_p_aa_p_1;
    for(i=0;i<10;i++) {  
      up2=up*up; up3=up2*up;
      dup=-(up3-2.0*x*up2+xx_p_aa_p_1*up-x)/(3.0*up2-4.0*x*up+xx_p_aa_p_1);
      up=up+dup;
    }
    upmx=up-x;

    //p=0.7;
    p=0.7+0.4/x;
    q=1.1;
    u2=p*x; 
    u3=q*x;
    u2mx=u2-x;
    u3mx=u3-x;

    exp_upsq=exp(up*up);  
    A=(u2mx*u2mx+aa)/xx*exp(2.0*u2/x);
    Acoef=A*xx/exp(1.0/xx);
    u1=invx+sqrt( - log(Acoef/exp_upsq/(upmx*upmx+aa))  );
    u1mx=u1-x;

    t2=atan(u2mx/a);
    t3=atan(u3mx/a);
    I1=I1t;
    I2=u1/exp_upsq/(upmx*upmx+aa);
    I3=halfsqrtpi/Acoef*(erf(u2-invx)-erf(u1-invx));
    I4=(t3-t2)/(a*exp(u2*u2));
    I5=halfsqrtpi*(1.0-erf(u3))/(u3mx*u3mx+aa);
    Itot=I1+I2+I3+I4+I5;
    C1=I1/Itot;
    C2=(I1+I2)/Itot;
    C3=(I1+I2+I3)/Itot;
    C4=(I1+I2+I3+I4)/Itot;

    /*
    printf("up u1 u2 u3 = %e %e %e %e t2 t3 =%f %f\n",up,u1,u2,u3,t2,t3);
    printf("I1 I2 I3 I4 I5 = %e %e %e %e %e\n",I1/Itot,I2/Itot,I3/Itot,I4/Itot,I5/Itot); 
    printf("C1 C2 C3 C4 = %e %e %e %e\n",C1,C2,C3,C4);
    exit(0);
    */

    do {
      Rc=ran_mt(vstate);
      if(Rc<=C1) {
        u=-thres_expu2(vstate,0.0);
        umx=u-x;
        ratio = xx_p_aa/(umx*umx+aa);
      }
      else if(Rc<=C2) {
        u=u1*ran_mt(vstate);
        umx=u-x;
        ratio = exp(-(u*u-up*up))*(upmx*upmx+aa)/(umx*umx+aa);
      }
      else if(Rc<=C3) {
        u=invx+bin2_expu2(vstate,u1-invx,u2-invx);
        umx=u-x;
        ratio = A*xx/exp(2.0*u/x)/(umx*umx+aa);
      }
      else if(Rc<=C4) {
        R=ran_mt(vstate);
        u = x + a*tan(R*t3+(1.0-R)*t2);
        ratio = exp(-(u*u-u2*u2));
      }
      else {
        u=thres_expu2(vstate,u3);
        umx=u-x;
        ratio = (u3mx*u3mx+aa)/(umx*umx+aa);
      }
      R=ran_mt(vstate);
    } while(R>ratio);
  }

  *uz=sign*u;
}

/** 
  optimized to choose the better performance (based on running time 
  tests for lgT=2,4,6 with different x's) 
  (ZZ20190927)
  See Appendix A of Li&Zheng(2026)
**/
void get_uz_opt(double *uz, double a, double xt, void *vstate)
{
  double halfsqrtpi=0.886226925452757940959713778284;
  double x,ratio,sign;
  double I1t,I1,I2,I3,I4,I5,Itot,C1,C2,C3,C4;
  double invx,xx,aa,p,q,u,u1,u2,u3,umx,u1mx,u2mx,u3mx;
  double xx_p_aa,xx_p_aa_p_1;
  double up,up2,up3,dup,upmx;
  double t0,t1,t2,t3;
  double exp_upsq,A,Acoef;
  double Rc,R;
  int    i;

  x=xt;
  sign=1.0;
  if(x<0) {x=-x; sign=-1.0;}

  /** based on the relative performance of get_uz0 and get_uz1 **/
  if(x<-0.63*log10(a)+1.4) get_uz0(uz,a,xt,vstate);
  else get_uz1(uz,a,xt,vstate);

}


/** Semelin+2007 method **/
void get_uz_Semelin(double *uz, double a, double x, void *vstate)
{
  double u,r,d,t,t1,emd2,P,ratio,sign;

  sign=1.0;
  if(x<0) {x=-x; sign=-1.0;}

  d=0.0;
//  if(x>3) d=1.85-log10(a)/6.73+log(log(x));
  if(x>3) d=1.85-log(a)/6.73+log(log(x));

  emd2=exp(-d*d);
  t1=atan((d-x)/a);
  P=(t1+0.5*M_PI)/(0.5*M_PI*(1+emd2)+t1*(1-emd2));

  do {
      if(ran_mt(vstate)<P)
        {
         t=-0.5*M_PI+ran_mt(vstate)*(t1+0.5*M_PI);
         u=a*tan(t)+x;
        }
       else
        {
         t=t1+ran_mt(vstate)*(0.5*M_PI-t1);
         u=a*tan(t)+x;
        }

       r=ran_mt(vstate);
       ratio=exp(-u*u);
       if(u>=d) ratio=ratio/emd2;
     }
  while(r > ratio);

  *uz=sign*u;
}


/** RASCAS method **/
void get_uz_RASCAS(double *uz, double a, double x, void *vstate)
{
  double u,r,d,t,t1,z,z2,emd2,P,ratio,sign;

  sign=1.0;
  if(x<0) {x=-x; sign=-1.0;}

  z=log10(a);
  z2=z*z;
  d =      2.648963 + 2.014446*z + 0.351479*z2
     + x*(-4.058673 - 3.675859*z - 0.640003*z2
     + x*( 3.017395 + 2.117133*z + 0.370294*z2
     + x*(-0.869789 - 0.565886*z - 0.096312*z2
     + x*( 0.110987 + 0.070103*z + 0.011557*z2
     + x*(-0.005200 - 0.003240*z - 0.000519*z2)))));

  emd2=exp(-d*d);
  t1=atan((d-x)/a);
  P=(t1+0.5*M_PI)/(0.5*M_PI*(1+emd2)+t1*(1-emd2));

  do {
      if(ran_mt(vstate)<P)
        {
         t=-0.5*M_PI+ran_mt(vstate)*(t1+0.5*M_PI);
         u=a*tan(t)+x;
        }
       else
        {
         t=t1+ran_mt(vstate)*(0.5*M_PI-t1);
         u=a*tan(t)+x;
        }

       r=ran_mt(vstate);
       ratio=exp(-u*u);
       if(u>=d) ratio=ratio/emd2;
     }
  while(r > ratio);

  *uz=sign*u;
}

/** Smith+2015 method **/
void get_uz_Smith(double *uz, double a, double x, void *vstate)
{
  double u,r,d,t,t1,xcw,emd2,P,ratio,sign;

  sign=1.0;
  if(x<0) {x=-x; sign=-1.0;}

  //xcw=6.9184721+81.766279/(log10(a)-14.651253);
  xcw=6.9184721+81.766279/(log(a)-14.651253); 
  /** based on reading the inset of their Fig.2 **/

  d=0.0;
  if(1<x && x<xcw) d=x-1/(x+exp(1-x*x)/a);
  if(xcw<x)        d=xcw-1/xcw+0.15*(x-xcw);

  emd2=exp(-d*d);
  t1=atan((d-x)/a);
  P=(t1+0.5*M_PI)/(0.5*M_PI*(1+emd2)+t1*(1-emd2));

  do {
      if(ran_mt(vstate)<P)
        {
         t=-0.5*M_PI+ran_mt(vstate)*(t1+0.5*M_PI);
         u=a*tan(t)+x;
        }
       else
        {
         t=t1+ran_mt(vstate)*(0.5*M_PI-t1);
         u=a*tan(t)+x;
        }

       r=ran_mt(vstate);
       ratio=exp(-u*u);
       if(u>=d) ratio=ratio/emd2;
     }
  while(r > ratio);

  *uz=sign*u;
}
