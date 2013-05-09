/*
 * SUMMARY:      CalcBlowingSnow.c - Calculate energy of sublimation from blowing snow
 * USAGE:        Part of VIC
 *
 * AUTHOR:       Laura Bowling
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              lbowling@u.washington.edu
 * ORIG-DATE:     3-Feb-2002 
 * LAST-MOD: 
 * DESCRIPTION:  Calculate blowing snow
 * DESCRIP-END.
 * FUNCTIONS:    CalcBlowingSnow()
 * COMMENTS:     
 * Modifications:
 *   2005-Aug-05 Merged with Laura Bowling's updated code to fix the following problems:
 *	         - Error in array declaration line 373 and 375
 *	         - Added calculation of blowing snow transport.  This is not really used
 *	           currently, but it is something Laura is experimenting with.
 *	         - Fixed RH profile.
 *	         - Fixed vertical integration functions.
 *										TJB
 *   2004-Oct-04 Merged with Laura Bowling's updated lake model code.		TJB
 *   2007-Apr-03 Module returns an ERROR value that can be trapped in main      GCT
 *   2011-Nov-04 Updated mtclim functions to MTCLIM 4.3.			TJB
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <mtclim_constants_vic.h>

static char vcid[] = "$Id$";

#define GRAMSPKG 1000.
#define CH_WATER 4186.8e3
#define JOULESPCAL     4.1868   /* Joules per calorie */
#define Ka  .0245187            /* thermal conductivity of air (W/mK) */
#define CSALT 0.68              /* saltation constant m/s */
#define UTHRESH 0.25            /* threshold shear velocity m/s */
#define KIN_VIS 1.3e-5          /* Kinemativ viscosity of air (m2/s) */
#define MAX_ITER 100             /* Max. iterations for numerical integration */
#define K 5
#define MACHEPS 1.0e-6          /* Accuracy tolerance for numerical integration */
#define SETTLING 0.3            /* Particle settling velocity m/s */
#define UPARTICLE 2.8*UTHRESH   /* Horizontal particle velocity m/s */
                                /* After Pomeroy and Gray (1990) */
#define NUMINCS 10              /* Number of prob intervals to solve for wind. */
#define LAPLACEK 1.                      /* Fit parameter of the laplace distribution. */ 
#define SIMPLE 0               /* SBSM (1) or Liston & Sturm (0) mass flux */
#define SPATIAL_WIND 1         /* Variable (1) or constant (0) wind distribution. */
#define VAR_THRESHOLD 1         /* Variable (1) or constant (0) threshold shear stress. */
#define FETCH 1               /* Include fetch dependence (1). */
#define CALC_PROB 1             /* Variable (1) or constant (0) probability of occurence. */

double qromb(double (*sub_with_height)(), double es, double Wind, double AirDens, double ZO, 
	     double EactAir, double F, double hsalt, double phi_r, double ushear, double Zrh, 
	     double a, double b);
double (*funcd)(double z,double es,  double Wind, double AirDens, double ZO,          
			  double EactAir,double F, double hsalt, double phi_r,         
			  double ushear, double Zrh);
double sub_with_height(double z,double es,  double Wind, double AirDens, double ZO,          
			  double EactAir,double F, double hsalt, double phi_r,         
			  double ushear, double Zrh);
double transport_with_height(double z,double es,  double Wind, double AirDens, double ZO,
				double EactAir,double F, double hsalt, double phi_r,         
				double ushear, double Zrh);
double rtnewt(double x1, double x2, double xacc, double Ur, double Zr);
void get_shear(double x, double *f, double *df, double Ur, double Zr);
double get_prob(double Tair, double Age, double SurfaceLiquidWater, double U10);
double get_thresh(double Tair, double SurfaceLiquidWater, double Zo_salt, int flag);
void shear_stress(double U10, double ZO,double *ushear, double *Zo_salt, double utshear);
double CalcSubFlux(double EactAir, double es, double Zrh, double AirDens, double utshear, 
		   double ushear, double fe, double Tsnow, double Tair, double U10, 
		   double Zo_salt, double F, double *Transport);

/*****************************************************************************
  Function name: CalcBlowingSnow()

  Purpose      : Calculate sublimation from blowing snow

  Required     :  double Dt;                     Model time step (hours) 
  double Tair;                    Air temperature (C) 
  int LastSnow;                   Time steps since last snowfall. 
  double SurfaceLiquidWater;      Liquid water in the surface layer (m) 
  double Wind;                    Wind speed (m/s), 2 m above snow 
  double Ls;                      Latent heat of sublimation (J/kg) 
  double AirDens;                 Density of air (kg/m3) 
  double Lv;                      Latent heat of vaporization (J/kg3) 
  double Press;                   Air pressure (Pa) 
  double EactAir;                 Actual vapor pressure of air (Pa) 
  double ZO;                      Snow roughness height (m)
  Returns      : BlowingMassFlux

  Modifies     : 

  Comments     : Called from SnowPackEnergyBalance
    Reference:  

*****************************************************************************/
double CalcBlowingSnow( double Dt, 
			double Tair,
			int LastSnow,
			double SurfaceLiquidWater,
			double Wind,
			double Ls,
			double AirDens,
			double Press,
			double EactAir,
			double ZO,
			double Zrh,
			double snowdepth,
			float lag_one,
			float sigma_slope,
			double Tsnow, 
			int iveg, 
			int Nveg, 
			float fe,
			double displacement,
			double roughness,
			double *TotalTransport)

{
  /* Local variables: */

  double Age;
  double U10, Uo, prob_occurence;
  double es, Ros, F;
  double SubFlux;
  double Diffusivity;
  double ushear, Qsalt, hsalt, phi_s, psi_s;
  double Tk;
  double Lv;
  double T, ztop;
  double ut10, utshear;
  int p;
  double upper, lower, Total;
  double area;
  double sigma_w;
  double undersat_2;
  double b, temp2; /* SBSM scaling parameter. */
  double temp, temp3;
  double Zo_salt;
  double ratio, wind10;
  double Uveg, hv, Nd;
  double Transport;
  int count=0;

  Lv = (2.501e6 - 0.002361e6 * Tsnow);
  /*******************************************************************/
  /* Calculate some general variables, that don't depend on wind speed. */

  /* Age in hours */
  Age = LastSnow*(Dt);
      
  /* Saturation density of water vapor, Liston A-8 */
  es = svp(Tair);

  Tk = Tair  + KELVIN;
    
  Ros = 0.622*es/(287*Tk);
  
  /* Diffusivity in m2/s, Liston eq. A-7 */
  Diffusivity = (2.06e-5) * pow(Tk/273.,1.75);

  // Essery et al. 1999, eq. 6 (m*s/kg)
  F = (Ls/(Ka*Tk))*(Ls*MW/(R*Tk) - 1.);
  F += 1./(Diffusivity*Ros);

  /* grid cell 10 m wind speed = 50th percentile wind */
  /* Wind speed at 2 m above snow was passed to this function. */

  wind10 = Wind*log(10./ZO)/log((2+ZO)/ZO);
  //  fprintf(stderr,"wind=%f, Uo=%f\n",Wind, Uo);

  /* Check for bare soil case. */
  if(iveg == Nveg) {
    fe = 1500;
    sigma_slope = .0002;
  }
  // sigma_w/uo:
  ratio = (2.44 - (0.43)*lag_one)*sigma_slope;
  //  sigma_w = wind10/(.69+(1/ratio));
  //  Uo = sigma_w/ratio;
  
  sigma_w = wind10*ratio;
  Uo = wind10;
  
  /*********** Parameters for roughness above snow. *****************/
    hv = (3./2.)*displacement;
    Nd = (4./3.)*(roughness/displacement);

  /*******************************************************************/
  /** Begin loop through wind probability function.                  */

  Total = 0.0;
  *TotalTransport = 0.0;
  area = 1./NUMINCS;
  
  if(snowdepth > 0.0) {
    if(SPATIAL_WIND && sigma_w != 0.) {
      for(p= 0; p< NUMINCS; p++) {
      
	SubFlux = lower = upper = 0.0;
      /* Find the limits of integration. */
      if(p==0) {
	lower = -9999;
	upper = Uo + sigma_w*log(2.*(p+1)*area);
      }
      else if(p > 0 && p < NUMINCS/2) {
	lower = Uo + sigma_w*log(2.*(p)*area);
	upper = Uo + sigma_w*log(2.*(p+1)*area);
      }
      else if(p < (NUMINCS-1) && p >= NUMINCS/2) {
	lower = Uo - sigma_w*log(2.-2.*(p*area));
	upper = Uo - sigma_w*log(2.-2.*((p+1.)*area));
      }
      else if(p == NUMINCS-1) {
	lower =  Uo - sigma_w*log(2.-2.*(p*area));
	upper = 9999;
      }

      if(lower > upper)  {/* Could happen if lower > Uo*2 */
	lower = upper;
	fprintf(stderr,"Warning: Error with probability boundaries in CalcBlowingSnow()\n");
      }


      /* Find expected value of wind speed for the interval. */
      U10 = Uo;
      if(lower >= Uo ) 
	U10 = -0.5*((upper+sigma_w)*exp((-1./sigma_w)*(upper - Uo))
		    - (lower+sigma_w)*exp((-1./sigma_w)*(lower - Uo)))/area;
      else if(upper <= Uo )
	U10 = 0.5*((upper-sigma_w)*exp((1./sigma_w)*(upper - Uo))
		   - (lower-sigma_w)*exp((1./sigma_w)*(lower - Uo)))/area;
      else {
	fprintf(stderr,"ERROR in CalcBlowingSnow.c: Problem with probability ranges\n");
	fprintf(stderr,"  Increment = %d, integration limits = %f - %f\n",p,upper, lower);
        return ( ERROR );
      }
   
      if(U10 < 0.4)
	U10 = .4;

      if(U10 > 25.) U10 = 25.;
      /*******************************************************************/
      /* Calculate parameters for probability of blowing snow occurence. */
      /* ( Li and Pomeroy 1997) */

      if(snowdepth < hv) {
		Uveg = U10/sqrt(1.+ 170*Nd*(hv - snowdepth));
      }
      else
	Uveg = U10;
  
      //  fprintf(stderr, "Uveg = %f, U10 = %f\n",Uveg, U10);

      prob_occurence = get_prob(Tair, Age, SurfaceLiquidWater, Uveg);

      //   printf("prob=%f\n",prob_occurence);
     
      /*******************************************************************/
      /* Calculate threshold shear stress. Send 0 for constant or  */
      /* 1 for variable threshold after Li and Pomeroy (1997)      */

      utshear = get_thresh(Tair, SurfaceLiquidWater, ZO, VAR_THRESHOLD);
   
      /* Iterate to find actual shear stress during saltation. */

      shear_stress(U10, ZO, &ushear, &Zo_salt, utshear);
      
      if(ushear > utshear) {
	
	SubFlux = CalcSubFlux(EactAir, es, Zrh, AirDens, utshear,ushear, fe, Tsnow, 
			      Tair, U10, Zo_salt, F, &Transport);
      }
      else {
	SubFlux=0.0;
	Transport = 0.0;
      }
 
      Total += (1./NUMINCS)*SubFlux*prob_occurence;
      *TotalTransport += (1./NUMINCS)*Transport*prob_occurence;

      }
   
    }
    else {
      U10=Uo;
       /*******************************************************************/
      /* Calculate parameters for probability of blowing snow occurence. */
      /* ( Li and Pomeroy 1997) */
  
      if(snowdepth < hv)
	Uveg = U10/sqrt(1.+ 170*Nd*(hv - snowdepth));
      else
	Uveg = U10;

      prob_occurence = get_prob(Tair, Age, SurfaceLiquidWater, Uveg);
    
      /*******************************************************************/
      /* Calculate threshold shear stress. Send 0 for constant or  */
      /* 1 for variable threshold after Li and Pomeroy (1997)      */

      utshear = get_thresh(Tair, SurfaceLiquidWater, ZO, VAR_THRESHOLD);

      /* Iterate to find actual shear stress during saltation. */
      
      shear_stress(Uo, ZO, &ushear, &Zo_salt, utshear);

      if(ushear > utshear) {
	SubFlux = CalcSubFlux(EactAir, es, Zrh, AirDens, utshear,ushear, fe, Tsnow, 
			      Tair, Uo, Zo_salt, F, &Transport);
      }
      else {
	SubFlux=0.0;
	Transport = 0.0;
      }
      Total = SubFlux*prob_occurence;
      *TotalTransport = Transport*prob_occurence;
    }
  }

  if(Total < -.00005)
    Total = -.00005;
 
  return Total;
  
}

double qromb(double (*funcd)(), double es, double Wind, double AirDens, double ZO, 
	     double EactAir, double F, double hsalt, double phi_r, double ushear, double Zrh, 
	     double a, double b)
     // Returns the integral of the function func from a to b.  Integration is performed 
     // by Romberg's method:  Numerical Recipes in C Section 4.3
{
  void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
  double trapzd(double (*funcd)(), double es, double Wind, double AirDens, 
		double ZO, double EactAir, double F, double hsalt, double phi_r, 
		double ushear, double Zrh, double a, double b, int n);

  double ss, dss;
  double s[MAX_ITER+1], h[MAX_ITER+2];
  int j;

  h[1] = 1.0;
  for(j=1; j<=MAX_ITER; j++) {
    s[j]=trapzd(funcd,es, Wind, AirDens, ZO, EactAir, F, hsalt, phi_r, 
		ushear, Zrh, a,b,j);
    if(j >= K) {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) <= MACHEPS*fabs(ss)) return ss;
    }
    h[j+1]=0.25*h[j];
  }
  nrerror("Too many steps in routine qromb");
  return 0.0;
}

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  int i, m, ns;
  double den, dif, dift, ho, hp, w;
  double *c,*d;

  ns=1;
  dif=fabs(x-xa[1]);
  c=(double *)malloc((size_t) ((n+1)*sizeof(double)));
  if(!c) nrerror("allocation failure in vector()");
  d=(double *)malloc((size_t) ((n+1)*sizeof(double)));
  if(!d) nrerror("allocation failure in vector()");

  for (i=1; i<=n; i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for(m=1;m<n;m++) {
    for(i=1; i<=n-m; i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
      den = w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free(d);
  free(c);
}



double trapzd(double (*funcd)(), double es, double Wind, double AirDens, double ZO, 
	      double EactAir, double F, double hsalt, double phi_r, double ushear, 
	      double Zrh, double a, double b, int n)
{
  double x, tnm, sum, del;
  static double s;
  int it, j;

  if (n==1) {
    return (s=0.5*(b-a)*((*funcd)(a, es, Wind, AirDens, ZO, EactAir, F, hsalt, 
					 phi_r, ushear, Zrh) + 
			 (*funcd)(b,es, Wind, AirDens, ZO, EactAir, F, hsalt, 
					 phi_r, ushear, Zrh)));
  }
  else {
    for (it=1, j=1; j<n-1; j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for(sum=0.0, j=1; j<=it; j++, x+=del ) sum += (*funcd)(x,es, Wind, AirDens, ZO, EactAir, 
							   F, hsalt, phi_r, ushear, Zrh);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}

double rtnewt(double x1, double x2, double acc, double Ur, double Zr)
{
  int j;
  double df, dx, dxold, f, fh, fl;
  double temp, xh, xl, rts;

    get_shear(x1,&fl,&df, Ur, Zr);
    get_shear(x2,&fh,&df, Ur, Zr);

    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
      fprintf(stderr, "Root must be bracketed in rtnewt.\n");
      exit(0);
    }

    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    if (fl < 0.0) {
      xl=x1;
      xh=x2; }
    else {
      xh=x1;
      xl=x2;
    }
    rts=0.5*(x1+x2);
    dxold=fabs(x2-x1);
    dx=dxold;
    get_shear(rts,&f,&df, Ur, Zr);
    for(j=1; j<=MAX_ITER; j++) {
      if((((rts-xh)*df-f)*((rts-x1)*df-f) > 0.0)
	 || (fabs(2.0*f) > fabs(dxold*df))) {
	dxold=dx;
	dx=0.5*(xh-xl);
	rts=xl+dx;
	if (xl == rts) return rts; }
      else {
	dxold=dx;
	dx=f/df;
	temp=rts;
	rts -= dx;
	if (temp == rts) return rts;
      }
      if(fabs(dx) < acc) return rts;
      // if(rts < .025) rts=.025;
      get_shear(rts,&f,&df, Ur, Zr);
       if(f<0.0)
	 xl=rts;
       else 
	 xh = rts;
    }
    fprintf(stderr, "Maximum number of iterations exceeded in rtnewt.\n");
    return 0.0;
}

void get_shear(double x, double *f, double *df, double Ur, double Zr) {

  *f = log(2.*G_STD*Zr/.12)+log(1/(x*x)) - von_K*Ur/x;
  *df = von_K*Ur/(x*x) - 2./x;
}

/*****************************************************************************
  Function name: sub_with_height()

  Purpose      : Calculate the sublimation rate for a given height above the boundary layer.

  Required     :
    double z               - Height of solution (m) 
    double Tair;           - Air temperature (C) 
    double Wind;           - Wind speed (m/s), 2 m above snow 
    double AirDens;        - Density of air (kg/m3) 
    double ZO;             - Snow roughness height (m)
    double EactAir;        - Actual vapor pressure of air (Pa) 
    double F;              - Denominator of dm/dt          
    double hsalt;          - Height of the saltation layer (m)
    double phi_r;          - Saltation layer mass concentration (kg/m3)
    double ushear;         - shear velocity (m/s) 
    double Zrh;            - Reference height of humidity measurements

  Returns      :
   double f(z)             - Sublimation rate in kg/m^3*s

  Modifies     : none

  Comments     :  Currently neglects radiation absorption of snow particles.
*****************************************************************************/
double sub_with_height(double z,
		       double es,          
		       double Wind,          
		       double AirDens,     
		       double ZO,          
		       double EactAir,        
		       double F,          
		       double hsalt,         
		       double phi_r,         
		       double ushear,        
		       double Zrh)
{

  /* Local variables */
  double Rrz, ALPHAz, Mz;
  double Rmean, terminal_v, fluctuat_v;
  double Vtz, Re, Nu;
  double sigz, dMdt;
  double temp;
  double psi_t, phi_t;

 
  //  Calculate sublimation loss rate (1/s)
  Rrz = 4.6e-5* pow (z, -.258);
  ALPHAz = 4.08 + 12.6 * z;
  Mz = (4./3.) * PI * ice_density * Rrz * Rrz * Rrz *(1. +(3./ALPHAz) + (2./(ALPHAz*ALPHAz)));

  Rmean = pow((3.*Mz)/(4.*PI*ice_density),1./3.);

  // Pomeroy and Male 1986
  terminal_v = 1.1e7 * pow(Rmean,1.8);

  // Pomeroy (1988)
  fluctuat_v = 0.005 * pow(Wind, 1.36);

  // Ventilation velocity for turbulent suspension Lee (1975)
  Vtz = terminal_v + 3.*fluctuat_v*cos(PI/4.);

  Re = 2. * Rmean * Vtz / KIN_VIS;
  Nu = 1.79 + 0.606 * pow(Re, 0.5);

  // LCB: found error in rh calc, 1/20/04, check impact
  // sigz = ((EactAir/es) - 1.) * (1 + .027*log(z) - .027*log(Zrh));
  sigz = ((EactAir/es) - 1.) * (1.019 + .027*log(z));

  dMdt = 2 * PI * Rmean * sigz * Nu / F;
  // sublimation loss rate coefficient (1/s)

  psi_t = dMdt/Mz;

  // Concentration of turbulent suspended snow Kind (1992)
  
  temp = (0.5*ushear*ushear)/(Wind*SETTLING);
  phi_t = phi_r* ( (temp + 1.) * pow((z/hsalt),(-1.*SETTLING)/(von_K*ushear)) - temp );
  
  return psi_t * phi_t;
}

/*******************************************************************/
/* Calculate parameters for probability of blowing snow occurence. */
/* ( Li and Pomeroy 1997) */
/*******************************************************************/

double get_prob(double Tair, double Age, double SurfaceLiquidWater, double U10) 
{
  double mean_u_occurence;
  double sigma_occurence;
  double prob_occurence;

  if(CALC_PROB) {
    if(SurfaceLiquidWater < 0.001) {
      mean_u_occurence = 11.2 + 0.365*Tair + 0.00706*Tair*Tair+0.9*log(Age);
      sigma_occurence = 4.3 + 0.145*Tair + 0.00196*Tair*Tair;
      
      prob_occurence = 1./(1.+exp(sqrt(PI)*(mean_u_occurence-U10)/sigma_occurence));
    }
    else {
      mean_u_occurence = 21.;
      sigma_occurence = 7.;
	
      prob_occurence = 1./(1.+exp(sqrt(PI)*(mean_u_occurence-U10)/sigma_occurence));
    }

    if(prob_occurence < 0.0)
      prob_occurence = 0.0;
    if(prob_occurence > 1.0)
      prob_occurence = 1.0;
  }
  else
    prob_occurence = 1.;

  return prob_occurence;
}

double get_thresh(double Tair, double SurfaceLiquidWater, double Zo_salt, int flag) 
{
  double ut10;
  double utshear;

  if(SurfaceLiquidWater < 0.001) {
    // Threshold wind speed after Li and Pomeroy (1997)
    ut10 = 9.43 + .18 * Tair + .0033 * Tair*Tair;
  }
  else {
    
    // Threshold wind speed after Li and Pomeroy (1997)
    ut10 = 9.9;
  }

  if(flag) {
    // Variable threshold, Li and Pomeroy 1997
    utshear = von_K * ut10 / log(10./Zo_salt);
  }
  // Constant threshold, i.e. Liston and Sturm
  else
    utshear = UTHRESH;
  
  return utshear;
}


void shear_stress(double U10, double ZO,double *ushear, double *Zo_salt, double utshear)
{
  double umin, umax, xacc;
  double fl, fh, df;

  /* Find min & max shear stress to bracket value. */
  umin = utshear;
  umax = von_K*U10;
  xacc = 0.10*umin;

  /* Check to see if value is bracketed. */
   get_shear(umin,&fl,&df, U10, 10.);
   get_shear(umax,&fh,&df, U10, 10.);

    if(fl < 0.0 && fh < 0.0) { 
      fprintf(stderr, "Solution in rtnewt surpasses upper boundary.\n");
      fprintf(stderr, "fl(%f)=%f, fh(%f)=%f\n",umin, fl, umax, fh);
      exit(0);
    }
    
    if(fl > 0.0 && fh > 0.0) {
      //      fprintf(stderr, "No solution possible that exceeds utshear.\n");
      //   fprintf(stderr, "utshear=%f, u10=%f\n",utshear, U10);
     
      *Zo_salt = ZO;
      *ushear = von_K * U10 / log(10./ZO);  
    }
    else {
      /* Iterate to find actual shear stress. */
      *ushear = rtnewt (umin, umax, xacc, U10, 10.);
      *Zo_salt = 0.12 *(*ushear) * (*ushear) / (2.* G_STD);
    }
    
}

double CalcSubFlux(double EactAir, double es, double Zrh, double AirDens, double utshear, 
		   double ushear, double fe, double Tsnow, double Tair, double U10, 
		   double Zo_salt, double F, double *Transport)

{
  float b, undersat_2;
  double SubFlux;
  double Qsalt, hsalt;
  double phi_s, psi_s;
  double T, ztop;
  double particle;
  double saltation_transport;
  double suspension_transport;

  SubFlux=0.0;
  particle = utshear*2.8;
  // SBSM:
  if(SIMPLE) {   
    b=.25;
    if(EactAir >= es)
      undersat_2 = 0.0;
    else
      undersat_2 = ((EactAir/es)-1.)*(1.-.027*log(Zrh)+0.027*log(2)); 
    //  fprintf(stderr,"RH=%f\n",EactAir/es);
    SubFlux =  b*undersat_2* pow(U10, 5.) / F;
  }

  else {
    //  Sublimation flux (kg/m2*s) = mass-concentration * sublimation rate * height
    //  for both the saltation layer and the suspension layer
	  
    //  Saltation layer is assumed constant with height
    // Maximum saltation transport rate (kg/m*s)
    // Liston and Sturm 1998, eq. 6
    Qsalt = ( CSALT * AirDens /  G_STD ) * (utshear / ushear) * (ushear*ushear - utshear*utshear); 
    if(FETCH)
      Qsalt *= (1.+(500./(3.*fe))*(exp(-3.*fe/500.)-1.));
    
    // Liston and Sturm (1998)
    // hsalt = 1.6 * ushear * ushear / ( 2. * G_STD );
	
    // Pomeroy and Male (1992)
    hsalt = 0.08436*pow(ushear,1.27);

    // Saltation layer mass concentration (kg/m3)
    phi_s = Qsalt / (hsalt * particle);

    T = 0.5*(ushear*ushear)/(U10*SETTLING);
    ztop = hsalt*pow(T/(T+1.), (von_K*ushear)/(-1.*SETTLING));

    if(EactAir >= es) {
      SubFlux = 0.0;
    }
    else 
      { 
	// Sublimation loss-rate for the saltation layer (s-1)
	psi_s =  sub_with_height(hsalt/2.,es, U10, AirDens, Zo_salt, EactAir, F, hsalt, 
				 phi_s, ushear, Zrh);
    
	// Sublimation from the saltation layer in kg/m2*s
	SubFlux = phi_s*psi_s*hsalt;
    
	//  Suspension layer must be integrated
	SubFlux += qromb(sub_with_height, es, U10, AirDens, Zo_salt, EactAir, F, hsalt,
			 phi_s, ushear, Zrh, hsalt, ztop);
      }

    // Transport out of the domain by saltation Qs(fe) (kg/m*s), eq 10 Liston and Sturm
    saltation_transport = Qsalt*(1-exp(-3.*fe/500.));

    // Transport in the suspension layer
    suspension_transport = qromb(transport_with_height, es, U10, AirDens, Zo_salt, 
    				 EactAir, F, hsalt, phi_s, ushear, Zrh, hsalt, ztop);

    // Transport at the downstream edge of the fetch in kg/m*s
    *Transport = (suspension_transport + saltation_transport);
    if(FETCH) 
      *Transport /= fe;
  }

  return SubFlux;
}

/*****************************************************************************
  Function name: transport_with_height()

  Purpose      : Calculate the transport rate for a given height above the boundary layer.

  Required     :
    double z               - Height of solution (m) 
    double Tair;           - Air temperature (C) 
    double Wind;           - Wind speed (m/s), 2 m above snow 
    double AirDens;        - Density of air (kg/m3) 
    double ZO;             - Snow roughness height (m)
    double EactAir;        - Actual vapor pressure of air (Pa) 
    double F;              - Denominator of dm/dt          
    double hsalt;          - Height of the saltation layer (m)
    double phi_r;          - Saltation layer mass concentration (kg/m3)
    double ushear;         - shear velocity (m/s) 
    double Zrh;            - Reference height of humidity measurements

  Returns      :
   double f(z)             - Transport rate in kg/m^2*s

  Modifies     : none

*****************************************************************************/
double transport_with_height(double z,
		       double es,          
		       double Wind,          
		       double AirDens,     
		       double ZO,          
		       double EactAir,        
		       double F,          
		       double hsalt,         
		       double phi_r,         
		       double ushear,        
		       double Zrh)
{

  /* Local variables */
  double u_z;
  double temp;
  double phi_t;
  
  // Find wind speed at current height

  u_z = ushear*log(z/ZO)/von_K;

  // Concentration of turbulent suspended snow Kind (1992)
  
  temp = (0.5*ushear*ushear)/(Wind*SETTLING);
  phi_t = phi_r* ( (temp + 1.) * pow((z/hsalt),(-1.*SETTLING)/(von_K*ushear)) - temp );
  
  return u_z * phi_t;
}
