/********************************************************************************
  filename  : aurad.c
  purpose   : Calculate daily net radiation
  interface : see below
  programmer: see below 
  date      : August 4, 1995
  changes   : changed include file from bgc.h to rad_and_vpd.h (08/04/1995)
  notes     : original version received August 4, 1995 from John Kimball,
              NTSG, School of Forestry, University of Montana

	      unexplained (in list below) variables:
	      tranday - daily total transmittance

	      references:
	      Buffo, J., L. Fritschen, and J. Murphy, Direct solar radiation on 
	      various slopes from 0 to 60 north latitude.  USDA Forest Service
	      Paper Research Paper PNW-142, 1972.
	      Swift, L. W. Jr., Algortihm for solar radiation on mountain slopes.
	      Water Resources Research, v.12, p.108-112, 1976.
	      
	      for daily total transmittance:
	      Bristow, K. L., and G. S. Campbell, On the relationship between 
	      incoming solar radiation and daily maximum and minimum temperature.
	      Agricultural and Forest Meteorology, v.31, p.159-166, 1984.
********************************************************************************/



/*    ****************************************************** 
c     *      AURAD.F                                       * 
c     ****************************************************** 

CM    CONTAINS: fltrad, shrad 
c 
c     ****************************************************** 
c     *      subroutine fltrad                             * 
c     *                                                    * 
C     ****************************************************** 
CF    FILE    : aurad.f 
CM    CONTAINS: fltrad 
CMA   AUTHOR  : jc coughlan 
CMS   SYSTEM  : autosys 
CMV   VERSION : #2 
CMX   STATUS  : 
CVF   VFILE   : 
CVI   VINPUT  : 
CVO   VOUTPUT : 
CR    RETURNS : 
CD    DESCRIP : given below 
CH    HISTORY : CREATED 5/2/1989 
ch              modified 5/ /1990 r nemani jp patterson    */
/* 
cd     flat surface radiation computation. 
cd     slope = 0    aspect = 0 
cd     this subroutine computes incident shortwave radiation and 
cd     net shortwave radiation for any given day based on 
cd     sun-earth geometry, transmissivity; 
cd     derived from: buffo et al 1972; swift 1976. 
cd 
cd     variable list: 
cd     sol=solar constant derived from solcon array 
cd     am=optical air mass for angles > 21 degrees 
cd     a=optical air mass array for angles between 
cd     0 and 21 degrees above horizon 
cd     decl=declination 
cd     jday=day of year 
cd     h=angle of sun from solar noon 
cd     tran=transmissivity constant 
cd     tram=transmissivity corrected for air mass 
cd     nnh=calculation interval in seconds 
cd     nc=seconds in one day (24 hours) 
cd     n=number of intervals of length nnh in one day 
cd     dt=direct solar perpendicular to sun on the 
cd     outside of atmosphere for interval (kj/m**2) 
cd     etf=direct solar on outside of atmosphere 
cd     parallel to earths surface for interval 
cd     grad=total daily radiation at given location (kj/m**2) 
cd     hrad=direct solar on earths surface (flat)       " 
cd     tdif=total daily diffuse radiation               " 
cd     difrad=diffuse on slope for interval             " 
cd     drad=direct on slope for interval                " 
cd     cza=cosine zenith angle 
cd     cbsa=cosine beam slope angle 
cd     globf=global radiation, determining diffuse 
cd     diffl=diffuse on flat for interval 
cd     dayl=daylength (hours) 
cd     albdo=albedo 					never used
cd     arad=absorbed radiation, w/m2 			removed
cd     radn=incident radiation  w/m2 
cd    (* . . . . . . . . . . . . . . . . . . . . . . . . . . . */ 
   
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>


static char vcid[] = "$Id$";

/* Since K&R C doesn't allow multiple assignments of arrays within functions
   this must be done as an external variable, but that's ok since this way
   both functions can use the data   */

    double  aa[22]  = 
	{0.0, 2.90,3.05,3.21,3.39, 3.69, 3.82, 4.07, 4.37, 4.72, 5.12,
	 5.60,6.18,6.88,7.77,8.90,10.39,12.44,15.36,19.79,26.96,30.00};

    double  dec[47] =
	{ 0.0,-23.0,-22.0,-21.0,-19.0,-17.0,-15.0,-12.0,-9.0,  -6.0,  -3.0,  
	  0.0,  3.0,  6.0,  9.0, 12.0, 14.0, 17.0, 19.0, 21.0, 22.0,  23.0,
	 23.5, 23.5, 23.0, 21.5, 20.0, 18.0, 16.0, 14.0, 12.0,  9.0,   6.0,
	  3.0,  0.0, -3.0, -6.0, -9.0,-12.0,-15.0,-17.0,-19.0, -21.0,-22.0,
	-23.0,-23.5,-23.5};

    double  solcon[13] =
	{ 0.0, 1.445,1.431,1.410,1.389,1.368,1.354,
	 1.354,1.375,1.403,1.424,1.438,1.445};


double fltrad(double  sehorz, 
	      double  swhorz, 
	      double  slat, 
	      int     jday, 
	      double  tranday,
	      double *dayl) 
{
  int    i, nnh, nc, n, nh, ml, mo, idec;
  
  double  decl,h, grad, radn; 
  double  xtran, sol, tdif, tdrad, dayl2;
  double  cosdec, cosdla, sindla, sindec, coshh;
  double  cza, am, tram, dt, etf, hrad, globf, diffl, cbsa, drad, se;
  double  difrad;
  
  xtran=tranday;
  nnh=3600;
  nc=86400; 
  n = nc/nnh + 1;
  dayl[0]=0.0; 
  

  mo = jday/30 + 1; 
  mo = (mo<12) ? mo:12;
  
  /*   solcon array is in units of kw/m**2 */

  sol=solcon[mo]; 
  idec= (int) (1.0 + (double) jday/8.0);
  decl = dec[idec];
  grad = 0.0; 
  tdif=0.0; 
  tdrad=0.0; 
  nh=0; 
  dayl2=0.0; 
  
  cosdec = cos(decl*DtoR); 
  cosdla = cos(slat*DtoR); 
  sindla = sin(slat*DtoR); 
  sindec = sin(decl*DtoR); 
  
  for (i=1; i<=n; i++) { 
    nh = nh+nnh;
    
    /*         determine angle from solar noon   */
    
    h=(double) (nh-43200)*0.0041667;
    coshh   = cos(h*DtoR); 
    cza=cosdec*cosdla*coshh+sindec*sindla;
    
    if (cza > 0.0) { 
      
      /* (* daylength based on solar elevation above a flat horizon *) */
      
      dayl2=dayl2+((double) (nnh)/3600.0);
      
      /*             (* determine optical air mass *) */
      
      am=1.0/(cza+0.0000001);
      
      if (am > 2.9) { 
	ml= ((int) (acos(cza)/0.0174533)) - 69;
	
	if(ml < 1)
	  ml=1;
	
	if(ml > 21)
	  ml=21; 
	
	am = aa[ml];

      }
      
      tram = pow(xtran,am);
      dt   = sol* (double) nnh; 
      etf=cza*dt; 
      hrad=etf*tram; 
      dt=dt*tram;
      globf=sqrt(hrad*etf); 
      diffl=globf*(1.0-globf/etf); 
   
      /*             (* compute flat surface angle */
      
      cbsa = cosdla*cosdec*coshh + sindla*sindec;
      
      /*     (* compute direct incident radiation if solar angle is positive */
      
      if(cbsa >= 0.0) { 
	drad=cbsa*dt;
	
	/*  (* the following three lines computes a topographic 
	    c                 reduction of direct radiation 
	    c                 ehe = east horizon elevation (degrees) 
	    c                 whe = west horizon elevation (degrees)  */
	
	se = 1.5715-(acos(cza));
	
	if (h < 0.0 && se < sehorz *DtoR) 
	  drad=0.0; 
	
	if(h > 0.0 && se < swhorz *DtoR)
	  drad=0.0; 
      }
      else 
	drad=0; 
      
      /*             the code below is a simplification of the 
		     c             commented line assuming a flat surface.   */
      
      difrad = diffl;
      
      /*             difrad=diffl*(cos(dslop*.5*DtoR))**2.)    */
   
      dayl[0]	=dayl[0]+((double) (nnh)/3600.0); 
      grad	=grad+drad+difrad; 
      tdif	=tdif+difrad; 
      tdrad	=tdrad+drad; 
    } 
  }
  
  /*    averad=(grad/(dayl*3600.0));   */
  
  radn = grad;
  /*    arad = averad;   */
  
  return((double) radn); 
}




/*    ****************************************************** 
c     *      subroutine shrad                              * 
c     *                                                    * 
C     ****************************************************** 
CF    FILE    : aurad.f 
CM    CONTAINS: shrad 
CMA   AUTHOR  : jc coughlan 
CMS   SYSTEM  : autosys 
CMV   VERSION : #1 
CMX   STATUS  : 
CVF   VFILE   : 
CVI   VINPUT  : 
CVO   VOUTPUT : 
CR    RETURNS :  
CD    DESCRIP : given below 
CH    HISTORY : CREATED 5/2/1989 
ch              modified 5/ /1990 r nemani and jp patterson 
c 

cd     this subroutine computes incident shortwave radiation and 
cd     net shortwave radiation for any given day based on surface 
cd     characteristics, sun-earth geometry, transmissivity; based on 
cd     buffo et al 1972; swift 1976. 
cd 
cd     variable list: 
cd     sol=solar constant derived from solcon array 
cd     am=optical air mass for angles > 21 degrees 
cd     a=optical air mass array for angles between 
cd     0 and 21 degrees above horizon 
cd     decl=declination 
cd     jday=day of year 
cd     asp=aspect in degrees 
cd     dslop=slope in degrees !!! 
cd     h=angle of sun from solar noon 
cd     tran=transmissivity constant 
cd     tram=transmissivity corrected for air mass 
cd     nnh=calculation interval in seconds 
cd     nc=seconds in one day (24 hours) 
cd     n=number of intervals of length nnh in one day 
cd     dt=direct solar perpendicular to sun on the 
cd     outside of atmosphere for interval (kj/m**2) 
cd     etf=direct solar on outside of atmosphere 
cd     parallel to earths surface for interval 
cd     grad=total daily radiation at given location (kj/m**2) 
cd     hrad=direct solar on earths surface (flat)       " 
cd     tdif=total daily diffuse radiation               " 
cd     difrad=diffuse on slope for interval             " 
cd     drad=direct on slope for interval                " 
cd     cza=cosine zenith angle 
cd     cbsa=cosine beam slope angle 
cd     globf=global radiation, determining diffuse 
cd     diffl=diffuse on flat for interval 
cd     dayl=daylength (hours) 
cd     albdo=albedo 					never used
cd     arad=absorbed radiation, w/m2 			removed
cd     radn=incident radiation  w/m2 
cd   */

double shrad (double asp, double dslop, double slat, double sehorz, double swhorz, int jday, double tranday)
{
  int   i, n, nh, nc, nnh, ml, mo, idec;
  
  double decl, h, grad, radn, tempf;
  double xtran, dayl, sol, tdif, tdrad, dayl2;
  double cosdec, sindsl, sinasp, cosdsl, cosdla, cosasp;
  double sindla, sindec, sinhh, coshh, cza, am, tram, dt, etf;
  double hrad, globf, diffl, cbsa, drad, se, difrad;
  
  /*     slope in degrees - no conversion! */

  xtran=tranday; 
  nnh =  3600; 
  nc  = 86400; 
  n = nc/nnh + 1; 
  dayl=0.0; 
  mo = jday/30+1; 
  mo = (mo<12) ? mo:12; 
  
  /*     solcon array is in units of kw/m**2   */
  
  sol  = solcon[mo]; 
  idec = (int) (1.0 + (double) jday/8.0); 
  decl = dec[idec]; 
  grad = 0.0; 
  tdif=0.0; 
  tdrad=0.0; 
  nh = 0; 
  dayl2=0.0; 
  
  tempf= (double) (decl*DtoR);
  cosdec = cos((double) tempf); 
  /*  cosdec = cos((double) decl);  */
  sindsl = sin((double) (dslop*DtoR)); 
  sinasp = sin((double) (asp*DtoR)); 
  cosdsl = cos((double) (dslop*DtoR)); 
  cosdla = cos((double) (slat*DtoR)); 
  cosasp = cos((double) (asp*DtoR)); 
  sindla = sin((double) (slat*DtoR)); 
  sindec = sin((double) (decl*DtoR)); 
  
  for (i=1;i<=n; i++) {
    nh=nh+nnh;
    
    /*  determine angle from solar noon   */
    
    h    = (double) (nh-43200)*0.0041667; 
    sinhh = sin(h*DtoR); 
    coshh = cos(h*DtoR); 
    cza  = cosdec*cosdla*coshh+sindec*sindla;
    
    if(cza > 0.0) { 
      
      /*  daylength based on solar elevation above a flat horizon */
      
      dayl2=dayl2+((double) (nnh)/3600.0);
      
      /*             next 6 lines, determine optical air mass */
      
      am = 1.0/(cza+0.0000001);
      
      if(am > 2.9) { 
	ml= ((int) (acos(cza)/0.0174533)) - 69;
	
	if(ml < 1)
	  ml=1;
	
	if(ml > 21)
	  ml=21; 
	
	am = aa[ml];
      } 
      
      tram = pow(xtran,am);
      dt=sol*(double) nnh;
      etf=cza*dt; 
      hrad=etf*tram; 
      dt=dt*tram; 
      globf=sqrt(hrad*etf); 
      diffl=globf*(1.0-globf/etf); 
      
/*             the following line is a simplification of the 
c             commented computation.  It is kept a documentation for 
c             the derivation of the simpler form.   */
   
      cbsa=(-sindsl)*sinasp*sinhh*cosdec 
	+ (-cosasp*sindsl*sindla+cosdsl*cosdla)*cosdec*coshh
	  +(cosasp*sindsl*cosdla+cosdsl*sindla)*sindec;
      
/*             cbsa=-sin(dslop*DtoR)*sin(asp*DtoR)* 
c             &         sin(h*DtoR)*cos(conv(decl))+(-cos(conv(asp)) 
c          &        *sin(conv(dslop))*sin(conv(dlat))+cos(conv(dslop))* 
c             &        cos(conv(dlat)))*cos(conv(decl)) 
c             &        *cos(conv(h))+(cos(conv(asp))*sin(conv(dslop))* 
c             &            cos(conv(dlat))+cos(conv(dslop))* 
c             &        sin(conv(dlat)))*sin(conv(decl))    */
   
      if (cbsa >= 0.0) { 
	drad=cbsa*dt;
	
/*            the following three lines computes a topographic reduction  of 
   
c                 direct radiation 
c                 ehe = east horizon elevation (degrees) 
c                 whe = west horizon elevation (degrees) */
	
	se = 1.5715-(acos(cza));
	
	if(h < 0.0 && se < sehorz*DtoR) 
	  drad=0.0;
	
	if(h > 0.0 && se < swhorz*DtoR) 
	  drad=0.0;
      }
      else 
	drad=0; 
      
      difrad = diffl * pow(cos(dslop*0.5*DtoR),2.0);
      
/*	    radt   = (drad+difrad)/(double) nnh;   */
   
      dayl=dayl+((double) (nnh)/3600.0); 
      grad=grad+drad+difrad; 
      tdif=tdif+difrad; 
      tdrad=tdrad+drad; 
    } 
  }
  
/*    averad=(grad/(dayl*3600.0));    */
  radn = grad; 
/*    arad = averad;   */
 
  return((double) radn);
}
 
/*cf   ***************************************************** 
cf     *            end of AURAD.F                         * 
cf     *****************************************************/

