#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

void read_snowband(FILE    *snowband,
		   int      gridcell,
		   double   elev,
		   double **Tfactor,
		   double **Pfactor,
		   double **AreaFract)
/**********************************************************************
  read_snowband		Keith Cherkauer		July 9, 1998

  This routine reads snow elevation band median elevaton, and 
  precipitation fraction for use with the snow model.

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  char    ErrStr[MAXSTRING];
  int     band;
  int     Nbands;
  int     cell;
  double  total;
  double  area_fract;
  double  band_elev;
  double  prec_frac;

  Nbands       = options.SNOW_BAND;
  Tfactor[0]   = (double *)calloc(Nbands,sizeof(double));
  Pfactor[0]   = (double *)calloc(Nbands,sizeof(double));
  AreaFract[0] = (double *)calloc(Nbands,sizeof(double));

  if(Nbands>1) {

    /** Find Current Grid Cell in SnowBand File **/
    rewind(snowband);
    fscanf(snowband, "%i", &cell);
    while(cell != gridcell && !feof(snowband)) {
      fgets(ErrStr,MAXSTRING,snowband);
      fscanf(snowband, "%i", &cell);
    }
    if(feof(snowband)) {
      sprintf(ErrStr,"Cannot find current gridcell (%i) in snow band file",
	      gridcell);
      nrerror(ErrStr);
    }

    /** Read Area Fraction **/
    total = 0.;
    for(band=0;band<Nbands;band++) {
      fscanf(snowband, "%lf", &area_fract);
      if(area_fract<0) {
	sprintf(ErrStr,"Negative snow band area fraction (%lf) read from file", 
		area_fract);
	nrerror(ErrStr);
      }
      AreaFract[0][band]  = area_fract;
      total              += area_fract;
    }
    if(total!=1.) {
      fprintf(stderr,"WARNING: Sum of the snow band area fractions does not equal 1 (%lf), dividing each fraction by the sum\n",
	      total);
      for(band=0;band<options.SNOW_BAND;band++) 
	AreaFract[0][band] /= total;
    }

    /** Read Band Elevation **/
    for(band=0;band<Nbands;band++) {
      fscanf(snowband, "%lf", &band_elev);
      if(band_elev<0) {
	sprintf(ErrStr,"WARNING: Negative snow band elevation (%lf) read from file", 
		band_elev);
      }
      Tfactor[0][band] = (elev - band_elev) / 1000. * T_lapse;
    }
    total = 0.;

    /** Read Precipitation Fraction **/
    for(band=0;band<options.SNOW_BAND;band++) {
      fscanf(snowband, "%lf", &prec_frac);
      if(prec_frac<0) {
	sprintf(ErrStr,"Snow band precipitation fraction (%lf) must be between 0 and 1", 
		prec_frac);
	nrerror(ErrStr);
      }
      if(prec_frac>0 && AreaFract[0][band]==0) {
	sprintf(ErrStr,"Snow band precipitation fraction (%lf) should be 0 when the area fraction is 0. (band = %i)", 
		prec_frac, band);
	nrerror(ErrStr);
      }
      Pfactor[0][band] = prec_frac;
      total += prec_frac;
    }
    if(total!=1.) {
      fprintf(stderr,"WARNING: Sum of the snow band precipitation fractions does not equal %i (%lf), dividing each fraction by the sum\n",
	      1, total);
      for(band=0;band<options.SNOW_BAND;band++) 
	Pfactor[0][band] /= total;
    }
    for(band=0;band<options.SNOW_BAND;band++) {
      if(AreaFract[0][band] > 0)
	Pfactor[0][band] /= AreaFract[0][band];
      else Pfactor[0][band]  = 0.;
    }
  }

  else if(Nbands==1) {
    /** If no snow bands, set factors to use unmodified forcing data **/
    AreaFract[0][0] = 1.;
    Tfactor[0][0]   = 0.;
    Pfactor[0][0]   = 1.;
  }

  else {
    sprintf(ErrStr,"Number of snow elevation bands must be > 0 (%i)",Nbands);
    nrerror(ErrStr);
  }

} 



