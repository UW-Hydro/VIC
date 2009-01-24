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
		   double **AreaFract,
		   char   **AboveTreeLine)
/**********************************************************************
  read_snowband		Keith Cherkauer		July 9, 1998

  This routine reads snow band median elevaton, and 
  precipitation fraction for use with the snow model.

  04-25-03 Modified to allocate treeline variable.            KAC
  2007-Sep-18 Modified to output a warning and set snow bands to 1
	      if no snowband information is present in the snowband
	      file for the current grid cell.				KAC via TJB
  2009-Jan-12 Changed wording of error messages from "snow elevation band"
	      to "snow band".						TJB

**********************************************************************/
{
  extern option_struct options;

  char    ErrStr[MAXSTRING];
  int     band;
  int     Nbands;
  int     cell;
  double  total;
  double  area_fract;
  double  band_elev;
  double  prec_frac;

  Nbands         = options.SNOW_BAND;
  *Tfactor       = (double *)calloc(Nbands,sizeof(double));
  *Pfactor       = (double *)calloc(Nbands,sizeof(double));
  *AreaFract     = (double *)calloc(Nbands,sizeof(double));
  *AboveTreeLine = (char *)calloc(Nbands,sizeof(char));

  if (*Tfactor == NULL || *Pfactor == NULL || *AreaFract == NULL) 
    nrerror("Memory allocation failure in read_snowband");

  if ( Nbands <= 0 ) {
    sprintf(ErrStr,"Number of snow bands must be > 0 (%d)",Nbands);
    nrerror(ErrStr);
  }

  /** Set default values for factors to use unmodified forcing data **/
  for (band = 0; band < Nbands; band++) {
    (*AreaFract)[band] = 0.;
    (*Tfactor)[band]   = 0.;
    (*Pfactor)[band]   = 1.;
  }
  (*AreaFract)[0] = 1.;

  if ( Nbands > 1 ) {

    /** Find Current Grid Cell in SnowBand File **/
#if !NO_REWIND
    rewind(snowband);
#endif

    fscanf(snowband, "%d", &cell);
    while ( cell != gridcell && !feof(snowband) ) {
      fgets(ErrStr,MAXSTRING,snowband);
      fscanf(snowband, "%d", &cell);
    }

    if ( feof(snowband) ) {
      fprintf(stderr, "WARNING: Cannot find current gridcell (%i) in snow band file; setting cell to have one elevation band.\n",
              gridcell);
      /** 1 band is the default; no action necessary **/
      return;
    }

    /** Read Area Fraction **/
    total = 0.;
    for( band = 0; band < Nbands; band++ ) {
      fscanf(snowband, "%lf", &area_fract);
      if(area_fract<0) {
	sprintf(ErrStr,"Negative snow band area fraction (%f) read from file", 
		area_fract);
	nrerror(ErrStr);
      }
      (*AreaFract)[band]  = area_fract;
      total              += area_fract;
    }
    if ( total != 1. ) {
      fprintf(stderr,"WARNING: Sum of the snow band area fractions does not equal 1 (%f), dividing each fraction by the sum\n",
	      total);
      for ( band = 0; band < options.SNOW_BAND; band++ ) 
	(*AreaFract)[band] /= total;
    }

    /** Read Band Elevation **/
    for ( band = 0; band < Nbands; band++ ) {
      fscanf(snowband, "%lf", &band_elev);
      if ( band_elev < 0 ) {
	fprintf(stderr,"Negative snow band elevation (%f) read from file\n", 
		band_elev);
      }
      (*Tfactor)[band] = ( elev - band_elev ) / 1000. * T_lapse;
    }
    total = 0.;

    /** Read Precipitation Fraction **/
    for ( band = 0; band < options.SNOW_BAND; band++ ) {
      fscanf(snowband, "%lf", &prec_frac);
      if(prec_frac<0) {
	sprintf(ErrStr,"Snow band precipitation fraction (%f) must be between 0 and 1", 
		prec_frac);
	nrerror(ErrStr);
      }
      if ( prec_frac > 0 && (*AreaFract)[band] == 0 ) {
	sprintf(ErrStr,"Snow band precipitation fraction (%f) should be 0 when the area fraction is 0. (band = %d)", 
		prec_frac, band);
	nrerror(ErrStr);
      }
      (*Pfactor)[band] = prec_frac;
      total += prec_frac;
    }
    if ( total != 1. ) {
      fprintf(stderr,"WARNING: Sum of the snow band precipitation fractions does not equal %d (%f), dividing each fraction by the sum\n",
	      1, total);
      for(band = 0; band < options.SNOW_BAND; band++) 
	(*Pfactor)[band] /= total;
    }
    for ( band = 0; band < options.SNOW_BAND; band++ ) {
      if ((*AreaFract)[band] > 0)
	(*Pfactor)[band] /= (*AreaFract)[band];
      else 
	(*Pfactor)[band]  = 0.;
    }

  }

} 
