#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

void read_initial_soil_thermal(FILE   *initsoil,
			       int     cellnum,
			       int     veg,
			       int     band,
			       int     Nnodes,
			       double  dp,
			       double *fdepth,
			       double *moist,
			       double *ice,
			       double *dz,
			       double *T)
/******************************************************************************
  read_initial_snow_thermal 	Keith Cherkauer	      January 12, 1999

  This routine reads parameters for initial soil thermal conditions - if defined

  Parameters Read from File:
  TYPE    NAME                    UNITS   DESCRIPTION
  int     cellnum                 N/A     number of current cell
  int     veg                     N/A     vegetation number within cell
  int     band                    N/A     snow band number within cell
  double *fdepth                  m       depth of freezing and thawing fronts
  double *moist                   mm/mm   moisture content of each soil layer
  double *ice                     mm/mm   ice content of each soil layer
  double *thermdepths             m       depth from surface of each thermal soln node
  double *T                       C       temperature at each thermal soln node depth

  Parameters Computed from those in the File:
  TYPE    NAME                    UNITS   DESCRIPTION
  double *dz                      m       thickness of soil layer between thermal nodes

  Modifications:

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  char             ErrStr[MAXSTRING];
  char             tmpstr[MAXSTRING];
  int              i, tmpint, tmpveg, tmpband, index, tmpnodes;
  double           *thermdepths;

  rewind(initsoil);

  fscanf(initsoil, "%s", tmpstr);
  if(tmpstr[0] != '#') {
    index = atoi(tmpstr);
    fscanf(initsoil, "%i %i", &tmpveg, &tmpband);
  }
  else index = tmpveg = tmpband = -999;
  while(index != cellnum && veg != tmpveg && band != tmpband) {
    fgets(tmpstr,MAXSTRING,initsoil);
    fscanf(initsoil, "%s", tmpstr);
    if(tmpstr[0] != '#') {
      index = atoi(tmpstr);
      fscanf(initsoil, "%i %i", &tmpveg, &tmpband);
    }
    else index = tmpveg = tmpband = -999;
  }

  if(!feof(initsoil) && options.FROZEN_SOIL) {

    /** Read Soil Thermal Initialization File for Activated Frozen Soil Model **/

    thermdepths = (double*)calloc(Nnodes,sizeof(double));
  
    fscanf(initsoil, "%i",  &tmpnodes);
    if(Nnodes!=tmpnodes) {
      sprintf(ErrStr,"Number of nodes defined in soil thermal initialization file (%i), not equal to number of nodes defined in model (%i).",tmpnodes,Nnodes);
    }
    for(i=0;i<2;i++)       fscanf(initsoil, "%lf",     &fdepth[i]);
    for(i=0;i<options.Nlayer;i++) fscanf(initsoil, "%lf",     &moist[i]);
    for(i=0;i<options.Nlayer;i++) fscanf(initsoil, "%lf",     &ice[i]);
    for(i=0;i<Nnodes;i++)  fscanf(initsoil, "%lf %lf", &thermdepths[i], &T[i]);

    compute_dz(dz, thermdepths, Nnodes, dp);
    
    free((char*)thermdepths);

  }

  else if(!feof(initsoil)) {

    /** Read Initial Soil Thermal Properties for Full Energy Balance Model **/

    thermdepths = (double*)calloc(Nnodes,sizeof(double));
  
    fscanf(initsoil, "%i",  &tmpnodes);
    if(Nnodes!=tmpnodes) {
      sprintf(ErrStr,"Number of nodes defined in soil thermal initialization file (%i), not equal to number of nodes defined in model (%i).",tmpnodes,Nnodes);
    }
    for(i=0;i<2;i++) {
      fscanf(initsoil, "%lf", &fdepth[i]);
      fdepth[i] = 0.;
    }
    for(i=0;i<options.Nlayer;i++) fscanf(initsoil, "%lf", &moist[i]);

    for(i=0;i<options.Nlayer;i++) {
      fscanf(initsoil, "%lf", &ice[i]);
      moist[i] += ice[i];
      ice[i] = 0.;
    }
    for(i=0;i<Nnodes;i++)  fscanf(initsoil, "%lf %lf", &thermdepths[i], &T[i]);

    compute_dz(dz, thermdepths, Nnodes, dp);

    free((char*)thermdepths);

  }
  else {

    /** No Initial Conditions Provided in File **/
    sprintf(ErrStr,"No initial soil energy balance conditions provided for cell %i, vegetation %i, and snow band %i.  Currently unable to handle this options, check input file.",cellnum,veg,band);
    nrerror(ErrStr);
  }

}

void compute_dz(double *dz, double *thermdepths, int Nnodes, double dp) {
 
  char ErrStr[MAXSTRING];
  int  i, j;
  double sum;

  if(thermdepths[Nnodes-1] != dp) {
    sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %lf, but is equal to %lf",
	    Nnodes-1,dp,thermdepths[Nnodes-1]);
    nrerror(ErrStr);
  }

  for(j=Nnodes-1;j>0;j--) {
    thermdepths[j] -= thermdepths[j-1];
    thermdepths[j] = rint(thermdepths[j] * 10000.) / 10000.;
  }

  sum = 0;
  dz[0] = rint(thermdepths[1] * 10000.) / 10000.;
  for(j=1;j<Nnodes;j++) {
    dz[j] = 2. * rint((thermdepths[j] - dz[j-1] / 2.) * 10000.) / 10000.;
    if(dz[j] < 0) {
      sprintf(ErrStr,"Check spacing between thermal layers %i and %i\n",
	      j-1,j);
      nrerror(ErrStr);
    }
    sum += (dz[j-1] + dz[j]) / 2.;
  }

  if(rint(sum*1000) != rint(dp*1000)) {
    sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %lf, but is equal to %lf",
	    Nnodes-1,dp,sum);
    nrerror(ErrStr);
  }

}
