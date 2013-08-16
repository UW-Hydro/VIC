#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id: compute_dz.c,v 4.2.2.2 2004/05/10 21:46:45 tbohn Exp $";

void compute_dz(double *dz, double *thermdepths, int Nnodes, double dp) {
/***********************************************************************
  compute_dz              Keith Cherkauer	 	May 16, 2000

  This routines computes the soil thermal node thicknesses using an
  array of node depths.

  07-May-04 Changed
		if(thermdepths[Nnodes-1] != dp) {
	    to
		if((int)(thermdepths[Nnodes-1]*1000 + 0.5) != (int)(dp*1000 + 0.5)) {
	    and changed
		rint(something)
	    to
		(float)(int)(something + 0.5)
	    to handle floating-point representation errors in the data
	    read from the state file without resorting to rint function.
								TJB
***********************************************************************/
 
  char ErrStr[MAXSTRING];
  int  j;
  double sum;

  if((int)(thermdepths[Nnodes-1]*1000 + 0.5) != (int)(dp*1000 + 0.5)) {
    sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %f, but is equal to %f",
	    Nnodes-1,dp,thermdepths[Nnodes-1]);
    nrerror(ErrStr);
  }

  for(j=Nnodes-1;j>0;j--) {
    thermdepths[j] -= thermdepths[j-1];
    thermdepths[j] = (float)(int)(thermdepths[j] * 10000. + 0.5) / 10000.;
  }

  sum = 0;
  dz[0] = (float)(int)(thermdepths[1] * 10000. + 0.5) / 10000.;
  for(j=1;j<Nnodes;j++) {
    dz[j] = 2. * (float)(int)((thermdepths[j] - dz[j-1] / 2.) * 10000. + 0.5) / 10000.;
    if(dz[j] < 0) {
      sprintf(ErrStr,"Check spacing between thermal layers %i and %i\n",
	      j-1,j);
      nrerror(ErrStr);
    }
    sum += (dz[j-1] + dz[j]) / 2.;
  }

  if((int)(sum*1000 + 0.5) != (int)(dp*1000 + 0.5)) {
    sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %f, but is equal to %f",
	    Nnodes-1,dp,sum);
    nrerror(ErrStr);
  }

}
