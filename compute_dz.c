#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id";

void compute_dz(double *dz, double *thermdepths, int Nnodes, double dp) {
 
  char ErrStr[MAXSTRING];
  int  j;
  double sum;

  if(thermdepths[Nnodes-1] != dp) {
    sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %f, but is equal to %f",
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
    sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %f, but is equal to %f",
	    Nnodes-1,dp,sum);
    nrerror(ErrStr);
  }

}
