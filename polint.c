#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

void polint(double *xa,double *ya,int n,double x,double *y,double *dy)
{
	int i,m,ns=0;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[0]);
	c=(double *)malloc(n*sizeof(double));
	d=(double *)malloc(n*sizeof(double));
	for (i=0;i<n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=0;i<=n-m-1;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine POLINT");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free((char *)d);
	free((char *)c);
}
