#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/** Number of Times To Run the Optimizer **/
#define N_TRY    15

/** Number of Random Parameter sets Used to Start the Optimizer **/
#define N_INIT  75

/** Number of Parameters to Optimize **/
#define N_PARAM   5

/** Final Optimized Results Tolerance (Stopping Criteria) **/
#define F_TOL     1.e-3

/** Maximum Number of iterations used by Amoeba to find the 
    function minimum **/
#define ITMAX 1000

#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0
#define INVALIDMAX 10
#define MINVAL 1
#define MAXVAL 0

/** Function Prototypes **/
void     nrerror(char*);
void     amoeba(float**,float*,int,float,int*,char*,FILE*);
float    solve_model(float*,char*,FILE*);
float    ran2(long *);
void     quick(float *, float **, int, int);
void     qs(float *, float **, int, int, int);
float  **set_param_limits(int);

main(int argc,char *argv[]) {
/**********************************************************************
  optimize_vic.c         Keith Cherkauer            January 20, 1999

  This program was written to optimize the VIC model using a random,
  autostart, with a simplex minimizer.  

  The program requires:
  (1) the name of a script which will change the parameters of the 
  model, run it, route the results, and compute an R^2 value (times 
  -1 so that an R^2 of -1 is a perfect fit) which is then returned 
  to this program and used to find optimized parameters.

  (2) the name of a log file to which it can write the results of
  all optimization attempts.

**********************************************************************/

  FILE         *fopti;
  float       **p, *y, ftol, **tmpp, *tmpy;
  unsigned int  tmpnum;
  int           Ninit, init, ndim, i, j, iter, try, Ntry, param;
  char          runstr[512];
  char          runoptstr[512];
  time_t        currtime, tmptime;
  float       **param_lim;
  long         *ran2seed;

  if(argc!=4) {
    fprintf(stderr,"Usage: %s <model start script> <model optimize script> <optimization log>\n",
	    argv[0]);
    fprintf(stderr,"\t<model start script> is the model run script to be used to generate the initial random population.\n");
    fprintf(stderr,"\t<model optimize script> is the model run script for the model that is to be optimized (can be the same as the <model start script>, or may run a more computationally intensive mode of the model.\n");
    exit(0);
  }
  
  strcpy(runstr,argv[1]);
  strcpy(runoptstr,argv[2]);
  if((fopti=fopen(argv[3],"w"))==NULL) {
    fprintf(stderr,"ERROR: Unable to open optimization log file %s.\n",
	    argv[2]);
    exit(0);
  }

  Ntry  = N_TRY;
  Ninit = N_INIT;
  ndim  = N_PARAM;
  ftol  = F_TOL;

  p    = (float**)calloc(ndim+1,sizeof(float*));
  y    = (float *)calloc(ndim+1,sizeof(float));
  tmpp = (float**)calloc(Ninit,sizeof(float*));
  tmpy = (float *)calloc(Ninit,sizeof(float));
  for(i=0;i<ndim+1;i++) 
    p[i] = (float*)calloc(ndim,sizeof(float));
  for(i=0;i<Ninit;i++) 
    tmpp[i] = (float*)calloc(ndim,sizeof(float));
  
  i=0;
  tmptime = time(&currtime);
  tmptime *= -1;
  ran2seed = &tmptime;
  ran2(ran2seed);
  ran2seed[0] = 1;

  param_lim = set_param_limits(ndim);

  for(try=0;try<Ntry;try++) {

    fprintf(fopti,"Determining starting parameters for %s trial %i\n",
	    runstr,try);

    for(init=0;init<Ninit;init++) {
    
      fprintf(fopti,"%i\t%i:\t",try,init);

      for(param=0;param<ndim;param++) {
	tmpp[init][param] = ((param_lim[MAXVAL][param] 
				- param_lim[MINVAL][param])
			       * (ran2(ran2seed))) 
	  + param_lim[MINVAL][param];
      }
      
      tmpy[init] = solve_model(tmpp[init],runstr,fopti);
      
    }
    
    quick(tmpy,tmpp,ndim,Ninit);

    for(i=0;i<ndim+1;i++) {
      for(j=0;j<ndim;j++)
	p[i][j] = tmpp[Ninit-i-1][j];
      if(strcmp(runstr,runoptstr)==0)
        y[i]    = tmpy[Ninit-i-1];
      else
        y[i]    = solve_model(p[i],runoptstr,fopti);
    }

    fprintf(fopti,"\nInitial Parameter Set\n");
    for(i=0;i<ndim+1;i++) {
      fprintf(fopti,"%i:\t",i);
      for(j=0;j<ndim;j++) {
	fprintf(fopti,"%f\t",p[i][j]);
      }
      fprintf(fopti,"=\t%f\n",y[i]);
    }
    fflush(fopti);

    amoeba(p,y,ndim,ftol,&iter,runoptstr,fopti);

    fprintf(fopti,"\nResults for Simplex Optimization %i\n",try);
    fprintf(fopti,"\tNumber of iterations needed = %i\n",iter);
    for(i=0;i<ndim+1;i++) {
      fprintf(fopti,"%i:\t",i);
      for(j=0;j<ndim;j++) {
	fprintf(fopti,"%.5g\t",p[i][j]);
      }
      fprintf(fopti,"=\t%.5g\n",y[i]);
    }
    fprintf(fopti,"\n");
  }
    
}

void nrerror(char *errstr) {
/********************************************************************
  Error handling routine
********************************************************************/
  fprintf(stderr,"NRERROR: %s\n",errstr);
  exit(0);
}

void amoeba(float **p,
	    float  *y,
	    int     ndim,
	    float   ftol,
	    int    *iter,
	    char   *runstr,
	    FILE   *fopti)
{
/***********************************************************************
  Simplex optimization routine from Numerical Recipies
***********************************************************************/
  int mpts,j,inhi,ilo,ihi,i;
  float yprr,ypr,rtol;
  float *pr,*prr,*pbar;
  int jk;
  
  pr=(float *)calloc(ndim,sizeof(float));
  prr=(float *)calloc(ndim,sizeof(float));
  pbar=(float *)calloc(ndim,sizeof(float));
  mpts=ndim+1;
  *iter=0;
  for (;;) {
    ilo=0;
    ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
    for (i=0;i<mpts;i++) {
      if (y[i] < y[ilo]) ilo=i;
      if (y[i] > y[ihi]) {
	inhi=ihi;
	ihi=i;
      } else if (y[i] > y[inhi])
	if (i != ihi) inhi=i;
    }
    if(fabs(y[ihi])+fabs(y[ilo]) > 0.)
      rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
    else rtol = 0.;
    if (rtol < ftol) break;
    if ((*iter)++ == ITMAX) nrerror("Too many iterations in AMOEBA");
    for (j=0;j<ndim;j++) pbar[j]=0.0;
    for (i=0;i<mpts;i++)
      if (i != ihi)
	for (j=0;j<ndim;j++) pbar[j] += p[i][j];
    for (j=0;j<ndim;j++) {
      pbar[j] /= ndim;
      pr[j]=(1.0+ALPHA)*pbar[j]-ALPHA*p[ihi][j];
    }
    ypr=solve_model(pr,runstr,fopti);
    if (ypr <= y[ilo]) {
      for (j=0;j<ndim;j++)
	prr[j]=GAMMA*pr[j]+(1.0-GAMMA)*pbar[j];
      yprr=solve_model(prr,runstr,fopti);
      if (yprr < y[ilo]) {
	for (j=0;j<ndim;j++) p[ihi][j]=prr[j];
	y[ihi]=yprr;
      } else {
	for (j=0;j<ndim;j++) p[ihi][j]=pr[j];
	y[ihi]=ypr;
      }
    } else if (ypr >= y[inhi]) {
      if (ypr < y[ihi]) {
	for (j=0;j<ndim;j++) p[ihi][j]=pr[j];
	y[ihi]=ypr;
      }
      for (j=0;j<ndim;j++)
	prr[j]=BETA*p[ihi][j]+(1.0-BETA)*pbar[j];
      yprr=solve_model(prr,runstr,fopti);
      if (yprr < y[ihi]) {
	for (j=0;j<ndim;j++) p[ihi][j]=prr[j];
	y[ihi]=yprr;
      } else {
	for (i=0;i<mpts;i++) {
	  if (i != ilo) {
	    for (j=0;j<ndim;j++) {
	      pr[j]=0.5*(p[i][j]+p[ilo][j]);
	      p[i][j]=pr[j];
	    }
	    y[i]=solve_model(pr,runstr,fopti);
	  }
	}
      }
    } else {
      for (j=0;j<ndim;j++) p[ihi][j]=pr[j];
      y[ihi]=ypr;
    }
  }
  free((char*)pbar);
  free((char*)prr);
  free((char*)pr);
}

float solve_model(float *p, char *runstr, FILE *fopti) {
/***********************************************************************
  solve_model

  This subroutine checks the parameters in *p to see if they are valid.
  If so, then it starts the script that runs the model with the new
  parameter set.  Before returning to the amoeba it opens and reads 
  the R squared value for the run just completed, so that it can be 
  returned.  Invalid parameters return a default value of INVALIDMAX.

  parameters:
       p[0] -> b_infilt
       p[1] -> Ds
       p[2] -> Ws
       p[3] -> D2
       p[4] -> D3

***********************************************************************/

  FILE *fin;
  char *cmdstr;
  char  filename[512];
  float Rsqr;
  int   i,j;

/*   Rsqr = 0.; */
/*   for(i=0;i<5;i++) { */
/*     Rsqr += (p[i]*p[i]); */
/*   } */

  if((p[0]>0) && (p[1]>0 && p[1]<1) && (p[2]>0 && p[2]<1) && (p[3]>0.10)
     && (p[4]>0.10)) {

    cmdstr = (char *)calloc(512,sizeof(char));
    
    sprintf(filename,"R2.simp.%s.txt",runstr);
    sprintf(cmdstr,"%s %f %f %f %f %f %s", 
	    runstr, p[0], p[1], p[2], p[3], p[4], filename);
    
    system(cmdstr);

    if((fin=fopen(filename,"r"))==NULL) {
      nrerror("Unable to open Rsquared file");
    }
    
    fscanf(fin,"%f",&Rsqr);
    fclose(fin);
    fprintf(fopti,"Rsquared (%f,%f,%f,%f,%f) = %f\n",
	    p[0],p[1],p[2],p[3],p[4],Rsqr);
    
    free((char *)cmdstr);

  }
  else {
    
    Rsqr = (float)INVALIDMAX;
    
    fprintf(fopti,"Invalid Parameter Set (%f,%f,%f,%f,%f) -> returning %f\n",
	    p[0],p[1],p[2],p[3],p[4],Rsqr);

  }
  fflush(fopti);

  return(Rsqr);

}

float **set_param_limits(int Nparam) {
/**************************************************************
  Sets limits for all parameters
**************************************************************/

  float **param_lim;
  int     i;

  param_lim = (float **)calloc(2,sizeof(float *));
  for(i=0;i<2;i++)
    param_lim[i] = (float *)calloc(Nparam,sizeof(float));

/*   for(i=0;i<Nparam;i++) { */
/*     param_lim[MAXVAL][i] = 10.; */
/*     param_lim[MINVAL][i] = -2.; */
/*   } */

  /* b_infilt */
  param_lim[MAXVAL][0] = 0.5000;
  param_lim[MINVAL][0] = 1.e-5;

  /* Ds */
  param_lim[MAXVAL][1] = 0.50000;
  param_lim[MINVAL][1] = 1.e-5;

  /* Ws */
  param_lim[MAXVAL][2] = 1.0000;
  param_lim[MINVAL][2] = 0.2000;

  /* D2 */
  param_lim[MAXVAL][3] = 2.50;
  param_lim[MINVAL][3] = 0.10;

  /* D3 */
  param_lim[MAXVAL][4] = 2.50;
  param_lim[MINVAL][4] = 0.10;

  return(param_lim);

}

void quick(float *item, float **param, int Nparam, int count)
/**********************************************************************
        this subroutine starts the quick sort
**********************************************************************/
{
  qs(item,param,Nparam,0,count-1);
}
 
void qs(float *item, float **param, int Nparam, int left, int right)
/**********************************************************************
        this is the quick sort subroutine - it returns the values in
	an array from high to low.
**********************************************************************/
{
  register int i,j;
  float x,y;
  float *xp, *yp;
  int p;

  xp = (float*)calloc(Nparam,sizeof(float));
  yp = (float*)calloc(Nparam,sizeof(float));
 
  i=left;
  j=right;
  x=item[(left+right)/2];
  for(p=0;p<Nparam;p++) xp[p] = param[(left+right)/2][p];
 
  do {
    while(item[i]>x && i<right) i++;
    while(x>item[j] && j>left) j--;
 
    if (i<=j) {
      y=item[i];
      for(p=0;p<Nparam;p++) yp[p] = param[i][p];
      item[i]=item[j];
      for(p=0;p<Nparam;p++) param[i][p] = param[j][p];
      item[j]=y;
      for(p=0;p<Nparam;p++) param[j][p] = yp[p];
      i++;
      j--;
    }
  } while (i<=j);
 
  if(left<j) qs(item,param,Nparam,left,j);
  if(i<right) qs(item,param,Nparam,i,right);

  free((char*)xp);
  free((char*)yp);

}

#define M 714025
#define IA 1366
#define IC 150889
 
float ran2(long *idum)
/******************************************************************
  Random number generator from Numerical Recipes
******************************************************************/
{
        static long iy,ir[98];
        static int iff=0;
        int j;
 
        if (*idum < 0 || iff == 0) {
                iff=1;
                if ((*idum=(IC-(*idum)) % M) < 0) *idum = -*idum;
                for (j=1;j<=97;j++) {
                        *idum=(IA*(*idum)+IC) % M;
                        ir[j]=(*idum);
                }
                *idum=(IA*(*idum)+IC) % M;
                iy=(*idum);
        }
        j=1 + 97.0*iy/M;
        if (j > 97 || j < 1) nrerror("RAN2: This cannot happen.");
        iy=ir[j];
        *idum=(IA*(*idum)+IC) % M;
        ir[j]=(*idum);
        return (float) iy/M;
}
 
#undef M
#undef IA
#undef IC
#undef ALPHA
#undef BETA
#undef GAMMA
#undef ITMAX
#undef INVALIDMAX

