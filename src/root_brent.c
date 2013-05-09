/*
 * SUMMARY:      RootBrent.c - Determine surface temperature iteratively
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * LAST-MOD: Mon Jan 24 12:06:38 2000 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Determine surface temperature iteratively using the Brent
 *               method.  
 * DESCRIP-END.
 * FUNCTIONS:    RootBrent()
 * COMMENTS:     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/*****
#include <stdarg.h>
#include "settings.h"
#include "brent.h"
#include "massenergy.h"
#include "DHSVMerror.h"
*****/

#define MAXTRIES 5
#define MAXITER 1000
#define MACHEPS 3e-8
#define TSTEP   10
#define T       1e-7   

/*****************************************************************************
  GENERAL DOCUMENTATION FOR THIS MODULE
  -------------------------------------

  Source: Brent, R. P., 1973, Algorithms for minimization without derivatives,
                        Prentice Hall, Inc., Englewood Cliffs, New Jersey
			Chapter 4
  This source includes an implementation of the algorithm in ALGOL-60, which
  was translated into C for this application.

  The method is also discussed in:
  Press, W. H., S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, 1992,
                Numerical Recipes in FORTRAN, The art of scientific computing,
		Second edition, Cambridge University Press
  (Be aware that this book discusses a Brent method for minimization (brent), 
  and one for root finding (zbrent).  The latter one is similar to the one 
  implemented here and is also copied from Brent [1973].)

  The function returns the surface temperature, TSurf, for which the sum
  of the energy balance terms is zero, with TSurf in the interval 
  [MinTSurf, MaxTSurf].  The surface temperature is calculated to within
  a tolerance (6 * MACHEPS * |TSurf| + 2 * T), where MACHEPS is the relative
  machine precision and T is a positive tolerance, as specified in brent.h.

  The function assures that f(MinTSurf) and f(MaxTSurf) have opposite signs.
  If this is not the case the program will abort.  In addition the program
  will perform not more than a certain number of iterations, as specified
  in brent.h, and will abort if more iterations are needed.
******************************************************************************/
  
/*****************************************************************************
  Function name: RootBrent()

  Purpose      : Calculate the surface temperature

  Required     :
    double LowerBound     - Lower bound for root
    double UpperBound     - Upper bound for root
    char *ErrorString     - For storing description of errors (if any)
    double (*Function)(double Estimate, va_list ap)
    ...                   - Variable arguments 
                            The number and order of arguments has to be
                            appropriate for the Function pointed to, since
                            the list of arguments after Nargs will be passed
                            on to Function.
                            See the appropriate Function for the correct
                            arguments. 

  Returns      :
    double b              - Effective surface temperature (C)

  Modifies     : 
    char *ErrorString     - Stores description of errors

  Comments     :
    04-Jun-04 Removed message announcing the dumping of variables,
	      since that doesn't always happen when root_brent fails.
	      Changed remaining message from ERROR to WARNING for
	      the same reason.						TJB
    21-Sep-04 No longer print warning to stderr from this routine;
	      instead store warning messages in parameter ErrorString.	TJB
  2007-Aug-31 Corrected handling of Function return value if Function returns
	      ERROR.  This can happen if Function is func_surf_energy_bal.	JCA
  2007-Sep-01 Removed the integer "eval" since it is never used for anything.	JCA
  2009-May-22 Modified root-bracketing scheme to handle case when one bound
	      yields garbage output from the target function.			TJB
*****************************************************************************/
double root_brent(double LowerBound, double UpperBound, char *ErrorString,
                double (*Function)(double Estimate, va_list ap), ...)
{
  const char *Routine = "RootBrent";
  va_list ap;                   /* Used in traversing variable argument list */ 
  double a;
  double b;
  double c;
  double d;
  double e;
  double fa;
  double fb;
  double fc;
  double m;
  double p;
  double q;
  double r;
  double s;
  double tol;
  double last_bad;
  double last_good;
  int which_err;
  int i;
  int j;

  /* initialize variable argument list */
  a = LowerBound;
  b = UpperBound;
  va_start(ap, Function);
  fa = Function(a, ap);
  va_start(ap, Function);
  fb = Function(b, ap);
 
  which_err = 0;

  // If Function returns values of ERROR for both bounds, give up
  if (fa == ERROR && fb == ERROR) {
    sprintf(ErrorString,"ERROR: %s: lower and upper bounds %f and %f failed to bracket the root because the given function was not defined at either point.\n",Routine,a,b);
    va_end(ap);
    return(ERROR);
  }      

  // If Function returns value of ERROR for one bound but not both bounds,
  // move the offending bound until the Function returns a valid value
  if(fa == ERROR || fb == ERROR) {

    if (fa == ERROR) {
      which_err = -1;
      last_bad = a;
      last_good = b;
    }
    else {
      which_err = 1;
      last_good = a;
      last_bad = b;
    }

    c = 0.5*(last_bad+last_good);
    va_start(ap, Function);
    fc = Function(c, ap);

    /* search for valid point via bisection */
    j = 0;
    while (fc == ERROR && j < MAXITER) {
      last_bad = c;
      c = 0.5*(last_bad+last_good);
      va_start(ap, Function);
      fc = Function(c, ap);
      j++;
    }

    if (fc == ERROR) {
      /* if we get here, we could not find a bound for which the function returns a valid value */
      sprintf(ErrorString,"ERROR: %s: the given function produced undefined values while attempting to bracket the root between %f and %f.\n",Routine,LowerBound,UpperBound);
      va_end(ap);
      return(ERROR);
    }
    else {
      if (which_err == -1) {
        a = c;
        fa = fc;
      }
      else {
        b = c;
        fb = fc;
      }
    }

  }

  // At this point, we have two bounds that yield valid values of the target function

  /*  if root not bracketed attempt to bracket the root */
  j = 0;
  while ((fa * fb) >= 0  && j < MAXTRIES) {
    /* Expansion of bounds depends on whether initial bounds encountered undefined function values */
    if (which_err == 0) { // No undefined values were encountered
      a -= TSTEP;
      b += TSTEP;
      va_start(ap, Function);
      fa = Function(a, ap);
      va_start(ap, Function);
      fb = Function(b, ap);
    }
    else { // Undefined values were encountered
      if (which_err == -1) { // Undefined values encountered in the lower direction
        b += TSTEP;
        va_start(ap, Function);
        fb = Function(b, ap);
        if (fb == ERROR) {
          /* Undefined function values in both directions - give up */
          sprintf(ErrorString,"ERROR: %s: the given function produced undefined values while attempting to bracket the root between %f and %f.\n",Routine,LowerBound,UpperBound);
          va_end(ap);
          return(ERROR);
        }
        last_good = a;
      }
      else { // Undefined values encountered in the upper direction
        a -= TSTEP;
        va_start(ap, Function);
        fa = Function(a, ap);
        if (fa == ERROR) {
          /* Undefined function values in both directions - give up */
          sprintf(ErrorString,"ERROR: %s: the given function produced undefined values while attempting to bracket the root between %f and %f.\n",Routine,LowerBound,UpperBound);
          va_end(ap);
          return(ERROR);
        }
        last_good = b;
      }

      /* search for valid point via bisection */
      c = 0.5*(last_good+last_bad);
      va_start(ap, Function);
      fc = Function(c, ap);
      i = 0;
      while (fc == ERROR && i < MAXITER) {
        last_bad = c;
        c = 0.5*(last_bad+last_good);
        va_start(ap, Function);
        fc = Function(c, ap);
        i++;
      }

      if (fc == ERROR) {
        /* if we get here, we could not find a bound for which the function returns a valid value */
        sprintf(ErrorString,"ERROR: %s: the given function produced undefined values while attempting to bracket the root between %f and %f.\n",Routine,LowerBound,UpperBound);
        va_end(ap);
        return(ERROR);
      }
      else {
        if (which_err == -1) {
          a = c;
          fa = fc;
        }
        else {
          b = c;
          fb = fc;
        }
      }

    }

    j++;
  }
  if ((fa * fb) >= 0) {
    /* if we get here, the lower and upper bounds did not bracket the root */
    sprintf(ErrorString,"WARNING: %s: lower and upper bounds %f and %f failed to bracket the root.\n",Routine,a,b);
    va_end(ap);
    return(ERROR);
  }

  // At this point, we have bracketed the root
 
  // Now search for the root

  fc = fb;

  for (i = 0; i < MAXITER; i++) {

    if (fb*fc > 0) {
      c = a;
      fc = fa;
      d = b - a;
      e = d;
    }
    
    if (fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    
    tol = 2 * MACHEPS * fabs(b) + T;
    m = 0.5 * (c - b);
    
    if (fabs(m) <= tol || fb == 0) {
      va_end(ap);
      return b;
    }
    
    else {
      if (fabs(e) < tol || fabs(fa) <= fabs(fb)) {
	d = m;
	e = d;
      }
      else {
	s = fb/fa;
	
	if (a == c) {
	  
	  /* linear interpolation */
          
	  p = 2 * m * s;
	  q = 1 - s;
	}
	
	else {
	  
	  /* inverse quadratic interpolation */
	  
	  q = fa/fc;
	  r = fb/fc;
	  p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
	  q = (q - 1) * (r - 1) * (s - 1);
	}
	
	if (p > 0)
	  q = -q;
	else
	  p = -p;
	s = e;
	e = d;
	if ((2 * p) < ( 3 * m * q - fabs(tol * q)) && p < fabs(0.5 * s * q))
	  d = p/q;
	else {
	  d = m;
	  e = d;
	}
      }
      a = b;
      fa = fb;
      b += (fabs(d) > tol) ? d : ((m > 0) ? tol : -tol);
      va_start(ap, Function);
      fb = Function(b, ap);

      // Catch ERROR values returned from Function
      if(fb == ERROR){
	sprintf(ErrorString,"ERROR returned to root_brent on iteration %d: temperature = %.4f\n",i+1,b);
	va_end(ap);
	return( ERROR );
      }      

    }
  }
  /* If we get here, there were too many iterations */
  sprintf(ErrorString,"WARNING: %s: too many iterations.\n",Routine);
  va_end(ap);
  return(ERROR);

}

#undef MAXTRIES
#undef MAXITER
#undef MACHEPS
#undef TSTEP
#undef T
