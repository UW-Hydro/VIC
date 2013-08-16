/*
 * Purpose: allocate and free memory for the atmos data struct
 * Usage  : Part of VIC
 * Author : Bart Nijssen
 * E-mail : nijssen@u.washington.edu
 * Created: Fri Aug 27 18:22:42 1999
 * Last Changed: Tue Feb  8 14:46:22 2000 by Keith Cherkauer <cherkaue@u.washington.edu>
 * Notes  :
 */

/******************************************************************************/
/*			    PREPROCESSOR DIRECTIVES                           */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <vicNl.h>

static char vcid[] = "$Id: alloc_atmos.c,v 4.1.2.2 2007/01/11 02:01:04 vicadmin Exp $";

/******************************************************************************/
/*				 alloc_atmos()                                */
/******************************************************************************/
void alloc_atmos(int nrecs, atmos_data_struct **atmos)
/*******************************************************************
  alloc_atmos    

  Modifications:
  01-11-00 Fixed allocation bug                             KAC
  2006-Dec-20 All atmos_data arrays are always dynamically allocated now.	TJB

*******************************************************************/
{
  extern param_set_struct param_set;

  int i;

  *atmos = (atmos_data_struct *) calloc(nrecs, sizeof(atmos_data_struct)); 
  if (*atmos == NULL)
    vicerror("Memory allocation error in alloc_atmos().");

   for (i = 0; i < nrecs; i++) {
     (*atmos)[i].prec = (double *) calloc(NR+1, sizeof(double));
     if ((*atmos)[i].prec == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
     (*atmos)[i].orig_prec = (double *) calloc(NR+1, sizeof(double)); //ingjerd jun 2009
     if ((*atmos)[i].orig_prec == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
     (*atmos)[i].air_temp = (double *) calloc(NR+1, sizeof(double));
     if ((*atmos)[i].air_temp == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
     (*atmos)[i].wind = (double *) calloc(NR+1, sizeof(double));
     if ((*atmos)[i].wind == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
     (*atmos)[i].vpd = (double *) calloc(NR+1, sizeof(double));
     if ((*atmos)[i].vpd == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
     (*atmos)[i].vp = (double *) calloc(NR+1, sizeof(double));
     if ((*atmos)[i].vp == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
     (*atmos)[i].pressure = (double *) calloc(NR+1, sizeof(double));
     if ((*atmos)[i].pressure == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
     (*atmos)[i].density = (double *) calloc(NR+1, sizeof(double));
     if ((*atmos)[i].density == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
     (*atmos)[i].shortwave = (double *) calloc(NR+1, sizeof(double));
     if ((*atmos)[i].shortwave == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
     (*atmos)[i].longwave = (double *) calloc(NR+1, sizeof(double));
     if ((*atmos)[i].longwave == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
     (*atmos)[i].snowflag = (char *) calloc(NR+1, sizeof(char));
     if ((*atmos)[i].snowflag == NULL)
       vicerror("Memory allocation error in alloc_atmos().");
   }
}

/******************************************************************************/
/*				  free_atmos()                                */
/******************************************************************************/
void free_atmos(int nrecs, atmos_data_struct **atmos)
/***************************************************************************
  Modifications:
  2006-Dec-20 All atmos_data arrays are always dynamically allocated now.	TJB
***************************************************************************/
{
  int i;

  if (*atmos == NULL)
    return;

   for (i = 0; i < nrecs; i++) {
     free((*atmos)[i].prec);
     free((*atmos)[i].air_temp);
     free((*atmos)[i].wind);
     free((*atmos)[i].vpd);
     free((*atmos)[i].vp);
     free((*atmos)[i].pressure);
     free((*atmos)[i].density);
     free((*atmos)[i].shortwave);
     free((*atmos)[i].longwave);
     free((*atmos)[i].snowflag);
   }
  free(*atmos);
}
