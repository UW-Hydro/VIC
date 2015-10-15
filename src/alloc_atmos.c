/*
 * Purpose: allocate and free memory for the atmos data struct
 * Usage  : Part of VIC
 * Author : Bart Nijssen
 * E-mail : nijssen@u.washington.edu
 * Created: Fri Aug 27 18:22:42 1999
 * Last Changed: Tue Sep  2 15:11:02 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
 * Notes  :
 */

/****************************************************************************/
/*			  PREPROCESSOR DIRECTIVES                           */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/****************************************************************************/
/*			       alloc_atmos()                                */
/****************************************************************************/
void alloc_atmos(int nrecs, atmos_data_struct **atmos)
/*******************************************************************
  alloc_atmos    

  Modifications:
  01-11-00 Fixed allocation bug                             KAC
  2006-Sep-23 Implemented flexible output configuration; removed
	      LDAS_OUTPUT and OPTIMIZE compile-time options.		TJB
  2006-Dec-20 All atmos_data arrays are always dynamically allocated
	      now.							TJB
  2010-Mar-31 Added runoff_in.						TJB
  2010-Sep-24 Renamed runoff_in to channel_in.				TJB
  2011-Nov-04 Added tskc.						TJB
  2013-Jul-25 Added Catm, coszen, fdir, and par.			TJB

*******************************************************************/
{
  extern param_set_struct param_set;

  int i;

  *atmos = (atmos_data_struct *) calloc(nrecs, sizeof(atmos_data_struct)); 
  if (*atmos == NULL)
    vicerror("Memory allocation error in alloc_atmos().");

  for (i = 0; i < nrecs; i++) {
    (*atmos)[i].air_temp = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].air_temp == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].Catm = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].Catm == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].channel_in = (double *) calloc(NR+1, sizeof(double));	
    if ((*atmos)[i].channel_in == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].coszen = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].coszen == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].density = (double *) calloc(NR+1, sizeof(double));	
    if ((*atmos)[i].density == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].fdir = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].fdir == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].longwave = (double *) calloc(NR+1, sizeof(double));	
    if ((*atmos)[i].longwave == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].par = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].par == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].prec = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].prec == NULL)
      vicerror("Memory allocation error in alloc_atmos().");      
    (*atmos)[i].pressure = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].pressure == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].shortwave = (double *) calloc(NR+1, sizeof(double));	
    if ((*atmos)[i].shortwave == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].snowflag = (char *) calloc(NR+1, sizeof(char));	
    if ((*atmos)[i].snowflag == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].tskc = (double *) calloc(NR+1, sizeof(double));	
    if ((*atmos)[i].tskc == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].vp = (double *) calloc(NR+1, sizeof(double));	
    if ((*atmos)[i].vp == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].vpd = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].vpd == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].wind = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].wind == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].irr_run = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].irr_run == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
    (*atmos)[i].irr_with = (double *) calloc(NR+1, sizeof(double));
    if ((*atmos)[i].irr_with == NULL)
      vicerror("Memory allocation error in alloc_atmos().");
  }    			

}

/****************************************************************************/
/*	      		  free_atmos()                                      */
/****************************************************************************/
void free_atmos(int nrecs, atmos_data_struct **atmos)
/***************************************************************************
  Modifications:
  09-02-2003 Added check for LINK_DEBUG global option.  If LINK_DEBUG is
             TRUE atmospheric data is not dynamically allocated, so it
             should not be freed.                                   KAC
  2006-Sep-23 (Port from 4.0.6) Implemented flexible output configuration;
	      removed LDAS_OUTPUT and OPTIMIZE compile-time options.	TJB
  2006-Dec-20 All atmos_data arrays are always dynamically allocated
	      now.							TJB
  2010-Mar-31 Added runoff_in.						TJB
  2010-Sep-24 Renamed runoff_in to channel_in.				TJB
  2011-Nov-04 Added tskc.						TJB
  2013-Jul-25 Added Catm, coszen, fdir, and par.			TJB
***************************************************************************/
{
  int i;

  if (*atmos == NULL)
    return;

  for (i = 0; i < nrecs; i++) {
    free((*atmos)[i].air_temp);
    free((*atmos)[i].Catm);
    free((*atmos)[i].channel_in);
    free((*atmos)[i].coszen);
    free((*atmos)[i].density);
    free((*atmos)[i].fdir);
    free((*atmos)[i].longwave);
    free((*atmos)[i].par);
    free((*atmos)[i].prec);
    free((*atmos)[i].pressure);
    free((*atmos)[i].shortwave);
    free((*atmos)[i].snowflag);
    free((*atmos)[i].tskc);
    free((*atmos)[i].vp);
    free((*atmos)[i].vpd);
    free((*atmos)[i].wind);
    free((*atmos)[i].irr_run);
    free((*atmos)[i].irr_with);
  }

  free(*atmos);
}
