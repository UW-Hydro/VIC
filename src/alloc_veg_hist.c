/*
 * Purpose: allocate and free memory for the veg_hist data struct
 * Usage  : Part of VIC
 * Author : Ted Bohn
 * E-mail : theodore.bohn@asu.edu
 * Created: 2013-Sep-03
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
/*			       alloc_veg_hist()                             */
/****************************************************************************/
void alloc_veg_hist(int nrecs, int nveg, veg_hist_struct ***veg_hist)
/*******************************************************************
  alloc_veg_hist    

  Modifications:
  2014-Apr-25 Added veg cover fraction.					TJB
*******************************************************************/
{
  int i,j;

  (*veg_hist) = (veg_hist_struct **) calloc(nrecs, sizeof(veg_hist_struct *)); 
  if ((*veg_hist) == NULL)
    vicerror("Memory allocation error in alloc_veg_hist().");

  for (i = 0; i < nrecs; i++) {
    (*veg_hist)[i] = (veg_hist_struct *) calloc(nveg, sizeof(veg_hist_struct));
    if ((*veg_hist)[i] == NULL)
      vicerror("Memory allocation error in alloc_veg_hist().");
    for (j = 0; j < nveg; j++) {
      (*veg_hist)[i][j].albedo = (double *) calloc(NR+1, sizeof(double));
      if ((*veg_hist)[i][j].albedo == NULL)
        vicerror("Memory allocation error in alloc_veg_hist().");
      (*veg_hist)[i][j].crop_frac = (double *) calloc(NR+1, sizeof(double));
      if ((*veg_hist)[i][j].crop_frac == NULL)
        vicerror("Memory allocation error in alloc_veg_hist().");
      (*veg_hist)[i][j].LAI = (double *) calloc(NR+1, sizeof(double));
      if ((*veg_hist)[i][j].LAI == NULL)
        vicerror("Memory allocation error in alloc_veg_hist().");
      (*veg_hist)[i][j].vegcover = (double *) calloc(NR+1, sizeof(double));
      if ((*veg_hist)[i][j].vegcover == NULL)
        vicerror("Memory allocation error in alloc_veg_hist().");
    }
  }

}

/****************************************************************************/
/*	      		  free_veg_hist()                                   */
/****************************************************************************/
void free_veg_hist(int nrecs, int nveg, veg_hist_struct ***veg_hist)
/***************************************************************************
  Modifications:
  2014-Apr-25 Added veg cover fraction.					TJB
***************************************************************************/
{
  int i,j;

  if (*veg_hist == NULL)
    return;

  for (i = 0; i < nrecs; i++) {
    for (j = 0; j < nveg; j++) {
      free((*veg_hist)[i][j].albedo);
      free((*veg_hist)[i][j].crop_frac);
      free((*veg_hist)[i][j].LAI);
      free((*veg_hist)[i][j].vegcover);
    }
    free((*veg_hist)[i]);
  }

  free(*veg_hist);
}
