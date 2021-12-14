/******************************************************************************
* \file
* \author  Keith Cherkauer <cherkaue@purdue.edu>
*
* \section DESCRIPTION
*
* This routine writes vegetation variables to stdout.  Used primarily
* for debugging purposes.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* \brief        This routine writes vegetation variables to stdout.
******************************************************************************/
void
write_vegvar(veg_var_struct *veg,
             int             n)
{
    printf("Vegetation Variables: vegtype %i\n", n);
    printf("\tcanopyevap  = %f\n", veg->canopyevap);
    printf("\tWdew        = %f\n", veg->Wdew);
    printf("\tthroughfall = %f\n", veg->throughfall);
}
