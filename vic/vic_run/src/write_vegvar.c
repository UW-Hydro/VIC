/******************************************************************************
* \file
* \author  Keith Cherkauer <cherkaue@purdue.edu>
*
* \section DESCRIPTION
*
* This routine writes vegetation variables to stdout.  Used primarily
* for debugging purposes.
*
* \section LICENSE
*
* The Variable Infiltration Capacity (VIC) macroscale hydrological model
* Copyright (C) 2014  The Land Surface Hydrology Group, Department of Civil
* and Environmental Engineering, University of Washington.
*
* The VIC model is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with
* this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
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
