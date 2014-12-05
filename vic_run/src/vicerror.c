/******************************************************************************
* @section DESCRIPTION
*
* This subroutine was written to handle numerical errors within the VIC model.
*
* This will flush all file buffers so that all records that have been run will
* be written to disk before the model is exited.
*
* @section LICENSE
*
* The Variable Infiltration Capacity (VIC) macroscale hydrological model
* Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
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

#include <vic_def.h>
#include <vic_run.h>

/******************************************************************************
* @brief        This subroutine handles numerical errors within the model.
******************************************************************************/
void
vicerror(char error_text[])
{
    extern option_struct options;
    extern Error_struct  Error;
    filenames_struct     fnames;
    void _exit();

    options.COMPRESS = false;   /* turn off compression of last set of files */

    fprintf(stderr, "VIC model run-time error...\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now writing output files...\n");
    close_files(&(Error.filep), Error.out_data_files, &fnames);
    fprintf(stderr, "...now exiting to system...\n");
    fflush(stdout);
    fflush(stderr);
    _exit(1);
}
