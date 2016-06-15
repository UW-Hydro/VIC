/******************************************************************************
 * @section DESCRIPTION
 *
 * MPI support routines for VIC's CESM driver
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
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
 *****************************************************************************/

#include <vic_driver_cesm.h>

/******************************************************************************
 * @brief   Wrapper around VIC's initialize mpi function, first translating the
            Fortran MPI Communicator to C.
 *****************************************************************************/
void
initialize_vic_cesm_mpi(MPI_Fint *MPI_COMM_VIC_F)
{
    extern MPI_Comm MPI_COMM_VIC;

    MPI_COMM_VIC = MPI_Comm_f2c(*MPI_COMM_VIC_F);

    initialize_mpi();
}
