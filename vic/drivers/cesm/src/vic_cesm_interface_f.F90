!>
!! @section DESCRIPTION
!!
!! Fortran interface for CESM driver.
!!
!! @section LICENSE
!!
!! The Variable Infiltration Capacity (VIC) macroscale hydrological model
!! Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
!! and Environmental Engineering, University of Washington.
!!
!! The VIC model is free software; you can redistribute it and/or
!! modify it under the terms of the GNU General Public License
!! as published by the Free Software Foundation; either version 2
!! of the License, or (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License along with
!! this program; if not, write to the Free Software Foundation, Inc.,
!! 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

MODULE vic_cesm_interface

  ! !PUBLIC TYPES:
  IMPLICIT NONE
  SAVE
  PRIVATE ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  PUBLIC :: vic_cesm_init_mpi
  PUBLIC :: vic_cesm_init
  PUBLIC :: vic_cesm_run
  PUBLIC :: vic_cesm_final

  ! Init MPI Interface
  INTERFACE
     INTEGER(C_INT) FUNCTION vic_cesm_init_mpi(MPI_COMM_VIC_F) BIND(C, name='vic_cesm_init_mpi')
       USE, INTRINSIC :: ISO_C_BINDING
       USE vic_cesm_def_mod
       IMPLICIT NONE
       INTEGER(C_INT), INTENT(inout) :: MPI_COMM_VIC_F
     END FUNCTION vic_cesm_init_mpi
  END INTERFACE

  ! Init Interface
  INTERFACE
     INTEGER(C_INT) FUNCTION vic_cesm_init(vic_global_param_file, caseid, &
                                           runtype, vclock) BIND(C, name='vic_cesm_init')
       USE, INTRINSIC :: ISO_C_BINDING
       USE vic_cesm_def_mod
       IMPLICIT NONE
       CHARACTER(len=1, kind=C_CHAR), DIMENSION(*), INTENT(in) :: vic_global_param_file
       CHARACTER(len=1, kind=C_CHAR), DIMENSION(*), INTENT(in) :: caseid
       CHARACTER(len=1, kind=C_CHAR), DIMENSION(*), INTENT(in) :: runtype
       TYPE(vic_clock), INTENT(in) :: vclock
       END FUNCTION vic_cesm_init
  END INTERFACE

  ! Run Interface
  INTERFACE
     INTEGER(C_INT) FUNCTION vic_cesm_run(vclock) BIND(C, name='vic_cesm_run')
       USE, INTRINSIC :: ISO_C_BINDING
       USE vic_cesm_def_mod
       IMPLICIT NONE
       TYPE(vic_clock), INTENT(inout) :: vclock
     END FUNCTION vic_cesm_run
  END INTERFACE

  ! Finalize Interface
  INTERFACE
     INTEGER(C_INT) FUNCTION vic_cesm_final() BIND(C, name='vic_cesm_final')
       USE, INTRINSIC :: ISO_C_BINDING
       USE vic_cesm_def_mod
       IMPLICIT NONE
     END FUNCTION vic_cesm_final
  END INTERFACE

END MODULE vic_cesm_interface
