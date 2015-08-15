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

module vic_cesm_interface

! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: vic_cesm_init
  public :: vic_cesm_run
  public :: vic_cesm_final

    ! Init Interface
    interface
        integer (C_INT) function vic_cesm_init(vic_clock) BIND(C, name='vic_cesm_init')
        use, intrinsic :: ISO_C_BINDING
        use vic_cesm_def_mod
        implicit none
          type(vic_clock_type), intent(inout) :: vic_clock
!         type (C_PTR), value :: sendbuf
!         integer (C_INT), value :: sendcount
!         type (C_PTR), value :: recvcounts
     end function vic_cesm_init
   end interface

    ! Run Interface
    interface
        integer (C_INT) function vic_cesm_run(vic_clock) BIND(C, name='vic_cesm_run')
        use, intrinsic :: ISO_C_BINDING
        use vic_cesm_def_mod
        implicit none
          type(vic_clock_type), intent(inout) :: vic_clock
!         type (C_PTR), value :: sendbuf
!         integer (C_INT), value :: sendcount
!         type (C_PTR), value :: recvcounts
     end function vic_cesm_run
   end interface

    ! Run Interface
    interface
        integer (C_INT) function vic_cesm_final() BIND(C, name='vic_cesm_final')
        use, intrinsic :: ISO_C_BINDING
        use vic_cesm_def_mod
        implicit none
     end function vic_cesm_final
   end interface

end module vic_cesm_interface
