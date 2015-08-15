!>
!! @section DESCRIPTION
!!
!! Fortran definitions of VIC C structures
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

module vic_cesm_def_mod

use, intrinsic :: iso_c_binding

! type declarations (see cesm_def.h for more on the c/fortran type defs)
    type, bind(C) :: vic_clock_type
        integer(C_INT) :: timestep;             ! timestep in seconds
        integer(C_INT16_T) :: start_year;       ! start year
        integer(C_INT16_T) :: start_month;      ! start month
        integer(C_INT16_T) :: start_day;        ! start day
        integer(C_INT) :: start_dayseconds;     ! start dayseconds
        integer(C_INT16_T) :: stop_year;        ! stop year
        integer(C_INT16_T) :: stop_month;       ! stop month
        integer(C_INT16_T) :: stop_day;         ! stop day
        integer(C_INT) :: stop_dayseconds;      ! stop dayseconds
        integer(C_INT16_T) :: current_year;     ! current year
        integer(C_INT16_T) :: current_month;    ! current month
        integer(C_INT16_T) :: current_day;      ! current day
        integer(C_INT) :: current_dayseconds;   ! current dayseconds
        integer(C_INT16_T) :: previous_year;    ! previous year
        integer(C_INT16_T) :: previous_month;   ! previous month
        integer(C_INT16_T) :: previous_day;     ! previous day
        integer(C_INT) :: previous_dayseconds;  ! previous dayseconds
        logical(C_BOOL) :: state_flag;          ! state flag
        logical(C_BOOL) :: stop_flag;           ! stop flag
    end type

end module vic_cesm_def_mod
