!>
!! @section DESCRIPTION
!!
!! Fortran print library for CESM VIC coupling
!!
!! @section LICENSE
!!
!! The Variable Infiltration Capacity (VIC) macroscale hydrological model
!! Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
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

MODULE vic_cesm_print_library

  USE vic_cesm_def_mod
  USE, INTRINSIC :: iso_c_binding

  IMPLICIT NONE
  SAVE
  PRIVATE ! except

  PUBLIC :: set_print_library_iulog
  PUBLIC :: print_vic_clock
  PUBLIC :: print_domain
  PUBLIC :: print_location
  PUBLIC :: print_x2l_data
  PUBLIC :: print_l2x_data

  INTEGER :: iulog = 6  ! "stdout" log file unit number, default is 6

CONTAINS


  !--------------------------------------------------------------------------
  !> @brief   Set print library iulog
  !--------------------------------------------------------------------------
  SUBROUTINE set_print_library_iulog(unit_num)

    INTEGER, INTENT(in) :: unit_num

    iulog = unit_num

  END SUBROUTINE set_print_library_iulog


  !--------------------------------------------------------------------------
  !> @brief   Print the VIC clock type
  !--------------------------------------------------------------------------
  SUBROUTINE print_vic_clock(vclock)

    TYPE(vic_clock), INTENT(in) :: vclock

    WRITE(iulog, *) 'vic_clock              :'
    WRITE(iulog, *) '    timestep           : ', vclock%timestep
    WRITE(iulog, *) '    current_year       : ', vclock%current_year
    WRITE(iulog, *) '    current_month      : ', vclock%current_month
    WRITE(iulog, *) '    current_day        : ', vclock%current_day
    WRITE(iulog, *) '    current_dayseconds : ', vclock%current_dayseconds
    WRITE(iulog, *) '    state_flag         : ', vclock%state_flag
    WRITE(iulog, *) '    stop_flag          : ', vclock%stop_flag
    WRITE(iulog, *) '    calendar           : ', TRIM(Copy_a2s(vclock%calendar))

  END SUBROUTINE print_vic_clock


  !--------------------------------------------------------------------------
  !> @brief   Print the VIC domain type
  !--------------------------------------------------------------------------
  SUBROUTINE print_domain(domain)

    TYPE(domain_struct), INTENT(in) :: domain

    INTEGER :: i
    TYPE(location_struct), DIMENSION(:), POINTER :: locations

    CALL c_f_pointer(domain%locations, locations, [domain%ncells_active])

    WRITE(iulog, *) 'domain         :'
    WRITE(iulog, *) '    ncells     : ', domain%ncells_active
    WRITE(iulog, *) '    n_nx       : ', domain%n_nx
    WRITE(iulog, *) '    n_ny       : ', domain%n_ny
    WRITE(iulog, *) '    locations  :'

    DO i = 1, domain%ncells_active
       WRITE(iulog, *) 'location ', i, '(Fortran indexing)'
       CALL print_location(locations(i))
    END DO

  END SUBROUTINE print_domain


  !--------------------------------------------------------------------------
  !> @brief   Print the VIC location struct
  !--------------------------------------------------------------------------
  SUBROUTINE print_location(location)

    TYPE(location_struct), INTENT(in) :: location

    WRITE(iulog, *) 'location             :'
    WRITE(iulog, *) '    latitude   : ', location%latitude
    WRITE(iulog, *) '    longitude  : ', location%longitude
    WRITE(iulog, *) '    area       : ', location%area
    WRITE(iulog, *) '    frac       : ', location%frac
    WRITE(iulog, *) '    nveg       : ', location%nveg
    WRITE(iulog, *) '    global_idx : ', location%global_idx
    WRITE(iulog, *) '    io_idx     : ', location%io_idx
    WRITE(iulog, *) '    local_idx  : ', location%local_idx

  END SUBROUTINE print_location


  !--------------------------------------------------------------------------
  !> @brief   Print the VIC x2l data type
  !--------------------------------------------------------------------------
  SUBROUTINE print_x2l_data(x2l_data)

    TYPE(x2l_data_struct), INTENT(in) :: x2l_data

    WRITE(iulog, *) 'x2l_data               :'
    WRITE(iulog, *) '    x2l_Sa_z           : ', x2l_data%x2l_Sa_z
    WRITE(iulog, *) '    x2l_Sa_u           : ', x2l_data%x2l_Sa_u
    WRITE(iulog, *) '    x2l_Sa_v           : ', x2l_data%x2l_Sa_v
    WRITE(iulog, *) '    x2l_Sa_ptem        : ', x2l_data%x2l_Sa_ptem
    WRITE(iulog, *) '    x2l_Sa_shum        : ', x2l_data%x2l_Sa_shum
    WRITE(iulog, *) '    x2l_Sa_pbot        : ', x2l_data%x2l_Sa_pbot
    WRITE(iulog, *) '    x2l_Sa_tbot        : ', x2l_data%x2l_Sa_tbot
    WRITE(iulog, *) '    x2l_Faxa_lwdn      : ', x2l_data%x2l_Faxa_lwdn
    WRITE(iulog, *) '    x2l_Faxa_rainc     : ', x2l_data%x2l_Faxa_rainc
    WRITE(iulog, *) '    x2l_Faxa_rainl     : ', x2l_data%x2l_Faxa_rainl
    WRITE(iulog, *) '    x2l_Faxa_snowc     : ', x2l_data%x2l_Faxa_snowc
    WRITE(iulog, *) '    x2l_Faxa_snowl     : ', x2l_data%x2l_Faxa_snowl
    WRITE(iulog, *) '    x2l_Faxa_swndr     : ', x2l_data%x2l_Faxa_swndr
    WRITE(iulog, *) '    x2l_Faxa_swvdr     : ', x2l_data%x2l_Faxa_swvdr
    WRITE(iulog, *) '    x2l_Faxa_swndf     : ', x2l_data%x2l_Faxa_swndf
    WRITE(iulog, *) '    x2l_Faxa_swvdf     : ', x2l_data%x2l_Faxa_swvdf
    WRITE(iulog, *) '    x2l_Sa_co2prog     : ', x2l_data%x2l_Sa_co2prog
    WRITE(iulog, *) '    x2l_Sa_co2diag     : ', x2l_data%x2l_Sa_co2diag
    WRITE(iulog, *) '    x2l_Faxa_bcphidry  : ', x2l_data%x2l_Faxa_bcphidry
    WRITE(iulog, *) '    x2l_Faxa_bcphodry  : ', x2l_data%x2l_Faxa_bcphodry
    WRITE(iulog, *) '    x2l_Faxa_bcphiwet  : ', x2l_data%x2l_Faxa_bcphiwet
    WRITE(iulog, *) '    x2l_Faxa_ocphidry  : ', x2l_data%x2l_Faxa_ocphidry
    WRITE(iulog, *) '    x2l_Faxa_ocphodry  : ', x2l_data%x2l_Faxa_ocphodry
    WRITE(iulog, *) '    x2l_Faxa_ocphiwet  : ', x2l_data%x2l_Faxa_ocphiwet
    WRITE(iulog, *) '    x2l_Faxa_dstwet1   : ', x2l_data%x2l_Faxa_dstwet1
    WRITE(iulog, *) '    x2l_Faxa_dstwet2   : ', x2l_data%x2l_Faxa_dstwet2
    WRITE(iulog, *) '    x2l_Faxa_dstwet3   : ', x2l_data%x2l_Faxa_dstwet3
    WRITE(iulog, *) '    x2l_Faxa_dstwet4   : ', x2l_data%x2l_Faxa_dstwet4
    WRITE(iulog, *) '    x2l_Faxa_dstdry1   : ', x2l_data%x2l_Faxa_dstdry1
    WRITE(iulog, *) '    x2l_Faxa_dstdry2   : ', x2l_data%x2l_Faxa_dstdry2
    WRITE(iulog, *) '    x2l_Faxa_dstdry3   : ', x2l_data%x2l_Faxa_dstdry3
    WRITE(iulog, *) '    x2l_Faxa_dstdry4   : ', x2l_data%x2l_Faxa_dstdry4
    WRITE(iulog, *) '    x2l_Flrr_flood     : ', x2l_data%x2l_Flrr_flood

  END SUBROUTINE print_x2l_data


  !--------------------------------------------------------------------------
  !> @brief   Print the VIC l2x data type
  !--------------------------------------------------------------------------
  SUBROUTINE print_l2x_data(l2x_data)

    TYPE(l2x_data_struct), INTENT(in) :: l2x_data

    WRITE(iulog, *) 'l2x_data               :'
    WRITE(iulog, *) '    l2x_Sl_t           : ', l2x_data%l2x_Sl_t
    WRITE(iulog, *) '    l2x_Sl_tref        : ', l2x_data%l2x_Sl_tref
    WRITE(iulog, *) '    l2x_Sl_qref        : ', l2x_data%l2x_Sl_qref
    WRITE(iulog, *) '    l2x_Sl_avsdr       : ', l2x_data%l2x_Sl_avsdr
    WRITE(iulog, *) '    l2x_Sl_anidr       : ', l2x_data%l2x_Sl_anidr
    WRITE(iulog, *) '    l2x_Sl_avsdf       : ', l2x_data%l2x_Sl_avsdf
    WRITE(iulog, *) '    l2x_Sl_anidf       : ', l2x_data%l2x_Sl_anidf
    WRITE(iulog, *) '    l2x_Sl_snowh       : ', l2x_data%l2x_Sl_snowh
    WRITE(iulog, *) '    l2x_Sl_u10         : ', l2x_data%l2x_Sl_u10
    WRITE(iulog, *) '    l2x_Sl_ddvel       : ', l2x_data%l2x_Sl_ddvel
    WRITE(iulog, *) '    l2x_Sl_fv          : ', l2x_data%l2x_Sl_fv
    WRITE(iulog, *) '    l2x_Sl_ram1        : ', l2x_data%l2x_Sl_ram1
    WRITE(iulog, *) '    l2x_Sl_logz0       : ', l2x_data%l2x_Sl_logz0
    WRITE(iulog, *) '    l2x_Fall_taux      : ', l2x_data%l2x_Fall_taux
    WRITE(iulog, *) '    l2x_Fall_tauy      : ', l2x_data%l2x_Fall_tauy
    WRITE(iulog, *) '    l2x_Fall_lat       : ', l2x_data%l2x_Fall_lat
    WRITE(iulog, *) '    l2x_Fall_sen       : ', l2x_data%l2x_Fall_sen
    WRITE(iulog, *) '    l2x_Fall_lwup      : ', l2x_data%l2x_Fall_lwup
    WRITE(iulog, *) '    l2x_Fall_evap      : ', l2x_data%l2x_Fall_evap
    WRITE(iulog, *) '    l2x_Fall_swnet     : ', l2x_data%l2x_Fall_swnet
    WRITE(iulog, *) '    l2x_Fall_fco2_lnd  : ', l2x_data%l2x_Fall_fco2_lnd
    WRITE(iulog, *) '    l2x_Fall_flxdst1   : ', l2x_data%l2x_Fall_flxdst1
    WRITE(iulog, *) '    l2x_Fall_flxdst2   : ', l2x_data%l2x_Fall_flxdst2
    WRITE(iulog, *) '    l2x_Fall_flxdst3   : ', l2x_data%l2x_Fall_flxdst3
    WRITE(iulog, *) '    l2x_Fall_flxdst4   : ', l2x_data%l2x_Fall_flxdst4
    WRITE(iulog, *) '    l2x_Fall_flxvoc    : ', l2x_data%l2x_Fall_flxvoc
    WRITE(iulog, *) '    l2x_Flrl_rofliq    : ', l2x_data%l2x_Flrl_rofliq
    WRITE(iulog, *) '    l2x_Flrl_rofice    : ', l2x_data%l2x_Flrl_rofice

  END SUBROUTINE print_l2x_data

  !--------------------------------------------------------------------------
  !> @brief   copy char array to string
  !--------------------------------------------------------------------------
  PURE FUNCTION Copy_a2s(a)  RESULT (s)    ! copy char array to string
    CHARACTER,INTENT(IN) :: a(:)
    CHARACTER(SIZE(a)) :: s
    INTEGER :: i
    DO i = 1,SIZE(a)
      s(i:i) = a(i)
    END DO
  END FUNCTION Copy_a2s

END MODULE vic_cesm_print_library
