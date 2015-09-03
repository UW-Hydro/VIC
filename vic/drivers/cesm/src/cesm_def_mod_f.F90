!>
!! @section DESCRIPTION
!!
!! Fortran definitions of VIC C structures see cesm_def.h for more on the
!! c/fortran type defs
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

MODULE vic_cesm_def_mod

  USE, INTRINSIC :: iso_c_binding

  INTEGER,PARAMETER :: VICMAXSTRING = 2048 ! same as MAXSTRING in vic_def.h

  !--------------------------------------------------------------------------
  !> @brief   This structure stores clock information. See also ESMF_Clock.
  !! @note    Order is important and any changes here must be echoed in
  !!          vic_cesm_def.h
  !--------------------------------------------------------------------------
  TYPE, bind(C) :: vic_clock
    INTEGER(C_INT)                                         :: timestep            !< timestep in seconds
    INTEGER(C_INT16_T)                                     :: current_year        !< current year
    INTEGER(C_INT16_T)                                     :: current_month       !< current month
    INTEGER(C_INT16_T)                                     :: current_day         !< current day
    INTEGER(C_INT)                                         :: current_dayseconds  !< current dayseconds
    LOGICAL(C_BOOL)                                        :: state_flag          !< state flag
    LOGICAL(C_BOOL)                                        :: stop_flag           !< stop flag
    CHARACTER(len=VICMAXSTRING, kind=C_CHAR), DIMENSION(1) :: calendar            !< calendar
  END TYPE vic_clock

  !--------------------------------------------------------------------------
  !> @brief   This type stores location information for individual grid cells.
  !!          See also vic_cesm_def.h.
  !! @note    Order is important and any changes here must be echoed in
  !!          vic_cesm_def.h
  !--------------------------------------------------------------------------
  TYPE, bind(C) :: location_struct
    REAL(C_DOUBLE)    :: latitude     !< latitude of grid cell center
    REAL(C_DOUBLE)    :: longitude    !< longitude of grid cell center
    REAL(C_DOUBLE)    :: area         !< area of grid cell
    REAL(C_DOUBLE)    :: frac         !< fraction of grid cell that is active
    INTEGER(C_SIZE_T) :: nveg         !< number of vegetation TYPE according to PARAMETER file
    INTEGER(C_SIZE_T) :: global_idx   !< index of grid cell in global list of grid cells
    INTEGER(C_SIZE_T) :: io_idx       !< index of cell in 1-D I/O arrays
    INTEGER(C_SIZE_T) :: local_idx    !< index of grid cell in local list of grid cells
  END TYPE location_struct

  !--------------------------------------------------------------------------
  !> @brief   This type stores domain information. See also
  !!          vic_cesm_def.h.
  !! @note    Order is important and any changes here must be echoed in
  !!          vic_cesm_def.h
  !--------------------------------------------------------------------------
  TYPE, bind(C) :: domain_struct
    INTEGER(C_SIZE_T) :: ncells     !< number of active grid cell in domain
    INTEGER(C_SIZE_T) :: n_nx       !< size of x-index
    INTEGER(C_SIZE_T) :: n_ny       !< size of y-index
    TYPE(C_PTR)       :: locations  !< locations structs for local domain
  END TYPE domain_struct

  !--------------------------------------------------------------------------
  !> @brief   This structure is a c type container for the x2l fields.
  !! @note    Order is important and any changes here must be echoed in
  !!          vic_cesm_def_mod_f.F90
  !--------------------------------------------------------------------------
  TYPE, bind(C) :: x2l_data_struct
    REAL(C_DOUBLE)  :: x2l_Sa_z           !<bottom atm level height
    REAL(C_DOUBLE)  :: x2l_Sa_u           !<bottom atm level zon wind
    REAL(C_DOUBLE)  :: x2l_Sa_v           !<bottom atm level mer wind
    REAL(C_DOUBLE)  :: x2l_Sa_ptem        !<bottom atm level pot temp
    REAL(C_DOUBLE)  :: x2l_Sa_shum        !<bottom atm level spec hum
    REAL(C_DOUBLE)  :: x2l_Sa_pbot        !<bottom atm level pressure
    REAL(C_DOUBLE)  :: x2l_Sa_tbot        !<bottom atm level temp
    REAL(C_DOUBLE)  :: x2l_Faxa_lwdn      !<downward lw heat flux
    REAL(C_DOUBLE)  :: x2l_Faxa_rainc     !<prec: liquid "convective"
    REAL(C_DOUBLE)  :: x2l_Faxa_rainl     !<prec: liquid "large scale"
    REAL(C_DOUBLE)  :: x2l_Faxa_snowc     !<prec: frozen "convective"
    REAL(C_DOUBLE)  :: x2l_Faxa_snowl     !<prec: frozen "large scale"
    REAL(C_DOUBLE)  :: x2l_Faxa_swndr     !<sw: nir direct  downward
    REAL(C_DOUBLE)  :: x2l_Faxa_swvdr     !<sw: vis direct  downward
    REAL(C_DOUBLE)  :: x2l_Faxa_swndf     !<sw: nir diffuse downward
    REAL(C_DOUBLE)  :: x2l_Faxa_swvdf     !<sw: vis diffuse downward
    REAL(C_DOUBLE)  :: x2l_Sa_co2prog     !<bottom atm level prognostic co2
    REAL(C_DOUBLE)  :: x2l_Sa_co2diag     !<bottom atm level diagnostic co2
    REAL(C_DOUBLE)  :: x2l_Faxa_bcphidry  !<flux: Black Carbon hydrophilic dry deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_bcphodry  !<flux: Black Carbon hydrophobic dry deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_bcphiwet  !<flux: Black Carbon hydrophilic wet deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_ocphidry  !<flux: Organic Carbon hydrophilic dry deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_ocphodry  !<flux: Organic Carbon hydrophobic dry deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_ocphiwet  !<flux: Organic Carbon hydrophilic dry deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_dstwet1   !<flux: Size 1 dust -- wet deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_dstwet2   !<flux: Size 2 dust -- wet deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_dstwet3   !<flux: Size 3 dust -- wet deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_dstwet4   !<flux: Size 4 dust -- wet deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_dstdry1   !<flux: Size 1 dust -- dry deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_dstdry2   !<flux: Size 2 dust -- dry deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_dstdry3   !<flux: Size 3 dust -- dry deposition
    REAL(C_DOUBLE)  :: x2l_Faxa_dstdry4   !<flux: Size 4 dust -- dry deposition
    REAL(C_DOUBLE)  :: x2l_Flrr_flood     !<rtm->lnd rof (flood) flux
    LOGICAL(C_BOOL) :: x2l_vars_set       !< x2l set flag
  END TYPE x2l_data_struct

  !--------------------------------------------------------------------------
  !> @brief   This structure is a c type container for the l2x fields.
  !! @note    Order is important and any changes here must be echoed in
  !!          vic_cesm_def_mod_f.F90
  !--------------------------------------------------------------------------
  TYPE, bind(C) :: l2x_data_struct
    REAL(C_DOUBLE)  :: l2x_Sl_t           !< temperature
    REAL(C_DOUBLE)  :: l2x_Sl_tref        !< 2m reference temperature
    REAL(C_DOUBLE)  :: l2x_Sl_qref        !< 2m reference specific humidity
    REAL(C_DOUBLE)  :: l2x_Sl_avsdr       !< albedo: direct , visible
    REAL(C_DOUBLE)  :: l2x_Sl_anidr       !< albedo: direct , near-ir
    REAL(C_DOUBLE)  :: l2x_Sl_avsdf       !< albedo: diffuse, visible
    REAL(C_DOUBLE)  :: l2x_Sl_anidf       !< albedo: diffuse, near-ir
    REAL(C_DOUBLE)  :: l2x_Sl_snowh       !< snow height
    REAL(C_DOUBLE)  :: l2x_Sl_u10         !< 10m wind
    REAL(C_DOUBLE)  :: l2x_Sl_ddvel       !< dry deposition velocities (optional)
    REAL(C_DOUBLE)  :: l2x_Sl_fv          !< friction velocity
    REAL(C_DOUBLE)  :: l2x_Sl_ram1        !< aerodynamical resistance
    REAL(C_DOUBLE)  :: l2x_Sl_logz0       !< log z0
    REAL(C_DOUBLE)  :: l2x_Fall_taux      !< wind stress, zonal
    REAL(C_DOUBLE)  :: l2x_Fall_tauy      !< wind stress, meridional
    REAL(C_DOUBLE)  :: l2x_Fall_lat       !< latent          heat flux
    REAL(C_DOUBLE)  :: l2x_Fall_sen       !< sensible        heat flux
    REAL(C_DOUBLE)  :: l2x_Fall_lwup      !< upward longwave heat flux
    REAL(C_DOUBLE)  :: l2x_Fall_evap      !< evaporation     water flux
    REAL(C_DOUBLE)  :: l2x_Fall_swnet     !< heat flux       shortwave net
    REAL(C_DOUBLE)  :: l2x_Fall_fco2_lnd  !< co2 flux **For testing set to 0
    REAL(C_DOUBLE)  :: l2x_Fall_flxdst1   !< dust flux size bin 1
    REAL(C_DOUBLE)  :: l2x_Fall_flxdst2   !< dust flux size bin 2
    REAL(C_DOUBLE)  :: l2x_Fall_flxdst3   !< dust flux size bin 3
    REAL(C_DOUBLE)  :: l2x_Fall_flxdst4   !< dust flux size bin 4
    REAL(C_DOUBLE)  :: l2x_Fall_flxvoc    !< MEGAN fluxes
    REAL(C_DOUBLE)  :: l2x_Flrl_rofliq    !< lnd->rtm input fluxes
    REAL(C_DOUBLE)  :: l2x_Flrl_rofice    !< lnd->rtm input fluxes
    LOGICAL(C_BOOL) :: l2x_vars_set       !< l2x set flag
  END TYPE l2x_data_struct

END MODULE vic_cesm_def_mod
