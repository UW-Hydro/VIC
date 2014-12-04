module CNPhenologyMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNPhenologyMod
!
! !DESCRIPTION:
! Module holding routines used in phenology model for coupled carbon
! nitrogen code.
!
! !USES:
  use clmtype
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varcon  , only: tfrz
  use clm_varpar  , only: numpft
  implicit none
  save
  private

! !PUBLIC MEMBER FUNCTIONS:
  public :: CNPhenologyInit      ! Initialization
  public :: CNPhenology          ! Update
!
! !REVISION HISTORY:
! 8/1/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: migrated all routines to vector data structures
! 2/4/08,  slevis: adding crop phenology from AgroIBIS
! 10/10/12: Adapted for use in VIC by Michael Brunke

! !PRIVATE DATA MEMBERS:

!  real(r8)           :: dt              ! radiation time step delta t (seconds)
  real(r8)           :: fracday         ! dtime as a fraction of day
  real(r8)           :: crit_dayl       ! critical daylength for offset (seconds)
  real(r8)           :: ndays_on        ! number of days to complete onset
  real(r8)           :: ndays_off       ! number of days to complete offset
  real(r8)           :: fstor2tran      ! fraction of storage to move to transfer on each onset
  real(r8)           :: crit_onset_fdd  ! critical number of freezing days
  real(r8)           :: crit_onset_swi  ! water stress days for offset trigger
  real(r8)           :: soilpsi_on      ! water potential for onset trigger (MPa)
  real(r8)           :: crit_offset_fdd ! critical number of freezing degree days
                                        ! to trigger offset
  real(r8)           :: crit_offset_swi ! water stress days for offset trigger
  real(r8)           :: soilpsi_off     ! water potential for offset trigger (MPa)
  real(r8)           :: lwtop           ! live wood turnover proportion (annual fraction)
  !
  ! CropPhenology variables and constants
  !
!  real(r8)           :: p1d, p1v            ! photoperiod factor constants for crop vernalization
!  real(r8)           :: hti                 ! cold hardening index threshold for vernalization
!  real(r8)           :: tbase               ! base temperature for vernalization
!  integer, parameter :: NOT_Planted   = 999 ! If not planted   yet in year
!  integer, parameter :: NOT_Harvested = 999 ! If not harvested yet in year
!  integer, parameter :: inNH       = 1      ! Northern Hemisphere
!  integer, parameter :: inSH       = 2      ! Southern Hemisphere
!  integer, pointer   :: inhemi(:)           ! Hemisphere that pft is in
!  integer            :: minplantjday(0:numpft,inSH) ! minimum planting julian day
!  integer            :: maxplantjday(0:numpft,inSH) ! maximum planting julian day
!  integer            :: jdayyrstart(inSH)           ! julian day of start of year

!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNPhenology
!
! !INTERFACE:
subroutine CNPhenology (yr, dt, num_soilc, num_soilp)
!
! !DESCRIPTION:
! Dynamic phenology routine for coupled carbon-nitrogen code (CN)
! 1. grass phenology
!
! !USES:
!
! !ARGUMENTS:
   integer, intent(in) :: yr              ! current year
   real(r8), intent(in) :: dt             ! timestep (s)
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 7/28/03: Created by Peter Thornton
! 9/05/03, Peter Thornton: moved from call with (p) to call with (c)
! 10/3/03, Peter Thornton: added subroutine calls for different phenology types
! 11/7/03, Peter Thornton: moved phenology type tests into phenology type
!    routines, and moved onset, offset, background litfall routines into
!    main phenology call.
! 10/10/12: Adapted for use in VIC by Michael Brunke
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
! local pointers to implicit in/out scalars
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   ! each of the following phenology type routines includes a filter
   ! to operate only on the relevant pfts

   call CNPhenologyClimate(yr, num_soilp)
   
   call CNEvergreenPhenology(yr, num_soilp)

   call CNSeasonDecidPhenology(dt, num_soilp)

   call CNStressDecidPhenology(yr, dt, num_soilp)

   ! the same onset and offset routines are called regardless of
   ! phenology type - they depend only on onset_flag, offset_flag, bglfr, and bgtr

   call CNOnsetGrowth(dt, num_soilp)

   call CNOffsetLitterfall(dt, num_soilp)

   call CNBackgroundLitterfall(num_soilp)

   call CNLivewoodTurnover(num_soilp)

   ! gather all pft-level litterfall fluxes to the column
   ! for litter C and N inputs

   call CNLitterToColumn(num_soilc)

end subroutine CNPhenology

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNPhenologyInit
!
! !INTERFACE:
subroutine CNPhenologyInit( dt, begp, endp )
!
! !DESCRIPTION:
! Initialization of CNPhenology. Must be called after time-manager is
! initialized, and after pftcon file is read in.
!
! !USES:
!   use surfrdMod       , only: crop_prog
   use clm_varcon      , only: secspday
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: dt         ! Time step
   integer, intent(IN) :: begp, endp ! Beginning and ending PFT index
! !CALLED FROM:
! subroutine clm_initialize2 in module clm_initializeMod.F90
!
! !REVISION HISTORY:
! 3/28/11: Created by Erik Kluzek
! 10/25/12: Adapted for use in VIC by Michael Brunke
!
! !LOCAL VARIABLES:
!EOP
!------------------------------------------------------------------------

    !
    ! Get time-step and what fraction of a day it is
    !
    fracday = dt/secspday

    ! set some local parameters - these will be moved into
    ! parameter file after testing

    ! -----------------------------------------
    ! Constants for CNSeasonDecidPhenology
    ! -----------------------------------------
    !
    ! critical daylength from Biome-BGC, v4.1.2
    crit_dayl = 39300._r8

    ! -----------------------------------------
    ! Constants for CNSeasonDecidPhenology and CNStressDecidPhenology
    ! -----------------------------------------
    ndays_on  = 30._r8
    ndays_off = 15._r8

    ! transfer parameters
    fstor2tran = 0.5_r8
    ! -----------------------------------------
    ! Constants for CNStressDecidPhenology
    ! -----------------------------------------

    ! onset parameters
    crit_onset_fdd = 15.0_r8
    ! critical onset gdd now being calculated as a function of annual
    ! average 2m temp.
    ! crit_onset_gdd = 150.0 ! c3 grass value
    ! crit_onset_gdd = 1000.0   ! c4 grass value
    crit_onset_swi = 15.0_r8
    soilpsi_on     = -2.0_r8

    ! offset parameters
    crit_offset_fdd = 15.0_r8
    crit_offset_swi = 15.0_r8
    soilpsi_off     = -2.0_r8

    ! -----------------------------------------
    ! Constants for CNLivewoodTurnover
    ! -----------------------------------------

    ! set the global parameter for livewood turnover rate
    ! define as an annual fraction (0.7), and convert to fraction per second
    lwtop = 0.7_r8 / 31536000.0_r8

    ! -----------------------------------------
    ! Call any subroutine specific initialization routines
    ! -----------------------------------------

end subroutine CNPhenologyInit
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNPhenologyClimate
!
! !INTERFACE:
subroutine CNPhenologyClimate (yr, num_soilp)
!
! !DESCRIPTION:
! For coupled carbon-nitrogen code (CN).
!
! !USES:
   use clm_time_manager, only: get_days_per_year
!   use CropRestMod     , only: CropRestYear, CropRestIncYear
!
! !ARGUMENTS:
   integer, intent(in) :: yr              ! current year
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 3/13/07: Created by Peter Thornton
! 10/10/12: Adapted for use in VIC by Michael Brunke
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)             ! pft vegetation type
   ! ecophysiological constants
   real(r8), pointer :: t_ref2m(:)         ! 2m air temperature (K)
   real(r8), pointer :: tempavg_t2m(:)     ! temp. avg 2m air temperature (K)
   real(r8), pointer :: gdd0(:)            ! growing deg. days base 0 deg C (ddays)
   real(r8), pointer :: gdd8(:)            !    "     "    "    "   8  "  "    "
   real(r8), pointer :: gdd10(:)           !    "     "    "    "  10  "  "    "
   real(r8), pointer :: gdd020(:)          ! 20-yr mean of gdd0 (ddays)
   real(r8), pointer :: gdd820(:)          ! 20-yr mean of gdd8 (ddays)
   real(r8), pointer :: gdd1020(:)         ! 20-yr mean of gdd10 (ddays)
   integer , pointer :: pgridcell(:)       ! pft's gridcell index
!
! local pointers to implicit in/out scalars
!
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: p                    ! indices
   integer, save :: nyrs = -999    ! number of years prognostic crop has run
   real(r8):: dayspyr              ! days per year (days)
   integer kyr                     ! current year
   integer kmo                     !         month of year  (1, ..., 12)
   integer kda                     !         day of month   (1, ..., 31)
   integer mcsec                   !         seconds of day (0, ..., seconds/day)
   real(r8), parameter :: yravg   = 20.0_r8      ! length of years to average for gdd
   real(r8), parameter :: yravgm1 = yravg-1.0_r8 ! minus 1 of above
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers to derived type arrays
   ivt                           => clm3%g%c%p%itype
   t_ref2m                       => clm3%g%c%p%pes%t_ref2m
   tempavg_t2m                   => clm3%g%c%p%pepv%tempavg_t2m

   gdd0                          => clm3%g%c%p%pps%gdd0
   gdd8                          => clm3%g%c%p%pps%gdd8
   gdd10                         => clm3%g%c%p%pps%gdd10
   gdd020                        => clm3%g%c%p%pps%gdd020
   gdd820                        => clm3%g%c%p%pps%gdd820
   gdd1020                       => clm3%g%c%p%pps%gdd1020
   pgridcell                     => clm3%g%c%p%gridcell

   ! set time steps

   dayspyr = get_days_per_year(yr)

   do p = 1,num_soilp
      tempavg_t2m(p) = tempavg_t2m(p) + t_ref2m(p) * (fracday/dayspyr)
   end do

end subroutine CNPhenologyClimate
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNEvergreenPhenology
!
! !INTERFACE:
subroutine CNEvergreenPhenology (yr, num_soilp)
!
! !DESCRIPTION:
! For coupled carbon-nitrogen code (CN).
!
! !USES:
   use clm_varcon      , only: secspday
   use clm_time_manager, only: get_days_per_year
!
! !ARGUMENTS:
   integer, intent(in) :: yr              ! current year
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/2/03: Created by Peter Thornton
! 10/10/12: Adapted for use in VIC by Michael Brunke
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)       ! pft vegetation type
   ! ecophysiological constants
   real(r8), pointer :: evergreen(:) ! binary flag for evergreen leaf habit (0 or 1)
   real(r8), pointer :: leaf_long(:) ! leaf longevity (yrs)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: bglfr(:)     ! background litterfall rate (1/s)
   real(r8), pointer :: bgtr(:)      ! background transfer growth rate (1/s)
   real(r8), pointer :: lgsf(:)      ! long growing season factor [0-1]
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   real(r8):: dayspyr                ! Days per year
   integer :: p                      ! indices
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers to derived type arrays
   ivt       => clm3%g%c%p%itype
   evergreen => pftcon%evergreen
   leaf_long => pftcon%leaf_long
   bglfr     => clm3%g%c%p%pepv%bglfr
   bgtr      => clm3%g%c%p%pepv%bgtr
   lgsf      => clm3%g%c%p%pepv%lgsf
   dayspyr   = get_days_per_year(yr)

   do p = 1,num_soilp
      if (evergreen(ivt(p)) == 1._r8) then
          bglfr(p) = 1._r8/(leaf_long(ivt(p))*dayspyr*secspday)
          bgtr(p)  = 0._r8
          lgsf(p)  = 0._r8
      end if
   end do

end subroutine CNEvergreenPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSeasonDecidPhenology
!
! !INTERFACE:
subroutine CNSeasonDecidPhenology (dt, num_soilp)
!
! !DESCRIPTION:
! For coupled carbon-nitrogen code (CN).
! This routine handles the seasonal deciduous phenology code (temperate
! deciduous vegetation that has only one growing season per year).
!
! !USES:
   use shr_const_mod   , only: SHR_CONST_TKFRZ, SHR_CONST_PI
   use clm_varcon      , only: secspday
!
! !ARGUMENTS:
   real(r8), intent(in) :: dt             ! timestep (s)
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/6/03: Created by Peter Thornton
! 10/24/03, Peter Thornton: migrated to vector data structures
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)                ! pft vegetation type
   integer , pointer :: pcolumn(:)            ! pft's column index
   integer , pointer :: pgridcell(:)          ! pft's gridcell index
   real(r8), pointer :: latdeg(:)             ! latitude (radians)
   real(r8), pointer :: decl(:)               ! solar declination (radians)
   real(r8), pointer :: t_soisno(:,:)         ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   real(r8), pointer :: soilpsi(:,:)          ! soil water potential in each soil layer (MPa)
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   ! ecophysiological constants
   real(r8), pointer :: season_decid(:) ! binary flag for seasonal-deciduous leaf habit (0 or 1)
   real(r8), pointer :: woody(:)        ! binary flag for woody lifeform (1=woody, 0=not woody)
!
! local pointers to implicit in/out scalars
   real(r8), pointer :: dormant_flag(:)    ! dormancy flag
   real(r8), pointer :: days_active(:)     ! number of days since last dormancy
   real(r8), pointer :: onset_flag(:)      ! onset flag
   real(r8), pointer :: onset_counter(:)   ! onset counter (seconds)
   real(r8), pointer :: onset_gddflag(:)   ! onset freeze flag
   real(r8), pointer :: onset_gdd(:)       ! onset growing degree days
   real(r8), pointer :: offset_flag(:)     ! offset flag
   real(r8), pointer :: offset_counter(:)  ! offset counter (seconds)
   real(r8), pointer :: dayl(:)            ! daylength (seconds)
   real(r8), pointer :: prev_dayl(:)       ! daylength from previous albedo timestep (seconds)
   real(r8), pointer :: annavg_t2m(:)      ! annual average 2m air temperature (K)
   real(r8), pointer :: prev_leafc_to_litter(:)  ! previous timestep leaf C litterfall flux (gC/m2/s)
   real(r8), pointer :: prev_frootc_to_litter(:) ! previous timestep froot C litterfall flux (gC/m2/s)
   real(r8), pointer :: lgsf(:)            ! long growing season factor [0-1]
   real(r8), pointer :: bglfr(:)           ! background litterfall rate (1/s)
   real(r8), pointer :: bgtr(:)            ! background transfer growth rate (1/s)
   real(r8), pointer :: leafc_xfer_to_leafc(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: leafn_xfer_to_leafn(:)
   real(r8), pointer :: frootn_xfer_to_frootn(:)
   real(r8), pointer :: livestemn_xfer_to_livestemn(:)
   real(r8), pointer :: deadstemn_xfer_to_deadstemn(:)
   real(r8), pointer :: livecrootn_xfer_to_livecrootn(:)
   real(r8), pointer :: deadcrootn_xfer_to_deadcrootn(:)
   real(r8), pointer :: leafc_xfer(:)      ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc_xfer(:)     ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc_xfer(:)  ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc_xfer(:)  ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc_xfer(:) ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc_xfer(:) ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: leafn_xfer(:)      ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn_xfer(:)     ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn_xfer(:)  ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn_xfer(:)  ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn_xfer(:) ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn_xfer(:) ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: leafc_storage_to_xfer(:)
   real(r8), pointer :: frootc_storage_to_xfer(:)
   real(r8), pointer :: livestemc_storage_to_xfer(:)
   real(r8), pointer :: deadstemc_storage_to_xfer(:)
   real(r8), pointer :: livecrootc_storage_to_xfer(:)
   real(r8), pointer :: deadcrootc_storage_to_xfer(:)
   real(r8), pointer :: gresp_storage_to_xfer(:)
   real(r8), pointer :: leafn_storage_to_xfer(:)
   real(r8), pointer :: frootn_storage_to_xfer(:)
   real(r8), pointer :: livestemn_storage_to_xfer(:)
   real(r8), pointer :: deadstemn_storage_to_xfer(:)
   real(r8), pointer :: livecrootn_storage_to_xfer(:)
   real(r8), pointer :: deadcrootn_storage_to_xfer(:)
!#if (defined CNDV)
!   logical , pointer :: pftmayexist(:)     ! exclude seasonal decid pfts from tropics
!#endif
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p            !indices
   real(r8):: ws_flag        !winter-summer solstice flag (0 or 1)
   real(r8):: crit_onset_gdd !critical onset growing degree-day sum
   real(r8):: soilt
   real(r8):: lat            !latitude (radians)
   real(r8):: temp           !temporary variable for daylength calculation

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
   ivt                           => clm3%g%c%p%itype
   pcolumn                       => clm3%g%c%p%column
   pgridcell                     => clm3%g%c%p%gridcell
   latdeg                        => clm3%g%latdeg
   decl                          => clm3%g%c%cps%decl
   t_soisno                      => clm3%g%c%ces%t_soisno
   leafc_storage                 => clm3%g%c%p%pcs%leafc_storage
   frootc_storage                => clm3%g%c%p%pcs%frootc_storage
   livestemc_storage             => clm3%g%c%p%pcs%livestemc_storage
   deadstemc_storage             => clm3%g%c%p%pcs%deadstemc_storage
   livecrootc_storage            => clm3%g%c%p%pcs%livecrootc_storage
   deadcrootc_storage            => clm3%g%c%p%pcs%deadcrootc_storage
   gresp_storage                 => clm3%g%c%p%pcs%gresp_storage
   leafn_storage                 => clm3%g%c%p%pns%leafn_storage
   frootn_storage                => clm3%g%c%p%pns%frootn_storage
   livestemn_storage             => clm3%g%c%p%pns%livestemn_storage
   deadstemn_storage             => clm3%g%c%p%pns%deadstemn_storage
   livecrootn_storage            => clm3%g%c%p%pns%livecrootn_storage
   deadcrootn_storage            => clm3%g%c%p%pns%deadcrootn_storage
   season_decid                  => pftcon%season_decid
   woody                         => pftcon%woody

   ! Assign local pointers to derived type arrays (out)
   dormant_flag                  => clm3%g%c%p%pepv%dormant_flag
   days_active                   => clm3%g%c%p%pepv%days_active
   onset_flag                    => clm3%g%c%p%pepv%onset_flag
   onset_counter                 => clm3%g%c%p%pepv%onset_counter
   onset_gddflag                 => clm3%g%c%p%pepv%onset_gddflag
   onset_gdd                     => clm3%g%c%p%pepv%onset_gdd
   offset_flag                   => clm3%g%c%p%pepv%offset_flag
   offset_counter                => clm3%g%c%p%pepv%offset_counter
   dayl                          => clm3%g%c%p%pepv%dayl
   prev_dayl                     => clm3%g%c%p%pepv%prev_dayl
   annavg_t2m                    => clm3%g%c%p%pepv%annavg_t2m
   prev_leafc_to_litter          => clm3%g%c%p%pepv%prev_leafc_to_litter
   prev_frootc_to_litter         => clm3%g%c%p%pepv%prev_frootc_to_litter
   bglfr                         => clm3%g%c%p%pepv%bglfr
   bgtr                          => clm3%g%c%p%pepv%bgtr
   lgsf                          => clm3%g%c%p%pepv%lgsf
   leafc_xfer_to_leafc           => clm3%g%c%p%pcf%leafc_xfer_to_leafc
   frootc_xfer_to_frootc         => clm3%g%c%p%pcf%frootc_xfer_to_frootc
   livestemc_xfer_to_livestemc   => clm3%g%c%p%pcf%livestemc_xfer_to_livestemc
   deadstemc_xfer_to_deadstemc   => clm3%g%c%p%pcf%deadstemc_xfer_to_deadstemc
   livecrootc_xfer_to_livecrootc => clm3%g%c%p%pcf%livecrootc_xfer_to_livecrootc
   deadcrootc_xfer_to_deadcrootc => clm3%g%c%p%pcf%deadcrootc_xfer_to_deadcrootc
   leafn_xfer_to_leafn           => clm3%g%c%p%pnf%leafn_xfer_to_leafn
   frootn_xfer_to_frootn         => clm3%g%c%p%pnf%frootn_xfer_to_frootn
   livestemn_xfer_to_livestemn   => clm3%g%c%p%pnf%livestemn_xfer_to_livestemn
   deadstemn_xfer_to_deadstemn   => clm3%g%c%p%pnf%deadstemn_xfer_to_deadstemn
   livecrootn_xfer_to_livecrootn => clm3%g%c%p%pnf%livecrootn_xfer_to_livecrootn
   deadcrootn_xfer_to_deadcrootn => clm3%g%c%p%pnf%deadcrootn_xfer_to_deadcrootn
   leafc_xfer                    => clm3%g%c%p%pcs%leafc_xfer
   frootc_xfer                   => clm3%g%c%p%pcs%frootc_xfer
   livestemc_xfer                => clm3%g%c%p%pcs%livestemc_xfer
   deadstemc_xfer                => clm3%g%c%p%pcs%deadstemc_xfer
   livecrootc_xfer               => clm3%g%c%p%pcs%livecrootc_xfer
   deadcrootc_xfer               => clm3%g%c%p%pcs%deadcrootc_xfer
   leafn_xfer                    => clm3%g%c%p%pns%leafn_xfer
   frootn_xfer                   => clm3%g%c%p%pns%frootn_xfer
   livestemn_xfer                => clm3%g%c%p%pns%livestemn_xfer
   deadstemn_xfer                => clm3%g%c%p%pns%deadstemn_xfer
   livecrootn_xfer               => clm3%g%c%p%pns%livecrootn_xfer
   deadcrootn_xfer               => clm3%g%c%p%pns%deadcrootn_xfer
   leafc_storage_to_xfer         => clm3%g%c%p%pcf%leafc_storage_to_xfer
   frootc_storage_to_xfer        => clm3%g%c%p%pcf%frootc_storage_to_xfer
   livestemc_storage_to_xfer     => clm3%g%c%p%pcf%livestemc_storage_to_xfer
   deadstemc_storage_to_xfer     => clm3%g%c%p%pcf%deadstemc_storage_to_xfer
   livecrootc_storage_to_xfer    => clm3%g%c%p%pcf%livecrootc_storage_to_xfer
   deadcrootc_storage_to_xfer    => clm3%g%c%p%pcf%deadcrootc_storage_to_xfer
   gresp_storage_to_xfer         => clm3%g%c%p%pcf%gresp_storage_to_xfer
   leafn_storage_to_xfer         => clm3%g%c%p%pnf%leafn_storage_to_xfer
   frootn_storage_to_xfer        => clm3%g%c%p%pnf%frootn_storage_to_xfer
   livestemn_storage_to_xfer     => clm3%g%c%p%pnf%livestemn_storage_to_xfer
   deadstemn_storage_to_xfer     => clm3%g%c%p%pnf%deadstemn_storage_to_xfer
   livecrootn_storage_to_xfer    => clm3%g%c%p%pnf%livecrootn_storage_to_xfer
   deadcrootn_storage_to_xfer    => clm3%g%c%p%pnf%deadcrootn_storage_to_xfer
!#if (defined CNDV)
!   pftmayexist                   => clm3%g%c%p%pdgvs%pftmayexist
!#endif

   ! start pft loop
   do p = 1,num_soilp
      c = pcolumn(p)

      if (season_decid(ivt(p)) == 1._r8) then

         ! set background litterfall rate, background transfer rate, and
         ! long growing season factor to 0 for seasonal deciduous types
         bglfr(p) = 0._r8
         bgtr(p) = 0._r8
         lgsf(p) = 0._r8

         ! onset gdd sum from Biome-BGC, v4.1.2
         crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_t2m(p) - SHR_CONST_TKFRZ))

         ! use solar declination information stored during Surface Albedo()
         ! and latitude from gps to calcluate daylength (convert latitude from degrees to radians)
         ! the constant 13750.9871 is the number of seconds per radian of hour-angle

         prev_dayl(p) = dayl(p)
         lat = (SHR_CONST_PI/180._r8)*latdeg(pgridcell(p))
         temp = -(sin(lat)*sin(decl(c)))/(cos(lat) * cos(decl(c)))
         temp = min(1._r8,max(-1._r8,temp))
         dayl(p) = 2.0_r8 * 13750.9871_r8 * acos(temp)

         ! set flag for solstice period (winter->summer = 1, summer->winter = 0)
         if (dayl(p) >= prev_dayl(p)) then
            ws_flag = 1._r8
         else
            ws_flag = 0._r8
         end if

         ! update offset_counter and test for the end of the offset period
         if (offset_flag(p) == 1.0_r8) then
            ! decrement counter for offset period
            offset_counter(p) = offset_counter(p) - dt

            ! if this is the end of the offset_period, reset phenology
            ! flags and indices
            if (offset_counter(p) == 0.0_r8) then
               ! this code block was originally handled by call cn_offset_cleanup(p)
               ! inlined during vectorization

               offset_flag(p) = 0._r8
               offset_counter(p) = 0._r8
               dormant_flag(p) = 1._r8
               days_active(p) = 0._r8
!#if (defined CNDV)
!               pftmayexist(p) = .true.
!#endif

               ! reset the previous timestep litterfall flux memory
               prev_leafc_to_litter(p) = 0._r8
               prev_frootc_to_litter(p) = 0._r8
            end if
         end if

         ! update onset_counter and test for the end of the onset period
         if (onset_flag(p) == 1.0_r8) then
            ! decrement counter for onset period
            onset_counter(p) = onset_counter(p) - dt

            ! if this is the end of the onset period, reset phenology
            ! flags and indices
            if (onset_counter(p) == 0.0_r8) then
               ! this code block was originally handled by call cn_onset_cleanup(p)
               ! inlined during vectorization

               onset_flag(p) = 0.0_r8
               onset_counter(p) = 0.0_r8
               ! set all transfer growth rates to 0.0
               leafc_xfer_to_leafc(p)   = 0.0_r8
               frootc_xfer_to_frootc(p) = 0.0_r8
               leafn_xfer_to_leafn(p)   = 0.0_r8
               frootn_xfer_to_frootn(p) = 0.0_r8
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_xfer_to_livestemc(p)   = 0.0_r8
                  deadstemc_xfer_to_deadstemc(p)   = 0.0_r8
                  livecrootc_xfer_to_livecrootc(p) = 0.0_r8
                  deadcrootc_xfer_to_deadcrootc(p) = 0.0_r8
                  livestemn_xfer_to_livestemn(p)   = 0.0_r8
                  deadstemn_xfer_to_deadstemn(p)   = 0.0_r8
                  livecrootn_xfer_to_livecrootn(p) = 0.0_r8
                  deadcrootn_xfer_to_deadcrootn(p) = 0.0_r8
               end if
               ! set transfer pools to 0.0
               leafc_xfer(p) = 0.0_r8
               leafn_xfer(p) = 0.0_r8
               frootc_xfer(p) = 0.0_r8
               frootn_xfer(p) = 0.0_r8
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_xfer(p) = 0.0_r8
                  livestemn_xfer(p) = 0.0_r8
                  deadstemc_xfer(p) = 0.0_r8
                  deadstemn_xfer(p) = 0.0_r8
                  livecrootc_xfer(p) = 0.0_r8
                  livecrootn_xfer(p) = 0.0_r8
                  deadcrootc_xfer(p) = 0.0_r8
                  deadcrootn_xfer(p) = 0.0_r8
               end if
            end if
         end if

         ! test for switching from dormant period to growth period
         if (dormant_flag(p) == 1.0_r8) then

            ! Test to turn on growing degree-day sum, if off.
            ! switch on the growing degree day sum on the winter solstice

            if (onset_gddflag(p) == 0._r8 .and. ws_flag == 1._r8) then
               onset_gddflag(p) = 1._r8
               onset_gdd(p) = 0._r8
            end if

            ! Test to turn off growing degree-day sum, if on.
            ! This test resets the growing degree day sum if it gets past
            ! the summer solstice without reaching the threshold value.
            ! In that case, it will take until the next winter solstice
            ! before the growing degree-day summation starts again.

            if (onset_gddflag(p) == 1._r8 .and. ws_flag == 0._r8) then
               onset_gddflag(p) = 0._r8
               onset_gdd(p) = 0._r8
            end if

            ! if the gdd flag is set, and if the soil is above freezing
            ! then accumulate growing degree days for onset trigger

            soilt = t_soisno(c,3)
            if (onset_gddflag(p) == 1.0_r8 .and. soilt > SHR_CONST_TKFRZ) then
               onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
            end if

            ! set onset_flag if critical growing degree-day sum is exceeded
            if (onset_gdd(p) > crit_onset_gdd) then
               onset_flag(p) = 1.0_r8
               dormant_flag(p) = 0.0_r8
               onset_gddflag(p) = 0.0_r8
               onset_gdd(p) = 0.0_r8
               onset_counter(p) = ndays_on * secspday

               ! move all the storage pools into transfer pools,
               ! where they will be transfered to displayed growth over the onset period.
               ! this code was originally handled with call cn_storage_to_xfer(p)
               ! inlined during vectorization

               ! set carbon fluxes for shifting storage pools to transfer pools
               leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
               frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
                  deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
                  livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
                  deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
                  gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
               end if

               ! set nitrogen fluxes for shifting storage pools to transfer pools
               leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
               frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
                  deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
                  livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
                  deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
               end if
            end if

         ! test for switching from growth period to offset period
         else if (offset_flag(p) == 0.0_r8) then
#if (defined CNDV)
            ! If days_active > 355, then remove pft in
            ! CNDVEstablishment at the end of the year.
            ! days_active > 355 is a symptom of seasonal decid. pfts occurring in
            ! gridcells where dayl never drops below crit_dayl.
            ! This results in TLAI>1e4 in a few gridcells.
            days_active(p) = days_active(p) + fracday
            if (days_active(p) > 355._r8) pftmayexist(p) = .false.
#endif

            ! only begin to test for offset daylength once past the summer sol
            if (ws_flag == 0._r8 .and. dayl(p) < crit_dayl) then
               offset_flag(p) = 1._r8
               offset_counter(p) = ndays_off * secspday
               prev_leafc_to_litter(p) = 0._r8
               prev_frootc_to_litter(p) = 0._r8
            end if
         end if

      end if ! end if seasonal deciduous

   end do ! end of pft loop

end subroutine CNSeasonDecidPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNStressDecidPhenology
!
! !INTERFACE:
subroutine CNStressDecidPhenology (yr, dt, num_soilp)
!
! !DESCRIPTION:
! This routine handles phenology for vegetation types, such as grasses and
! tropical drought deciduous trees, that respond to cold and drought stress
! signals and that can have multiple growing seasons in a given year.
! This routine allows for the possibility that leaves might persist year-round
! in the absence of a suitable stress trigger, by switching to an essentially
! evergreen habit, but maintaining a deciduous leaf longevity, while waiting
! for the next stress trigger.  This is in contrast to the seasonal deciduous
! algorithm (for temperate deciduous trees) that forces a single growing season
! per year.
!
! !USES:
   use clm_time_manager, only: get_days_per_year
   use clm_varcon      , only: secspday
   use shr_const_mod   , only: SHR_CONST_TKFRZ, SHR_CONST_PI
!
! !ARGUMENTS:
   integer, intent(in) :: yr              ! current year
   real(r8), intent(in) :: dt             ! timestep (s)
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/27/03: Created by Peter Thornton
! 01/29/04: Made onset_gdd critical sum a function of temperature, as in
!           seasonal deciduous algorithm.
! 10/10/12: Modified for use in VIC by Michael Brunke
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)                ! pft vegetation type
   integer , pointer :: pcolumn(:)            ! pft's column index
   integer , pointer :: pgridcell(:)          ! pft's gridcell index
   real(r8), pointer :: latdeg(:)             ! latitude (radians)
   real(r8), pointer :: decl(:)               ! solar declination (radians)
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: t_soisno(:,:)         ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   real(r8), pointer :: soilpsi(:,:)          ! soil water potential in each soil layer (MPa)
   real(r8), pointer :: leaf_long(:)          ! leaf longevity (yrs)
   real(r8), pointer :: stress_decid(:)       ! binary flag for stress-deciduous leaf habit (0 or 1)
   real(r8), pointer :: woody(:)              ! binary flag for woody lifeform (1=woody, 0=not woody)

!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: dormant_flag(:)    ! dormancy flag
   real(r8), pointer :: days_active(:)     ! number of days since last dormancy
   real(r8), pointer :: onset_flag(:)      ! onset flag
   real(r8), pointer :: onset_counter(:)   ! onset counter (seconds)
   real(r8), pointer :: onset_gddflag(:)   ! onset freeze flag
   real(r8), pointer :: onset_fdd(:)       ! onset freezing degree days counter
   real(r8), pointer :: onset_gdd(:)       ! onset growing degree days
   real(r8), pointer :: onset_swi(:)       ! onset soil water index
   real(r8), pointer :: offset_flag(:)     ! offset flag
   real(r8), pointer :: offset_counter(:)  ! offset counter (seconds)
   real(r8), pointer :: dayl(:)            ! daylength (seconds)
   real(r8), pointer :: offset_fdd(:)      ! offset freezing degree days counter
   real(r8), pointer :: offset_swi(:)      ! offset soil water index
   real(r8), pointer :: annavg_t2m(:)      ! annual average 2m air temperature (K)
   real(r8), pointer :: lgsf(:)            ! long growing season factor [0-1]
   real(r8), pointer :: bglfr(:)           ! background litterfall rate (1/s)
   real(r8), pointer :: bgtr(:)            ! background transfer growth rate (1/s)
   real(r8), pointer :: prev_leafc_to_litter(:)  ! previous timestep leaf C litterfall flux (gC/m2/s)
   real(r8), pointer :: prev_frootc_to_litter(:) ! previous timestep froot C litterfall flux (gC/m2/s)
   real(r8), pointer :: leafc_xfer_to_leafc(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: leafn_xfer_to_leafn(:)
   real(r8), pointer :: frootn_xfer_to_frootn(:)
   real(r8), pointer :: livestemn_xfer_to_livestemn(:)
   real(r8), pointer :: deadstemn_xfer_to_deadstemn(:)
   real(r8), pointer :: livecrootn_xfer_to_livecrootn(:)
   real(r8), pointer :: deadcrootn_xfer_to_deadcrootn(:)
   real(r8), pointer :: leafc_xfer(:)      ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc_xfer(:)     ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc_xfer(:)  ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc_xfer(:)  ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc_xfer(:) ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc_xfer(:) ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: leafn_xfer(:)      ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn_xfer(:)     ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn_xfer(:)  ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn_xfer(:)  ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn_xfer(:) ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn_xfer(:) ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: leafc_storage_to_xfer(:)
   real(r8), pointer :: frootc_storage_to_xfer(:)
   real(r8), pointer :: livestemc_storage_to_xfer(:)
   real(r8), pointer :: deadstemc_storage_to_xfer(:)
   real(r8), pointer :: livecrootc_storage_to_xfer(:)
   real(r8), pointer :: deadcrootc_storage_to_xfer(:)
   real(r8), pointer :: gresp_storage_to_xfer(:)
   real(r8), pointer :: leafn_storage_to_xfer(:)
   real(r8), pointer :: frootn_storage_to_xfer(:)
   real(r8), pointer :: livestemn_storage_to_xfer(:)
   real(r8), pointer :: deadstemn_storage_to_xfer(:)
   real(r8), pointer :: livecrootn_storage_to_xfer(:)
   real(r8), pointer :: deadcrootn_storage_to_xfer(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   real(r8),parameter :: secspqtrday = secspday / 4  ! seconds per quarter day
   integer :: c,p             ! indices
   real(r8):: dayspyr         ! days per year
   real(r8):: crit_onset_gdd  ! degree days for onset trigger
   real(r8):: soilt           ! temperature of top soil layer
   real(r8):: psi             ! water stress of top soil layer
   real(r8):: lat             !latitude (radians)
   real(r8):: temp            !temporary variable for daylength calculation
!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%c%p%itype
    pcolumn                        => clm3%g%c%p%column
    pgridcell                      => clm3%g%c%p%gridcell
    latdeg                         => clm3%g%latdeg
    decl                           => clm3%g%c%cps%decl
    leafc_storage                  => clm3%g%c%p%pcs%leafc_storage
    frootc_storage                 => clm3%g%c%p%pcs%frootc_storage
    livestemc_storage              => clm3%g%c%p%pcs%livestemc_storage
    deadstemc_storage              => clm3%g%c%p%pcs%deadstemc_storage
    livecrootc_storage             => clm3%g%c%p%pcs%livecrootc_storage
    deadcrootc_storage             => clm3%g%c%p%pcs%deadcrootc_storage
    gresp_storage                  => clm3%g%c%p%pcs%gresp_storage
    leafn_storage                  => clm3%g%c%p%pns%leafn_storage
    frootn_storage                 => clm3%g%c%p%pns%frootn_storage
    livestemn_storage              => clm3%g%c%p%pns%livestemn_storage
    deadstemn_storage              => clm3%g%c%p%pns%deadstemn_storage
    livecrootn_storage             => clm3%g%c%p%pns%livecrootn_storage
    deadcrootn_storage             => clm3%g%c%p%pns%deadcrootn_storage
    soilpsi                        => clm3%g%c%cps%soilpsi
    t_soisno                       => clm3%g%c%ces%t_soisno
    leaf_long                      => pftcon%leaf_long
    woody                          => pftcon%woody
    stress_decid                   => pftcon%stress_decid

   ! Assign local pointers to derived type arrays (out)
    dormant_flag                   => clm3%g%c%p%pepv%dormant_flag
    days_active                    => clm3%g%c%p%pepv%days_active
    onset_flag                     => clm3%g%c%p%pepv%onset_flag
    onset_counter                  => clm3%g%c%p%pepv%onset_counter
    onset_gddflag                  => clm3%g%c%p%pepv%onset_gddflag
    onset_fdd                      => clm3%g%c%p%pepv%onset_fdd
    onset_gdd                      => clm3%g%c%p%pepv%onset_gdd
    onset_swi                      => clm3%g%c%p%pepv%onset_swi
    offset_flag                    => clm3%g%c%p%pepv%offset_flag
    offset_counter                 => clm3%g%c%p%pepv%offset_counter
    dayl                           => clm3%g%c%p%pepv%dayl
    offset_fdd                     => clm3%g%c%p%pepv%offset_fdd
    offset_swi                     => clm3%g%c%p%pepv%offset_swi
    annavg_t2m                     => clm3%g%c%p%pepv%annavg_t2m
    prev_leafc_to_litter           => clm3%g%c%p%pepv%prev_leafc_to_litter
    prev_frootc_to_litter          => clm3%g%c%p%pepv%prev_frootc_to_litter
    lgsf                           => clm3%g%c%p%pepv%lgsf
    bglfr                          => clm3%g%c%p%pepv%bglfr
    bgtr                           => clm3%g%c%p%pepv%bgtr
    leafc_xfer_to_leafc            => clm3%g%c%p%pcf%leafc_xfer_to_leafc
    frootc_xfer_to_frootc          => clm3%g%c%p%pcf%frootc_xfer_to_frootc
    livestemc_xfer_to_livestemc    => clm3%g%c%p%pcf%livestemc_xfer_to_livestemc
    deadstemc_xfer_to_deadstemc    => clm3%g%c%p%pcf%deadstemc_xfer_to_deadstemc
    livecrootc_xfer_to_livecrootc  => clm3%g%c%p%pcf%livecrootc_xfer_to_livecrootc
    deadcrootc_xfer_to_deadcrootc  => clm3%g%c%p%pcf%deadcrootc_xfer_to_deadcrootc
    leafn_xfer_to_leafn            => clm3%g%c%p%pnf%leafn_xfer_to_leafn
    frootn_xfer_to_frootn          => clm3%g%c%p%pnf%frootn_xfer_to_frootn
    livestemn_xfer_to_livestemn    => clm3%g%c%p%pnf%livestemn_xfer_to_livestemn
    deadstemn_xfer_to_deadstemn    => clm3%g%c%p%pnf%deadstemn_xfer_to_deadstemn
    livecrootn_xfer_to_livecrootn  => clm3%g%c%p%pnf%livecrootn_xfer_to_livecrootn
    deadcrootn_xfer_to_deadcrootn  => clm3%g%c%p%pnf%deadcrootn_xfer_to_deadcrootn
    leafc_xfer                     => clm3%g%c%p%pcs%leafc_xfer
    frootc_xfer                    => clm3%g%c%p%pcs%frootc_xfer
    livestemc_xfer                 => clm3%g%c%p%pcs%livestemc_xfer
    deadstemc_xfer                 => clm3%g%c%p%pcs%deadstemc_xfer
    livecrootc_xfer                => clm3%g%c%p%pcs%livecrootc_xfer
    deadcrootc_xfer                => clm3%g%c%p%pcs%deadcrootc_xfer
    leafn_xfer                     => clm3%g%c%p%pns%leafn_xfer
    frootn_xfer                    => clm3%g%c%p%pns%frootn_xfer
    livestemn_xfer                 => clm3%g%c%p%pns%livestemn_xfer
    deadstemn_xfer                 => clm3%g%c%p%pns%deadstemn_xfer
    livecrootn_xfer                => clm3%g%c%p%pns%livecrootn_xfer
    deadcrootn_xfer                => clm3%g%c%p%pns%deadcrootn_xfer
    leafc_storage_to_xfer          => clm3%g%c%p%pcf%leafc_storage_to_xfer
    frootc_storage_to_xfer         => clm3%g%c%p%pcf%frootc_storage_to_xfer
    livestemc_storage_to_xfer      => clm3%g%c%p%pcf%livestemc_storage_to_xfer
    deadstemc_storage_to_xfer      => clm3%g%c%p%pcf%deadstemc_storage_to_xfer
    livecrootc_storage_to_xfer     => clm3%g%c%p%pcf%livecrootc_storage_to_xfer
    deadcrootc_storage_to_xfer     => clm3%g%c%p%pcf%deadcrootc_storage_to_xfer
    gresp_storage_to_xfer          => clm3%g%c%p%pcf%gresp_storage_to_xfer
    leafn_storage_to_xfer          => clm3%g%c%p%pnf%leafn_storage_to_xfer
    frootn_storage_to_xfer         => clm3%g%c%p%pnf%frootn_storage_to_xfer
    livestemn_storage_to_xfer      => clm3%g%c%p%pnf%livestemn_storage_to_xfer
    deadstemn_storage_to_xfer      => clm3%g%c%p%pnf%deadstemn_storage_to_xfer
    livecrootn_storage_to_xfer     => clm3%g%c%p%pnf%livecrootn_storage_to_xfer
    deadcrootn_storage_to_xfer     => clm3%g%c%p%pnf%deadcrootn_storage_to_xfer

   ! set time steps
   dayspyr = get_days_per_year(yr)

   do p = 1,num_soilp
      c = pcolumn(p)

      if (stress_decid(ivt(p)) == 1._r8) then
         soilt = t_soisno(c,3)
         psi = soilpsi(c,3)

         ! use solar declination information stored during Surface Albedo()
         ! and latitude from gps to calcluate daylength (convert latitude from degrees to radians)
         ! the constant 13750.9871 is the number of seconds per radian of hour-angle

         lat = (SHR_CONST_PI/180._r8)*latdeg(pgridcell(p))
         temp = -(sin(lat)*sin(decl(c)))/(cos(lat) * cos(decl(c)))
         temp = min(1._r8,max(-1._r8,temp))
         dayl(p) = 2.0_r8 * 13750.9871_r8 * acos(temp)

         ! onset gdd sum from Biome-BGC, v4.1.2
         crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_t2m(p) - SHR_CONST_TKFRZ))


         ! update offset_counter and test for the end of the offset period
         if (offset_flag(p) == 1._r8) then
            ! decrement counter for offset period
            offset_counter(p) = offset_counter(p) - dt

            ! if this is the end of the offset_period, reset phenology
            ! flags and indices
            if (offset_counter(p) == 0._r8) then
               ! this code block was originally handled by call cn_offset_cleanup(p)
               ! inlined during vectorization
               offset_flag(p) = 0._r8
               offset_counter(p) = 0._r8
               dormant_flag(p) = 1._r8
               days_active(p) = 0._r8

               ! reset the previous timestep litterfall flux memory
               prev_leafc_to_litter(p) = 0._r8
               prev_frootc_to_litter(p) = 0._r8
            end if
         end if

         ! update onset_counter and test for the end of the onset period
         if (onset_flag(p) == 1.0_r8) then
            ! decrement counter for onset period
            onset_counter(p) = onset_counter(p) - dt

            ! if this is the end of the onset period, reset phenology
            ! flags and indices
            if (onset_counter(p) == 0.0_r8) then
               ! this code block was originally handled by call cn_onset_cleanup(p)
               ! inlined during vectorization
               onset_flag(p) = 0._r8
               onset_counter(p) = 0._r8
               ! set all transfer growth rates to 0.0
               leafc_xfer_to_leafc(p)   = 0._r8
               frootc_xfer_to_frootc(p) = 0._r8
               leafn_xfer_to_leafn(p)   = 0._r8
               frootn_xfer_to_frootn(p) = 0._r8
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_xfer_to_livestemc(p)   = 0._r8
                  deadstemc_xfer_to_deadstemc(p)   = 0._r8
                  livecrootc_xfer_to_livecrootc(p) = 0._r8
                  deadcrootc_xfer_to_deadcrootc(p) = 0._r8
                  livestemn_xfer_to_livestemn(p)   = 0._r8
                  deadstemn_xfer_to_deadstemn(p)   = 0._r8
                  livecrootn_xfer_to_livecrootn(p) = 0._r8
                  deadcrootn_xfer_to_deadcrootn(p) = 0._r8
               end if
               ! set transfer pools to 0.0
               leafc_xfer(p) = 0._r8
               leafn_xfer(p) = 0._r8
               frootc_xfer(p) = 0._r8
               frootn_xfer(p) = 0._r8
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_xfer(p) = 0._r8
                  livestemn_xfer(p) = 0._r8
                  deadstemc_xfer(p) = 0._r8
                  deadstemn_xfer(p) = 0._r8
                  livecrootc_xfer(p) = 0._r8
                  livecrootn_xfer(p) = 0._r8
                  deadcrootc_xfer(p) = 0._r8
                  deadcrootn_xfer(p) = 0._r8
               end if
            end if
         end if

         ! test for switching from dormant period to growth period
         if (dormant_flag(p) == 1._r8) then

            ! keep track of the number of freezing degree days in this
            ! dormancy period (only if the freeze flag has not previously been set
            ! for this dormancy period

            if (onset_gddflag(p) == 0._r8 .and. soilt < SHR_CONST_TKFRZ) onset_fdd(p) = onset_fdd(p) + fracday

            ! if the number of freezing degree days exceeds a critical value,
            ! then onset will require both wet soils and a critical soil
            ! temperature sum.  If this case is triggered, reset any previously
            ! accumulated value in onset_swi, so that onset now depends on
            ! the accumulated soil water index following the freeze trigger

            if (onset_fdd(p) > crit_onset_fdd) then
                onset_gddflag(p) = 1._r8
                onset_fdd(p) = 0._r8
                onset_swi(p) = 0._r8
            end if

            ! if the freeze flag is set, and if the soil is above freezing
            ! then accumulate growing degree days for onset trigger

            if (onset_gddflag(p) == 1._r8 .and. soilt > SHR_CONST_TKFRZ) then
               onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
            end if

            ! if soils are wet, accumulate soil water index for onset trigger
            if (psi >= soilpsi_on) onset_swi(p) = onset_swi(p) + fracday

            ! if critical soil water index is exceeded, set onset_flag, and
            ! then test for soil temperature criteria

            if (onset_swi(p) > crit_onset_swi) then
                onset_flag(p) = 1._r8

                ! only check soil temperature criteria if freeze flag set since
                ! beginning of last dormancy.  If freeze flag set and growing
                ! degree day sum (since freeze trigger) is lower than critical
                ! value, then override the onset_flag set from soil water.

                if (onset_gddflag(p) == 1._r8 .and. onset_gdd(p) < crit_onset_gdd) onset_flag(p) = 0._r8
            end if
            
            ! only allow onset if dayl > 6hrs
            if (onset_flag(p) == 1._r8 .and. dayl(p) <= secspqtrday) then
                onset_flag(p) = 0._r8
            end if

            ! if this is the beginning of the onset period
            ! then reset the phenology flags and indices

            if (onset_flag(p) == 1._r8) then
               dormant_flag(p) = 0._r8
               days_active(p) = 0._r8
               onset_gddflag(p) = 0._r8
               onset_fdd(p) = 0._r8
               onset_gdd(p) = 0._r8
               onset_swi(p) = 0._r8
               onset_counter(p) = ndays_on * secspday

               ! call subroutine to move all the storage pools into transfer pools,
               ! where they will be transfered to displayed growth over the onset period.
               ! this code was originally handled with call cn_storage_to_xfer(p)
               ! inlined during vectorization

               ! set carbon fluxes for shifting storage pools to transfer pools
               leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
               frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
                  deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
                  livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
                  deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
                  gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
               end if

               ! set nitrogen fluxes for shifting storage pools to transfer pools
               leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
               frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
                  deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
                  livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
                  deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
               end if
            end if

         ! test for switching from growth period to offset period
         else if (offset_flag(p) == 0._r8) then

            ! if soil water potential lower than critical value, accumulate
            ! as stress in offset soil water index

            if (psi <= soilpsi_off) then
               offset_swi(p) = offset_swi(p) + fracday

               ! if the offset soil water index exceeds critical value, and
               ! if this is not the middle of a previously initiated onset period,
               ! then set flag to start the offset period and reset index variables

               if (offset_swi(p) >= crit_offset_swi .and. onset_flag(p) == 0._r8) offset_flag(p) = 1._r8

            ! if soil water potential higher than critical value, reduce the
            ! offset water stress index.  By this mechanism, there must be a
            ! sustained period of water stress to initiate offset.

            else if (psi >= soilpsi_on) then
               offset_swi(p) = offset_swi(p) - fracday
               offset_swi(p) = max(offset_swi(p),0._r8)
            end if

            ! decrease freezing day accumulator for warm soil
            if (offset_fdd(p) > 0._r8 .and. soilt > SHR_CONST_TKFRZ) then
                offset_fdd(p) = offset_fdd(p) - fracday
                offset_fdd(p) = max(0._r8, offset_fdd(p))
            end if

            ! increase freezing day accumulator for cold soil
            if (soilt <= SHR_CONST_TKFRZ) then
               offset_fdd(p) = offset_fdd(p) + fracday

               ! if freezing degree day sum is greater than critical value, initiate offset
               if (offset_fdd(p) > crit_offset_fdd .and. onset_flag(p) == 0._r8) offset_flag(p) = 1._r8
            end if
            
            ! force offset if daylength is < 6 hrs
            if (dayl(p) <= secspqtrday) then
               offset_flag(p) = 1._r8
            end if

            ! if this is the beginning of the offset period
            ! then reset flags and indices
            if (offset_flag(p) == 1._r8) then
               offset_fdd(p) = 0._r8
               offset_swi(p) = 0._r8
               offset_counter(p) = ndays_off * secspday
               prev_leafc_to_litter(p) = 0._r8
               prev_frootc_to_litter(p) = 0._r8
            end if
         end if

         ! keep track of number of days since last dormancy for control on
         ! fraction of new growth to send to storage for next growing season

         if (dormant_flag(p) == 0.0_r8) then
             days_active(p) = days_active(p) + fracday
         end if

         ! calculate long growing season factor (lgsf)
         ! only begin to calculate a lgsf greater than 0.0 once the number
         ! of days active exceeds days/year.
         lgsf(p) = max(min((days_active(p)-dayspyr)/dayspyr, 1._r8),0._r8)

         ! set background litterfall rate, when not in the phenological offset period
         if (offset_flag(p) == 1._r8) then
            bglfr(p) = 0._r8
         else
            ! calculate the background litterfall rate (bglfr)
            ! in units 1/s, based on leaf longevity (yrs) and correction for long growing season

            bglfr(p) = (1._r8/(leaf_long(ivt(p))*dayspyr*secspday))*lgsf(p)
         end if

         ! set background transfer rate when active but not in the phenological onset period
         if (onset_flag(p) == 1._r8) then
            bgtr(p) = 0._r8
         else
            ! the background transfer rate is calculated as the rate that would result
            ! in complete turnover of the storage pools in one year at steady state,
            ! once lgsf has reached 1.0 (after 730 days active).

            bgtr(p) = (1._r8/(dayspyr*secspday))*lgsf(p)

            ! set carbon fluxes for shifting storage pools to transfer pools

            leafc_storage_to_xfer(p)  = leafc_storage(p) * bgtr(p)
            frootc_storage_to_xfer(p) = frootc_storage(p) * bgtr(p)
            if (woody(ivt(p)) == 1.0_r8) then
               livestemc_storage_to_xfer(p)  = livestemc_storage(p) * bgtr(p)
               deadstemc_storage_to_xfer(p)  = deadstemc_storage(p) * bgtr(p)
               livecrootc_storage_to_xfer(p) = livecrootc_storage(p) * bgtr(p)
               deadcrootc_storage_to_xfer(p) = deadcrootc_storage(p) * bgtr(p)
               gresp_storage_to_xfer(p)      = gresp_storage(p) * bgtr(p)
            end if

            ! set nitrogen fluxes for shifting storage pools to transfer pools
            leafn_storage_to_xfer(p)  = leafn_storage(p) * bgtr(p)
            frootn_storage_to_xfer(p) = frootn_storage(p) * bgtr(p)
            if (woody(ivt(p)) == 1.0_r8) then
               livestemn_storage_to_xfer(p)  = livestemn_storage(p) * bgtr(p)
               deadstemn_storage_to_xfer(p)  = deadstemn_storage(p) * bgtr(p)
               livecrootn_storage_to_xfer(p) = livecrootn_storage(p) * bgtr(p)
               deadcrootn_storage_to_xfer(p) = deadcrootn_storage(p) * bgtr(p)
            end if
         end if

      end if ! end if stress deciduous

   end do ! end of pft loop

end subroutine CNStressDecidPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNOnsetGrowth
!
! !INTERFACE:
subroutine CNOnsetGrowth (dt, num_soilp)
!
! !DESCRIPTION:
! Determines the flux of stored C and N from transfer pools to display
! pools during the phenological onset period.
!
! !USES:
    use pftvarcon , only: npcropmin
!!
! !ARGUMENTS:
   real(r8), intent(in) :: dt             ! timestep (s)
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/27/03: Created by Peter Thornton
! 10/10/12: Adapted for use in VIC by Michael Brunke
! 02/07/14: Added ignoring crops in PFT loops (Michael Brunke).
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)             ! pft vegetation type
   real(r8), pointer :: onset_flag(:)      ! onset flag
   real(r8), pointer :: onset_counter(:)   ! onset days counter
   real(r8), pointer :: leafc_xfer(:)      ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc_xfer(:)     ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc_xfer(:)  ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc_xfer(:)  ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc_xfer(:) ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc_xfer(:) ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: leafn_xfer(:)      ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn_xfer(:)     ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn_xfer(:)  ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn_xfer(:)  ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn_xfer(:) ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn_xfer(:) ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: woody(:)           ! binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: bgtr(:)            ! background transfer growth rate (1/s)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: leafc_xfer_to_leafc(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: leafn_xfer_to_leafn(:)
   real(r8), pointer :: frootn_xfer_to_frootn(:)
   real(r8), pointer :: livestemn_xfer_to_livestemn(:)
   real(r8), pointer :: deadstemn_xfer_to_deadstemn(:)
   real(r8), pointer :: livecrootn_xfer_to_livecrootn(:)
   real(r8), pointer :: deadcrootn_xfer_to_deadcrootn(:)
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: p            ! indices
   real(r8):: t1           ! temporary variable

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%c%p%itype
    onset_flag                     => clm3%g%c%p%pepv%onset_flag
    onset_counter                  => clm3%g%c%p%pepv%onset_counter
    leafc_xfer                     => clm3%g%c%p%pcs%leafc_xfer
    frootc_xfer                    => clm3%g%c%p%pcs%frootc_xfer
    livestemc_xfer                 => clm3%g%c%p%pcs%livestemc_xfer
    deadstemc_xfer                 => clm3%g%c%p%pcs%deadstemc_xfer
    livecrootc_xfer                => clm3%g%c%p%pcs%livecrootc_xfer
    deadcrootc_xfer                => clm3%g%c%p%pcs%deadcrootc_xfer
    leafn_xfer                     => clm3%g%c%p%pns%leafn_xfer
    frootn_xfer                    => clm3%g%c%p%pns%frootn_xfer
    livestemn_xfer                 => clm3%g%c%p%pns%livestemn_xfer
    deadstemn_xfer                 => clm3%g%c%p%pns%deadstemn_xfer
    livecrootn_xfer                => clm3%g%c%p%pns%livecrootn_xfer
    deadcrootn_xfer                => clm3%g%c%p%pns%deadcrootn_xfer
    bgtr                           => clm3%g%c%p%pepv%bgtr
    woody                          => pftcon%woody

   ! assign local pointers to derived type arrays (out)
    leafc_xfer_to_leafc            => clm3%g%c%p%pcf%leafc_xfer_to_leafc
    frootc_xfer_to_frootc          => clm3%g%c%p%pcf%frootc_xfer_to_frootc
    livestemc_xfer_to_livestemc    => clm3%g%c%p%pcf%livestemc_xfer_to_livestemc
    deadstemc_xfer_to_deadstemc    => clm3%g%c%p%pcf%deadstemc_xfer_to_deadstemc
    livecrootc_xfer_to_livecrootc  => clm3%g%c%p%pcf%livecrootc_xfer_to_livecrootc
    deadcrootc_xfer_to_deadcrootc  => clm3%g%c%p%pcf%deadcrootc_xfer_to_deadcrootc
    leafn_xfer_to_leafn            => clm3%g%c%p%pnf%leafn_xfer_to_leafn
    frootn_xfer_to_frootn          => clm3%g%c%p%pnf%frootn_xfer_to_frootn
    livestemn_xfer_to_livestemn    => clm3%g%c%p%pnf%livestemn_xfer_to_livestemn
    deadstemn_xfer_to_deadstemn    => clm3%g%c%p%pnf%deadstemn_xfer_to_deadstemn
    livecrootn_xfer_to_livecrootn  => clm3%g%c%p%pnf%livecrootn_xfer_to_livecrootn
    deadcrootn_xfer_to_deadcrootn  => clm3%g%c%p%pnf%deadcrootn_xfer_to_deadcrootn

   ! pft loop
   do p = 1,num_soilp
     if(ivt(p) < npcropmin) then

      ! only calculate these fluxes during onset period
      if (onset_flag(p) == 1._r8) then

         ! The transfer rate is a linearly decreasing function of time,
         ! going to zero on the last timestep of the onset period

         if (onset_counter(p) == dt) then
             t1 = 1.0_r8 / dt
         else
             t1 = 2.0_r8 / (onset_counter(p))
         end if
         leafc_xfer_to_leafc(p)   = t1 * leafc_xfer(p)
         frootc_xfer_to_frootc(p) = t1 * frootc_xfer(p)
         leafn_xfer_to_leafn(p)   = t1 * leafn_xfer(p)
         frootn_xfer_to_frootn(p) = t1 * frootn_xfer(p)
         if (woody(ivt(p)) == 1.0_r8) then
             livestemc_xfer_to_livestemc(p)   = t1 * livestemc_xfer(p)
             deadstemc_xfer_to_deadstemc(p)   = t1 * deadstemc_xfer(p)
             livecrootc_xfer_to_livecrootc(p) = t1 * livecrootc_xfer(p)
             deadcrootc_xfer_to_deadcrootc(p) = t1 * deadcrootc_xfer(p)
             livestemn_xfer_to_livestemn(p)   = t1 * livestemn_xfer(p)
             deadstemn_xfer_to_deadstemn(p)   = t1 * deadstemn_xfer(p)
             livecrootn_xfer_to_livecrootn(p) = t1 * livecrootn_xfer(p)
             deadcrootn_xfer_to_deadcrootn(p) = t1 * deadcrootn_xfer(p)
         end if

      end if ! end if onset period

      ! calculate the background rate of transfer growth (used for stress
      ! deciduous algorithm). In this case, all of the mass in the transfer
      ! pools should be moved to displayed growth in each timestep.

      if (bgtr(p) > 0._r8) then
         leafc_xfer_to_leafc(p)   = leafc_xfer(p) / dt
         frootc_xfer_to_frootc(p) = frootc_xfer(p) / dt
         leafn_xfer_to_leafn(p)   = leafn_xfer(p) / dt
         frootn_xfer_to_frootn(p) = frootn_xfer(p) / dt
         if (woody(ivt(p)) == 1.0_r8) then
             livestemc_xfer_to_livestemc(p)   = livestemc_xfer(p) / dt
             deadstemc_xfer_to_deadstemc(p)   = deadstemc_xfer(p) / dt
             livecrootc_xfer_to_livecrootc(p) = livecrootc_xfer(p) / dt
             deadcrootc_xfer_to_deadcrootc(p) = deadcrootc_xfer(p) / dt
             livestemn_xfer_to_livestemn(p)   = livestemn_xfer(p) / dt
             deadstemn_xfer_to_deadstemn(p)   = deadstemn_xfer(p) / dt
             livecrootn_xfer_to_livecrootn(p) = livecrootn_xfer(p) / dt
             deadcrootn_xfer_to_deadcrootn(p) = deadcrootn_xfer(p) / dt
         end if
      end if ! end if bgtr

     end if
   end do ! end pft loop

end subroutine CNOnsetGrowth
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNOffsetLitterfall
!
! !INTERFACE:
subroutine CNOffsetLitterfall (dt, num_soilp)
!
! !DESCRIPTION:
! Determines the flux of C and N from displayed pools to litter
! pools during the phenological offset period.
!
! !USES:
   use pftvarcon       , only: npcropmin
!
! !ARGUMENTS:
   real(r8), intent(in) :: dt             ! timestep (s)
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/27/03: Created by Peter Thornton
! 10/10/12: Adapted for use in VIC by Michael Brunke
! 02/07/14: Added ignoring crops in PFT loops (Michael Brunke).
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)                   ! pft vegetation type
   real(r8), pointer :: offset_flag(:)           ! offset flag
   real(r8), pointer :: offset_counter(:)        ! offset days counter
   real(r8), pointer :: leafc(:)                 ! (gC/m2) leaf C
   real(r8), pointer :: frootc(:)                ! (gC/m2) fine root C
   real(r8), pointer :: cpool_to_leafc(:)        ! allocation to leaf C (gC/m2/s)
   real(r8), pointer :: cpool_to_frootc(:)       ! allocation to fine root C (gC/m2/s)
!  integer , pointer :: pcolumn(:)               ! pft's column index
   real(r8), pointer :: grainc(:)                ! (gC/m2) grain C
   real(r8), pointer :: livestemc(:)             ! (gC/m2) livestem C
   real(r8), pointer :: cpool_to_grainc(:)       ! allocation to grain C (gC/m2/s)
   real(r8), pointer :: cpool_to_livestemc(:)    ! allocation to live stem C (gC/m2/s)
   real(r8), pointer :: livewdcn(:)              ! live wood C:N (gC/gN)
   real(r8), pointer :: graincn(:)               ! grain C:N (gC/gN)
   real(r8), pointer :: leafcn(:)                ! leaf C:N (gC/gN)
   real(r8), pointer :: lflitcn(:)               ! leaf litter C:N (gC/gN)
   real(r8), pointer :: frootcn(:)               ! fine root C:N (gC/gN)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: prev_leafc_to_litter(:)  ! previous timestep leaf C litterfall flux (gC/m2/s)
   real(r8), pointer :: prev_frootc_to_litter(:) ! previous timestep froot C litterfall flux (gC/m2/s)
   real(r8), pointer :: leafc_to_litter(:)       ! leaf C litterfall (gC/m2/s)
   real(r8), pointer :: frootc_to_litter(:)      ! fine root C litterfall (gC/m2/s)
   real(r8), pointer :: leafn_to_litter(:)       ! leaf N litterfall (gN/m2/s)
   real(r8), pointer :: leafn_to_retransn(:)     ! leaf N to retranslocated N pool (gN/m2/s)
   real(r8), pointer :: frootn_to_litter(:)      ! fine root N litterfall (gN/m2/s)
   real(r8), pointer :: livestemc_to_litter(:)   ! live stem C litterfall (gC/m2/s)
   real(r8), pointer :: grainc_to_food(:)        ! grain C to food (gC/m2/s)
   real(r8), pointer :: livestemn_to_litter(:)   ! livestem N to litter (gN/m2/s)
   real(r8), pointer :: grainn_to_food(:)        ! grain N to food (gN/m2/s)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p, c         ! indices
   real(r8):: t1           ! temporary variable

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%c%p%itype
    offset_flag                    => clm3%g%c%p%pepv%offset_flag
    offset_counter                 => clm3%g%c%p%pepv%offset_counter
    leafc                          => clm3%g%c%p%pcs%leafc
    frootc                         => clm3%g%c%p%pcs%frootc
    grainc                         => clm3%g%c%p%pcs%grainc
    livestemc                      => clm3%g%c%p%pcs%livestemc
    cpool_to_grainc                => clm3%g%c%p%pcf%cpool_to_grainc
    cpool_to_livestemc             => clm3%g%c%p%pcf%cpool_to_livestemc
    cpool_to_leafc                 => clm3%g%c%p%pcf%cpool_to_leafc
    cpool_to_frootc                => clm3%g%c%p%pcf%cpool_to_frootc
    leafcn                         => pftcon%leafcn
    lflitcn                        => pftcon%lflitcn
    frootcn                        => pftcon%frootcn
    livewdcn                       => pftcon%livewdcn
    graincn                        => pftcon%graincn

   ! assign local pointers to derived type arrays (out)
    prev_leafc_to_litter           => clm3%g%c%p%pepv%prev_leafc_to_litter
    prev_frootc_to_litter          => clm3%g%c%p%pepv%prev_frootc_to_litter
    leafc_to_litter                => clm3%g%c%p%pcf%leafc_to_litter
    frootc_to_litter               => clm3%g%c%p%pcf%frootc_to_litter
    livestemc_to_litter            => clm3%g%c%p%pcf%livestemc_to_litter
    grainc_to_food                 => clm3%g%c%p%pcf%grainc_to_food
    livestemn_to_litter            => clm3%g%c%p%pnf%livestemn_to_litter
    grainn_to_food                 => clm3%g%c%p%pnf%grainn_to_food
    leafn_to_litter                => clm3%g%c%p%pnf%leafn_to_litter
    leafn_to_retransn              => clm3%g%c%p%pnf%leafn_to_retransn
    frootn_to_litter               => clm3%g%c%p%pnf%frootn_to_litter

   ! The litterfall transfer rate starts at 0.0 and increases linearly
   ! over time, with displayed growth going to 0.0 on the last day of litterfall

   do p = 1,num_soilp
     if(ivt(p) < npcropmin) then

      ! only calculate fluxes during offset period
      if (offset_flag(p) == 1._r8) then

         if (offset_counter(p) == dt) then
             t1 = 1.0_r8 / dt
             leafc_to_litter(p)  = t1 * leafc(p)  + cpool_to_leafc(p)
             frootc_to_litter(p) = t1 * frootc(p) + cpool_to_frootc(p)
             ! this assumes that offset_counter == dt for crops
             ! if this were ever changed, we'd need to add code to the "else"
!             if (ivt(p) >= npcropmin) then
!                grainc_to_food(p) = t1 * grainc(p)  + cpool_to_grainc(p) 
!                livestemc_to_litter(p) = t1 * livestemc(p)  + cpool_to_livestemc(p)
!             end if
         else
             t1 = dt * 2.0_r8 / (offset_counter(p) * offset_counter(p))
             leafc_to_litter(p)  = prev_leafc_to_litter(p)  + t1*(leafc(p)  - prev_leafc_to_litter(p)*offset_counter(p))
             frootc_to_litter(p) = prev_frootc_to_litter(p) + t1*(frootc(p) - prev_frootc_to_litter(p)*offset_counter(p))
         end if

         ! calculate the leaf N litterfall and retranslocation
         leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
         leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

         ! calculate fine root N litterfall (no retranslocation of fine root N)
         frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))

!         if (ivt(p) >= npcropmin) then
!            livestemn_to_litter(p) = livestemc_to_litter(p) / livewdcn(ivt(p))
!            grainn_to_food(p) = grainc_to_food(p) / graincn(ivt(p))
!         end if

         ! save the current litterfall fluxes
         prev_leafc_to_litter(p)  = leafc_to_litter(p)
         prev_frootc_to_litter(p) = frootc_to_litter(p)

      end if ! end if offset period

     end if
   end do ! end pft loop

end subroutine CNOffsetLitterfall
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNBackgroundLitterfall
!
! !INTERFACE:
subroutine CNBackgroundLitterfall (num_soilp)
!
! !DESCRIPTION:
! Determines the flux of C and N from displayed pools to litter
! pools as the result of background litter fall.
!
! !USES:
   use pftvarcon       , only: npcropmin
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/2/03: Created by Peter Thornton
! 10/24/03, Peter Thornton: migrated to vector data structures
! 10/10/12: Adapted for use in VIC by Michael Brunke
! 02/07/14: Added ignoring crops in PFT loops (Michael Brunke).
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! pft level
   integer , pointer :: ivt(:)       ! pft vegetation type
   real(r8), pointer :: bglfr(:)     ! background litterfall rate (1/s)
   real(r8), pointer :: leafc(:)     ! (gC/m2) leaf C
   real(r8), pointer :: frootc(:)    ! (gC/m2) fine root C
   ! ecophysiological constants
   real(r8), pointer :: leafcn(:)    ! leaf C:N (gC/gN)
   real(r8), pointer :: lflitcn(:)   ! leaf litter C:N (gC/gN)
   real(r8), pointer :: frootcn(:)   ! fine root C:N (gC/gN)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: leafc_to_litter(:)
   real(r8), pointer :: frootc_to_litter(:)
   real(r8), pointer :: leafn_to_litter(:)
   real(r8), pointer :: leafn_to_retransn(:)
   real(r8), pointer :: frootn_to_litter(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p            ! indices

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%c%p%itype
    bglfr                          => clm3%g%c%p%pepv%bglfr
    leafc                          => clm3%g%c%p%pcs%leafc
    frootc                         => clm3%g%c%p%pcs%frootc
    leafcn                         => pftcon%leafcn
    lflitcn                        => pftcon%lflitcn
    frootcn                        => pftcon%frootcn

   ! assign local pointers to derived type arrays (out)
    leafc_to_litter                => clm3%g%c%p%pcf%leafc_to_litter
    frootc_to_litter               => clm3%g%c%p%pcf%frootc_to_litter
    leafn_to_litter                => clm3%g%c%p%pnf%leafn_to_litter
    leafn_to_retransn              => clm3%g%c%p%pnf%leafn_to_retransn
    frootn_to_litter               => clm3%g%c%p%pnf%frootn_to_litter

   ! pft loop
   do p = 1,num_soilp
     if(ivt(p) < npcropmin) then

      ! only calculate these fluxes if the background litterfall rate is non-zero
      if (bglfr(p) > 0._r8) then
         ! units for bglfr are already 1/s
         leafc_to_litter(p)  = bglfr(p) * leafc(p)
         frootc_to_litter(p) = bglfr(p) * frootc(p)

         ! calculate the leaf N litterfall and retranslocation
         leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
         leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

         ! calculate fine root N litterfall (no retranslocation of fine root N)
         frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))

      end if

     end if
   end do

end subroutine CNBackgroundLitterfall
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNLivewoodTurnover
!
! !INTERFACE:
subroutine CNLivewoodTurnover (num_soilp)
!
! !DESCRIPTION:
! Determines the flux of C and N from live wood to
! dead wood pools, for stem and coarse root.
!
! !USES:
   use pftvarcon       , only: npcropmin
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 12/5/03: created by Peter Thornton
! 10/10/12: Adapted for use in VIC by Michael Brunke
! 02/07/14: Added ignoring crops in PFT loops (Michael Brunke).
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! pft level
   integer , pointer :: ivt(:)         ! pft vegetation type
   real(r8), pointer :: livestemc(:)   ! (gC/m2) live stem C
   real(r8), pointer :: livecrootc(:)  ! (gC/m2) live coarse root C
   real(r8), pointer :: livestemn(:)   ! (gN/m2) live stem N
   real(r8), pointer :: livecrootn(:)  ! (gN/m2) live coarse root N
   ! ecophysiological constants
   real(r8), pointer :: woody(:)       ! binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: livewdcn(:)    ! live wood (phloem and ray parenchyma) C:N (gC/gN)
   real(r8), pointer :: deadwdcn(:)    ! dead wood (xylem and heartwood) C:N (gC/gN)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: livestemc_to_deadstemc(:)
   real(r8), pointer :: livecrootc_to_deadcrootc(:)
   real(r8), pointer :: livestemn_to_deadstemn(:)
   real(r8), pointer :: livestemn_to_retransn(:)
   real(r8), pointer :: livecrootn_to_deadcrootn(:)
   real(r8), pointer :: livecrootn_to_retransn(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p            ! indices
   real(r8):: ctovr        ! temporary variable for carbon turnover
   real(r8):: ntovr        ! temporary variable for nitrogen turnover

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%c%p%itype
    livestemc                      => clm3%g%c%p%pcs%livestemc
    livecrootc                     => clm3%g%c%p%pcs%livecrootc
    livestemn                      => clm3%g%c%p%pns%livestemn
    livecrootn                     => clm3%g%c%p%pns%livecrootn
    woody                          => pftcon%woody
    livewdcn                       => pftcon%livewdcn
    deadwdcn                       => pftcon%deadwdcn

   ! assign local pointers to derived type arrays (out)
    livestemc_to_deadstemc         => clm3%g%c%p%pcf%livestemc_to_deadstemc
    livecrootc_to_deadcrootc       => clm3%g%c%p%pcf%livecrootc_to_deadcrootc
    livestemn_to_deadstemn         => clm3%g%c%p%pnf%livestemn_to_deadstemn
    livestemn_to_retransn          => clm3%g%c%p%pnf%livestemn_to_retransn
    livecrootn_to_deadcrootn       => clm3%g%c%p%pnf%livecrootn_to_deadcrootn
    livecrootn_to_retransn         => clm3%g%c%p%pnf%livecrootn_to_retransn

   ! pft loop
   do p = 1,num_soilp
     if(ivt(p) < npcropmin) then

      ! only calculate these fluxes for woody types
      if (woody(ivt(p)) > 0._r8) then

         ! live stem to dead stem turnover

         ctovr = livestemc(p) * lwtop
         ntovr = ctovr / livewdcn(ivt(p))
         livestemc_to_deadstemc(p) = ctovr
         livestemn_to_deadstemn(p) = ctovr / deadwdcn(ivt(p))
         livestemn_to_retransn(p)  = ntovr - livestemn_to_deadstemn(p)

         ! live coarse root to dead coarse root turnover

         ctovr = livecrootc(p) * lwtop
         ntovr = ctovr / livewdcn(ivt(p))
         livecrootc_to_deadcrootc(p) = ctovr
         livecrootn_to_deadcrootn(p) = ctovr / deadwdcn(ivt(p))
         livecrootn_to_retransn(p)  = ntovr - livecrootn_to_deadcrootn(p)

      end if

     end if
   end do

end subroutine CNLivewoodTurnover
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNLitterToColumn
!
! !INTERFACE:
subroutine CNLitterToColumn (num_soilc)
!
! !DESCRIPTION:
! called at the end of cn_phenology to gather all pft-level litterfall fluxes
! to the column level and assign them to the three litter pools
!
! !USES:
  use clm_varpar, only : max_pft_per_col
  use pftvarcon , only : npcropmin
!
! !ARGUMENTS:
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
! 10/10/12: Adapted for use in VIC by Michael Brunke
! 02/07/14: Added ignoring crops in PFT loops (Michael Brunke).
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)          ! pft vegetation type
   real(r8), pointer :: wtcol(:)        ! weight (relative to column) for this pft (0-1)
   real(r8), pointer :: pwtgcell(:)     ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: leafc_to_litter(:)     ! leaf C litterfall (gC/m2/s)
   real(r8), pointer :: frootc_to_litter(:)    ! fine root N litterfall (gN/m2/s)
   real(r8), pointer :: livestemc_to_litter(:) ! live stem C litterfall (gC/m2/s)
   real(r8), pointer :: grainc_to_food(:)      ! grain C to food (gC/m2/s)
   real(r8), pointer :: livestemn_to_litter(:) ! livestem N to litter (gN/m2/s)
   real(r8), pointer :: grainn_to_food(:)      ! grain N to food (gN/m2/s)
   real(r8), pointer :: leafn_to_litter(:)     ! leaf N litterfall (gN/m2/s)
   real(r8), pointer :: frootn_to_litter(:)    ! fine root N litterfall (gN/m2/s)
   real(r8), pointer :: lf_flab(:)      ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)      ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)      ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)      ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)      ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)      ! fine root litter lignin fraction
   integer , pointer :: npfts(:)        ! number of pfts for each column
   integer , pointer :: pfti(:)         ! beginning pft index for each column
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: leafc_to_litr1c(:)     ! leaf C litterfall to litter 1 C (gC/m2/s)
   real(r8), pointer :: leafc_to_litr2c(:)     ! leaf C litterfall to litter 2 C (gC/m2/s)
   real(r8), pointer :: leafc_to_litr3c(:)     ! leaf C litterfall to litter 3 C (gC/m2/s)
   real(r8), pointer :: frootc_to_litr1c(:)    ! fine root C litterfall to litter 1 C (gC/m2/s)
   real(r8), pointer :: frootc_to_litr2c(:)    ! fine root C litterfall to litter 2 C (gC/m2/s)
   real(r8), pointer :: frootc_to_litr3c(:)    ! fine root C litterfall to litter 3 C (gC/m2/s)
   real(r8), pointer :: livestemc_to_litr1c(:) ! livestem C litterfall to litter 1 C (gC/m2/s)
   real(r8), pointer :: livestemc_to_litr2c(:) ! livestem C litterfall to litter 2 C (gC/m2/s)
   real(r8), pointer :: livestemc_to_litr3c(:) ! livestem C litterfall to litter 3 C (gC/m2/s)
   real(r8), pointer :: livestemn_to_litr1n(:) ! livestem N litterfall to litter 1 N (gN/m2/s)
   real(r8), pointer :: livestemn_to_litr2n(:) ! livestem N litterfall to litter 2 N (gN/m2/s)
   real(r8), pointer :: livestemn_to_litr3n(:) ! livestem N litterfall to litter 2 N (gN/m2/s)
   real(r8), pointer :: grainc_to_litr1c(:)    ! grain C litterfall to litter 1 C (gC/m2/s)
   real(r8), pointer :: grainc_to_litr2c(:)    ! grain C litterfall to litter 2 C (gC/m2/s)
   real(r8), pointer :: grainc_to_litr3c(:)    ! grain C litterfall to litter 3 C (gC/m2/s)
   real(r8), pointer :: grainn_to_litr1n(:)    ! grain N litterfall to litter 1 N (gN/m2/s)
   real(r8), pointer :: grainn_to_litr2n(:)    ! grain N litterfall to litter 2 N (gN/m2/s)
   real(r8), pointer :: grainn_to_litr3n(:)    ! grain N litterfall to litter 3 N (gN/m2/s)
   real(r8), pointer :: leafn_to_litr1n(:)     ! leaf N litterfall to litter 1 N (gN/m2/s)
   real(r8), pointer :: leafn_to_litr2n(:)     ! leaf N litterfall to litter 2 N (gN/m2/s)
   real(r8), pointer :: leafn_to_litr3n(:)     ! leaf N litterfall to litter 3 N (gN/m2/s)
   real(r8), pointer :: frootn_to_litr1n(:)    ! fine root N litterfall to litter 1 N (gN/m2/s)
   real(r8), pointer :: frootn_to_litr2n(:)    ! fine root N litterfall to litter 2 N (gN/m2/s)
   real(r8), pointer :: frootn_to_litr3n(:)    ! fine root N litterfall to litter 3 N (gN/m2/s)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
    integer :: c,pi,p        ! indices
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%c%p%itype
    wtcol                          => clm3%g%c%p%wtcol
    pwtgcell                       => clm3%g%c%p%wtgcell  
    leafc_to_litter                => clm3%g%c%p%pcf%leafc_to_litter
    frootc_to_litter               => clm3%g%c%p%pcf%frootc_to_litter
    livestemc_to_litter            => clm3%g%c%p%pcf%livestemc_to_litter
    grainc_to_food                 => clm3%g%c%p%pcf%grainc_to_food
    livestemn_to_litter            => clm3%g%c%p%pnf%livestemn_to_litter
    grainn_to_food                 => clm3%g%c%p%pnf%grainn_to_food
    leafn_to_litter                => clm3%g%c%p%pnf%leafn_to_litter
    frootn_to_litter               => clm3%g%c%p%pnf%frootn_to_litter
    npfts                          => clm3%g%c%npfts
    pfti                           => clm3%g%c%pfti
    lf_flab                        => pftcon%lf_flab
    lf_fcel                        => pftcon%lf_fcel
    lf_flig                        => pftcon%lf_flig
    fr_flab                        => pftcon%fr_flab
    fr_fcel                        => pftcon%fr_fcel
    fr_flig                        => pftcon%fr_flig

   ! assign local pointers to derived type arrays (out)
    leafc_to_litr1c                => clm3%g%c%ccf%leafc_to_litr1c
    leafc_to_litr2c                => clm3%g%c%ccf%leafc_to_litr2c
    leafc_to_litr3c                => clm3%g%c%ccf%leafc_to_litr3c
    frootc_to_litr1c               => clm3%g%c%ccf%frootc_to_litr1c
    frootc_to_litr2c               => clm3%g%c%ccf%frootc_to_litr2c
    frootc_to_litr3c               => clm3%g%c%ccf%frootc_to_litr3c
    grainc_to_litr1c               => clm3%g%c%ccf%grainc_to_litr1c
    grainc_to_litr2c               => clm3%g%c%ccf%grainc_to_litr2c
    grainc_to_litr3c               => clm3%g%c%ccf%grainc_to_litr3c
    livestemc_to_litr1c            => clm3%g%c%ccf%livestemc_to_litr1c
    livestemc_to_litr2c            => clm3%g%c%ccf%livestemc_to_litr2c
    livestemc_to_litr3c            => clm3%g%c%ccf%livestemc_to_litr3c
    livestemn_to_litr1n            => clm3%g%c%cnf%livestemn_to_litr1n
    livestemn_to_litr2n            => clm3%g%c%cnf%livestemn_to_litr2n
    livestemn_to_litr3n            => clm3%g%c%cnf%livestemn_to_litr3n
    grainn_to_litr1n               => clm3%g%c%cnf%grainn_to_litr1n
    grainn_to_litr2n               => clm3%g%c%cnf%grainn_to_litr2n
    grainn_to_litr3n               => clm3%g%c%cnf%grainn_to_litr3n
    leafn_to_litr1n                => clm3%g%c%cnf%leafn_to_litr1n
    leafn_to_litr2n                => clm3%g%c%cnf%leafn_to_litr2n
    leafn_to_litr3n                => clm3%g%c%cnf%leafn_to_litr3n
    frootn_to_litr1n               => clm3%g%c%cnf%frootn_to_litr1n
    frootn_to_litr2n               => clm3%g%c%cnf%frootn_to_litr2n
    frootn_to_litr3n               => clm3%g%c%cnf%frootn_to_litr3n

   do pi = 1,max_pft_per_col
      do c = 1,num_soilc

         if ( pi <=  npfts(c) ) then
            p = pfti(c) + pi - 1
            if (pwtgcell(p)>0._r8 .and. ivt(p) < npcropmin) then

               ! leaf litter carbon fluxes
               leafc_to_litr1c(c) = leafc_to_litr1c(c) + leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               leafc_to_litr2c(c) = leafc_to_litr2c(c) + leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               leafc_to_litr3c(c) = leafc_to_litr3c(c) + leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! leaf litter nitrogen fluxes
               leafn_to_litr1n(c) = leafn_to_litr1n(c) + leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               leafn_to_litr2n(c) = leafn_to_litr2n(c) + leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               leafn_to_litr3n(c) = leafn_to_litr3n(c) + leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root litter carbon fluxes
               frootc_to_litr1c(c) = frootc_to_litr1c(c) + frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               frootc_to_litr2c(c) = frootc_to_litr2c(c) + frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               frootc_to_litr3c(c) = frootc_to_litr3c(c) + frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! fine root litter nitrogen fluxes
               frootn_to_litr1n(c) = frootn_to_litr1n(c) + frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               frootn_to_litr2n(c) = frootn_to_litr2n(c) + frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               frootn_to_litr3n(c) = frootn_to_litr3n(c) + frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)


               ! agroibis puts crop stem litter together with leaf litter
               ! so I've used the leaf lf_f* parameters instead of making
               ! new ones for now (slevis)
               ! also for simplicity I've put "food" into the litter pools
!               if (ivt(p) >= npcropmin) then ! add livestemc to litter
                  ! stem litter carbon fluxes
!                  livestemc_to_litr1c(c) = livestemc_to_litr1c(c) + livestemc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
!                  livestemc_to_litr2c(c) = livestemc_to_litr2c(c) + livestemc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
!                  livestemc_to_litr3c(c) = livestemc_to_litr3c(c) + livestemc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

                  ! stem litter nitrogen fluxes
!                  livestemn_to_litr1n(c) = livestemn_to_litr1n(c) + livestemn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
!                  livestemn_to_litr2n(c) = livestemn_to_litr2n(c) + livestemn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
!                  livestemn_to_litr3n(c) = livestemn_to_litr3n(c) + livestemn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

                  ! grain litter carbon fluxes
!                  grainc_to_litr1c(c) = grainc_to_litr1c(c) + grainc_to_food(p) * lf_flab(ivt(p)) * wtcol(p)
!                  grainc_to_litr2c(c) = grainc_to_litr2c(c) + grainc_to_food(p) * lf_fcel(ivt(p)) * wtcol(p)
!                  grainc_to_litr3c(c) = grainc_to_litr3c(c) + grainc_to_food(p) * lf_flig(ivt(p)) * wtcol(p)

                  ! grain litter nitrogen fluxes
!                  grainn_to_litr1n(c) = grainn_to_litr1n(c) + grainn_to_food(p) * lf_flab(ivt(p)) * wtcol(p)
!                  grainn_to_litr2n(c) = grainn_to_litr2n(c) + grainn_to_food(p) * lf_fcel(ivt(p)) * wtcol(p)
!                  grainn_to_litr3n(c) = grainn_to_litr3n(c) + grainn_to_food(p) * lf_flig(ivt(p)) * wtcol(p)
!               end if
            end if
         end if

      end do

   end do

end subroutine CNLitterToColumn
!-----------------------------------------------------------------------

end module CNPhenologyMod
