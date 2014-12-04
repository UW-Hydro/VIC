!#include <misc.h>
!#include <preproc.h>

module CNiniTimeVarMod

  implicit none
  save
  private
  public :: CNiniTimeVar

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNiniTimeVar
!
! !INTERFACE:
subroutine CNiniTimeVar(begg, endg, begc, endc, begp, endp)

!#ifdef CN
!
! !DESCRIPTION:
! Initializes time varying variables used only in
! coupled carbon-nitrogen mode (CN):
!
! !USES:
   use clmtype
!   use clm_atmlnd  , only: clm_a2l
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pftvarcon   , only: noveg
!   use decompMod   , only: get_proc_bounds
!
! !ARGUMENTS:
   implicit none
   integer :: begp, endp   ! per-clump/proc beginning and ending pft indices
   integer :: begc, endc   ! per-clump/proc beginning and ending column indices
   integer :: begg, endg   ! per-clump/proc gridcell ending gridcell indices
!
! !CALLED FROM:
! subroutine clm_initialize2 in file clm_initializeMod.F90
!
! !REVISION HISTORY:
! 10/21/03: Created by Peter Thornton
! 10/22/12: Adapted for use in VIC by Michael Brunke
!
!
! local pointers to implicit in arguments
!
   real(r8), pointer :: evergreen(:) ! binary flag for evergreen leaf habit (0 or 1)
   real(r8), pointer :: woody(:)     ! binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: leafcn(:)    ! leaf C:N (gC/gN)
   real(r8), pointer :: deadwdcn(:)  ! dead wood (xylem and heartwood) C:N (gC/gN)
   integer , pointer :: ivt(:)       ! pft vegetation type
   logical , pointer :: lakpoi(:)    ! true => landunit is a lake point
   integer , pointer :: plandunit(:) ! landunit index associated with each pft
   integer , pointer :: clandunit(:) ! landunit index associated with each column
   integer , pointer :: itypelun(:)  ! landunit type
!
! local pointers to implicit out arguments
!
   real(r8), pointer :: forc_hgt_u_pft(:)    !observational height of wind at pft-level [m]
   real(r8), pointer :: annsum_counter(:) ! seconds since last annual accumulator turnover
   real(r8), pointer :: cannsum_npp(:)    ! annual sum of NPP, averaged from pft-level (gC/m2/yr)
   real(r8), pointer :: cannavg_t2m(:)    !annual average of 2m air temperature, averaged from pft-level (K)
   real(r8), pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real(r8), pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real(r8), pointer :: soil1c(:)             ! (gC/m2) soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)             ! (gC/m2) soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)             ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)             ! (gC/m2) soil organic matter C (slowest pool)
   real(r8), pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real(r8), pointer :: litr1n(:)             ! (gN/m2) litter labile N
   real(r8), pointer :: litr2n(:)             ! (gN/m2) litter cellulose N
   real(r8), pointer :: litr3n(:)             ! (gN/m2) litter lignin N
   real(r8), pointer :: soil1n(:)             ! (gN/m2) soil organic matter N (fast pool)
   real(r8), pointer :: soil2n(:)             ! (gN/m2) soil organic matter N (medium pool)
   real(r8), pointer :: soil3n(:)             ! (gN/m2) soil orgainc matter N (slow pool)
   real(r8), pointer :: soil4n(:)             ! (gN/m2) soil orgainc matter N (slowest pool)
   real(r8), pointer :: sminn(:)              ! (gN/m2) soil mineral N
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: xsmrpool(:)           ! (gC/m2) abstract C pool to meet excess MR demand
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real(r8), pointer :: npool(:)              ! (gN/m2) temporary plant N pool
   real(r8), pointer :: psnsun(:)             ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha(:)             ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: laisun(:)             ! sunlit projected leaf area index
   real(r8), pointer :: laisha(:)             ! shaded projected leaf area index
   real(r8), pointer :: dormant_flag(:)       ! dormancy flag
   real(r8), pointer :: days_active(:)        ! number of days since last dormancy
   real(r8), pointer :: onset_flag(:)         ! onset flag
   real(r8), pointer :: onset_counter(:)      ! onset days counter
   real(r8), pointer :: onset_gddflag(:)      ! onset flag for growing degree day sum
   real(r8), pointer :: onset_fdd(:)          ! onset freezing degree days counter
   real(r8), pointer :: onset_gdd(:)          ! onset growing degree days
   real(r8), pointer :: onset_swi(:)          ! onset soil water index
   real(r8), pointer :: offset_flag(:)        ! offset flag
   real(r8), pointer :: offset_counter(:)     ! offset days counter
   real(r8), pointer :: offset_fdd(:)         ! offset freezing degree days counter
   real(r8), pointer :: offset_swi(:)         ! offset soil water index
   real(r8), pointer :: lgsf(:)               ! long growing season factor [0-1]
   real(r8), pointer :: bglfr(:)              ! background litterfall rate (1/s)
   real(r8), pointer :: bgtr(:)               ! background transfer rate (1/s)
   real(r8), pointer :: dayl(:)               ! daylength (seconds)
   real(r8), pointer :: prev_dayl(:)          ! daylength from previous timestep (seconds)
   real(r8), pointer :: annavg_t2m(:)         ! annual average 2m air temperature (K)
   real(r8), pointer :: tempavg_t2m(:)        ! temporary average 2m air temperature (K)
   real(r8), pointer :: gpp(:)                ! GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: availc(:)             ! C flux available for allocation (gC/m2/s)
   real(r8), pointer :: xsmrpool_recover(:)   ! C flux assigned to recovery of negative cpool (gC/m2/s)
   real(r8), pointer :: alloc_pnow(:)         ! fraction of current allocation to display as new growth (DIM)
   real(r8), pointer :: c_allometry(:)        ! C allocation index (DIM)
   real(r8), pointer :: n_allometry(:)        ! N allocation index (DIM)
   real(r8), pointer :: plant_ndemand(:)      ! N flux required to support initial GPP (gN/m2/s)
   real(r8), pointer :: tempsum_potential_gpp(:) ! temporary annual sum of plant_ndemand
   real(r8), pointer :: annsum_potential_gpp(:)  ! annual sum of plant_ndemand
   real(r8), pointer :: tempmax_retransn(:)   ! temporary max of retranslocated N pool (gN/m2)
   real(r8), pointer :: annmax_retransn(:)    ! annual max of retranslocated N pool (gN/m2)
   real(r8), pointer :: avail_retransn(:)     ! N flux available from retranslocation pool (gN/m2/s)
   real(r8), pointer :: plant_nalloc(:)       ! total allocated N flux (gN/m2/s)
   real(r8), pointer :: plant_calloc(:)       ! total allocated C flux (gC/m2/s)
   real(r8), pointer :: excess_cflux(:)       ! C flux not allocated due to downregulation (gC/m2/s)
   real(r8), pointer :: downreg(:)            ! fractional reduction in GPP due to N limitation (DIM)
   real(r8), pointer :: tempsum_npp(:)        ! temporary annual sum of NPP
   real(r8), pointer :: annsum_npp(:)         ! annual sum of NPP
!#if (defined CNDV)
!   real(r8), pointer :: tempsum_litfall(:)    ! temporary annual sum of litfall
!   real(r8), pointer :: annsum_litfall(:)     ! annual sum of litfall
!#endif
   real(r8), pointer :: qflx_drain(:)         ! sub-surface runoff (mm H2O /s)
   ! new variables for fire
   real(r8), pointer :: wf(:)                 ! soil moisture in top 0.5 m
   real(r8), pointer :: me(:)                 ! moisture of extinction (proportion)
   real(r8), pointer :: fire_prob(:)          ! daily fire probability (0-1)
   real(r8), pointer :: mean_fire_prob(:)     ! e-folding mean of daily fire probability (0-1)
   real(r8), pointer :: fireseasonl(:)        ! annual fire season length (days, <= 365)
   real(r8), pointer :: farea_burned(:)       ! timestep fractional area burned (proportion)
   real(r8), pointer :: ann_farea_burned(:)   ! annual total fractional area burned (proportion)
   real(r8), pointer :: col_ctrunc(:)         ! (gC/m2) column-level sink for C truncation
   real(r8), pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
   real(r8), pointer :: totecosysc(:)         ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
   real(r8), pointer :: totlitc(:)            ! (gC/m2) total litter carbon
   real(r8), pointer :: totsomc(:)            ! (gC/m2) total soil organic matter carbon

   real(r8), pointer :: col_ntrunc(:)         ! (gN/m2) column-level sink for N truncation
   real(r8), pointer :: totcoln(:)            ! (gN/m2) total column nitrogen, incl veg
   real(r8), pointer :: totecosysn(:)         ! (gN/m2) total ecosystem nitrogen, incl veg
   real(r8), pointer :: totlitn(:)            ! (gN/m2) total litter nitrogen
   real(r8), pointer :: totsomn(:)            ! (gN/m2) total soil organic matter nitrogen
   real(r8), pointer :: dispvegc(:)           ! (gC/m2) displayed veg carbon, excluding storage and cpool
   real(r8), pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
   real(r8), pointer :: storvegc(:)           ! (gC/m2) stored vegetation carbon, excluding cpool
   real(r8), pointer :: totpftc(:)            ! (gC/m2) total pft-level carbon, including cpool
   real(r8), pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool
   real(r8), pointer :: prev_frootc_to_litter(:)!previous timestep froot C litterfall flux (gC/m2/s)
   real(r8), pointer :: prev_leafc_to_litter(:) !previous timestep leaf C litterfall flux (gC/m2/s)
   real(r8), pointer :: dispvegn(:)           ! (gN/m2) displayed veg nitrogen, excluding storage
   real(r8), pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
   real(r8), pointer :: storvegn(:)           ! (gN/m2) stored vegetation nitrogen
   real(r8), pointer :: totpftn(:)            ! (gN/m2) total pft-level nitrogen
   real(r8), pointer :: totvegn(:)            ! (gN/m2) total vegetation nitrogen
   real(r8), pointer :: lncsha(:)             ! leaf N concentration per unit projected LAI (gN leaf/m^2)
   real(r8), pointer :: lncsun(:)             ! leaf N concentration per unit projected LAI (gN leaf/m^2)
   real(r8), pointer :: vcmxsha(:)            ! shaded leaf Vcmax (umolCO2/m^2/s)
   real(r8), pointer :: vcmxsun(:)            ! sunlit leaf Vcmax (umolCO2/m^2/s)
   ! dynamic landuse variables
   real(r8), pointer :: seedc(:)              ! (gC/m2) column-level pool for seeding new PFTs
   real(r8), pointer :: prod10c(:)            ! (gC/m2) wood product C pool, 10-year lifespan
   real(r8), pointer :: prod100c(:)           ! (gC/m2) wood product C pool, 100-year lifespan
   real(r8), pointer :: totprodc(:)           ! (gC/m2) total wood product C
   real(r8), pointer :: seedn(:)              ! (gN/m2) column-level pool for seeding new PFTs
   real(r8), pointer :: prod10n(:)            ! (gN/m2) wood product N pool, 10-year lifespan
   real(r8), pointer :: prod100n(:)           ! (gN/m2) wood product N pool, 100-year lifespan
   real(r8), pointer :: totprodn(:)           ! (gN/m2) total wood product N
!
! !LOCAL VARIABLES:
   integer :: g,l,c,p      ! indices
!EOP
!-----------------------------------------------------------------------

    ! assign local pointers at the gridcell level

    ! assign local pointers at the column level
    clandunit                      => clm3%g%c%landunit
    annsum_counter                 => clm3%g%c%cps%annsum_counter
    cannsum_npp                    => clm3%g%c%cps%cannsum_npp
    cannavg_t2m                    => clm3%g%c%cps%cannavg_t2m
    wf                             => clm3%g%c%cps%wf
    me                             => clm3%g%c%cps%me
    fire_prob                      => clm3%g%c%cps%fire_prob
    mean_fire_prob                 => clm3%g%c%cps%mean_fire_prob
    fireseasonl                    => clm3%g%c%cps%fireseasonl
    farea_burned                   => clm3%g%c%cps%farea_burned
    ann_farea_burned               => clm3%g%c%cps%ann_farea_burned
    qflx_drain                     => clm3%g%c%cwf%qflx_drain
    cwdc                           => clm3%g%c%ccs%cwdc
    litr1c                         => clm3%g%c%ccs%litr1c
    litr2c                         => clm3%g%c%ccs%litr2c
    litr3c                         => clm3%g%c%ccs%litr3c
    soil1c                         => clm3%g%c%ccs%soil1c
    soil2c                         => clm3%g%c%ccs%soil2c
    soil3c                         => clm3%g%c%ccs%soil3c
    soil4c                         => clm3%g%c%ccs%soil4c
    
    ! dynamic landuse variables
    seedc                          => clm3%g%c%ccs%seedc
    prod10c                        => clm3%g%c%ccs%prod10c
    prod100c                       => clm3%g%c%ccs%prod100c
    totprodc                       => clm3%g%c%ccs%totprodc
    seedn                          => clm3%g%c%cns%seedn
    prod10n                        => clm3%g%c%cns%prod10n
    prod100n                       => clm3%g%c%cns%prod100n
    totprodn                       => clm3%g%c%cns%totprodn
    
    cwdn                           => clm3%g%c%cns%cwdn
    litr1n                         => clm3%g%c%cns%litr1n
    litr2n                         => clm3%g%c%cns%litr2n
    litr3n                         => clm3%g%c%cns%litr3n
    soil1n                         => clm3%g%c%cns%soil1n
    soil2n                         => clm3%g%c%cns%soil2n
    soil3n                         => clm3%g%c%cns%soil3n
    soil4n                         => clm3%g%c%cns%soil4n
    sminn                          => clm3%g%c%cns%sminn
    col_ctrunc                     => clm3%g%c%ccs%col_ctrunc
    totcolc                        => clm3%g%c%ccs%totcolc
    totecosysc                     => clm3%g%c%ccs%totecosysc
    totlitc                        => clm3%g%c%ccs%totlitc
    totsomc                        => clm3%g%c%ccs%totsomc

    col_ntrunc                     => clm3%g%c%cns%col_ntrunc
    totcoln                        => clm3%g%c%cns%totcoln
    totecosysn                     => clm3%g%c%cns%totecosysn
    totlitn                        => clm3%g%c%cns%totlitn
    totsomn                        => clm3%g%c%cns%totsomn

    ! assign local pointers at the pft level
    ivt                            => clm3%g%c%p%itype
    plandunit                      => clm3%g%c%p%landunit
    leafc                          => clm3%g%c%p%pcs%leafc
    leafc_storage                  => clm3%g%c%p%pcs%leafc_storage
    leafc_xfer                     => clm3%g%c%p%pcs%leafc_xfer
    frootc                         => clm3%g%c%p%pcs%frootc
    frootc_storage                 => clm3%g%c%p%pcs%frootc_storage
    frootc_xfer                    => clm3%g%c%p%pcs%frootc_xfer
    livestemc                      => clm3%g%c%p%pcs%livestemc
    livestemc_storage              => clm3%g%c%p%pcs%livestemc_storage
    livestemc_xfer                 => clm3%g%c%p%pcs%livestemc_xfer
    deadstemc                      => clm3%g%c%p%pcs%deadstemc
    deadstemc_storage              => clm3%g%c%p%pcs%deadstemc_storage
    deadstemc_xfer                 => clm3%g%c%p%pcs%deadstemc_xfer
    livecrootc                     => clm3%g%c%p%pcs%livecrootc
    livecrootc_storage             => clm3%g%c%p%pcs%livecrootc_storage
    livecrootc_xfer                => clm3%g%c%p%pcs%livecrootc_xfer
    deadcrootc                     => clm3%g%c%p%pcs%deadcrootc
    deadcrootc_storage             => clm3%g%c%p%pcs%deadcrootc_storage
    deadcrootc_xfer                => clm3%g%c%p%pcs%deadcrootc_xfer
    gresp_storage                  => clm3%g%c%p%pcs%gresp_storage
    gresp_xfer                     => clm3%g%c%p%pcs%gresp_xfer
    cpool                          => clm3%g%c%p%pcs%cpool
    xsmrpool                       => clm3%g%c%p%pcs%xsmrpool
    forc_hgt_u_pft                 => clm3%g%c%p%pps%forc_hgt_u_pft

    leafn                          => clm3%g%c%p%pns%leafn
    leafn_storage                  => clm3%g%c%p%pns%leafn_storage
    leafn_xfer                     => clm3%g%c%p%pns%leafn_xfer
    frootn                         => clm3%g%c%p%pns%frootn
    frootn_storage                 => clm3%g%c%p%pns%frootn_storage
    frootn_xfer                    => clm3%g%c%p%pns%frootn_xfer
    livestemn                      => clm3%g%c%p%pns%livestemn
    livestemn_storage              => clm3%g%c%p%pns%livestemn_storage
    livestemn_xfer                 => clm3%g%c%p%pns%livestemn_xfer
    deadstemn                      => clm3%g%c%p%pns%deadstemn
    deadstemn_storage              => clm3%g%c%p%pns%deadstemn_storage
    deadstemn_xfer                 => clm3%g%c%p%pns%deadstemn_xfer
    livecrootn                     => clm3%g%c%p%pns%livecrootn
    livecrootn_storage             => clm3%g%c%p%pns%livecrootn_storage
    livecrootn_xfer                => clm3%g%c%p%pns%livecrootn_xfer
    deadcrootn                     => clm3%g%c%p%pns%deadcrootn
    deadcrootn_storage             => clm3%g%c%p%pns%deadcrootn_storage
    deadcrootn_xfer                => clm3%g%c%p%pns%deadcrootn_xfer
    retransn                       => clm3%g%c%p%pns%retransn
    npool                          => clm3%g%c%p%pns%npool
    psnsun                         => clm3%g%c%p%pcf%psnsun
    psnsha                         => clm3%g%c%p%pcf%psnsha
    laisun                         => clm3%g%c%p%pps%laisun
    laisha                         => clm3%g%c%p%pps%laisha
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
    offset_fdd                     => clm3%g%c%p%pepv%offset_fdd
    offset_swi                     => clm3%g%c%p%pepv%offset_swi
    lgsf                           => clm3%g%c%p%pepv%lgsf
    bglfr                          => clm3%g%c%p%pepv%bglfr
    bgtr                           => clm3%g%c%p%pepv%bgtr
    dayl                           => clm3%g%c%p%pepv%dayl
    prev_dayl                      => clm3%g%c%p%pepv%prev_dayl
    annavg_t2m                     => clm3%g%c%p%pepv%annavg_t2m
    tempavg_t2m                    => clm3%g%c%p%pepv%tempavg_t2m
    gpp                            => clm3%g%c%p%pepv%gpp
    availc                         => clm3%g%c%p%pepv%availc
    xsmrpool_recover                  => clm3%g%c%p%pepv%xsmrpool_recover
    alloc_pnow                     => clm3%g%c%p%pepv%alloc_pnow
    c_allometry                    => clm3%g%c%p%pepv%c_allometry
    n_allometry                    => clm3%g%c%p%pepv%n_allometry
    plant_ndemand                  => clm3%g%c%p%pepv%plant_ndemand
    tempsum_potential_gpp          => clm3%g%c%p%pepv%tempsum_potential_gpp
    annsum_potential_gpp           => clm3%g%c%p%pepv%annsum_potential_gpp
    tempmax_retransn               => clm3%g%c%p%pepv%tempmax_retransn
    annmax_retransn                => clm3%g%c%p%pepv%annmax_retransn
    avail_retransn                 => clm3%g%c%p%pepv%avail_retransn
    plant_nalloc                   => clm3%g%c%p%pepv%plant_nalloc
    plant_calloc                   => clm3%g%c%p%pepv%plant_calloc
    excess_cflux                   => clm3%g%c%p%pepv%excess_cflux
    downreg                        => clm3%g%c%p%pepv%downreg
    tempsum_npp                    => clm3%g%c%p%pepv%tempsum_npp
    annsum_npp                     => clm3%g%c%p%pepv%annsum_npp
    pft_ctrunc                     => clm3%g%c%p%pcs%pft_ctrunc
    dispvegc                       => clm3%g%c%p%pcs%dispvegc
    storvegc                       => clm3%g%c%p%pcs%storvegc
    totpftc                        => clm3%g%c%p%pcs%totpftc
    totvegc                        => clm3%g%c%p%pcs%totvegc
    prev_frootc_to_litter          => clm3%g%c%p%pepv%prev_frootc_to_litter
    prev_leafc_to_litter           => clm3%g%c%p%pepv%prev_leafc_to_litter
    dispvegn                       => clm3%g%c%p%pns%dispvegn
    pft_ntrunc                     => clm3%g%c%p%pns%pft_ntrunc
    storvegn                       => clm3%g%c%p%pns%storvegn
    totpftn                        => clm3%g%c%p%pns%totpftn
    totvegn                        => clm3%g%c%p%pns%totvegn
    lncsha                         => clm3%g%c%p%pps%lncsha
    lncsun                         => clm3%g%c%p%pps%lncsun
    vcmxsha                        => clm3%g%c%p%pps%vcmxsha
    vcmxsun                        => clm3%g%c%p%pps%vcmxsun
    
    ! assign local pointers for ecophysiological constants
    evergreen                      => pftcon%evergreen
    woody                          => pftcon%woody
    leafcn                         => pftcon%leafcn
    deadwdcn                       => pftcon%deadwdcn

   ! Added 5/4/04, PET: initialize forc_hgt_u (gridcell-level),
   ! since this is not initialized before first call to CNVegStructUpdate,
   ! and it is required to set the upper bound for canopy top height.
   ! Changed 3/21/08, KO: still needed but don't have sufficient information 
   ! to set this properly (e.g., pft-level displacement height and roughness 
   ! length). So leave at 30m.
!dir$ concurrent
!cdir nodep
   do p = begp, endp
      forc_hgt_u_pft(p) = 30._r8
   end do

   ! initialize column-level variables
!dir$ concurrent
!cdir nodep
   do c = begc, endc

         ! column physical state variables
         annsum_counter(c) = 0._r8
         cannsum_npp(c)    = 0._r8
         cannavg_t2m(c)    = 280._r8
         wf(c) = 1.0_r8  ! it needs to be non zero so the first time step has no fires
         me(c) = 0._r8
         fire_prob(c) = 0._r8
         mean_fire_prob(c) = 0._r8
         fireseasonl(c) = 0._r8
         farea_burned(c) = 0._r8
         ann_farea_burned(c) = 0._r8
         
         ! needed for CNNLeaching
         qflx_drain(c) = 0._r8

         ! column carbon state variable initialization
         cwdc(c)   = 0._r8
         litr1c(c) = 0._r8
         litr2c(c) = 0._r8
         litr3c(c) = 0._r8
         soil1c(c) = 0._r8
         soil2c(c) = 0._r8
         soil3c(c) = 0._r8
         soil4c(c) = 10._r8
         col_ctrunc(c) = 0._r8
         totlitc(c)    = 0._r8
         totsomc(c)    = 0._r8
         totecosysc(c) = 0._r8
         totcolc(c)    = 0._r8

         ! column nitrogen state variables
         cwdn(c)   = cwdc(c) / 500._r8
         litr1n(c) = litr1c(c) / 90._r8
         litr2n(c) = litr2c(c) / 90._r8
         litr3n(c) = litr3c(c) / 90._r8
         soil1n(c) = soil1c(c) / 12._r8
         soil2n(c) = soil2c(c) / 12._r8
         soil3n(c) = soil3c(c) / 10._r8
         soil4n(c) = soil4c(c) / 10._r8
         sminn(c) = 0._r8
         col_ntrunc(c) = 0._r8
         totlitn(c)    = 0._r8
         totsomn(c)    = 0._r8
         totecosysn(c) = 0._r8
         totcoln(c)    = 0._r8

	 ! dynamic landcover state variables
         seedc(c)  = 0._r8
	 prod10c(c)    = 0._r8
	 prod100c(c)   = 0._r8
	 totprodc(c)   = 0._r8
	 seedn(c)      = 0._r8
	 prod10n(c)    = 0._r8
	 prod100n(c)   = 0._r8
	 totprodn(c)   = 0._r8
	 
	 ! also initialize dynamic landcover fluxes so that they have
	 ! real values on first timestep, prior to calling pftdyn_cnbal
	 clm3%g%c%ccf%dwt_seedc_to_leaf(c) = 0._r8
	 clm3%g%c%ccf%dwt_seedc_to_deadstem(c) = 0._r8
	 clm3%g%c%ccf%dwt_conv_cflux(c) = 0._r8
	 clm3%g%c%ccf%dwt_prod10c_gain(c) = 0._r8
	 clm3%g%c%ccf%prod10c_loss(c) = 0._r8
	 clm3%g%c%ccf%dwt_prod100c_gain(c) = 0._r8
	 clm3%g%c%ccf%prod100c_loss(c) = 0._r8
	 clm3%g%c%ccf%dwt_frootc_to_litr1c(c) = 0._r8
	 clm3%g%c%ccf%dwt_frootc_to_litr2c(c) = 0._r8
	 clm3%g%c%ccf%dwt_frootc_to_litr3c(c) = 0._r8
	 clm3%g%c%ccf%dwt_livecrootc_to_cwdc(c) = 0._r8
	 clm3%g%c%ccf%dwt_deadcrootc_to_cwdc(c) = 0._r8
	 clm3%g%c%ccf%dwt_closs(c) = 0._r8
	 clm3%g%c%cnf%dwt_seedn_to_leaf(c) = 0._r8
	 clm3%g%c%cnf%dwt_seedn_to_deadstem(c) = 0._r8
	 clm3%g%c%cnf%dwt_conv_nflux(c) = 0._r8
	 clm3%g%c%cnf%dwt_prod10n_gain(c) = 0._r8
	 clm3%g%c%cnf%prod10n_loss(c) = 0._r8
	 clm3%g%c%cnf%dwt_prod100n_gain(c) = 0._r8
	 clm3%g%c%cnf%prod100n_loss(c) = 0._r8
	 clm3%g%c%cnf%dwt_frootn_to_litr1n(c) = 0._r8
	 clm3%g%c%cnf%dwt_frootn_to_litr2n(c) = 0._r8
	 clm3%g%c%cnf%dwt_frootn_to_litr3n(c) = 0._r8
	 clm3%g%c%cnf%dwt_livecrootn_to_cwdn(c) = 0._r8
	 clm3%g%c%cnf%dwt_deadcrootn_to_cwdn(c) = 0._r8
	 clm3%g%c%cnf%dwt_nloss(c) = 0._r8
   end do

   ! initialize pft-level variables
!dir$ concurrent
!cdir nodep
   do p = begp, endp

         ! carbon state variables
         if (ivt(p) == noveg) then
            leafc(p) = 0._r8
            leafc_storage(p) = 0._r8
         else
            if (evergreen(ivt(p)) == 1._r8) then
               leafc(p) = 1._r8
               leafc_storage(p) = 0._r8
            else
               leafc(p) = 0._r8
               leafc_storage(p) = 1._r8
            end if
         end if
         
         leafc_xfer(p) = 0._r8
         frootc(p) = 0._r8
         frootc_storage(p) = 0._r8
         frootc_xfer(p) = 0._r8
         livestemc(p) = 0._r8
         livestemc_storage(p) = 0._r8
         livestemc_xfer(p) = 0._r8

         ! tree types need to be initialized with some stem mass so that
         ! roughness length is not zero in canopy flux calculation

         if (woody(ivt(p)) == 1._r8) then
            deadstemc(p) = 0.1_r8
         else
            deadstemc(p) = 0._r8
         end if

         deadstemc_storage(p) = 0._r8
         deadstemc_xfer(p) = 0._r8
         livecrootc(p) = 0._r8
         livecrootc_storage(p) = 0._r8
         livecrootc_xfer(p) = 0._r8
         deadcrootc(p) = 0._r8
         deadcrootc_storage(p) = 0._r8
         deadcrootc_xfer(p) = 0._r8
         gresp_storage(p) = 0._r8
         gresp_xfer(p) = 0._r8
         cpool(p) = 0._r8
         xsmrpool(p) = 0._r8
         pft_ctrunc(p) = 0._r8
         dispvegc(p) = 0._r8
         storvegc(p) = 0._r8
         totpftc(p)  = 0._r8

         ! calculate totvegc explicitly so that it is available for the isotope 	
         ! code on the first time step.
         totvegc(p)  = leafc(p) + leafc_storage(p) + leafc_xfer(p) + frootc(p) +  &
            frootc_storage(p) + frootc_xfer(p) + livestemc(p) + livestemc_storage(p) +  &
            livestemc_xfer(p) + deadstemc(p) + deadstemc_storage(p) + deadstemc_xfer(p) +  &
            livecrootc(p) + livecrootc_storage(p) + livecrootc_xfer(p) + deadcrootc(p) +  &
            deadcrootc_storage(p) + deadcrootc_xfer(p) + gresp_storage(p) +  &
            gresp_xfer(p) + cpool(p)

         ! nitrogen state variables
         if (ivt(p) == noveg) then
            leafn(p) = 0._r8
            leafn_storage(p) = 0._r8
         else
            leafn(p) = leafc(p) / leafcn(ivt(p))
            leafn_storage(p) = leafc_storage(p) / leafcn(ivt(p))
         end if

         leafn_xfer(p) = 0._r8
         frootn(p) = 0._r8
         frootn_storage(p) = 0._r8
         frootn_xfer(p) = 0._r8
         livestemn(p) = 0._r8
         livestemn_storage(p) = 0._r8
         livestemn_xfer(p) = 0._r8

         ! tree types need to be initialized with some stem mass so that
         ! roughness length is not zero in canopy flux calculation

         if (woody(ivt(p)) == 1._r8) then
            deadstemn(p) = deadstemc(p) / deadwdcn(ivt(p))
         else
            deadstemn(p) = 0._r8
         end if

         deadstemn_storage(p) = 0._r8
         deadstemn_xfer(p) = 0._r8
         livecrootn(p) = 0._r8
         livecrootn_storage(p) = 0._r8
         livecrootn_xfer(p) = 0._r8
         deadcrootn(p) = 0._r8
         deadcrootn_storage(p) = 0._r8
         deadcrootn_xfer(p) = 0._r8
         retransn(p) = 0._r8
         npool(p) = 0._r8
         pft_ntrunc(p) = 0._r8
         dispvegn(p) = 0._r8
         storvegn(p) = 0._r8
         totvegn(p)  = 0._r8
         totpftn(p)  = 0._r8

         ! initialization for psnsun and psnsha required for
         ! proper arbitrary initialization of allocation routine
         ! in initial ecosysdyn call

         psnsun(p) = 0._r8
         psnsha(p) = 0._r8
         laisun(p) = 0._r8
         laisha(p) = 0._r8
         lncsun(p) = 0._r8
         lncsha(p) = 0._r8
         vcmxsun(p) = 0._r8
         vcmxsha(p) = 0._r8

         ! ecophysiological variables
         ! phenology variables
         dormant_flag(p) = 1._r8
         days_active(p) = 0._r8
         onset_flag(p) = 0._r8
         onset_counter(p) = 0._r8
         onset_gddflag(p) = 0._r8
         onset_fdd(p) = 0._r8
         onset_gdd(p) = 0._r8
         onset_swi(p) = 0.0_r8
         offset_flag(p) = 0._r8
         offset_counter(p) = 0._r8
         offset_fdd(p) = 0._r8
         offset_swi(p) = 0._r8
         lgsf(p) = 0._r8
         bglfr(p) = 0._r8
         bgtr(p) = 0._r8
         dayl(p) = 0._r8                  ! Daylength initialized, MAB 3/21/14
         prev_dayl(p) = 0._r8
         annavg_t2m(p) = 280._r8
         tempavg_t2m(p) = 0._r8

         ! non-phenology variables
         gpp(p) = 0._r8
         availc(p) = 0._r8
         xsmrpool_recover(p) = 0._r8
         alloc_pnow(p) = 1._r8
         c_allometry(p) = 0._r8
         n_allometry(p) = 0._r8
         plant_ndemand(p) = 0._r8
         tempsum_potential_gpp(p) = 0._r8
         annsum_potential_gpp(p) = 0._r8
         tempmax_retransn(p) = 0._r8
         annmax_retransn(p) = 0._r8
         avail_retransn(p) = 0._r8
         plant_nalloc(p) = 0._r8
         plant_calloc(p) = 0._r8
         excess_cflux(p) = 0._r8
         downreg(p) = 0._r8
         prev_leafc_to_litter(p) = 0._r8
         prev_frootc_to_litter(p) = 0._r8
         tempsum_npp(p) = 0._r8
         annsum_npp(p) = 0._r8
!#if (defined CNDV)
!         tempsum_litfall(p) = 0._r8
!         annsum_litfall(p) = 0._r8
!#endif
		 
		 

   end do   ! end of loop over pfts  
!#endif

end subroutine CNiniTimeVar

end module CNiniTimeVarMod
