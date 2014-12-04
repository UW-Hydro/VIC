module CNEcosystemDynMod
!#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNEcosystemDynMod
!
! !DESCRIPTION:
! Ecosystem dynamics: phenology, vegetation
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNEcosystemDynInit   ! Ecosystem dynamics initialization
  public :: CNEcosystemDyn       ! Ecosystem dynamics: phenology, vegetation
!
! !REVISION HISTORY:
! Created by Peter Thornton
! 19 May 2009: PET - modified to include call to harvest routine
! 10/10/12: Adapted for use in VIC by Michael Brunke
!
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !PRIVATE TYPES:
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNEcosystemDynInit
!
! !INTERFACE:
  subroutine CNEcosystemDynInit(lbc, ubc, lbp, ubp, dt )
!
! !DESCRIPTION:
! Initialzation of the CN Ecosystem dynamics.
!
! !USES:
    use CNAllocationMod, only : CNAllocationInit
    use CNPhenologyMod , only : CNPhenologyInit
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc        ! column bounds
    integer, intent(in) :: lbp, ubp        ! pft bounds
    real(r8), intent(in) :: dt             ! timestep
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 04/05/11, Erik Kluzek creation
! 10/10/12: Adapted for use in VIC by Michael Brunke
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
     call CNAllocationInit ( lbc, ubc, lbp, ubp, dt )
     call CNPhenologyInit  ( dt, lbp, ubp )

  end subroutine CNEcosystemDynInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNEcosystemDyn
!
! !INTERFACE:
  subroutine CNEcosystemDyn(adspinup, nlevgrnd, nstep, yr, jdy, dt, lbg, ubg, &
	lbc, ubc, lbp, ubp, num_soilc, num_soilp)
!
! !DESCRIPTION:
! The core CN code is executed here. Calculates fluxes for maintenance
! respiration, decomposition, allocation, phenology, and growth respiration.
! These routines happen on the radiation time step so that canopy structure
! stays synchronized with albedo calculations.
!
! !USES:
    use clmtype
    use CNSetValueMod        , only: CNZeroFluxes
    use CNNDynamicsMod       , only: CNNDeposition,CNNFixation, CNNLeaching
    use CNMRespMod           , only: CNMResp
    use CNDecompMod          , only: CNDecompAlloc
    use CNPhenologyMod       , only: CNPhenology
    use CNGRespMod           , only: CNGResp
    use CNCStateUpdate1Mod   , only: CStateUpdate1,CStateUpdate0
    use CNNStateUpdate1Mod   , only: NStateUpdate1
    use CNGapMortalityMod    , only: CNGapMortality
    use CNCStateUpdate2Mod   , only: CStateUpdate2, CStateUpdate2h
    use CNNStateUpdate2Mod   , only: NStateUpdate2, NStateUpdate2h
    use CNFireMod            , only: CNFireArea, CNFireFluxes
    use CNCStateUpdate3Mod   , only: CStateUpdate3
    use CNNStateUpdate3Mod   , only: NStateUpdate3
    use CNBalanceCheckMod    , only: CBalanceCheck, NBalanceCheck
    use CNPrecisionControlMod, only: CNPrecisionControl
    use CNVegStructUpdateMod , only: CNVegStructUpdate
    use CNAnnualUpdateMod    , only: CNAnnualUpdate
    use CNSummaryMod         , only: CSummary, NSummary
    use CNWoodProductsMod    , only: CNWoodProducts
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: adspinup        ! flag for AD spin-up
    real(r8), intent(in) :: dt             ! time step
    integer, intent(in) :: nstep           ! current time step number
    integer, intent(in) :: yr		   ! current year
    integer, intent(in) :: jdy             ! day in year
    integer, intent(in) :: nlevgrnd        ! number of soil layers
    integer, intent(in) :: lbg, ubg        ! gridcell bounds
    integer, intent(in) :: lbc, ubc        ! column bounds
    integer, intent(in) :: lbp, ubp        ! pft bounds
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 10/22/03, Peter Thornton: created from EcosystemDyn during migration to
!                           new vector code.
! 11/3/03, Peter Thornton: removed update of elai, esai, frac_veg_nosno_alb.
!     These are now done in CNVegStructUpdate(), which is called
!     prior to SurfaceAlbedo().
! 11/13/03, Peter Thornton: switched from nolake to soil filtering.
! 10/10/12: Adapted for use in VIC by Michael Brunke
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
! local pointers to implicit out arguments
    real(r8), pointer :: leafc(:)
!
! !OTHER LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

 !   if (doalb) then

    leafc => clm3%g%c%p%pcs%leafc

       ! Call the main CN routines
       call CNZeroFluxes(num_soilc, num_soilp)

       call CNNDeposition(yr, jdy, lbg, ubg, lbc, ubc)

       call CNNFixation(yr, num_soilc)

       call CNMResp(nlevgrnd, lbc, ubc, num_soilc, num_soilp)

       call CNDecompAlloc(adspinup, dt, lbp, ubp, lbc, ubc, num_soilc, &
	num_soilp)

       ! CNphenology needs to be called after CNdecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call CNPhenology(yr, dt, num_soilc, num_soilp)

       call CNGResp(num_soilp)
       
       call CStateUpdate0(dt, num_soilp)

       call CStateUpdate1(dt, num_soilc, num_soilp)

       call NStateUpdate1(dt, num_soilc, num_soilp)

       call CNGapMortality(yr, num_soilc, num_soilp)

       call CStateUpdate2(dt, num_soilc, num_soilp)

       call NStateUpdate2(dt, num_soilc, num_soilp)
       
!       if (fpftdyn /= ' ') then
!          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp)
!       end if 

       call CStateUpdate2h(dt, num_soilc, num_soilp)

       call NStateUpdate2h(dt, num_soilc, num_soilp)
       
       call CNWoodProducts(dt, num_soilc)
       
       call CNFireArea(dt, nstep, yr, num_soilc)

       call CNFireFluxes(dt, lbp, ubp, lbc, ubc, num_soilc, num_soilp)

       ! see full_energy for soil moisture
       call CNNLeaching(nlevgrnd, dt, lbc, ubc, num_soilc)

       call CStateUpdate3(dt, num_soilc, num_soilp)

       call NStateUpdate3(dt, num_soilc, num_soilp)

       call CNPrecisionControl(num_soilc, num_soilp)

       call CNVegStructUpdate(dt, num_soilp)

       call CSummary(lbp, ubp, lbc, ubc, num_soilc, num_soilp)
       
       call NSummary(lbp, ubp, lbc, ubc, num_soilc, num_soilp)

!    end if  !end of if-doalb block

  end subroutine CNEcosystemDyn
!#endif
!-----------------------------------------------------------------------
end  module CNEcosystemDynMod
