module CNAnnualUpdateMod
!#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNAnnualUpdateMod
!
! !DESCRIPTION:
! Module for updating annual summation variables
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: CNAnnualUpdate
!
! !REVISION HISTORY:
! 4/23/2004: Created by Peter Thornton
! 1/1/2012: Adapted for use in VIC by Michael Brunke
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNAnnualUpdate
!
! !INTERFACE:
subroutine CNAnnualUpdate(dt, yr, lbc, ubc, lbp, ubp, num_soilc, num_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update annual summation variables
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_days_per_year
   use clm_varcon      , only: secspday
   use pftvarcon       , only: npcropmin
   use subgridAveMod   , only: p2c_1d
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: dt           ! radiation time step (seconds)
   integer, intent(in) :: yr            ! current year
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: lbp, ubp        ! pft bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 10/1/03: Created by Peter Thornton
! 1/1/2012: Adapted for use in VIC by Michael Brunke
! 02/07/14: Added ignoring crops in PFT loops (Michael Brunke).
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: pcolumn(:)               ! index into column level
                                                 ! quantities
   integer , pointer :: ivt(:)                   ! pft vegetation index
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: annsum_counter(:)        ! seconds since last annual accumulator turnover
   real(r8), pointer :: tempsum_potential_gpp(:) ! temporary annual sum of potential GPP
   real(r8), pointer :: annsum_potential_gpp(:)  ! annual sum of potential GPP
   real(r8), pointer :: tempmax_retransn(:)      ! temporary annual max of retranslocated N pool (gN/m2)
   real(r8), pointer :: annmax_retransn(:)       ! annual max of retranslocated N pool (gN/m2)
   real(r8), pointer :: tempavg_t2m(:)           ! temporary average 2m air temperature (K)
   real(r8), pointer :: annavg_t2m(:)            ! annual average 2m air temperature (K)
   real(r8), pointer :: tempsum_npp(:)           ! temporary sum NPP (gC/m2/yr)
   real(r8), pointer :: annsum_npp(:)            ! annual sum NPP (gC/m2/yr)
   real(r8), pointer :: cannsum_npp(:)           ! column annual sum NPP (gC/m2/yr)
   real(r8), pointer :: cannavg_t2m(:)    !annual average of 2m air temperature, averaged from pft-level (K)
!#if (defined CNDV)
!   real(r8), pointer :: tempsum_litfall(:)       ! temporary sum litfall (gC/m2/yr)
!   real(r8), pointer :: annsum_litfall(:)        ! annual sum litfall (gC/m2/yr)
!#endif
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p          ! indices

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays
   annsum_counter        => clm3%g%c%cps%annsum_counter
   tempsum_potential_gpp => clm3%g%c%p%pepv%tempsum_potential_gpp
   annsum_potential_gpp  => clm3%g%c%p%pepv%annsum_potential_gpp
   tempmax_retransn      => clm3%g%c%p%pepv%tempmax_retransn
   annmax_retransn       => clm3%g%c%p%pepv%annmax_retransn
   tempavg_t2m           => clm3%g%c%p%pepv%tempavg_t2m
   annavg_t2m            => clm3%g%c%p%pepv%annavg_t2m
   tempsum_npp           => clm3%g%c%p%pepv%tempsum_npp
   annsum_npp            => clm3%g%c%p%pepv%annsum_npp
   cannsum_npp           => clm3%g%c%cps%cannsum_npp
   cannavg_t2m           => clm3%g%c%cps%cannavg_t2m
!#if (defined CNDV)
!   tempsum_litfall       => clm3%g%c%p%pepv%tempsum_litfall
!   annsum_litfall        => clm3%g%c%p%pepv%annsum_litfall
!#endif
   pcolumn               => clm3%g%c%p%column
   ivt                   => clm3%g%c%p%itype

   ! column loop
   do c = 1,num_soilc
      annsum_counter(c) = annsum_counter(c) + dt
   end do

!#if (defined CNDV) || (defined CROP)
   ! In the future -- ONLY use this code and remove the similar part below
   ! So the #ifdef on CNDV and CROP would be removed
!   if (annsum_counter(filter_soilc(1)) >= get_days_per_year() * secspday) then ! new (slevis)
!#endif
   ! pft loop
   do p = 1,num_soilp
     if(ivt(p) < npcropmin) then
!#if (!defined CNDV)
      ! In the future -- REMOVE this code and use the equivalent code above always
!      c = pcolumn(p)                                                ! old (slevis)
!      if (annsum_counter(c) >= get_days_per_year() * secspday) then ! old (slevis)
!#endif
         ! update annual plant ndemand accumulator
         annsum_potential_gpp(p)  = tempsum_potential_gpp(p)
         tempsum_potential_gpp(p) = 0._r8

         ! update annual total N retranslocation accumulator
         annmax_retransn(p)  = tempmax_retransn(p)
         tempmax_retransn(p) = 0._r8

         ! update annual average 2m air temperature accumulator
         annavg_t2m(p)  = tempavg_t2m(p)
         tempavg_t2m(p) = 0._r8

         ! update annual NPP accumulator, convert to annual total
         annsum_npp(p) = tempsum_npp(p) * dt
         tempsum_npp(p) = 0._r8

!#if (defined CNDV)
         ! update annual litfall accumulator, convert to annual total
!         annsum_litfall(p) = tempsum_litfall(p) * dt
!         tempsum_litfall(p) = 0._r8
!#endif
!#if (!defined CNDV) && (!defined CROP)
!      end if ! old (slevis)
!#endif
     end if
   end do

   ! use p2c routine to get selected column-average pft-level fluxes and states
   call p2c_1d(lbp, ubp, lbc, ubc, annsum_npp, cannsum_npp, 'unity')
   call p2c_1d(lbp, ubp, lbc, ubc, annavg_t2m, cannavg_t2m, 'unity')
!#if (defined CNDV) || (defined CROP)
!   end if ! new (slevis)
!#endif

   ! column loop
   do c = 1,num_soilc
      if (annsum_counter(c) >= get_days_per_year(yr) * secspday) annsum_counter(c) = 0._r8
   end do

end subroutine CNAnnualUpdate
!-----------------------------------------------------------------------
!#endif

end module CNAnnualUpdateMod
