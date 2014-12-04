module CNWoodProductsMod
!#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNWoodProductsMod
!
! !DESCRIPTION:
! Calculate loss fluxes from wood products pools, and update product pool state variables
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: CNWoodProducts
!
! !REVISION HISTORY:
! 5/20/2009: Created by Peter Thornton
! 1/12/2013: Revised for use in VIC by Michael Brunke
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNWoodProducts
!
! !INTERFACE:
subroutine CNWoodProducts(dt, num_soilc)
!
! !DESCRIPTION:
! Update all loss fluxes from wood product pools, and update product pool state variables
! for both loss and gain terms.  Gain terms are calculated in pftdyn_cnbal() for gains associated
! with changes in landcover, and in CNHarvest(), for gains associated with wood harvest.
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   real(r8):: dt        ! time step (seconds)
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 5/21/09: Created by Peter Thornton
! 1/12/2013: Revised for use in VIC by Michael Brunke
!
! !LOCAL VARIABLES:

   integer :: c         ! indices
   type(column_type),   pointer :: cptr         ! pointer to column derived subtype
   real(r8) :: kprod10       ! decay constant for 10-year product pool
   real(r8) :: kprod100      ! decay constant for 100-year product pool

!EOP
!-----------------------------------------------------------------------

   cptr => clm3%g%c
	
   ! calculate column-level losses from product pools
	! the following (1/s) rate constants result in ~90% loss of initial state over 10 and 100 years,
	! respectively, using a discrete-time fractional decay algorithm.
	kprod10 = 7.2e-9
	kprod100 = 7.2e-10

   do c = 1,num_soilc

		! calculate fluxes (1/sec)
		cptr%ccf%prod10c_loss(c)    = cptr%ccs%prod10c(c)    * kprod10
		cptr%ccf%prod100c_loss(c)   = cptr%ccs%prod100c(c)   * kprod100
		cptr%cnf%prod10n_loss(c)    = cptr%cns%prod10n(c)    * kprod10
		cptr%cnf%prod100n_loss(c)   = cptr%cns%prod100n(c)   * kprod100
	end do

   ! update wood product state variables
   ! column loop
   do c = 1,num_soilc

      ! column-level fluxes

      ! fluxes into wood product pools, from landcover change
      cptr%ccs%prod10c(c)    = cptr%ccs%prod10c(c)    + cptr%ccf%dwt_prod10c_gain(c)*dt
      cptr%ccs%prod100c(c)   = cptr%ccs%prod100c(c)   + cptr%ccf%dwt_prod100c_gain(c)*dt
      cptr%cns%prod10n(c)    = cptr%cns%prod10n(c)    + cptr%cnf%dwt_prod10n_gain(c)*dt
      cptr%cns%prod100n(c)   = cptr%cns%prod100n(c)   + cptr%cnf%dwt_prod100n_gain(c)*dt

      ! fluxes into wood product pools, from harvest
      cptr%ccs%prod10c(c)    = cptr%ccs%prod10c(c)    + cptr%ccf%hrv_deadstemc_to_prod10c(c)*dt
      cptr%ccs%prod100c(c)   = cptr%ccs%prod100c(c)   + cptr%ccf%hrv_deadstemc_to_prod100c(c)*dt
#if (defined C13)
      cptr%cc13s%prod10c(c)  = cptr%cc13s%prod10c(c)  + cptr%cc13f%hrv_deadstemc_to_prod10c(c)*dt
      cptr%cc13s%prod100c(c) = cptr%cc13s%prod100c(c) + cptr%cc13f%hrv_deadstemc_to_prod100c(c)*dt
#endif
      cptr%cns%prod10n(c)    = cptr%cns%prod10n(c)    + cptr%cnf%hrv_deadstemn_to_prod10n(c)*dt
      cptr%cns%prod100n(c)   = cptr%cns%prod100n(c)   + cptr%cnf%hrv_deadstemn_to_prod100n(c)*dt
     
      ! fluxes out of wood product pools, from decomposition
      cptr%ccs%prod10c(c)    = cptr%ccs%prod10c(c)    - cptr%ccf%prod10c_loss(c)*dt
      cptr%ccs%prod100c(c)   = cptr%ccs%prod100c(c)   - cptr%ccf%prod100c_loss(c)*dt
      cptr%cns%prod10n(c)    = cptr%cns%prod10n(c)    - cptr%cnf%prod10n_loss(c)*dt
      cptr%cns%prod100n(c)   = cptr%cns%prod100n(c)   - cptr%cnf%prod100n_loss(c)*dt
 
   end do ! end of column loop

end subroutine CNWoodProducts
!-----------------------------------------------------------------------

!#endif

end module CNWoodProductsMod
