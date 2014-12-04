!#include <misc.h>
!#include <preproc.h>

module CNiniSpecialMod

  implicit none
  save
  private
  public :: CNiniSpecial

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNiniSpecial
!
! !INTERFACE:
subroutine CNiniSpecial (nlevgrnd, begg, endg, begc, endc, begp, endp)

!#ifdef CN
!
! !DESCRIPTION:
! One-time initialization of CN variables for special landunits
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use clm_varcon  , only: spval
   use clmtype
   use CNSetValueMod, only : CNSetCps, CNSetCcs, CNSetCns, CNSetCcf, &
	CNSetCnf, CNSetPps, CNSetPcs, CNSetPns, CNSetPcf, CNSetPnf, &
	CNSetPepv
!
! !ARGUMENTS:
   implicit none

   integer :: nlevgrnd     ! number of soil layers
   integer :: begp, endp   ! per-clump/proc beginning and ending pft indices
   integer :: begc, endc   ! per-clump/proc beginning and ending column indices
   integer :: begg, endg   ! per-clump/proc gridcell ending gridcell indices
!
! !CALLED FROM:
! subroutine iniTimeConst in file iniTimeConst.F90
!
! !REVISION HISTORY:
! 11/13/03: Created by Peter Thornton
! 10/22/12: Adapted for use in VIC by Michael Brunke
!
!
! local pointers to implicit out arguments
!
! !LOCAL VARIABLES:
!EOP
   integer :: fc,fp,l,c,p  ! indices
!-----------------------------------------------------------------------
   ! Determine subgrid bounds on this processor
!   call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

   ! initialize column-level fields
   call CNSetCps(nlevgrnd, endc, spval, clm3%g%c%cps)
   call CNSetCcs(endc, 0._r8, clm3%g%c%ccs)
   call CNSetCns(endc, 0._r8, clm3%g%c%cns)
   call CNSetCcf(endc, 0._r8, clm3%g%c%ccf)
   call CNSetCnf(endc, 0._r8, clm3%g%c%cnf)

   ! initialize column-average pft fields
   call CNSetPps(endc, spval, clm3%g%c%cps%pps_a)
   call CNSetPcs(endc, 0._r8, clm3%g%c%ccs%pcs_a)
   call CNSetPns(endc, 0._r8, clm3%g%c%cns%pns_a)
   call CNSetPcf(endc, 0._r8, clm3%g%c%ccf%pcf_a)
   call CNSetPnf(endc, 0._r8, clm3%g%c%cnf%pnf_a)

   ! initialize pft-level fields
   call CNSetPepv(endp, spval, clm3%g%c%p%pepv)
   call CNSetPps(endp, spval, clm3%g%c%p%pps)
   call CNSetPcs(endp, 0._r8, clm3%g%c%p%pcs)
   call CNSetPns(endp, 0._r8, clm3%g%c%p%pns)
   call CNSetPcf(endp, 0._r8, clm3%g%c%p%pcf)
   call CNSetPnf(endp, 0._r8, clm3%g%c%p%pnf)

   ! now loop through special filters and explicitly set the variables that
   ! have to be in place for SurfaceAlbedo and biogeophysics
   ! also set pcf%psnsun and pcf%psnsha to 0 (not included in CNSetPcf())

!dir$ concurrent
!cdir nodep
   do p = 1,endp
      clm3%g%c%p%pps%tlai(p) = 0._r8
      clm3%g%c%p%pps%tsai(p) = 0._r8
      clm3%g%c%p%pps%elai(p) = 0._r8
      clm3%g%c%p%pps%esai(p) = 0._r8
      clm3%g%c%p%pps%htop(p) = 0._r8
      clm3%g%c%p%pps%hbot(p) = 0._r8
      clm3%g%c%p%pps%fwet(p) = 0._r8
      clm3%g%c%p%pps%fdry(p) = 0._r8
      clm3%g%c%p%pps%frac_veg_nosno_alb(p) = 0._r8
      clm3%g%c%p%pps%frac_veg_nosno(p) = 0._r8
      clm3%g%c%p%pcf%psnsun(p) = 0._r8
      clm3%g%c%p%pcf%psnsha(p) = 0._r8
      
   end do

!dir$ concurrent
!cdir nodep
   do c = 1,endc
      clm3%g%c%ccf%pcf_a%psnsun(c) = 0._r8
      clm3%g%c%ccf%pcf_a%psnsha(c) = 0._r8
      
	  ! adding dynpft code
	  clm3%g%c%ccs%seedc(c) = 0._r8
	  clm3%g%c%ccs%prod10c(c) = 0._r8	  
	  clm3%g%c%ccs%prod100c(c) = 0._r8	  
	  clm3%g%c%ccs%totprodc(c) = 0._r8	  
	  clm3%g%c%cns%seedn(c) = 0._r8
	  clm3%g%c%cns%prod10n(c) = 0._r8	  
	  clm3%g%c%cns%prod100n(c) = 0._r8	  
	  clm3%g%c%cns%totprodn(c) = 0._r8	  
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
	  clm3%g%c%ccf%landuseflux(c) = 0._r8
	  clm3%g%c%ccf%landuptake(c) = 0._r8
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

!#endif

end subroutine CNiniSpecial

end module CNiniSpecialMod
