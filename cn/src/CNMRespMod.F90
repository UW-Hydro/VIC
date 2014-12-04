module CNMRespMod
!#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNMRespMod
!
! !DESCRIPTION:
! Module holding maintenance respiration routines for coupled carbon
! nitrogen code.
!
! !USES:
   use shr_kind_mod , only: r8 => shr_kind_r8
   use shr_const_mod, only: SHR_CONST_TKFRZ
   implicit none
   save
   private
! !PUBLIC MEMBER FUNCTIONS:
   public :: CNMResp
!
! !REVISION HISTORY:
! 8/14/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: Migrated all subroutines to vector data structures.
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNMResp
!
! !INTERFACE:
subroutine CNMResp(nlevgrnd, lbc, ubc, num_soilc, num_soilp)
!
! !DESCRIPTION:
!
! !USES:
   use clmtype
   use pftvarcon    , only : npcropmin
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: nlevgrnd                  ! number of soil layers
   integer, intent(in) :: lbc, ubc                    ! column-index bounds
   integer, intent(in) :: num_soilc                 ! number of soil points in column filter
   integer, intent(in) :: num_soilp                 ! number of soil points in pft filter
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 8/14/03: Created by Peter Thornton
! 10/10/12: Adapted for use in VIC by Michael Brunke
! 02/06/14: Added ignoring crops in PFT loop (Michael Brunke).
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   ! column level
   real(r8), pointer :: t_soisno(:,:) ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   ! pft level
   real(r8), pointer :: t_ref2m(:)    ! 2 m height surface air temperature (Kelvin)
   real(r8), pointer :: leafn(:)      ! (gN/m2) leaf N
   real(r8), pointer :: frootn(:)     ! (gN/m2) fine root N
   real(r8), pointer :: livestemn(:)  ! (gN/m2) live stem N
   real(r8), pointer :: livecrootn(:) ! (gN/m2) live coarse root N
   real(r8), pointer :: rootfr(:,:)   ! fraction of roots in each soil layer  (nlevgrnd)
   integer , pointer :: ivt(:)        ! pft vegetation type
   integer , pointer :: pcolumn(:)    ! index into column level quantities
   integer , pointer :: plandunit(:)  ! index into landunit level quantities
   integer , pointer :: clandunit(:)  ! index into landunit level quantities
   integer , pointer :: itypelun(:)   ! landunit type
   ! ecophysiological constants
   real(r8), pointer :: woody(:)      ! binary flag for woody lifeform (1=woody, 0=not woody)
   logical , pointer :: croplive(:)   ! Flag, true if planted, not harvested
!
! local pointers to implicit in/out arrays
!
   ! pft level
   real(r8), pointer :: leaf_mr(:)
   real(r8), pointer :: froot_mr(:)
   real(r8), pointer :: livestem_mr(:)
   real(r8), pointer :: livecroot_mr(:)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p,j          ! indices
   real(r8):: mr             ! maintenance respiration (gC/m2/s)
   real(r8):: br             ! base rate (gC/gN/s)
   real(r8):: q10            ! temperature dependence
   real(r8):: tc             ! temperature correction, 2m air temp (unitless)
   real(r8):: tcsoi(lbc:ubc,nlevgrnd) ! temperature correction by soil layer (unitless)
!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays
   t_soisno       => clm3%g%c%ces%t_soisno
   t_ref2m        => clm3%g%c%p%pes%t_ref2m
   leafn          => clm3%g%c%p%pns%leafn
   frootn         => clm3%g%c%p%pns%frootn
   livestemn      => clm3%g%c%p%pns%livestemn
   livecrootn     => clm3%g%c%p%pns%livecrootn
   rootfr         => clm3%g%c%p%pps%rootfr
   leaf_mr        => clm3%g%c%p%pcf%leaf_mr
   froot_mr       => clm3%g%c%p%pcf%froot_mr
   livestem_mr    => clm3%g%c%p%pcf%livestem_mr
   livecroot_mr   => clm3%g%c%p%pcf%livecroot_mr
   ivt            => clm3%g%c%p%itype
   woody          => pftcon%woody
   croplive       => clm3%g%c%p%pps%croplive

   ! base rate for maintenance respiration is from:
   ! M. Ryan, 1991. Effects of climate change on plant respiration.
   ! Ecological Applications, 1(2), 157-167.
   ! Original expression is br = 0.0106 molC/(molN h)
   ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
   br = 2.525e-6_r8
   ! Peter Thornton: 3/13/09 
   ! Q10 was originally set to 2.0, an arbitrary choice, but reduced to 1.5 as part of the tuning
   ! to improve seasonal cycle of atmospheric CO2 concentration in global
   ! simulatoins
   q10 = 1.5_r8

   ! column loop to calculate temperature factors in each soil layer
   do j=1,nlevgrnd
      do c = 1, num_soilc

         ! calculate temperature corrections for each soil layer, for use in
         ! estimating fine root maintenance respiration with depth

         tcsoi(c,j) = q10**((t_soisno(c,j)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
      end do
   end do

   ! pft loop for leaves and live wood
   do p = 1, num_soilp
     if(ivt(p) < npcropmin) then

      ! calculate maintenance respiration fluxes in
      ! gC/m2/s for each of the live plant tissues.
      ! Leaf and live wood MR

      tc = q10**((t_ref2m(p)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
      leaf_mr(p) = leafn(p)*br*tc
      if (woody(ivt(p)) == 1) then
         livestem_mr(p) = livestemn(p)*br*tc
         livecroot_mr(p) = livecrootn(p)*br*tc
      else
         livestem_mr(p) = livestemn(p)*br*tc
      end if
     end if
   end do

   ! soil and pft loop for fine root
   do j = 1,nlevgrnd
      do p = 1,num_soilp
       if(ivt(p) < npcropmin) then
        do c = 1, num_soilc

         ! Fine root MR
         ! rootfr(j) sums to 1.0 over all soil layers, and
         ! describes the fraction of root mass that is in each
         ! layer.  This is used with the layer temperature correction
         ! to estimate the total fine root maintenance respiration as a
         ! function of temperature and N content.

         froot_mr(p) = froot_mr(p) + frootn(p)*br*tcsoi(c,j)*rootfr(p,j)
        end do
       end if
      end do
   end do

end subroutine CNMResp

!#endif

end module CNMRespMod
