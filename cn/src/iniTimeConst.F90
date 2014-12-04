module iniTimeConstMod

  implicit none
  save
  private
  public :: iniTimeConst

contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: iniTimeConst
!
! !INTERFACE:
subroutine iniTimeConst(nlevgrnd, begg, endg, begc, endc, begp, endp)
!
! !DESCRIPTION:
! Initialize time invariant clm variables
! 1) removed references to shallow lake - since it is not used
! 2) ***Make c%z, c%zi and c%dz allocatable depending on if you
!    have lake or soil
! 3) rootfr only initialized for soil points
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan
  use clmtype
  use clm_varpar  , only : nlevsoi, nlevlak, lsmlon, lsmlat, numpft, numrad, nlevurb
  use clm_varcon  , only : spval
!  use clm_varcon  , only : albsat, albdry
  use pftvarcon   , only : noveg, ntree, roota_par, rootb_par,  &
                           smpso, smpsc, fnitr, nbrdlf_dcd_brl_shrub, &
                           z0mr, displar, dleaf, rhol, rhos, taul, taus, xl, &
                           qe25, mp, c3psn, slatop, dsladlai, leafcn, flnr, woody, &
                           lflitcn, frootcn, livewdcn, deadwdcn, froot_leaf, stem_leaf, croot_stem, &
                           flivewd, fcur, lf_flab, lf_fcel, lf_flig, fr_flab, fr_fcel, fr_flig, &
                           leaf_long, evergreen, stress_decid, season_decid, &
                           resist, pftpar20, pftpar28, pftpar29, pftpar30, pftpar31, &
                           allom1s, allom2s, &
                           allom1 , allom2 , allom3  , reinickerp, dwood
  use CNiniSpecialMod  , only : CNiniSpecial
!  use organicFileMod  , only : organicrd 
!  use clm_varctl      , only : fsnowoptics, fsnowaging
!  use SNICARMod       , only : SnowAge_init, SnowOptics_init

!
! !ARGUMENTS:
  implicit none
  integer  :: nlevgrnd         ! number of soil layers
  integer  :: begp, endp       ! per-proc beginning and ending pft indices
  integer  :: begc, endc       ! per-proc beginning and ending column indices
  integer  :: begl, endl       ! per-proc beginning and ending landunit indices
  integer  :: begg, endg       ! per-proc gridcell ending gridcell indices
!
! !CALLED FROM:
! subroutine initialize in module initializeMod.
!
! !REVISION HISTORY:
! Created by Gordon Bonan.
! Updated to clm2.1 data structrues by Mariana Vertenstein
! 4/26/05, Peter Thornton: Eliminated exponential decrease in saturated hydraulic
!   conductivity (hksat) with depth. 
! Updated: Colette L. Heald (05/06) reading in VOC emission factors
! 27 February 2008: Keith Oleson; Qing Liu (2004) saturated hydraulic conductivity 
! and matric potential
! 29 February 2008: David Lawrence; modified soil thermal and hydraulic properties to
! account for organic matter
! 18 March 2008: David Lawrence; nlevgrnd changes
! 03/28/08 Mark Flanner, read in netcdf files for SNICAR parameters
! 10/23/12 Adapted for use in VIC by Michael Brunke
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: ivt(:)             !  vegetation type index
  integer , pointer :: pcolumn(:)         ! column index of corresponding pft
  integer , pointer :: pgridcell(:)       ! gridcell index of corresponding pft
  integer , pointer :: cgridcell(:)       ! gridcell index of column
  integer , pointer :: ctype(:)           ! column type index
  real(r8), pointer :: lat(:)             ! gridcell latitude (radians)
!
! local pointers to implicit out arguments
!
  real(r8), pointer :: rootfr(:,:)        ! fraction of roots in each soil layer
  real(r8), pointer :: rresis(:,:)        !root resistance by layer (0-1)  (nlevgrnd)	
  real(r8), pointer :: z(:,:)             ! layer depth (m)
  real(r8), pointer :: zi(:,:)            ! interface level below a "z" level (m)
  real(r8), pointer :: dz(:,:)            ! layer thickness depth (m)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
  integer  :: ncid             ! netCDF file id 
  integer  :: n,j,ib,lev,bottom! indices
  integer  :: g,  c,p          ! indices
  integer  :: m                ! vegetation type index
  real(r8) :: bd               ! bulk density of dry soil material [kg/m^3]
  real(r8) :: tkm              ! mineral conductivity
  real(r8) :: xksat            ! maximum hydraulic conductivity of soil [mm/s]
  real(r8) :: scalez = 0.025_r8   ! Soil layer thickness discretization (m)
  real(r8) :: clay,sand        ! temporaries
  real(r8) :: slope,intercept        ! temporary, for rooting distribution
  real(r8) :: temp, max_decl   ! temporary, for calculation of max_dayl

  real(r8),pointer :: arrayl(:)      ! generic global array
  integer ,pointer :: irrayg(:)      ! generic global array
  integer  :: start(3),count(3)      ! netcdf start/count arrays
  integer  :: varid                  ! netCDF id's
  integer  :: ret

  integer  :: ier                                ! error status
  character(len=256) :: locfn                    ! local filename
  character(len= 32) :: subname = 'iniTimeConst' ! subroutine name
!------------------------------------------------------------------------

  integer :: closelatidx,closelonidx
  real(r8):: closelat,closelon
  integer :: iostat

!------------------------------------------------------------------------

  ! Assign local pointers to derived subtypes components (gridcell-level)
  lat             => clm3%g%lat
     
  ! Assign local pointers to derived subtypes components (column-level)

  ctype           => clm3%g%c%itype
  cgridcell       => clm3%g%c%gridcell
  z               => clm3%g%c%cps%z
  dz              => clm3%g%c%cps%dz
  zi              => clm3%g%c%cps%zi

  ! Assign local pointers to derived subtypes components (pft-level)

  ivt             => clm3%g%c%p%itype
  pgridcell       => clm3%g%c%p%gridcell
  pcolumn         => clm3%g%c%p%column
  rootfr          => clm3%g%c%p%pps%rootfr
  rresis          => clm3%g%c%p%pps%rresis

  ! --------------------------------------------------------------------
  ! Initialize time constant arrays of ecophysiological constants and
  ! arrays of dgvm ecophysiological constants
  ! --------------------------------------------------------------------

   do m = 0,numpft
      if (m <= ntree) then
         pftcon%tree(m) = 1
      else
         pftcon%tree(m) = 0
      end if
      pftcon%z0mr(m) = z0mr(m)
      pftcon%displar(m) = displar(m)
      pftcon%dleaf(m) = dleaf(m)
      pftcon%xl(m) = xl(m)
      do ib = 1,numrad
         pftcon%rhol(m,ib) = rhol(m,ib)
         pftcon%rhos(m,ib) = rhos(m,ib)
         pftcon%taul(m,ib) = taul(m,ib)
         pftcon%taus(m,ib) = taus(m,ib)
      end do
      pftcon%qe25(m) = qe25(m)
!      pftcon%vcmx25(m) = vcmx25(m)
      pftcon%mp(m) = mp(m)
      pftcon%c3psn(m) = c3psn(m)
      pftcon%slatop(m) = slatop(m)
      pftcon%dsladlai(m) = dsladlai(m)
      pftcon%leafcn(m) = leafcn(m)
      pftcon%flnr(m) = flnr(m)
      pftcon%smpso(m) = smpso(m)
      pftcon%smpsc(m) = smpsc(m)
      pftcon%fnitr(m) = fnitr(m)
      pftcon%woody(m) = woody(m)
      pftcon%lflitcn(m) = lflitcn(m)
      pftcon%frootcn(m) = frootcn(m)
      pftcon%livewdcn(m) = livewdcn(m)
      pftcon%deadwdcn(m) = deadwdcn(m)
      pftcon%froot_leaf(m) = froot_leaf(m)
      pftcon%stem_leaf(m) = stem_leaf(m)
      pftcon%croot_stem(m) = croot_stem(m)
      pftcon%flivewd(m) = flivewd(m)
      pftcon%fcur(m) = fcur(m)
      pftcon%lf_flab(m) = lf_flab(m)
      pftcon%lf_fcel(m) = lf_fcel(m)
      pftcon%lf_flig(m) = lf_flig(m)
      pftcon%fr_flab(m) = fr_flab(m)
      pftcon%fr_fcel(m) = fr_fcel(m)
      pftcon%fr_flig(m) = fr_flig(m)
!      pftcon%dw_fcel(m) = dw_fcel(m)
!      pftcon%dw_flig(m) = dw_flig(m)
      pftcon%leaf_long(m) = leaf_long(m)
      pftcon%evergreen(m) = evergreen(m)
      pftcon%stress_decid(m) = stress_decid(m)
      pftcon%season_decid(m) = season_decid(m)
      pftcon%resist(m) = resist(m)
      pftcon%dwood(m) = dwood
   end do

!#ifdef CNDV
!   do m = 0,numpft
!      dgv_pftcon%crownarea_max(m) = pftpar20(m)
!      dgv_pftcon%tcmin(m) = pftpar28(m)
!      dgv_pftcon%tcmax(m) = pftpar29(m)
!      dgv_pftcon%gddmin(m) = pftpar30(m)
!      dgv_pftcon%twmax(m) = pftpar31(m)
!      dgv_pftcon%reinickerp(m) = reinickerp
!      dgv_pftcon%allom1(m) = allom1
!      dgv_pftcon%allom2(m) = allom2
!      dgv_pftcon%allom3(m) = allom3
      ! modification for shrubs by X.D.Z
!      if (m > ntree .and. m <= nbrdlf_dcd_brl_shrub ) then 
!         dgv_pftcon%allom1(m) = allom1s
!         dgv_pftcon%allom2(m) = allom2s
!      end if
!   end do
!#endif

   ! --------------------------------------------------------------------
   ! Initialize nitrogen deposition values 
   ! for now these are constants by gridcell, eventually they
   ! will be variables from the atmosphere, and at some point in between
   ! they will be specified time varying fields.
   ! --------------------------------------------------------------------

   ! Grid level initialization
!   do g = begg, endg

      ! VOC emission factors
      ! Set gridcell and landunit indices
!      efisop(:,g)=efisop2d(:,g)

!   end do


   ! pft level initialization
   do p = begp, endp

      ! Initialize root fraction (computing from surface, d is depth in meter):
      ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
      ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with
      ! beta & d_obs given in Zeng et al. (1998).

      c = pcolumn(p)
!      if (ivt(p) /= noveg) then
!         do lev = 1, nlevgrnd
!            rootfr(p,lev) = 0._r8
!         enddo
!         do lev = 1, nlevsoi-1
!            rootfr(p,lev) = .5_r8*( exp(-roota_par(ivt(p)) * zi(c,lev-1))  &
!                               + exp(-rootb_par(ivt(p)) * zi(c,lev-1))  &
!                               - exp(-roota_par(ivt(p)) * zi(c,lev  ))  &
!                               - exp(-rootb_par(ivt(p)) * zi(c,lev  )) )
!         end do
!         rootfr(p,nlevsoi) = .5_r8*( exp(-roota_par(ivt(p)) * zi(c,nlevsoi-1))  &
!                                + exp(-rootb_par(ivt(p)) * zi(c,nlevsoi-1)) )
!         rootfr(p,nlevsoi+1:nlevgrnd) =  0.0_r8

! Already commented out.
!#if (defined CN)
!        ! replacing the exponential rooting distribution
!        ! with a linear decrease, going to zero at the bottom of the lowest
!        ! soil layer for woody pfts, but going to zero at the bottom of
!        ! layer 8 for non-woody pfts.  This corresponds to 3.43 m for woody
!        ! bottom, vs 1.38 m for non-woody bottom.
!        if (woody(ivt(p)) == 1) then
!           bottom = nlevsoi
!           slope = -2._r8/(zi(c,bottom)*zi(c,bottom))
!           intercept   = 2._r8/zi(c,bottom)
!           do lev = 1, bottom
!              rootfr(p,lev) = dz(c,lev) * 0.5_r8 * ((intercept+slope*zi(c,lev-1)) + (intercept+slope*zi(c,lev)))
!           end do
!           if (bottom < nlevsoi) then
!              do lev=bottom+1,nlevgrnd
!                 rootfr(p,lev) = 0._r8
!              end do
!           end if
!        else
!           bottom = 8
!           slope = -2._r8/(zi(c,bottom)*zi(c,bottom))
!           intercept   = 2._r8/zi(c,bottom)
!           do lev=1,bottom
!              rootfr(p,lev) = dz(c,lev) * 0.5_r8 * ((intercept+slope*zi(c,lev-1)) + (intercept+slope*zi(c,lev)))
!           end do
!           if (bottom < nlevsoi) then
!              do lev=bottom+1,nlevgrnd
!                 rootfr(p,lev) = 0._r8
!              end do
!           end if
!        end if
!#endif
!      else
!         rootfr(p,1:nlevsoi) = 0._r8
!      endif
      
      ! initialize rresis, for use in ecosystemdyn
      do lev = 1,nlevgrnd
         rresis(p,lev) = 0._r8
      end do

   end do ! end pft level initialization
   
   ! initialize the CN variables for special landunits, including lake points
   call CNiniSpecial(nlevgrnd, begg, endg, begc, endc, begp, endp)

end subroutine iniTimeConst

end module iniTimeConstMod
