module pftvarcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: pftvarcon
!
! !DESCRIPTION:
! Module containing vegetation constants and method to
! read and initialize vegetation (PFT) constants.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar  , only : mxpft, numpft, numrad, ivis, inir
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Vegetation type constants
!
  character(len=40) pftname(0:mxpft) !PFT description

  integer :: noveg                  !value for not vegetated 
  integer :: ndllf_evr_tmp_tree     !value for Needleleaf evergreen temperate tree
  integer :: ndllf_evr_brl_tree     !value for Needleleaf evergreen boreal tree
  integer :: ndllf_dcd_brl_tree     !value for Needleleaf deciduous boreal tree
  integer :: nbrdlf_evr_trp_tree    !value for Broadleaf evergreen tropical tree
  integer :: nbrdlf_evr_tmp_tree    !value for Broadleaf evergreen temperate tree
  integer :: nbrdlf_dcd_trp_tree    !value for Broadleaf deciduous tropical tree
  integer :: nbrdlf_dcd_tmp_tree    !value for Broadleaf deciduous temperate tree
  integer :: nbrdlf_dcd_brl_tree    !value for Broadleaf deciduous boreal tree
  integer :: ntree                  !value for last type of tree
  integer :: nbrdlf_evr_shrub       !value for Broadleaf evergreen shrub
  integer :: nbrdlf_dcd_tmp_shrub   !value for Broadleaf deciduous temperate shrub
  integer :: nbrdlf_dcd_brl_shrub   !value for Broadleaf deciduous boreal shrub
  integer :: nc3_arctic_grass       !value for C3 arctic grass
  integer :: nc3_nonarctic_grass    !value for C3 non-arctic grass
  integer :: nc4_grass              !value for C4 grass
  integer :: npcropmin              !value for first crop
  integer :: ncorn                  !value for corn
  integer :: nscereal               !value for spring temperate cereal
  integer :: nwcereal               !value for winter temperate cereal
  integer :: nsoybean               !value for soybean
  integer :: npcropmax              !value for last prognostic crop in list
  integer :: nc3crop                !value for generic crop
  integer :: nirrig                 !value for irrigated generic crop

  real(r8):: dleaf(0:mxpft)       !characteristic leaf dimension (m)
  real(r8):: c3psn(0:mxpft)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(r8):: mp(0:mxpft)          !slope of conductance-to-photosynthesis relationship
  real(r8):: qe25(0:mxpft)        !quantum efficiency at 25C (umol CO2 / umol photon)
  real(r8):: xl(0:mxpft)          !leaf/stem orientation index
  real(r8):: rhol(0:mxpft,numrad) !leaf reflectance: 1=vis, 2=nir
  real(r8):: rhos(0:mxpft,numrad) !stem reflectance: 1=vis, 2=nir
  real(r8):: taul(0:mxpft,numrad) !leaf transmittance: 1=vis, 2=nir
  real(r8):: taus(0:mxpft,numrad) !stem transmittance: 1=vis, 2=nir
  real(r8):: z0mr(0:mxpft)        !ratio of momentum roughness length to canopy top height (-)
  real(r8):: displar(0:mxpft)     !ratio of displacement height to canopy top height (-)
  real(r8):: roota_par(0:mxpft)   !CLM rooting distribution parameter [1/m]
  real(r8):: rootb_par(0:mxpft)   !CLM rooting distribution parameter [1/m]
  real(r8):: crop(0:mxpft)        ! crop pft: 0. = not crop, 1. = crop pft
  real(r8):: irrigated(0:mxpft)   ! irrigated pft: 0. = not, 1. = irrigated
  real(r8):: smpso(0:mxpft)       !soil water potential at full stomatal opening (mm)
  real(r8):: smpsc(0:mxpft)       !soil water potential at full stomatal closure (mm)
  real(r8):: fnitr(0:mxpft)       !foliage nitrogen limitation factor (-)
  ! begin new pft parameters for CN code
  real(r8):: slatop(0:mxpft)      !SLA at top of canopy [m^2/gC]
  real(r8):: dsladlai(0:mxpft)    !dSLA/dLAI [m^2/gC]
  real(r8):: leafcn(0:mxpft)      !leaf C:N [gC/gN]
  real(r8):: flnr(0:mxpft)        !fraction of leaf N in Rubisco [no units]
  real(r8):: woody(0:mxpft)       !woody lifeform flag (0 or 1)
  real(r8):: lflitcn(0:mxpft)      !leaf litter C:N (gC/gN)
  real(r8):: frootcn(0:mxpft)      !fine root C:N (gC/gN)
  real(r8):: livewdcn(0:mxpft)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(r8):: deadwdcn(0:mxpft)     !dead wood (xylem and heartwood) C:N (gC/gN)
  real(r8):: grperc(0:mxpft)       !growth respiration parameter
  real(r8):: grpnow(0:mxpft)       !growth respiration parameter

! for crop
  real(r8):: graincn(0:mxpft)      !grain C:N (gC/gN)
  real(r8):: mxtmp(0:mxpft)        !parameter used in accFlds
  real(r8):: baset(0:mxpft)        !parameter used in accFlds
  real(r8):: declfact(0:mxpft)     !parameter used in CNAllocation
  real(r8):: bfact(0:mxpft)        !parameter used in CNAllocation
  real(r8):: aleaff(0:mxpft)       !parameter used in CNAllocation
  real(r8):: arootf(0:mxpft)       !parameter used in CNAllocation
  real(r8):: astemf(0:mxpft)       !parameter used in CNAllocation
  real(r8):: arooti(0:mxpft)       !parameter used in CNAllocation
  real(r8):: fleafi(0:mxpft)       !parameter used in CNAllocation
  real(r8):: allconsl(0:mxpft)     !parameter used in CNAllocation
  real(r8):: allconss(0:mxpft)     !parameter used in CNAllocation
  real(r8):: ztopmx(0:mxpft)       !parameter used in CNVegStructUpdate
  real(r8):: laimx(0:mxpft)        !parameter used in CNVegStructUpdate
  real(r8):: gddmin(0:mxpft)       !parameter used in CNPhenology
  real(r8):: hybgdd(0:mxpft)       !parameter used in CNPhenology
  real(r8):: lfemerg(0:mxpft)      !parameter used in CNPhenology
  real(r8):: grnfill(0:mxpft)      !parameter used in CNPhenology
  integer :: mxmat(0:mxpft)        !parameter used in CNPhenology
  integer :: mnNHplantdate(0:mxpft)!minimum planting date for NorthHemisphere (YYYYMMDD)
  integer :: mxNHplantdate(0:mxpft)!maximum planting date for NorthHemisphere (YYYYMMDD)
  integer :: mnSHplantdate(0:mxpft)!minimum planting date for SouthHemisphere (YYYYMMDD)
  integer :: mxSHplantdate(0:mxpft)!maximum planting date for SouthHemisphere (YYYYMMDD)
  real(r8):: planttemp(0:mxpft)    !planting temperature used in CNPhenology (K)
  real(r8):: minplanttemp(0:mxpft) !mininum planting temperature used in CNPhenology (K)
  real(r8):: froot_leaf(0:mxpft)   !allocation parameter: new fine root C per new leaf C (gC/gC) 
  real(r8):: stem_leaf(0:mxpft)    !allocation parameter: new stem c per new leaf C (gC/gC)
  real(r8):: croot_stem(0:mxpft)   !allocation parameter: new coarse root C per new stem C (gC/gC)
  real(r8):: flivewd(0:mxpft)      !allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
  real(r8):: fcur(0:mxpft)         !allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
  real(r8):: fcurdv(0:mxpft)       !alternate fcur for use with cndv
  real(r8):: lf_flab(0:mxpft)      !leaf litter labile fraction
  real(r8):: lf_fcel(0:mxpft)      !leaf litter cellulose fraction
  real(r8):: lf_flig(0:mxpft)      !leaf litter lignin fraction
  real(r8):: fr_flab(0:mxpft)      !fine root litter labile fraction
  real(r8):: fr_fcel(0:mxpft)      !fine root litter cellulose fraction
  real(r8):: fr_flig(0:mxpft)      !fine root litter lignin fraction
  real(r8):: leaf_long(0:mxpft)    !leaf longevity (yrs)
  real(r8):: evergreen(0:mxpft)    !binary flag for evergreen leaf habit (0 or 1)
  real(r8):: stress_decid(0:mxpft) !binary flag for stress-deciduous leaf habit (0 or 1)
  real(r8):: season_decid(0:mxpft) !binary flag for seasonal-deciduous leaf habit (0 or 1)
  real(r8):: pconv(0:mxpft)        !proportion of deadstem to conversion flux
  real(r8):: pprod10(0:mxpft)      !proportion of deadstem to 10-yr product pool
  real(r8):: pprod100(0:mxpft)     !proportion of deadstem to 100-yr product pool
  real(r8):: pprodharv10(0:mxpft)  !harvest mortality proportion of deadstem to 10-yr pool

  ! new pft parameters for CN-fire code
  real(r8):: resist(0:mxpft)       !resistance to fire (no units)

  ! pft parameters for CNDV code
  ! from LPJ subroutine pftparameters
  real(r8) pftpar20(0:mxpft)       !tree maximum crown area (m2)
  real(r8) pftpar28(0:mxpft)       !min coldest monthly mean temperature
  real(r8) pftpar29(0:mxpft)       !max coldest monthly mean temperature
  real(r8) pftpar30(0:mxpft)       !min growing degree days (>= 5 deg C)
  real(r8) pftpar31(0:mxpft)       !upper limit of temperature of the warmest month (twmax)
  real(r8), parameter :: reinickerp = 1.6_r8 !parameter in allometric equation
  real(r8), parameter :: dwood  = 2.5e5_r8   !cn wood density (gC/m3); lpj:2.0e5
  real(r8), parameter :: allom1 = 100.0_r8   !parameters in
  real(r8), parameter :: allom2 =  40.0_r8   !...allometric
  real(r8), parameter :: allom3 =   0.5_r8   !...equations
  real(r8), parameter :: allom1s = 250.0_r8  !modified for shrubs by
  real(r8), parameter :: allom2s =   8.0_r8  !X.D.Z
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: pftconrd ! Read and initialize vegetation (PFT) constants
!
! !REVISION HISTORY:
! Created by Sam Levis (put into module form by Mariana Vertenstein)
! 10/21/03, Peter Thornton: Added new variables for CN code
! 06/24/09, Erik Kluzek: Add indices for all pft types, and add expected_pftnames array and comparision
! 09/17/10, David Lawrence: Modified code to read in netCDF pft physiology file
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pftconrd
!
! !INTERFACE:
  subroutine pftconrd
!
! !DESCRIPTION:
! Read and initialize vegetation (PFT) constants
!
! !USES:
    use clm_varcon, only : tfrz
    use nanMod    , only : nan
    use clmtype   , only : pftcon
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Gordon Bonan
! 3/26/13: Adapted for use in VIC by Michael Brunke
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256) :: locfn ! local file name
    integer :: i,n              ! loop indices
    integer :: npft             ! number of pfts on pft-physiology file
    integer :: pftnum           ! pft number
    character(len=32) :: subname = 'pftconrd'              ! subroutine name
    !
    ! Expected PFT names: The names expected on the fpftcon file and the order they are expected to be in.
    ! NOTE: similar types are assumed to be together, first trees (ending with broadleaf_deciduous_boreal_tree
    !       then shrubs, ending with broadleaf_deciduous_boreal_shrub, then grasses starting with c3_arctic_grass
    !       and finally crops, ending with soybean
    ! DO NOT CHANGE THE ORDER -- WITHOUT MODIFYING OTHER PARTS OF THE CODE WHERE THE ORDER MATTERS!
    !
    character(len=40), parameter :: expected_pftnames(0:mxpft) = (/ &
                 'not_vegetated                      '  &
               , 'needleleaf_evergreen_temperate_tree'  &
               , 'needleleaf_evergreen_boreal_tree   '  &
               , 'needleleaf_deciduous_boreal_tree   '  &
               , 'broadleaf_evergreen_tropical_tree  '  &
               , 'broadleaf_evergreen_temperate_tree '  &
               , 'broadleaf_deciduous_tropical_tree  '  &
               , 'broadleaf_deciduous_temperate_tree '  &
               , 'broadleaf_deciduous_boreal_tree    '  &
               , 'broadleaf_evergreen_shrub          '  &
               , 'broadleaf_deciduous_temperate_shrub'  &
               , 'broadleaf_deciduous_boreal_shrub   '  &
               , 'c3_arctic_grass                    '  &
               , 'c3_non-arctic_grass                '  &
               , 'c4_grass                           '  &
               , 'c3_crop                            '  &
               , 'c3_irrigated                       '  &
               , 'corn                               '  &
               , 'spring_temperate_cereal            '  &
               , 'winter_temperate_cereal            '  &
               , 'soybean                            '  &
    /)
!-----------------------------------------------------------------------

    ! Set specific vegetation type values

    locfn = 'pft-physiology.c110425'
    open(10, file=locfn)
    read(10, '(I2)') npft

    do n = 1, npft
      read(10, '(I2,1X,A40,80F12.3,4I5)') pftnum, pftname(n-1), z0mr(n-1), &
	displar(n-1), dleaf(n-1), c3psn(n-1), mp(n-1), qe25(n-1), &
	rhol(n-1,ivis), rhol(n-1,inir), rhos(n-1,ivis), rhos(n-1,inir), &
	taul(n-1,ivis), taul(n-1,inir), taus(n-1,ivis), taus(n-1,inir), &
	xl(n-1), roota_par(n-1), rootb_par(n-1), slatop(n-1), dsladlai(n-1), &
	leafcn(n-1), flnr(n-1), smpso(n-1), smpsc(n-1), fnitr(n-1), &
	woody(n-1), lflitcn(n-1), frootcn(n-1), livewdcn(n-1), deadwdcn(n-1), &
	grperc(n-1), grpnow(n-1), froot_leaf(n-1), stem_leaf(n-1), &
	croot_stem(n-1), flivewd(n-1), fcur(n-1), fcurdv(n-1), lf_flab(n-1), &
	lf_fcel(n-1), lf_flig(n-1), fr_flab(n-1), fr_fcel(n-1), fr_flig(n-1), &
	leaf_long(n-1), evergreen(n-1), stress_decid(n-1), season_decid(n-1), &
	resist(n-1), pftpar20(n-1), pftpar28(n-1), pftpar29(n-1), &
	pftpar30(n-1), pftpar31(n-1), pconv(n-1), pprod10(n-1), &
	pprodharv10(n-1), pprod100(n-1), graincn(n-1), mxtmp(n-1), &
	baset(n-1), declfact(n-1), bfact(n-1), aleaff(n-1), arootf(n-1), &
	astemf(n-1), arooti(n-1), fleafi(n-1), allconsl(n-1), allconss(n-1), &
	crop(n-1), irrigated(n-1), ztopmx(n-1), laimx(n-1), gddmin(n-1), &
	hybgdd(n-1), lfemerg(n-1), grnfill(n-1), mxmat(n-1), planttemp(n-1), &
	minplanttemp(n-1), mnNHplantdate(n-1), mnSHplantdate(n-1), &
	mxNHplantdate(n-1), mxSHplantdate(n-1)

!      pftcon%z0mr(n-1) = z0mr(n-1)
!      pftcon%displar(n-1) = displar(n-1)
!      pftcon%dleaf(n-1) = dleaf(n-1)
!      pftcon%c3psn(n-1) = c3psn(n-1)
!      pftcon%mp(n-1) = mp(n-1)
!      pftcon%qe25(n-1) = qe25(n-1)
!      do i = 1, 2
!        pftcon%rhol(n-1,i) = rhol(n-1,i)
!        pftcon%rhos(n-1,i) = rhos(n-1,i)
!        pftcon%taul(n-1,i) = taul(n-1,i)
!        pftcon%taus(n-1,i) = taus(n-1,i)
!      end do
!      pftcon%xl(n-1) = xl(n-1)
!      pftcon%roota_par(n-1) = roota_par(n-1)
!      pftcon%rootb_par(n-1) = rootb_par(n-1)
!      pftcon%slatop(n-1) = slatop(n-1)
!      pftcon%dsladlai(n-1) = dsladlai(n-1)
!      pftcon%leafcn(n-1) = leafcn(n-1)
!      pftcon%flnr(n-1) = flnr(n-1)
!      pftcon%smpso(n-1) = smpso(n-1)
!      pftcon%smpsc(n-1) = smpsc(n-1)
!      pftcon%fnitr(n-1) = fnitr(n-1)
!      pftcon%woody(n-1) = woody(n-1)
!      pftcon%lflitcn(n-1) = lflitcn(n-1)
!      pftcon%frootcn(n-1) = frootcn(n-1)
!      pftcon%livewdcn(n-1) = livewdcn(n-1)
!      pftcon%deadwdcn(n-1) = deadwdcn(n-1)
!      pftcon%froot_leaf(n-1) = froot_leaf(n-1)
!      pftcon%stem_leaf(n-1) = stem_leaf(n-1)
!      pftcon%croot_stem(n-1) = croot_stem(n-1)
!      pftcon%flivewd(n-1) = flivewd(n-1)
!      pftcon%fcur(n-1) = fcur(n-1)
!      pftcon%lf_flab(n-1) = lf_flab(n-1)
!      pftcon%lf_fcel(n-1) = lf_fcel(n-1)
!      pftcon%lf_flig(n-1) = lf_flig(n-1)
!      pftcon%fr_flab(n-1) = fr_flab(n-1)
!      pftcon%fr_fcel(n-1) = fr_fcel(n-1)
!      pftcon%fr_flig(n-1) = fr_flig(n-1)
!      pftcon%leaf_long(n-1) = leaf_long(n-1)
!      pftcon%evergreen(n-1) = evergreen(n-1)
!      pftcon%stress_decid(n-1) = stress_decid(n-1)
!      pftcon%season_decid(n-1) = season_decid(n-1)
!      pftcon%resist(n-1) = resist(n-1) 
    enddo
    close(10)

    do i = 0, mxpft
       if (pftname(i) /= expected_pftnames(i)) then
          write(*,*)'pftconrd: pftname is NOT what is expected, name = ', &
                        pftname(i), ', expected name = ', expected_pftnames(i)
          return
       end if
       if (pftname(i) == 'not_vegetated'                       ) noveg               = i
       if (pftname(i) == 'needleleaf_evergreen_temperate_tree' ) ndllf_evr_tmp_tree  = i
       if (pftname(i) == 'needleleaf_evergreen_boreal_tree'    ) ndllf_evr_brl_tree  = i
       if (pftname(i) == 'needleleaf_deciduous_boreal_tree'    ) ndllf_dcd_brl_tree  = i
       if (pftname(i) == 'broadleaf_evergreen_tropical_tree'   ) nbrdlf_evr_trp_tree  = i
       if (pftname(i) == 'broadleaf_evergreen_temperate_tree'  ) nbrdlf_evr_tmp_tree  = i
       if (pftname(i) == 'broadleaf_deciduous_tropical_tree'   ) nbrdlf_dcd_trp_tree  = i
       if (pftname(i) == 'broadleaf_deciduous_temperate_tree'  ) nbrdlf_dcd_tmp_tree  = i
       if (pftname(i) == 'broadleaf_deciduous_boreal_tree'     ) nbrdlf_dcd_brl_tree  = i
       if (pftname(i) == 'broadleaf_evergreen_shrub'           ) nbrdlf_evr_shrub     = i
       if (pftname(i) == 'broadleaf_deciduous_temperate_shrub' ) nbrdlf_dcd_tmp_shrub = i
       if (pftname(i) == 'broadleaf_deciduous_boreal_shrub'    ) nbrdlf_dcd_brl_shrub = i
       if (pftname(i) == 'c3_arctic_grass'                     ) nc3_arctic_grass     = i
       if (pftname(i) == 'c3_non-arctic_grass'                 ) nc3_nonarctic_grass  = i
       if (pftname(i) == 'c4_grass'                            ) nc4_grass            = i
       if (pftname(i) == 'c3_crop'                             ) nc3crop              = i
       if (pftname(i) == 'c3_irrigated'                        ) nirrig               = i
       if (pftname(i) == 'corn'                                ) ncorn                = i
       if (pftname(i) == 'spring_temperate_cereal'             ) nscereal             = i
       if (pftname(i) == 'winter_temperate_cereal'             ) nwcereal             = i
       if (pftname(i) == 'soybean'                             ) nsoybean             = i
    end do

    ntree                = nbrdlf_dcd_brl_tree  ! value for last type of tree
    npcropmin            = ncorn                ! first prognostic crop
    npcropmax            = nsoybean             ! last prognostic crop in list

#if (defined CNDV)
    fcur(:) = fcurdv(:)
#endif

    !
    ! Do some error checking
    !
    if ( npcropmax /= mxpft )then
       write(*,*) 'pftconrd ERROR: npcropmax is NOT the last value'
       return
    end if
    do i = 0, mxpft
       if (     (irrigated(i) == 1.0_r8) .and. i == nirrig )then
          ! correct
       else if ( irrigated(i) == 0.0_r8 )then
          ! correct
       else
          write (*,*) 'pftconrd ERROR: irrigated has wrong values'
          return
       end if
       if (      crop(i) == 1.0_r8 .and. (i >= nc3crop .and. i <= npcropmax) )then
          ! correct
       else if ( crop(i) == 0.0_r8 )then
          ! correct
       else
          write(*,*) 'pftconrd ERROR: crop has wrong values'
          return
       end if
       if ( (i /= noveg) .and. (i < npcropmin) .and. &
            abs(pconv(i)+pprod10(i)+pprod100(i) - 1.0_r8) > 1.e-7_r8 )then
          write(*,*) 'pftconrd ERROR: pconv+pprod10+pprod100 do NOT sum to one.'
          return
       end if
       if ( pprodharv10(i) > 1.0_r8 .or. pprodharv10(i) < 0.0_r8 )then
          write(*,*) 'pftconrd ERROR: pprodharv10 outside of range.'
          return
       end if
    end do

  end subroutine pftconrd

end module pftvarcon

