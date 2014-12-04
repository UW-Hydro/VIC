  subroutine vic2clmtype(nlevgrnd, rec, nrec, adspinup, yr, mo, day, &
	secs, jday, yr1, jday1, dt, lat, lon, begg, endg, begc, endc, begp, &
	endp, num_soilc, num_soilp, psfc, tair, vp, vpd, lwrad, swrad, swrd, &
	swri, alb, dep, thick, baseflow, moist, ice, tsoisno, t2m, tveg, &
	snodep, fvegwet, rootf, satpsi, soipsi, coeff, rveg, zo, zos, zov, &
	displ, lai, soilcfast, soilcmid, soilcslo1, soilcslo2, litrlabc, &
	litrcellc, litrligc, cwoodc, leafcc, finrtc, livstemc, deadstemcc, &
	livcorsrtc, deadcorsrtc, woodcc, soilnfast, soilnmid, soilnslo1, &
	soilnslo2, soilminn, litrlabn, litrcelln, litrlign, cwoodn, leafnn, &
	finrtn, livstemn, deadstemnn, livcorsrtn, deadcorsrtn, vegctot, &
	litctot, somctot, prod1gb4, prod1g, prod1n, darkresp, maintresp, groresp, autoresp, heteroresp, &
	litresp, ecoexchn, prodecon, dormancy, ndaysact, onsetflg, onsetcount,&
	onsetgddflg, onsetfdd, onsetgdd, onsetswi, offsetflg, offsetcount, &
	offsetfdd, offsetswi, lgsfact, backlfr, backtgr, daylen, &
	prevdaylen, annavgt2m, tempavgt2m, cavail, cflxrecov, allocpnow, &
	callom, nallom, plantndem, tempsumpotgpp, annsumpotgpp, &
	tempmxretransn, annmxretransn, availretransn, plantnalloc, &
	plantcalloc, cflxex, dwnreg, prevleafc2litr, prevfinrtc2litr, &
	tempsumnpp, annsumnpp, leafcstor, leafctrans, finrtcstor, &
	finrtctrans, livstemcstor, livstemctrans, deadstemcstor, &
	deadstemctrans, livcorsrtcstor, livcorsrtctrans, deadcorsrtcstor, &
	deadcorsrtctrans, grorespstor, groresptrans, photocpool, mrcpool, &
	ctruncpft, leafnstor, leafntrans, finrtnstor, finrtntrans, &
	livstemnstor, livstemntrans, deadstemnstor, deadstemntrans, &
	livcorsrtnstor, livcorsrtntrans, deadcorsrtnstor, deadcorsrtntrans, &
	nretrans, photonpool, ntruncpft, declin, potimmobfract, potgppfract, &
	annsumcount, colannsumnpp, colannsumt2m, fldcapwat, extinctmoist, &
	fireprob, meanfireprob, fireseasnlen, areafractburn, &
	annareafractburn, seedcc, colctrunc, colctot, woodprodc10, &
	woodprodc100, seednn, colntrunc, colntot, woodprodn10, woodprodn100, &
	litrfall, photosynth, intco2, stomresist, abspar, init_state)

! 03/06/2014  Added sending of CN variables for inclusion in VIC 
!             state files.                                          MAB
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_const_mod, only: SHR_CONST_PI, SHR_CONST_TKFRZ
   use clmtype
   use clm_varpar,        only : max_bands, max_layers, max_nodes, max_veg, mxpft, numrad
   use shr_orb_mod
   use SurfaceAlbedoMod, only : SurfaceAlbedo
   use SurfaceRadiationMod, only : SurfaceRadiation
   use CanopyFluxesMod, only : CanopyFluxes
   use CNEcosystemDynMod, only : CNEcosystemDyn
   use CNAnnualUpdateMod, only : CNAnnualUpdate
   use CNBalanceCheckMod, only : CBalanceCheck, NBalanceCheck
!
! !ARGUMENTS
   implicit none

   integer, intent(in) :: nlevgrnd ! Number of below ground layers
   integer, intent(in) :: yr       ! Current year
   integer, intent(in) :: mo       ! Current month
   integer, intent(in) :: day      ! Current day
   integer, intent(in) :: secs     ! Current time of the day in seconds
   real(r8), intent(in) :: jday     ! Julian day
   integer, intent(in) :: yr1      ! Next time step year
   real(r8), intent(in) :: jday1   ! Next time step Julian day
   real(r8), intent(in) :: dt      ! Time step in seconds
   integer, intent(in) :: rec   ! Record number
   integer, intent(in) :: nrec  ! Total number of records
!   integer, intent(in) :: nspinup ! Year in spin-up
   integer, intent(in) :: adspinup ! Flag for AD spin-up
   real(r8), intent(in) :: lat  ! Latitude (degrees N)
   real(r8), intent(in) :: lon  ! Longitude (degrees E)
   integer, intent(in) :: begg	! Beginning grid cell index = 1
   integer, intent(in) :: endg  ! Ending grid cell index = 1
   integer, intent(inout) :: begc  ! Beginning precip. dist. index = 1
   integer, intent(inout) :: endc  ! Ending precip. dist. index = Ndist
   integer, intent(inout) :: begp  ! Beginning veg. type index = 1
   integer, intent(inout) :: endp  ! Ending veg. type index = Nveg
   integer, intent(in) :: num_soilc ! number of soil columns in filter
   integer, intent(in) :: num_soilp ! number of soil pfts in filter
   real(r8), intent(in) :: psfc     ! surface pressure (Pa)
   real(r8), intent(in) :: tair     ! atmospheric temperature (C)
   real(r8), intent(in) :: vp       ! atmos. vapor pressure (Pa)
   real(r8), intent(in) :: vpd      ! vapor pressure deficit (Pa)
   real(r8), intent(in) :: lwrad    ! longwave radiation (W/m^2)
   real(r8), intent(in) :: swrad    ! shortwave radiation (W/m^2)
   real(r8), intent(in) :: swrd(numrad)     ! direct SW radiation (W/m^2)
   real(r8), intent(in) :: swri(numrad)     ! diffuse SW radiation (W/m^2)
   real(r8), intent(in) :: alb(max_bands)    ! surface albedo
   real(r8), intent(in) :: dep(max_nodes+2,max_bands)  ! layer depth (m)
   real(r8), intent(in) :: thick(max_nodes+2,max_bands)! layer thickness (m)
   real(r8), intent(in) :: baseflow(max_bands)         ! baseflow
   real(r8), intent(in) :: moist(max_nodes+2,max_bands)! soil moisture
   real(r8), intent(in) :: ice(max_nodes+2,max_bands)  ! soil ice content
   real(r8), intent(in) :: tsoisno(max_nodes+2,max_bands) ! layer temperature
   real(r8), intent(in) :: t2m(mxpft)                 ! 2-m air temperature
   real(r8), intent(in) :: tveg(0:mxpft,max_bands)     ! veg temperature
   real(r8), intent(in) :: snodep(max_bands)           ! snow depth
   real(r8), intent(in) :: fvegwet(0:mxpft,max_bands)  ! fract of wet veg
   real(r8), intent(in) :: rootf(max_nodes,mxpft,max_bands) ! root fraction
   real(r8), intent(in) :: satpsi(max_nodes,max_bands,2) ! sat matric potential
   real(r8), intent(in) :: soipsi(max_nodes,max_bands,2) ! soil matric potential
   real(r8), intent(in) :: coeff(max_nodes,max_bands,2)  ! Clapp-Hornberger coef
   real(r8), intent(in) :: rveg(0:mxpft,max_bands)       ! leaf resistance
   real(r8), intent(in) :: zo        ! soil roughness length
   real(r8), intent(in) :: zos       ! snow roughness length
   real(r8), intent(in) :: zov(0:mxpft,max_bands)   ! veg. roughness length
   real(r8), intent(in) :: displ(0:mxpft,max_bands) ! displacement height
   real(r8), intent(inout) :: lai(0:mxpft,max_bands)! leaf area index
   real(r8), intent(inout) :: soilcfast(max_bands)  ! fast soil C pool
   real(r8), intent(inout) :: soilcmid(max_bands)   ! medium soil C pool
   real(r8), intent(inout) :: soilcslo1(max_bands)  ! slow soil C pool
   real(r8), intent(inout) :: soilcslo2(max_bands)  ! slowest soil C pool
   real(r8), intent(inout) :: litrlabc(max_bands)   ! litter labile C
   real(r8), intent(inout) :: litrcellc(max_bands)  ! litter cellulose C
   real(r8), intent(inout) :: litrligc(max_bands)   ! litter lignin C
   real(r8), intent(inout) :: cwoodc(max_bands)     ! course woody debris C
   real(r8), intent(inout) :: leafcc(0:mxpft,max_bands) ! leaf C
   real(r8), intent(inout) :: finrtc(0:mxpft,max_bands)       ! fine root C
   real(r8), intent(inout) :: livstemc(0:mxpft,max_bands)    ! live stem C
   real(r8), intent(inout) :: deadstemcc(0:mxpft,max_bands)    ! dead stem C
   real(r8), intent(inout) :: livcorsrtc(0:mxpft,max_bands) ! live coarse root C
   real(r8), intent(inout) :: deadcorsrtc(0:mxpft,max_bands)! dead coarse root C
   real(r8), intent(inout) :: woodcc(0:mxpft,max_bands)        ! wood C
   real(r8), intent(inout) :: soilnfast(max_bands)  ! fast soil N
   real(r8), intent(inout) :: soilnmid(max_bands)   ! medium soil N
   real(r8), intent(inout) :: soilnslo1(max_bands)  ! slow soil N
   real(r8), intent(inout) :: soilnslo2(max_bands)  ! slowest soil N
   real(r8), intent(inout) :: soilminn(max_bands)   ! soil mineral N
   real(r8), intent(inout) :: litrlabn(max_bands)   ! litter labile N
   real(r8), intent(inout) :: litrcelln(max_bands)  ! litter cellulose N
   real(r8), intent(inout) :: litrlign(max_bands)   ! litter lignin N
   real(r8), intent(inout) :: cwoodn(max_bands)     ! course woody debris N
   real(r8), intent(inout) :: leafnn(0:mxpft,max_bands) ! leaf N
   real(r8), intent(inout) :: finrtn(0:mxpft,max_bands)       ! fine root N
   real(r8), intent(inout) :: livstemn(0:mxpft,max_bands)    ! live stem N
   real(r8), intent(inout) :: deadstemnn(0:mxpft,max_bands)    ! dead stem N
   real(r8), intent(inout) :: livcorsrtn(0:mxpft,max_bands) ! live coarse root N
   real(r8), intent(inout) :: deadcorsrtn(0:mxpft,max_bands)! dead coarse root N
   real(r8), intent(inout) :: vegctot(0:mxpft,max_bands)    ! total veg C
   real(r8), intent(inout) :: litctot(max_bands)            ! total litter C
   real(r8), intent(inout) :: somctot(max_bands)            ! total SOM C
   real(r8), intent(inout) :: prod1gb4(0:mxpft,max_bands)   ! GPP before downreg
   real(r8), intent(inout) :: prod1g(0:mxpft,max_bands)     ! GPP
   real(r8), intent(inout) :: prod1n(0:mxpft,max_bands)     ! NPP
   real(r8), intent(inout) :: darkresp(0:mxpft,max_bands)   ! dark respiration
   real(r8), intent(inout) :: maintresp(0:mxpft,max_bands)  ! maintenance resp
   real(r8), intent(inout) :: groresp(0:mxpft,max_bands)    ! growth respiration
   real(r8), intent(inout) :: autoresp(0:mxpft,max_bands)   ! auto respiration
   real(r8), intent(inout) :: heteroresp(max_bands)         ! hetero respiration
   real(r8), intent(inout) :: litresp(max_bands)         ! litter respiration
   real(r8), intent(inout) :: ecoexchn(max_bands)                     ! NEE
   real(r8), intent(inout) :: prodecon(max_bands)                     ! NEP
   real(r8), intent(inout) :: dormancy(0:mxpft,max_bands)   ! dormancy
   real(r8), intent(inout) :: ndaysact(0:mxpft,max_bands)   ! # days since dormancy
   real(r8), intent(inout) :: onsetflg(0:mxpft,max_bands)   ! onset flag
   real(r8), intent(inout) :: onsetcount(0:mxpft,max_bands) ! onset counter
   real(r8), intent(inout) :: onsetgddflg(0:mxpft,max_bands) ! onset flag for grow deg sum
   real(r8), intent(inout) :: onsetfdd(0:mxpft,max_bands)   ! onset freeze days
   real(r8), intent(inout) :: onsetgdd(0:mxpft,max_bands)   ! onset grow deg days
   real(r8), intent(inout) :: onsetswi(0:mxpft,max_bands)   ! onset soil water index
   real(r8), intent(inout) :: offsetflg(0:mxpft,max_bands)   ! # days since dormancy
   real(r8), intent(inout) :: offsetcount(0:mxpft,max_bands) ! onset counter
   real(r8), intent(inout) :: offsetfdd(0:mxpft,max_bands)   ! onset freeze days
   real(r8), intent(inout) :: offsetswi(0:mxpft,max_bands)   ! onset soil water index
   real(r8), intent(inout) :: lgsfact(0:mxpft,max_bands)    ! long grow season factor
   real(r8), intent(inout) :: backlfr(0:mxpft,max_bands)    ! background litterfall
   real(r8), intent(inout) :: backtgr(0:mxpft,max_bands)    ! background transfer growth
   real(r8), intent(inout) :: daylen(0:mxpft,max_bands)     ! daylength
   real(r8), intent(inout) :: prevdaylen(0:mxpft,max_bands) ! previous daylength
   real(r8), intent(inout) :: annavgt2m(0:mxpft,max_bands)  ! annual avg 2-m air temp.
   real(r8), intent(inout) :: tempavgt2m(0:mxpft,max_bands) ! temp avg 2-m air temp.
   real(r8), intent(inout) :: cavail(0:mxpft,max_bands)     ! C flx avail for allocation
   real(r8), intent(inout) :: cflxrecov(0:mxpft,max_bands)  ! C flx for recovery
   real(r8), intent(inout) :: allocpnow(0:mxpft,max_bands)  ! Alloc fract for new growth
   real(r8), intent(inout) :: callom(0:mxpft,max_bands)     ! C allocation index
   real(r8), intent(inout) :: nallom(0:mxpft,max_bands)     ! N allocation index
   real(r8), intent(inout) :: plantndem(0:mxpft,max_bands)  ! plant N demand
   real(r8), intent(inout) :: tempsumpotgpp(0:mxpft,max_bands) ! temp sum pot GPP
   real(r8), intent(inout) :: annsumpotgpp(0:mxpft,max_bands) ! ann sum pot GPP
   real(r8), intent(inout) :: tempmxretransn(0:mxpft,max_bands) ! temp max retranslocated N
   real(r8), intent(inout) :: annmxretransn(0:mxpft,max_bands) ! ann max retranslocated N
   real(r8), intent(inout) :: availretransn(0:mxpft,max_bands) ! avail retranslocated N
   real(r8), intent(inout) :: plantnalloc(0:mxpft,max_bands) ! allocated N flux
   real(r8), intent(inout) :: plantcalloc(0:mxpft,max_bands) ! allocated C flux
   real(r8), intent(inout) :: cflxex(0:mxpft,max_bands)      ! C flux not alloc
   real(r8), intent(inout) :: dwnreg(0:mxpft,max_bands)      ! GPP reduct from N limit
   real(r8), intent(inout) :: prevleafc2litr(0:mxpft,max_bands) ! prev leaf C to litter
   real(r8), intent(inout) :: prevfinrtc2litr(0:mxpft,max_bands) ! prev fine root C to litter
   real(r8), intent(inout) :: tempsumnpp(0:mxpft,max_bands)  ! temp ann sum NPP
   real(r8), intent(inout) :: annsumnpp(0:mxpft,max_bands)   ! ann sum NPP
   real(r8), intent(inout) :: leafcstor(0:mxpft,max_bands)   ! leaf C storage
   real(r8), intent(inout) :: leafctrans(0:mxpft,max_bands)  ! leaf C transfer
   real(r8), intent(inout) :: finrtcstor(0:mxpft,max_bands)  ! fine root C storage
   real(r8), intent(inout) :: finrtctrans(0:mxpft,max_bands) ! fine root C transfer
   real(r8), intent(inout) :: livstemcstor(0:mxpft,max_bands) ! fine root C storage
   real(r8), intent(inout) :: livstemctrans(0:mxpft,max_bands) ! fine root C transfer
   real(r8), intent(inout) :: deadstemcstor(0:mxpft,max_bands) ! fine root C storage
   real(r8), intent(inout) :: deadstemctrans(0:mxpft,max_bands) ! fine root C transfer
   real(r8), intent(inout) :: livcorsrtcstor(0:mxpft,max_bands) ! live coarse root C storage
   real(r8), intent(inout) :: livcorsrtctrans(0:mxpft,max_bands) ! live coarse root C transfer
   real(r8), intent(inout) :: deadcorsrtcstor(0:mxpft,max_bands) ! dead coarse root C storage
   real(r8), intent(inout) :: deadcorsrtctrans(0:mxpft,max_bands) ! dead coarse root C transfer
   real(r8), intent(inout) :: grorespstor(0:mxpft,max_bands)  ! growth respiration storage
   real(r8), intent(inout) :: groresptrans(0:mxpft,max_bands) ! growth respiration transfer
   real(r8), intent(inout) :: photocpool(0:mxpft,max_bands)   ! photosynthate C pool
   real(r8), intent(inout) :: mrcpool(0:mxpft,max_bands)      ! C pool for excess MR demand
   real(r8), intent(inout) :: ctruncpft(0:mxpft,max_bands)    ! PFT sink for C truncation
   real(r8), intent(inout) :: leafnstor(0:mxpft,max_bands)   ! leaf N storage
   real(r8), intent(inout) :: leafntrans(0:mxpft,max_bands)  ! leaf N transfer
   real(r8), intent(inout) :: finrtnstor(0:mxpft,max_bands)  ! fine root N storage
   real(r8), intent(inout) :: finrtntrans(0:mxpft,max_bands) ! fine root N transfer
   real(r8), intent(inout) :: livstemnstor(0:mxpft,max_bands) ! fine root N storage
   real(r8), intent(inout) :: livstemntrans(0:mxpft,max_bands) ! fine root N transfer
   real(r8), intent(inout) :: deadstemnstor(0:mxpft,max_bands) ! dead stem N storage
   real(r8), intent(inout) :: deadstemntrans(0:mxpft,max_bands) ! dead stem N transfer
   real(r8), intent(inout) :: livcorsrtnstor(0:mxpft,max_bands) ! live coarse root N storage
   real(r8), intent(inout) :: livcorsrtntrans(0:mxpft,max_bands) ! live coarse root N transfer
   real(r8), intent(inout) :: deadcorsrtnstor(0:mxpft,max_bands) ! dead coarse root N storage
   real(r8), intent(inout) :: deadcorsrtntrans(0:mxpft,max_bands) ! dead coarse root N transfer
   real(r8), intent(inout) :: nretrans(0:mxpft,max_bands)     ! retranslocated N
   real(r8), intent(inout) :: photonpool(0:mxpft,max_bands)   ! photosynthate N
   real(r8), intent(inout) :: ntruncpft(0:mxpft,max_bands)    ! PFT sink for N truncation
   real(r8), intent(inout) :: declin(max_bands)              ! solar declination
   real(r8), intent(inout) :: potimmobfract(max_bands)        ! pot. immobilization fract
   real(r8), intent(inout) :: potgppfract(max_bands)          ! pot GPP fraction
   real(r8), intent(inout) :: annsumcount(max_bands)          ! ann sum counter
   real(r8), intent(inout) :: colannsumnpp(max_bands)         ! ann sum of NPP
   real(r8), intent(inout) :: colannsumt2m(max_bands)         ! ann sum of 2-m air temp
   real(r8), intent(inout) :: fldcapwat(max_nodes,max_bands)  ! vol soil wat @ field capacity
   real(r8), intent(inout) :: extinctmoist(max_bands)         ! moisture of extinction
   real(r8), intent(inout) :: fireprob(max_bands)             ! fire probability
   real(r8), intent(inout) :: meanfireprob(max_bands)         ! mean fire prob
   real(r8), intent(inout) :: fireseasnlen(max_bands)         ! fire season length
   real(r8), intent(inout) :: areafractburn(max_bands)       ! fract area burned
   real(r8), intent(inout) :: annareafractburn(max_bands)     ! ann fract area burned
   real(r8), intent(inout) :: seedcc(max_bands)               ! C pool for seeding new PFTs
   real(r8), intent(inout) :: colctrunc(max_bands)            ! col-lev C sink for truncation
   real(r8), intent(inout) :: colctot(max_bands)              ! total column C
   real(r8), intent(inout) :: woodprodc10(max_bands)          ! 10-yr lifespan wood product C
   real(r8), intent(inout) :: woodprodc100(max_bands)         ! 100-yr lifespan wood product C
   real(r8), intent(inout) :: seednn(max_bands)               ! N pool for seeding new PFTs
   real(r8), intent(inout) :: colntrunc(max_bands)            ! col-lev N sink for truncation
   real(r8), intent(inout) :: colntot(max_bands)              ! total column N
   real(r8), intent(inout) :: woodprodn10(max_bands)          ! 10-yr lifespan wood product N
   real(r8), intent(inout) :: woodprodn100(max_bands)         ! 100-yr lifespan wood product N
   real(r8), intent(inout) :: photosynth(0:mxpft,max_bands) ! photosynthesis
   real(r8), intent(inout) :: litrfall(0:mxpft,max_bands)   ! litterfall
   real(r8), intent(inout) :: intco2(0:mxpft,max_bands)     ! intracellular CO2
   real(r8), intent(inout) :: stomresist(0:mxpft,max_bands) ! stomatal resistance
   real(r8), intent(inout) :: abspar(0:mxpft,max_bands)     ! absorbed PAR
   integer , intent(in)    :: init_state                     ! flag to determine where initial state comes from

   real(r8), pointer :: z(:,:)          ! layer depth (m)
   real(r8), pointer :: dz(:,:)         ! layer thickness (m)
   real(r8), pointer :: qflx_drain(:)   ! sub-surface runoff (mm H20/s)
   real(r8), pointer :: h2osoi_liq(:,:) ! liquid water (kg / m2)
   real(r8), pointer :: h2osoi_ice(:,:) ! ice lens (kg / m2)
   real(r8), pointer :: t_soisno(:,:)     ! soil temperature (K)
   real(r8), pointer :: t_ref2m(:)      ! 2 m height temperature (K) 
   real(r8), pointer :: t_veg(:)        ! vegetation temperature (K)
   real(r8), pointer :: fwet(:)       ! wet vegetation fraction
   real(r8), pointer :: snowdp(:)       ! snow depth (m)
   real(r8), pointer :: rootfr(:,:)     ! fraction of roots in ea. soil layer
   real(r8), pointer :: psisat(:,:)     ! soil water at saturation (MPa) for CN
   real(r8), pointer :: soilpsi(:,:)    ! soil water potential (MPa) for CN
   real(r8), pointer :: bsw2(:,:)        ! Clapp-Hornberger coefficient for CN
   real(r8), pointer :: bsw(:,:)        ! Clapp-Hornberger coefficient
   real(r8), pointer :: sucsat(:,:)     ! minimum soil suction (mm)
   real(r8), pointer :: forc_pbot(:)    ! atmospheric pressure (Pa)
   real(r8), pointer :: forc_t(:)       ! atmospheric temperatuer (K)
   real(r8), pointer :: forc_q(:)       ! atmos. spec. hum. (kg/kg)
   real(r8), pointer :: forc_vp(:)      ! atmos. vapor pressure (Pa)
   real(r8), pointer :: qg(:)           ! surface spec. hum. (kg/kg)
   real(r8), pointer :: thm(:)          ! forc_t + 0.0098 * forc_hgt_t_pft
   real(r8), pointer :: forc_hgt_u(:)   ! wind forcing height (m)
   real(r8), pointer :: forc_hgt_t(:)   ! temperature forcing height (m)
   real(r8), pointer :: forc_hgt_q(:)   ! humidity forcing height (m)
   real(r8), pointer :: forc_hgt_u_pft(:) ! forc_hgt_u + z0m + d
   real(r8), pointer :: forc_hgt_t_pft(:) ! forc_hgt_t + z0m + d
   real(r8), pointer :: forc_hgt_q_pft(:) ! forc_hgt_q + z0m + d
   real(r8), pointer :: forc_pco2(:)    ! CO2 partial pressure (Pa)
   real(r8), pointer :: forc_po2(:)     ! O2 partial pressure (Pa)
   real(r8), pointer :: forc_lwrad(:)   ! longwave radiation (W/m^2)
   real(r8), pointer :: forc_solar(:)   ! shortwave radiation (W/m^2)
   real(r8), pointer :: forc_solad(:,:) ! direct SW radiation (W/m^2)
   real(r8), pointer :: forc_solai(:,:) ! diffuse SW radiation (W/m^2)
   real(r8), pointer :: soil1c(:)       ! soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)       ! soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)       ! soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)       ! soil organic matter C (slowest pool)
   real(r8), pointer :: litr1c(:)       ! litter labile C
   real(r8), pointer :: litr2c(:)       ! litter cellulose C
   real(r8), pointer :: litr3c(:)       ! litter lignin C
   real(r8), pointer :: cwdc(:)         ! coarse woody debris C
   real(r8), pointer :: leafc(:)        ! leaf C
   real(r8), pointer :: frootc(:)       ! fine root C
   real(r8), pointer :: livestemc(:)    ! live stem C
   real(r8), pointer :: deadstemc(:)    ! dead stem C
   real(r8), pointer :: livecrootc(:)   ! live coarse root C
   real(r8), pointer :: deadcrootc(:)   ! dead coarse root C
   real(r8), pointer :: woodc(:)        ! wood C
   real(r8), pointer :: totvegc(:)      ! total vegetation C (g C/m^2)
   real(r8), pointer :: totlitc(:)      ! total litter C (g C/m^2)
   real(r8), pointer :: totsomc(:)      ! total SOM C (g C/m^2)
   real(r8), pointer :: soil1n(:)       ! soil organic matter N (fast pool)
   real(r8), pointer :: soil2n(:)       ! soil organic matter N (medium pool)
   real(r8), pointer :: soil3n(:)       ! soil organic matter N (slow pool)
   real(r8), pointer :: soil4n(:)       ! soil organic matter N (slowest pool)
   real(r8), pointer :: sminn(:)        ! soil mineral N
   real(r8), pointer :: litr1n(:)       ! litter labile N
   real(r8), pointer :: litr2n(:)       ! litter cellulose N
   real(r8), pointer :: litr3n(:)       ! litter lignin N
   real(r8), pointer :: cwdn(:)         ! coarse woody debris N
   real(r8), pointer :: leafn(:)        ! leaf N
   real(r8), pointer :: frootn(:)       ! fine root N
   real(r8), pointer :: livestemn(:)    ! live stem N
   real(r8), pointer :: deadstemn(:)    ! dead stem N
   real(r8), pointer :: livecrootn(:)   ! live coarse root N
   real(r8), pointer :: deadcrootn(:)   ! dead coarse root N
   real(r8), pointer :: cpool_to_leafc(:)
   real(r8), pointer :: gpp2(:)         ! GPP before downregulation
   real(r8), pointer :: gpp(:)          ! GPP
   real(r8), pointer :: npp(:)          ! NPP
   real(r8), pointer :: fpsn(:)         ! photosynthesis
   real(r8), pointer :: leaf_mr(:)      ! leaf maintenance respiration
   real(r8), pointer :: mr(:)           ! maintenance respiration
   real(r8), pointer :: gr(:)           ! growth respiration
   real(r8), pointer :: ar(:)           ! autotrophic respiration
   real(r8), pointer :: hr(:)           ! heterotrophic respiration
   real(r8), pointer :: lithr(:)        ! litter hetero respiration
   real(r8), pointer :: nee(:)          ! NEE
   real(r8), pointer :: nep(:)          ! NEP
   real(r8), pointer :: dormant_flag(:) ! dormant flag
   real(r8), pointer :: days_active(:)  ! # days since dormancy
   real(r8), pointer :: onset_flag(:)   ! onset flag
   real(r8), pointer :: onset_counter(:) ! onset counter
   real(r8), pointer :: onset_gddflag(:) ! onset flag for growing deg days
   real(r8), pointer :: onset_fdd(:)    ! onset freezing deg days counter
   real(r8), pointer :: onset_gdd(:)    ! onset growing deg days counter
   real(r8), pointer :: onset_swi(:)    ! onset soil water index
   real(r8), pointer :: offset_flag(:)  ! offset flag
   real(r8), pointer :: offset_counter(:) ! onset counter
   real(r8), pointer :: offset_fdd(:)   ! onset freezing deg days counter
   real(r8), pointer :: offset_swi(:)   ! onset soil water index
   real(r8), pointer :: lgsf(:)         ! long growing season factor
   real(r8), pointer :: bglfr(:)        ! background litterfall rate
   real(r8), pointer :: bgtr(:)         ! background transfer growth rate
   real(r8), pointer :: dayl(:)         ! daylength
   real(r8), pointer :: prev_dayl(:)    ! previous timestep daylength
   real(r8), pointer :: annavg_t2m(:)   ! annual avg 2-m air temperature
   real(r8), pointer :: tempavg_t2m(:)  ! temp avg 2-ma air temperature
   real(r8), pointer :: availc(:)       ! C flux for allocation
   real(r8), pointer :: xsmrpool_recover(:) ! C flux for recovery
   real(r8), pointer :: alloc_pnow(:)   ! fract of alloc for new growth
   real(r8), pointer :: c_allometry(:)  ! C allocation index
   real(r8), pointer :: n_allometry(:)  ! N allocation index
   real(r8), pointer :: plant_ndemand(:) ! N flux for initial GPP
   real(r8), pointer :: tempsum_potential_gpp(:) ! temp ann sum of pot GPP
   real(r8), pointer :: annsum_potential_gpp(:) ! ann sum of potential GPP
   real(r8), pointer :: tempmax_retransn(:) ! temp ann max retransloc N
   real(r8), pointer :: annmax_retransn(:) ! ann max of retranslocated N
   real(r8), pointer :: avail_retransn(:) ! N flux avail for retranslocation
   real(r8), pointer :: plant_nalloc(:) ! total allocated N flux
   real(r8), pointer :: plant_calloc(:) ! total allocated C flux
   real(r8), pointer :: excess_cflux(:) ! C flux not allocated
   real(r8), pointer :: downreg(:)      ! fract reduct in GPP from N limit
   real(r8), pointer :: prev_leafc_to_litter(:) ! prev leaf C to litterfall
   real(r8), pointer :: prev_frootc_to_litter(:) ! prev fine root C to litter
   real(r8), pointer :: tempsum_npp(:)  ! temp ann sum of NPP
   real(r8), pointer :: annsum_npp(:)   ! annual sum of NPP
   real(r8), pointer :: leafc_storage(:) ! leaf C storage
   real(r8), pointer :: leafc_xfer(:)   ! leaf C transfer
   real(r8), pointer :: frootc_storage(:) ! fine root C storage
   real(r8), pointer :: frootc_xfer(:)  ! fine root C transfer
   real(r8), pointer :: livestemc_storage(:) ! live stem C storage
   real(r8), pointer :: livestemc_xfer(:) ! live stem C transfer
   real(r8), pointer :: deadstemc_storage(:) ! dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:) ! dead stem C transfer
   real(r8), pointer :: livecrootc_storage(:) ! live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:) ! live coarse root C transfer
   real(r8), pointer :: deadcrootc_storage(:) ! dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:) ! dead coarse root C transfer
   real(r8), pointer :: gresp_storage(:) ! growth respiration storage
   real(r8), pointer :: gresp_xfer(:)   ! growth respiration transfer
   real(r8), pointer :: cpool(:)        ! temp photosynthate C pool
   real(r8), pointer :: xsmrpool(:)     ! abstract C pool for excess MR demand
   real(r8), pointer :: pft_ctrunc(:)   ! PFT-lev sink for C truncation
   real(r8), pointer :: leafn_storage(:) ! leaf N storage
   real(r8), pointer :: leafn_xfer(:)   ! leaf N transfer
   real(r8), pointer :: frootn_storage(:) ! fine root N storage
   real(r8), pointer :: frootn_xfer(:)  ! fine root N transfer
   real(r8), pointer :: livestemn_storage(:) ! live stem N storage
   real(r8), pointer :: livestemn_xfer(:) ! live stem N transfer
   real(r8), pointer :: deadstemn_storage(:) ! dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:) ! dead stem N transfer
   real(r8), pointer :: livecrootn_storage(:) ! live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:) ! live coarse root N transfer
   real(r8), pointer :: deadcrootn_storage(:) ! dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:) ! dead coarse root N transfer
   real(r8), pointer :: retransn(:)     ! retranslocated N
   real(r8), pointer :: npool(:)        ! temp photosynthate N pool
   real(r8), pointer :: pft_ntrunc(:)   ! PFT-lev sink for N truncation
   real(r8), pointer :: fpi(:)          ! fract of potential immobilization
   real(r8), pointer :: fpg(:)          ! fract of potential GPP
   real(r8), pointer :: annsum_counter(:) ! ann sum turnover counter
   real(r8), pointer :: cannsum_npp(:)  ! col-avg ann sum of NPP
   real(r8), pointer :: cannsum_t2m(:)  ! col-avg ann sum of 2-m air temp
   real(r8), pointer :: watfc(:,:)      ! vol soil water at field capacity
   real(r8), pointer :: me(:)           ! moisture of extinction
   real(r8), pointer :: fire_prob(:)    ! daily fire probability
   real(r8), pointer :: mean_fire_prob(:) ! mean of daily fire probability
   real(r8), pointer :: fireseasonl(:)  ! fire season length
   real(r8), pointer :: farea_burned(:) ! fract area burned
   real(r8), pointer :: ann_farea_burned(:) ! ann sum of fract area burned
   real(r8), pointer :: seedc(:)        ! C pool for seeding new PFTs
   real(r8), pointer :: col_ctrunc(:)   ! col-lev sink for C truncation
   real(r8), pointer :: totcolc(:)      ! total column C
   real(r8), pointer :: prod10c(:)      ! 10-yr wood product C
   real(r8), pointer :: prod100c(:)     ! 100-yr wood product C
   real(r8), pointer :: seedn(:)        ! N pool for seeding new PFTs
   real(r8), pointer :: col_ntrunc(:)   ! col-lev sink for N truncation
   real(r8), pointer :: totcoln(:)      ! total column N
   real(r8), pointer :: prod10n(:)      ! 10-yr wood product N
   real(r8), pointer :: prod100n(:)     ! 100-yr wood product N
   real(r8), pointer :: tlai(:)         ! total leaf area index
   real(r8), pointer :: elai(:)         ! effective leaf area index
   real(r8), pointer :: litfall(:)      ! litterfall
   integer, pointer :: ivt(:)          ! PFT index
   integer, pointer :: ncol(:)          ! column index
   real(r8), pointer :: latrad(:)       ! latitude (radians)
   real(r8), pointer :: lonrad(:)       ! longitude (radians)
   real(r8), pointer :: latdeg(:)       ! latitude (degrees)
   real(r8), pointer :: londeg(:)       ! longitude (degrees)
   real(r8), pointer :: lat_a(:)        ! "atm" latitude (radians)
   real(r8), pointer :: lon_a(:)        ! "atm" longitude (radians)
   real(r8), pointer :: latdeg_a(:)     ! "atm" latitude (degrees)
   real(r8), pointer :: londeg_a(:)     ! "atm" longitude (degrees)
   real(r8), pointer :: coszen(:)       ! cosine of zenith angle
   real(r8), pointer :: decl(:)       ! declination angle (radians)
   real(r8), pointer :: albgrd(:,:)     ! direct surface reflectance
   real(r8), pointer :: albgri(:,:)     ! diffuse surface reflectance
   real(r8), pointer :: rb(:)           ! canopy resistance
   real(r8), pointer :: rssun(:)        ! sunlit stomatal resistance
   real(r8), pointer :: rssha(:)        ! shaded stomatal resistance
   real(r8), pointer :: cisun(:)        ! sunlit intracellular CO2
   real(r8), pointer :: cisha(:)        ! shaded intracellular CO2
   real(r8), pointer :: parsun(:)       ! sunlit absorbed PAR
   real(r8), pointer :: parsha(:)       ! shaded absorbed PAR
   real(r8), pointer :: laisun(:)       ! sunlit LAI
   real(r8), pointer :: laisha(:)       ! shaded LAI
   real(r8), pointer :: displa(:)       ! displacement height
   real(r8), pointer :: z0mv(:)         ! roughness length over vegetation
   real(r8), pointer :: z0hv(:)
   real(r8), pointer :: z0qv(:)
   real(r8), pointer :: z0mg(:)         ! roughness length over ground
   real(r8), pointer :: z0hg(:)
   real(r8), pointer :: z0qg(:)
   integer, pointer :: frac_veg_nosno(:)

   integer g, c, p, k, r, doy
   integer, parameter :: nlevsno = 2
   real(r8) :: vg                       ! surface vapor pressure (Pa)
   real(r8) :: eccen                    ! orbital eccentricity
   real(r8) :: obliq                    ! obliquity (degrees)
   real(r8) :: mvelp                    ! moving vernal equinox longitude
   real(r8) :: obliqr                   ! obliquity (radians)
   real(r8) :: lambm0                   ! mean lon of perihelion at vernal equinox (radians)
   real(r8) :: mvelpp                   ! moving vernal equinox lon + pi (rads)
   real(r8) :: delta                    ! declination angle (radians)
   real(r8) :: eccf                     ! Earth-sun distance factor
   real(r8) :: cosz                     ! cosine of zenith angle
   real(r8) :: eccen1                    ! orbital eccentricity
   real(r8) :: obliq1                    ! obliquity (degrees)
   real(r8) :: mvelp1                    ! moving vernal equinox longitude
   real(r8) :: obliqr1                   ! obliquity (radians)
   real(r8) :: lambm01                   ! mean lon of perihelion at vernal equinox (radians)
   real(r8) :: mvelpp1                   ! moving vernal equinox lon + pi (rads)
   real(r8) :: delta1                    ! declination angle (radians)
   real(r8) :: eccf1                     ! Earth-sun distance factor

   ! Assign local pointers to derived type arrays
   ncol       => clm3%g%c%p%column
   ivt        => clm3%g%c%p%itype
   z          => clm3%g%c%cps%z
   dz         => clm3%g%c%cps%dz
   frac_veg_nosno => clm3%g%c%p%pps%frac_veg_nosno
   qflx_drain => clm3%g%c%cwf%qflx_drain
   h2osoi_liq => clm3%g%c%cws%h2osoi_liq
   h2osoi_ice => clm3%g%c%cws%h2osoi_ice
   t_soisno   => clm3%g%c%ces%t_soisno
   t_ref2m    => clm3%g%c%p%pes%t_ref2m
   t_veg      => clm3%g%c%p%pes%t_veg
   fwet       => clm3%g%c%p%pps%fwet
   snowdp     => clm3%g%c%cps%snowdp
   rootfr     => clm3%g%c%p%pps%rootfr
   psisat     => clm3%g%c%cps%psisat
   soilpsi    => clm3%g%c%cps%soilpsi
   bsw2        => clm3%g%c%cps%bsw2
   sucsat     => clm3%g%c%cps%sucsat
   bsw        => clm3%g%c%cps%bsw
   watfc      => clm3%g%c%cps%watfc
   forc_t     => clm_a2l%forc_t
   forc_q     => clm_a2l%forc_q
   forc_pbot  => clm_a2l%forc_pbot
   forc_vp    => clm_a2l%forc_vp
   qg         => clm3%g%c%cws%qg
   thm        => clm3%g%c%p%pes%thm
   forc_hgt_u => clm_a2l%forc_hgt_u
   forc_hgt_t => clm_a2l%forc_hgt_t
   forc_hgt_q => clm_a2l%forc_hgt_q
   forc_hgt_u_pft => clm3%g%c%p%pps%forc_hgt_u_pft
   forc_hgt_t_pft => clm3%g%c%p%pps%forc_hgt_t_pft
   forc_hgt_q_pft => clm3%g%c%p%pps%forc_hgt_q_pft
   rb         => clm3%g%c%p%pps%rb1
   tlai       => clm3%g%c%p%pps%tlai
   elai       => clm3%g%c%p%pps%elai
   displa     => clm3%g%c%p%pps%displa
   z0mv       => clm3%g%c%p%pps%z0mv
   z0hv       => clm3%g%c%p%pps%z0hv
   z0qv       => clm3%g%c%p%pps%z0qv
   z0mg       => clm3%g%c%cps%z0mg
   z0hg       => clm3%g%c%cps%z0hg
   z0qg       => clm3%g%c%cps%z0qg
   forc_pco2  => clm_a2l%forc_pco2
   forc_po2   => clm_a2l%forc_po2
   forc_lwrad => clm_a2l%forc_lwrad
   forc_solar => clm_a2l%forc_solar
   forc_solad => clm_a2l%forc_solad
   forc_solai => clm_a2l%forc_solai
   latrad     => clm3%g%lat
   lonrad     => clm3%g%lon
   latdeg     => clm3%g%latdeg
   londeg     => clm3%g%londeg
   lat_a      => clm3%g%lat_a
   lon_a      => clm3%g%lon_a
   latdeg_a   => clm3%g%latdeg_a
   londeg_a   => clm3%g%londeg_a
   coszen     => clm3%g%c%cps%coszen
   albgrd     => clm3%g%c%cps%albgrd
   albgri     => clm3%g%c%cps%albgri
   decl       => clm3%g%c%cps%decl

   ! Assign VIC data to local pointers

   do g = begg, endg

   ! Added by MAB, 8/27/13

     latdeg(g) = lat
     latdeg_a(g) = lat
     if(lon .lt. 0.) then
       londeg(g) = 360._r8 + lon
     else
       londeg(g) = lon
     endif
     londeg_a(g) = londeg(g)

     latrad(g) = latdeg(g) * SHR_CONST_PI / 180.
     lonrad(g) = londeg(g) * SHR_CONST_PI / 180.
     lat_a(g) = latrad(g)
     lon_a(g) = lonrad(g)

     call shr_orb_params(yr, eccen, obliq, mvelp, obliqr, lambm0, mvelpp)
     call shr_orb_decl(jday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
     cosz = shr_orb_cosz(jday, latrad(g), lonrad(g), delta)
     do c = begc, endc
       coszen((g - 1) * endc + c) = cosz
       decl((g - 1) * endc + c) = delta
     end do

   ! Added by MAB, 10/8/13

     if(yr1 .gt. 0 .and. jday1 .gt. 0.) then
       call shr_orb_params(yr1, eccen1, obliq1, mvelp1, obliqr1, lambm01, &
	mvelpp1)
       call shr_orb_decl(jday1, eccen1, mvelpp1, lambm01, obliqr1, delta1, &
	eccf1)
     end if

   ! Added by MAB, 8/29/13

     forc_t(g) = tair + 273.16
     forc_q(g) = 0.622_r8 * vp / psfc
     forc_vp(g) = vp
     vg = vp + vpd
     do c = begc, endc
       qg((g - 1) * endc + c) = 0.622_r8 * vg / psfc
     end do

   ! Added by MAB, 10/11/13

     forc_hgt_u(g) = 10._r8
     forc_hgt_t(g) = 2._r8
     forc_hgt_q(g) = 2._r8

   ! Added by MAB, 8/14/13

     forc_pbot(g) = psfc
     forc_pco2(g) = psfc * 3.8e-4
     forc_po2(g) = psfc * 0.2095
     forc_lwrad(g) = lwrad
     forc_solar(g) = swrad
     do r = 1, numrad
       forc_solad(g,r) = swrd(r)
       forc_solai(g,r) = swri(r)
     end do

   do c = begc, endc

     qflx_drain((g - 1) * endc + c) = baseflow(c)
     snowdp((g - 1) * endc + c) = snodep(c)

     ! Added by MAB, 10/11/13
     if(snowdp((g - 1) * endc + c) .gt. 0.) then
       z0mg((g - 1) * endc + c) = zos
     else
       z0mg((g - 1) * endc + c) = zo
     end if
     z0hg((g - 1) * endc + c) = z0mv(c)
     z0qg((g - 1) * endc + c) = z0mv(c)

     do r = 1, numrad
       albgrd((g - 1) * endc + c,r) = alb(c)
       albgri((g - 1) * endc + c,r) = alb(c)
     end do

     do k = -nlevsno + 1, nlevgrnd  ! -nlevsno+1:nlevgrnd
         z((g - 1) * endc + c,k) = dep(k + nlevsno, c)
         dz((g - 1) * endc + c,k) = thick(k + nlevsno, c)
         h2osoi_liq((g - 1) * endc + c,k) = moist(k + nlevsno,c)
         h2osoi_ice((g - 1) * endc + c,k) = ice(k + nlevsno, c)
         t_soisno((g - 1) * endc + c,k) = tsoisno(k + nlevsno,c) + &
		SHR_CONST_TKFRZ
     end do

     do k = 1, nlevgrnd
         bsw((g - 1) * endc + c,k) = coeff(k,c,1)
         sucsat((g - 1) * endc + c,k) = satpsi(k,c,1)
         bsw2((g - 1) * endc + c,k) = coeff(k,c,2)
         psisat((g - 1) * endc + c,k) = satpsi(k,c,2)
         soilpsi((g - 1) * endc + c,k) = soipsi(k,c,2)
         watfc((g - 1) * endc + c,k) = fldcapwat(k,c)
     end do

     ! How to find the right p's for g != 1?  (g - 1) * endc * endp?

     do p = begp, endp
         t_ref2m((g - 1) * endc * endp + (c - 1) * endp + p) = t2m(p) + 273.16
         do k = 1, nlevgrnd
           rootfr((g - 1) * endc * endp + (c - 1) * endp + p,k) = rootf(k,p,c)
         end do
         t_veg((g - 1) * endc * endp + (c - 1) * endp + p) = tveg(ivt(p),c) + &
		273.16
         fwet((g - 1) * endc * endp + (c - 1) * endp + p) = fvegwet(ivt(p),c)
         rb((g - 1) * endc * endp + (c - 1) * endp + p) = rveg(ivt(p),c)

         ! Added by MAB, 10/11/13
!         if(rec == 0 .and. nspinup == 0) then
         if(rec == 0) then
           tlai((g - 1) * endc * endp + (c - 1) * endp + p) = lai(ivt(p),c)
         end if
         displa((g - 1) * endc * endp + (c - 1) * endp + p) = displ(ivt(p),c)
         z0mv((g - 1) * endc * endp + (c - 1) * endp + p) = zov(ivt(p),c)
         z0hv((g - 1) * endc * endp + (c - 1) * endp + p) = zov(ivt(p),c)
         z0qv((g - 1) * endc * endp + (c - 1) * endp + p) = zov(ivt(p),c)
         if(frac_veg_nosno((g - 1) * endc * endp + (c - 1) * endp + p) == 0) &
		then
           forc_hgt_u_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_u(1) + z0mg(c) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
           forc_hgt_t_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_t(1) + z0mg(c) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
           forc_hgt_q_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_q(1) + z0mg(c) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
         else
           forc_hgt_u_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_u(1) + &
		z0mv((g - 1) * endc * endp + (c - 1) * endp + p) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
           forc_hgt_t_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_t(1) + &
		z0mv((g - 1) * endc * endp + (c - 1) * endp + p) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
           forc_hgt_q_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_q(1) + &
		z0mv((g - 1) * endc * endp + (c - 1) * endp + p) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
         end if
         thm((g - 1) * endc * endp + (c - 1) * endp + p) = forc_t(1) + &
		0.0098_r8 * &
		forc_hgt_t_pft((g - 1) * endc * endp + (c - 1) * endp + p)
       end do
     end do
   end do

   begc = 1
   endc = num_soilc
   begp = 1
   endp = num_soilc * num_soilp

   soil1c => clm3%g%c%ccs%soil1c
   soil2c => clm3%g%c%ccs%soil2c
   soil3c => clm3%g%c%ccs%soil3c
   soil4c => clm3%g%c%ccs%soil4c
   litr1c => clm3%g%c%ccs%litr1c
   litr2c => clm3%g%c%ccs%litr2c
   litr3c => clm3%g%c%ccs%litr3c
   cwdc => clm3%g%c%ccs%cwdc
   leafc => clm3%g%c%p%pcs%leafc
   frootc => clm3%g%c%p%pcs%frootc
   livestemc => clm3%g%c%p%pcs%livestemc
   deadstemc => clm3%g%c%p%pcs%deadstemc
   livecrootc => clm3%g%c%p%pcs%livecrootc
   deadcrootc => clm3%g%c%p%pcs%deadcrootc
   woodc => clm3%g%c%p%pcs%woodc
   totvegc => clm3%g%c%p%pcs%totvegc
   totlitc => clm3%g%c%ccs%totlitc
   totsomc => clm3%g%c%ccs%totsomc
   soil1n => clm3%g%c%cns%soil1n
   soil2n => clm3%g%c%cns%soil2n
   soil3n => clm3%g%c%cns%soil3n
   soil4n => clm3%g%c%cns%soil4n
   sminn => clm3%g%c%cns%sminn
   litr1n => clm3%g%c%cns%litr1n
   litr2n => clm3%g%c%cns%litr2n
   litr3n => clm3%g%c%cns%litr3n
   cwdn => clm3%g%c%cns%cwdn
   leafn => clm3%g%c%p%pns%leafn
   frootn => clm3%g%c%p%pns%frootn
   livestemn => clm3%g%c%p%pns%livestemn
   deadstemn => clm3%g%c%p%pns%deadstemn
   livecrootn => clm3%g%c%p%pns%livecrootn
   deadcrootn => clm3%g%c%p%pns%deadcrootn
   cpool_to_leafc => clm3%g%c%p%pcf%cpool_to_leafc
   gpp2 => clm3%g%c%p%pepv%gpp
   gpp => clm3%g%c%p%pcf%gpp
   npp => clm3%g%c%p%pcf%npp
   ar => clm3%g%c%p%pcf%ar
   hr => clm3%g%c%ccf%hr
   nee => clm3%g%c%ccf%nee
   nep => clm3%g%c%ccf%nep
   tlai       => clm3%g%c%p%pps%tlai
   elai       => clm3%g%c%p%pps%elai
   dormant_flag => clm3%g%c%p%pepv%dormant_flag
   days_active => clm3%g%c%p%pepv%days_active
   onset_flag => clm3%g%c%p%pepv%onset_flag
   onset_counter => clm3%g%c%p%pepv%onset_counter
   onset_gddflag => clm3%g%c%p%pepv%onset_gddflag
   onset_fdd => clm3%g%c%p%pepv%onset_fdd
   onset_gdd => clm3%g%c%p%pepv%onset_gdd
   onset_swi => clm3%g%c%p%pepv%onset_swi
   offset_flag => clm3%g%c%p%pepv%offset_flag
   offset_counter => clm3%g%c%p%pepv%offset_counter
   offset_fdd => clm3%g%c%p%pepv%offset_fdd
   offset_swi => clm3%g%c%p%pepv%offset_swi
   lgsf => clm3%g%c%p%pepv%lgsf
   bglfr => clm3%g%c%p%pepv%bglfr
   bgtr => clm3%g%c%p%pepv%bgtr
   dayl => clm3%g%c%p%pepv%dayl
   prev_dayl => clm3%g%c%p%pepv%prev_dayl
   annavg_t2m => clm3%g%c%p%pepv%annavg_t2m
   tempavg_t2m => clm3%g%c%p%pepv%tempavg_t2m
   availc => clm3%g%c%p%pepv%availc
   xsmrpool_recover => clm3%g%c%p%pepv%xsmrpool_recover
   alloc_pnow => clm3%g%c%p%pepv%alloc_pnow
   c_allometry => clm3%g%c%p%pepv%c_allometry
   n_allometry => clm3%g%c%p%pepv%n_allometry
   plant_ndemand => clm3%g%c%p%pepv%plant_ndemand
   tempsum_potential_gpp => clm3%g%c%p%pepv%tempsum_potential_gpp
   annsum_potential_gpp => clm3%g%c%p%pepv%annsum_potential_gpp
   tempmax_retransn => clm3%g%c%p%pepv%tempmax_retransn
   annmax_retransn => clm3%g%c%p%pepv%annmax_retransn
   avail_retransn => clm3%g%c%p%pepv%avail_retransn
   plant_nalloc => clm3%g%c%p%pepv%plant_nalloc
   plant_calloc => clm3%g%c%p%pepv%plant_calloc
   excess_cflux => clm3%g%c%p%pepv%excess_cflux
   downreg => clm3%g%c%p%pepv%downreg
   prev_leafc_to_litter => clm3%g%c%p%pepv%prev_leafc_to_litter
   prev_frootc_to_litter => clm3%g%c%p%pepv%prev_frootc_to_litter
   tempsum_npp => clm3%g%c%p%pepv%tempsum_npp
   annsum_npp => clm3%g%c%p%pepv%annsum_npp
   leafc_storage => clm3%g%c%p%pcs%leafc_storage
   leafc_xfer => clm3%g%c%p%pcs%leafc_xfer
   frootc_storage => clm3%g%c%p%pcs%frootc_storage
   frootc_xfer => clm3%g%c%p%pcs%frootc_xfer
   livestemc_storage => clm3%g%c%p%pcs%livestemc_storage
   livestemc_xfer => clm3%g%c%p%pcs%livestemc_xfer
   deadstemc_storage => clm3%g%c%p%pcs%deadstemc_storage
   deadstemc_xfer => clm3%g%c%p%pcs%deadstemc_xfer
   livecrootc_storage => clm3%g%c%p%pcs%livecrootc_storage
   livecrootc_xfer => clm3%g%c%p%pcs%livecrootc_xfer
   deadcrootc_storage => clm3%g%c%p%pcs%deadcrootc_storage
   deadcrootc_xfer => clm3%g%c%p%pcs%deadcrootc_xfer
   gresp_storage => clm3%g%c%p%pcs%gresp_storage
   gresp_xfer => clm3%g%c%p%pcs%gresp_xfer
   cpool => clm3%g%c%p%pcs%cpool
   xsmrpool => clm3%g%c%p%pcs%xsmrpool
   pft_ctrunc => clm3%g%c%p%pcs%pft_ctrunc
   leafn_storage => clm3%g%c%p%pns%leafn_storage
   leafn_xfer => clm3%g%c%p%pns%leafn_xfer
   frootn_storage => clm3%g%c%p%pns%frootn_storage
   frootn_xfer => clm3%g%c%p%pns%frootn_xfer
   livestemn_storage => clm3%g%c%p%pns%livestemn_storage
   livestemn_xfer => clm3%g%c%p%pns%livestemn_xfer
   deadstemn_storage => clm3%g%c%p%pns%deadstemn_storage
   deadstemn_xfer => clm3%g%c%p%pns%deadstemn_xfer
   livecrootn_storage => clm3%g%c%p%pns%livecrootn_storage
   livecrootn_xfer => clm3%g%c%p%pns%livecrootn_xfer
   deadcrootn_storage => clm3%g%c%p%pns%deadcrootn_storage
   deadcrootn_xfer => clm3%g%c%p%pns%deadcrootn_xfer
   retransn => clm3%g%c%p%pns%retransn
   npool => clm3%g%c%p%pns%npool
   pft_ntrunc => clm3%g%c%p%pns%pft_ntrunc
   decl => clm3%g%c%cps%decl
   fpi => clm3%g%c%cps%fpi
   fpg => clm3%g%c%cps%fpg
   annsum_counter => clm3%g%c%cps%annsum_counter
   cannsum_npp => clm3%g%c%cps%cannsum_npp
   cannsum_t2m => clm3%g%c%cps%cannavg_t2m
   me => clm3%g%c%cps%me
   fire_prob => clm3%g%c%cps%fire_prob
   mean_fire_prob => clm3%g%c%cps%mean_fire_prob
   fireseasonl => clm3%g%c%cps%fireseasonl
   farea_burned => clm3%g%c%cps%farea_burned
   ann_farea_burned => clm3%g%c%cps%ann_farea_burned
   seedc => clm3%g%c%ccs%seedc
   col_ctrunc => clm3%g%c%ccs%col_ctrunc
   totcolc => clm3%g%c%ccs%totcolc
   prod10c => clm3%g%c%ccs%prod10c
   prod100c => clm3%g%c%ccs%prod100c
   seedn => clm3%g%c%cns%seedn
   col_ntrunc => clm3%g%c%cns%col_ntrunc
   totcoln => clm3%g%c%cns%totcoln
   prod10n => clm3%g%c%cns%prod10n
   prod100n => clm3%g%c%cns%prod100n

   if(rec .eq. 0 .and. init_state .eq. 1) then

    do c = begc, endc
     soil1c(c) = soilcfast(c)
     soil2c(c) = soilcmid(c)
     soil3c(c) = soilcslo1(c)
     soil4c(c) = soilcslo2(c)
     litr1c(c) = litrlabc(c)
     litr2c(c) = litrcellc(c)
     litr3c(c) = litrligc(c)
     cwdc(c) = cwoodc(c)
     totlitc(c) = litctot(c)
     totsomc(c) = somctot(c)
     soil1n(c) = soilnfast(c)
     soil2n(c) = soilnmid(c)
     soil3n(c) = soilnslo1(c)
     soil4n(c) = soilnslo2(c)
     sminn(c) = soilminn(c)
     litr1n(c) = litrlabn(c)
     litr2n(c) = litrcelln(c)
     litr3n(c) = litrlign(c)
     cwdn(c) = cwoodn(c)
     hr(c) = heteroresp(c)
     nee(c) = ecoexchn(c)
     nep(c) = prodecon(c)
     decl(c) = declin(c)
     fpi(c) = potimmobfract(c)
     fpg(c) = potgppfract(c)
     annsum_counter(c) = annsumcount(c)
     cannsum_npp(c) = colannsumnpp(c)
     cannsum_t2m(c) = colannsumt2m(c)
     do k = 1, nlevgrnd
       watfc(c,k) = fldcapwat(k,c)
     end do
     me(c) = extinctmoist(c)
     fire_prob(c) = fireprob(c)
     mean_fire_prob(c) = meanfireprob(c)
     fireseasonl(c) = fireseasnlen(c)
     farea_burned(c) = areafractburn(c)
     ann_farea_burned(c) = annareafractburn(c)
     seedc(c) = seedcc(c)
     col_ctrunc(c) = colctrunc(c)
     totcolc(c) = colctot(c)
     prod10c(c) = woodprodc10(c)
     prod100c(c) = woodprodc100(c)
     seedn(c) = seednn(c)
     col_ntrunc(c) = colntrunc(c)
     totcoln(c) = colntot(c)
     prod10n(c) = woodprodn10(c)
     prod100n(c) = woodprodn100(c)
    end do

    do p = begp, endp
     c = ncol(p)
       elai(p) = lai(ivt(p),c)
       leafc(p) = leafcc(ivt(p),c)
       frootc(p) = finrtc(ivt(p),c)
       livestemc(p) = livstemc(ivt(p),c)
       deadstemc(p) = deadstemcc(ivt(p),c)
       livecrootc(p) = livcorsrtc(ivt(p),c)
       deadcrootc(p) = deadcorsrtc(ivt(p),c)
       woodc(p) = woodcc(ivt(p),c)
       totvegc(p) = vegctot(ivt(p),c)
       leafn(p) = leafnn(ivt(p),c)
       frootn(p) = finrtn(ivt(p),c)
       livestemn(p) = livstemn(ivt(p),c)
       deadstemn(p) = deadstemnn(ivt(p),c)
       livecrootn(p) = livcorsrtn(ivt(p),c)
       deadcrootn(p) = deadcorsrtn(ivt(p),c)
       gpp2(p) = prod1gb4(ivt(p),c)
       gpp(p) = prod1g(ivt(p),c)
       npp(p) = prod1n(ivt(p),c)
       ar(p) = autoresp(ivt(p),c)
       dormant_flag(p) = dormancy(ivt(p),c)
       days_active(p) = ndaysact(ivt(p),c)
       onset_flag(p) = onsetflg(ivt(p),c)
       onset_counter(p) = onsetcount(ivt(p),c)
       onset_gddflag(p) = onsetgddflg(ivt(p),c)
       onset_fdd(p) = onsetfdd(ivt(p),c)
       onset_gdd(p) = onsetgdd(ivt(p),c)
       onset_swi(p) = onsetswi(ivt(p),c)
       offset_flag(p) = offsetflg(ivt(p),c)
       offset_counter(p) = offsetcount(ivt(p),c)
       offset_fdd(p) = offsetfdd(ivt(p),c)
       offset_swi(p) = offsetswi(ivt(p),c)
       lgsf(p) = lgsfact(ivt(p),c)
       bglfr(p) = backlfr(ivt(p),c)
       bgtr(p) = backtgr(ivt(p),c)
       dayl(p) = daylen(ivt(p),c)
       prev_dayl(p) = prevdaylen(ivt(p),c)
       annavg_t2m(p) = annavgt2m(ivt(p),c)
       tempavg_t2m(p) = tempavgt2m(ivt(p),c)
       availc(p) = cavail(ivt(p),c)
       xsmrpool_recover(p) = cflxrecov(ivt(p),c)
       alloc_pnow(p) = allocpnow(ivt(p),c)
       c_allometry(p) = callom(ivt(p),c)
       n_allometry(p) = nallom(ivt(p),c)
       plant_ndemand(p) = plantndem(ivt(p),c)
       tempsum_potential_gpp(p) = tempsumpotgpp(ivt(p),c)
       annsum_potential_gpp(p) = annsumpotgpp(ivt(p),c)
       tempmax_retransn(p) = tempmxretransn(ivt(p),c)
       annmax_retransn(p) = annmxretransn(ivt(p),c)
       avail_retransn(p) = availretransn(ivt(p),c)
       plant_nalloc(p) = plantnalloc(ivt(p),c)
       plant_calloc(p) = plantcalloc(ivt(p),c)
       excess_cflux(p) = cflxex(ivt(p),c)
       downreg(p) = dwnreg(ivt(p),c)
       prev_leafc_to_litter(p) = prevleafc2litr(ivt(p),c)
       prev_frootc_to_litter(p) = prevfinrtc2litr(ivt(p),c)
       tempsum_npp(p) = tempsumnpp(ivt(p),c)
       annsum_npp(p) = annsumnpp(ivt(p),c)
       leafc_storage(p) = leafcstor(ivt(p),c)
       leafc_xfer(p) = leafctrans(ivt(p),c)
       frootc_storage(p) = finrtcstor(ivt(p),c)
       frootc_xfer(p) = finrtctrans(ivt(p),c)
       livestemc_storage(p) = livstemcstor(ivt(p),c)
       livestemc_xfer(p) = livstemctrans(ivt(p),c)
       deadstemc_storage(p) = deadstemcstor(ivt(p),c)
       deadstemc_xfer(p) = deadstemctrans(ivt(p),c)
       livecrootc_storage(p) = livcorsrtcstor(ivt(p),c)
       livecrootc_xfer(p) = livcorsrtctrans(ivt(p),c)
       deadcrootc_storage(p) = deadcorsrtcstor(ivt(p),c)
       deadcrootc_xfer(p) = deadcorsrtctrans(ivt(p),c)
       gresp_storage(p) = grorespstor(ivt(p),c)
       gresp_xfer(p) = groresptrans(ivt(p),c)
       cpool(p) = photocpool(ivt(p),c)
       xsmrpool(p) = mrcpool(ivt(p),c)
       pft_ctrunc(p) = ctruncpft(ivt(p),c)
       leafn_storage(p) = leafnstor(ivt(p),c)
       leafn_xfer(p) = leafntrans(ivt(p),c)
       frootn_storage(p) = finrtnstor(ivt(p),c)
       frootn_xfer(p) = finrtntrans(ivt(p),c)
       livestemn_storage(p) = livstemnstor(ivt(p),c)
       livestemn_xfer(p) = livstemntrans(ivt(p),c)
       deadstemn_storage(p) = deadstemnstor(ivt(p),c)
       deadstemn_xfer(p) = deadstemntrans(ivt(p),c)
       livecrootn_storage(p) = livcorsrtnstor(ivt(p),c)
       livecrootn_xfer(p) = livcorsrtntrans(ivt(p),c)
       deadcrootn_storage(p) = deadcorsrtnstor(ivt(p),c)
       deadcrootn_xfer(p) = deadcorsrtntrans(ivt(p),c)
       retransn(p) = nretrans(ivt(p),c)
       npool(p) = photonpool(ivt(p),c)
       pft_ntrunc(p) = ntruncpft(ivt(p),c)
    end do
   end if

   if(jday1 .gt. 0.) then
     call SurfaceAlbedo(begg, endg, begc, endc, begp, endp, jday1, delta1)
   end if
   call SurfaceRadiation(yr, mo, day, secs, dt, begp, endp)
   call CanopyFluxes(dt, nlevgrnd, begg, endg, begc, endc, begp, endp)

   doy = int(jday)
   call CNEcosystemDyn(adspinup, nlevgrnd, rec, yr, doy, dt, begg, endg, &
	begc, endc, begp, endp, num_soilc, num_soilp)
   call CNAnnualUpdate(dt, yr, begc, endc, begp, endp, num_soilc, num_soilp)

   ! Check the carbon and nitrogen balance
   call CBalanceCheck(dt, begc, endc, num_soilc)
   call NBalanceCheck(dt, begc, endc, num_soilc)

   soil1c => clm3%g%c%ccs%soil1c
   soil2c => clm3%g%c%ccs%soil2c
   soil3c => clm3%g%c%ccs%soil3c
   soil4c => clm3%g%c%ccs%soil4c
   litr1c => clm3%g%c%ccs%litr1c
   litr2c => clm3%g%c%ccs%litr2c
   litr3c => clm3%g%c%ccs%litr3c
   cwdc => clm3%g%c%ccs%cwdc
   leafc => clm3%g%c%p%pcs%leafc
   frootc => clm3%g%c%p%pcs%frootc
   livestemc => clm3%g%c%p%pcs%livestemc
   deadstemc => clm3%g%c%p%pcs%deadstemc
   livecrootc => clm3%g%c%p%pcs%livecrootc
   deadcrootc => clm3%g%c%p%pcs%deadcrootc
   woodc => clm3%g%c%p%pcs%woodc
   totvegc => clm3%g%c%p%pcs%totvegc
   totlitc => clm3%g%c%ccs%totlitc
   totsomc => clm3%g%c%ccs%totsomc
   soil1n => clm3%g%c%cns%soil1n
   soil2n => clm3%g%c%cns%soil2n
   soil3n => clm3%g%c%cns%soil3n
   soil4n => clm3%g%c%cns%soil4n
   sminn => clm3%g%c%cns%sminn
   litr1n => clm3%g%c%cns%litr1n
   litr2n => clm3%g%c%cns%litr2n
   litr3n => clm3%g%c%cns%litr3n
   cwdn => clm3%g%c%cns%cwdn
   leafn => clm3%g%c%p%pns%leafn
   frootn => clm3%g%c%p%pns%frootn
   livestemn => clm3%g%c%p%pns%livestemn
   deadstemn => clm3%g%c%p%pns%deadstemn
   livecrootn => clm3%g%c%p%pns%livecrootn
   deadcrootn => clm3%g%c%p%pns%deadcrootn
   cpool_to_leafc => clm3%g%c%p%pcf%cpool_to_leafc
   gpp2 => clm3%g%c%p%pepv%gpp
   gpp => clm3%g%c%p%pcf%gpp
   npp => clm3%g%c%p%pcf%npp
   mr => clm3%g%c%p%pcf%mr
   leaf_mr => clm3%g%c%p%pcf%leaf_mr
   gr => clm3%g%c%p%pcf%gr
   ar => clm3%g%c%p%pcf%ar
   hr => clm3%g%c%ccf%hr
   lithr => clm3%g%c%ccf%lithr
   nee => clm3%g%c%ccf%nee
   nep => clm3%g%c%ccf%nep
   tlai       => clm3%g%c%p%pps%tlai
   elai       => clm3%g%c%p%pps%elai
   dormant_flag => clm3%g%c%p%pepv%dormant_flag
   days_active => clm3%g%c%p%pepv%days_active
   onset_flag => clm3%g%c%p%pepv%onset_flag
   onset_counter => clm3%g%c%p%pepv%onset_counter
   onset_gddflag => clm3%g%c%p%pepv%onset_gddflag
   onset_fdd => clm3%g%c%p%pepv%onset_fdd
   onset_gdd => clm3%g%c%p%pepv%onset_gdd
   onset_swi => clm3%g%c%p%pepv%onset_swi
   offset_flag => clm3%g%c%p%pepv%offset_flag
   offset_counter => clm3%g%c%p%pepv%offset_counter
   offset_fdd => clm3%g%c%p%pepv%offset_fdd
   offset_swi => clm3%g%c%p%pepv%offset_swi
   lgsf => clm3%g%c%p%pepv%lgsf
   bglfr => clm3%g%c%p%pepv%bglfr
   bgtr => clm3%g%c%p%pepv%bgtr
   dayl => clm3%g%c%p%pepv%dayl
   prev_dayl => clm3%g%c%p%pepv%prev_dayl
   annavg_t2m => clm3%g%c%p%pepv%annavg_t2m
   tempavg_t2m => clm3%g%c%p%pepv%tempavg_t2m
   availc => clm3%g%c%p%pepv%availc
   xsmrpool_recover => clm3%g%c%p%pepv%xsmrpool_recover
   alloc_pnow => clm3%g%c%p%pepv%alloc_pnow
   c_allometry => clm3%g%c%p%pepv%c_allometry
   n_allometry => clm3%g%c%p%pepv%n_allometry
   plant_ndemand => clm3%g%c%p%pepv%plant_ndemand
   tempsum_potential_gpp => clm3%g%c%p%pepv%tempsum_potential_gpp
   annsum_potential_gpp => clm3%g%c%p%pepv%annsum_potential_gpp
   tempmax_retransn => clm3%g%c%p%pepv%tempmax_retransn
   annmax_retransn => clm3%g%c%p%pepv%annmax_retransn
   avail_retransn => clm3%g%c%p%pepv%avail_retransn
   plant_nalloc => clm3%g%c%p%pepv%plant_nalloc
   plant_calloc => clm3%g%c%p%pepv%plant_calloc
   excess_cflux => clm3%g%c%p%pepv%excess_cflux
   downreg => clm3%g%c%p%pepv%downreg
   prev_leafc_to_litter => clm3%g%c%p%pepv%prev_leafc_to_litter
   prev_frootc_to_litter => clm3%g%c%p%pepv%prev_frootc_to_litter
   tempsum_npp => clm3%g%c%p%pepv%tempsum_npp
   annsum_npp => clm3%g%c%p%pepv%annsum_npp
   leafc_storage => clm3%g%c%p%pcs%leafc_storage
   leafc_xfer => clm3%g%c%p%pcs%leafc_xfer
   frootc_storage => clm3%g%c%p%pcs%frootc_storage
   frootc_xfer => clm3%g%c%p%pcs%frootc_xfer
   livestemc_storage => clm3%g%c%p%pcs%livestemc_storage
   livestemc_xfer => clm3%g%c%p%pcs%livestemc_xfer
   deadstemc_storage => clm3%g%c%p%pcs%deadstemc_storage
   deadstemc_xfer => clm3%g%c%p%pcs%deadstemc_xfer
   livecrootc_storage => clm3%g%c%p%pcs%livecrootc_storage
   livecrootc_xfer => clm3%g%c%p%pcs%livecrootc_xfer
   deadcrootc_storage => clm3%g%c%p%pcs%deadcrootc_storage
   deadcrootc_xfer => clm3%g%c%p%pcs%deadcrootc_xfer
   gresp_storage => clm3%g%c%p%pcs%gresp_storage
   gresp_xfer => clm3%g%c%p%pcs%gresp_xfer
   cpool => clm3%g%c%p%pcs%cpool
   xsmrpool => clm3%g%c%p%pcs%xsmrpool
   pft_ctrunc => clm3%g%c%p%pcs%pft_ctrunc
   leafn_storage => clm3%g%c%p%pns%leafn_storage
   leafn_xfer => clm3%g%c%p%pns%leafn_xfer
   frootn_storage => clm3%g%c%p%pns%frootn_storage
   frootn_xfer => clm3%g%c%p%pns%frootn_xfer
   livestemn_storage => clm3%g%c%p%pns%livestemn_storage
   livestemn_xfer => clm3%g%c%p%pns%livestemn_xfer
   deadstemn_storage => clm3%g%c%p%pns%deadstemn_storage
   deadstemn_xfer => clm3%g%c%p%pns%deadstemn_xfer
   livecrootn_storage => clm3%g%c%p%pns%livecrootn_storage
   livecrootn_xfer => clm3%g%c%p%pns%livecrootn_xfer
   deadcrootn_storage => clm3%g%c%p%pns%deadcrootn_storage
   deadcrootn_xfer => clm3%g%c%p%pns%deadcrootn_xfer
   retransn => clm3%g%c%p%pns%retransn
   npool => clm3%g%c%p%pns%npool
   pft_ntrunc => clm3%g%c%p%pns%pft_ntrunc
   decl => clm3%g%c%cps%decl
   fpi => clm3%g%c%cps%fpi
   fpg => clm3%g%c%cps%fpg
   annsum_counter => clm3%g%c%cps%annsum_counter
   cannsum_npp => clm3%g%c%cps%cannsum_npp
   cannsum_t2m => clm3%g%c%cps%cannavg_t2m
   me => clm3%g%c%cps%me
   fire_prob => clm3%g%c%cps%fire_prob
   mean_fire_prob => clm3%g%c%cps%mean_fire_prob
   fireseasonl => clm3%g%c%cps%fireseasonl
   farea_burned => clm3%g%c%cps%farea_burned
   ann_farea_burned => clm3%g%c%cps%ann_farea_burned
   seedc => clm3%g%c%ccs%seedc
   col_ctrunc => clm3%g%c%ccs%col_ctrunc
   totcolc => clm3%g%c%ccs%totcolc
   prod10c => clm3%g%c%ccs%prod10c
   prod100c => clm3%g%c%ccs%prod100c
   seedn => clm3%g%c%cns%seedn
   col_ntrunc => clm3%g%c%cns%col_ntrunc
   totcoln => clm3%g%c%cns%totcoln
   prod10n => clm3%g%c%cns%prod10n
   prod100n => clm3%g%c%cns%prod100n
   litfall => clm3%g%c%p%pcf%litfall
   fpsn => clm3%g%c%p%pcf%fpsn
   cisun => clm3%g%c%p%pps%cisun
   cisha => clm3%g%c%p%pps%cisha
   rssun => clm3%g%c%p%pps%rssun
   rssha => clm3%g%c%p%pps%rssha
   parsun => clm3%g%c%p%pef%parsun
   parsha => clm3%g%c%p%pef%parsha
   laisun => clm3%g%c%p%pps%laisun
   laisha => clm3%g%c%p%pps%laisha

   do c = begc, endc
     soilcfast(c) = soil1c(c)
     soilcmid(c) = soil2c(c)
     soilcslo1(c) = soil3c(c)
     soilcslo2(c) = soil4c(c)
     litrlabc(c) = litr1c(c)
     litrcellc(c) = litr2c(c)
     litrligc(c) = litr3c(c)
     cwoodc(c) = cwdc(c)
     litctot(c) = totlitc(c)
     somctot(c) = totsomc(c)
     soilnfast(c) = soil1n(c)
     soilnmid(c) = soil2n(c)
     soilnslo1(c) = soil3n(c)
     soilnslo2(c) = soil4n(c)
     soilminn(c) = sminn(c)
     litrlabn(c) = litr1n(c)
     litrcelln(c) = litr2n(c)
     litrlign(c) = litr3n(c)
     cwoodn(c) = cwdn(c)
     heteroresp(c) = hr(c)
     litresp(c) = lithr(c)
     ecoexchn(c) = nee(c)
     prodecon(c) = nep(c)
     declin(c) = decl(c)
     potimmobfract(c) = fpi(c)
     potgppfract(c) = fpg(c)
     annsumcount(c) = annsum_counter(c)
     colannsumnpp(c) = cannsum_npp(c)
     colannsumt2m(c) = cannsum_t2m(c)
     do k = 1, nlevgrnd
       fldcapwat(k,c) = watfc(c,k)
     end do
     extinctmoist(c) = me(c)
     fireprob(c) = fire_prob(c)
     meanfireprob(c) = mean_fire_prob(c)
     fireseasnlen(c) = fireseasonl(c)
     areafractburn(c) = farea_burned(c)
     annareafractburn(c) = ann_farea_burned(c)
     seedcc(c) = seedc(c)
     colctrunc(c) = col_ctrunc(c)
     colctot(c) = totcolc(c)
     woodprodc10(c) = prod10c(c)
     woodprodc100(c) = prod100c(c)
     seednn(c) = seedn(c)
     colntrunc(c) = col_ntrunc(c)
     colntot(c) = totcoln(c)
     woodprodn10(c) = prod10n(c)
     woodprodn100(c) = prod100n(c)
   end do

   do p = begp, endp
     c = ncol(p)
       lai(ivt(p),c) = elai(p)
       leafcc(ivt(p),c) = leafc(p)
       finrtc(ivt(p),c) = frootc(p)
       livstemc(ivt(p),c) = livestemc(p)
       deadstemcc(ivt(p),c) = deadstemc(p)
       livcorsrtc(ivt(p),c) = livecrootc(p)
       deadcorsrtc(ivt(p),c) = deadcrootc(p)
       woodcc(ivt(p),c) = woodc(p)
       vegctot(ivt(p),c) = totvegc(p)
       leafnn(ivt(p),c) = leafn(p)
       finrtn(ivt(p),c) = frootn(p)
       livstemn(ivt(p),c) = livestemn(p)
       deadstemnn(ivt(p),c) = deadstemn(p)
       livcorsrtn(ivt(p),c) = livecrootn(p)
       deadcorsrtn(ivt(p),c) = deadcrootn(p)
       prod1gb4(ivt(p),c) = gpp2(p)
       prod1g(ivt(p),c) = gpp(p)
       prod1n(ivt(p),c) = npp(p)
       darkresp(ivt(p),c) = leaf_mr(p)
       maintresp(ivt(p),c) = mr(p)
       groresp(ivt(p),c) = gr(p)
       autoresp(ivt(p),c) = ar(p)
       dormancy(ivt(p),c) = dormant_flag(p)
       ndaysact(ivt(p),c) = days_active(p)
       onsetflg(ivt(p),c) = onset_flag(p)
       onsetcount(ivt(p),c) = onset_counter(p)
       onsetgddflg(ivt(p),c) = onset_gddflag(p)
       onsetfdd(ivt(p),c) = onset_fdd(p)
       onsetgdd(ivt(p),c) = onset_gdd(p)
       onsetswi(ivt(p),c) = onset_swi(p)
       offsetflg(ivt(p),c) = offset_flag(p)
       offsetcount(ivt(p),c) = offset_counter(p)
       offsetfdd(ivt(p),c) = offset_fdd(p)
       offsetswi(ivt(p),c) = offset_swi(p)
       lgsfact(ivt(p),c) = lgsf(p)
       backlfr(ivt(p),c) = bglfr(p)
       backtgr(ivt(p),c) = bgtr(p)
       daylen(ivt(p),c) = dayl(p)
       prevdaylen(ivt(p),c) = prev_dayl(p)
       annavgt2m(ivt(p),c) = annavg_t2m(p)
       tempavgt2m(ivt(p),c) = tempavg_t2m(p)
       cavail(ivt(p),c) = availc(p)
       cflxrecov(ivt(p),c) = xsmrpool_recover(p)
       allocpnow(ivt(p),c) = alloc_pnow(p)
       callom(ivt(p),c) = c_allometry(p)
       nallom(ivt(p),c) = n_allometry(p)
       plantndem(ivt(p),c) = plant_ndemand(p)
       tempsumpotgpp(ivt(p),c) = tempsum_potential_gpp(p)
       annsumpotgpp(ivt(p),c) = annsum_potential_gpp(p)
       tempmxretransn(ivt(p),c) = tempmax_retransn(p)
       annmxretransn(ivt(p),c) = annmax_retransn(p)
       availretransn(ivt(p),c) = avail_retransn(p)
       plantnalloc(ivt(p),c) = plant_nalloc(p)
       plantcalloc(ivt(p),c) = plant_calloc(p)
       cflxex(ivt(p),c) = excess_cflux(p)
       dwnreg(ivt(p),c) = downreg(p)
       prevleafc2litr(ivt(p),c) = prev_leafc_to_litter(p)
       prevfinrtc2litr(ivt(p),c) = prev_frootc_to_litter(p)
       tempsumnpp(ivt(p),c) = tempsum_npp(p)
       annsumnpp(ivt(p),c) = annsum_npp(p)
       leafcstor(ivt(p),c) = leafc_storage(p)
       leafctrans(ivt(p),c) = leafc_xfer(p)
       finrtcstor(ivt(p),c) = frootc_storage(p)
       finrtctrans(ivt(p),c) = frootc_xfer(p)
       livstemcstor(ivt(p),c) = livestemc_storage(p)
       livstemctrans(ivt(p),c) = livestemc_xfer(p)
       deadstemcstor(ivt(p),c) = deadstemc_storage(p)
       deadstemctrans(ivt(p),c) = deadstemc_xfer(p)
       livcorsrtcstor(ivt(p),c) = livecrootc_storage(p)
       livcorsrtctrans(ivt(p),c) = livecrootc_xfer(p)
       deadcorsrtcstor(ivt(p),c) = deadcrootc_storage(p)
       deadcorsrtctrans(ivt(p),c) = deadcrootc_xfer(p)
       grorespstor(ivt(p),c) = gresp_storage(p)
       groresptrans(ivt(p),c) = gresp_xfer(p)
       photocpool(ivt(p),c) = cpool(p)
       mrcpool(ivt(p),c) = xsmrpool(p)
       ctruncpft(ivt(p),c) = pft_ctrunc(p)
       leafnstor(ivt(p),c) = leafn_storage(p)
       leafntrans(ivt(p),c) = leafn_xfer(p)
       finrtnstor(ivt(p),c) = frootn_storage(p)
       finrtntrans(ivt(p),c) = frootn_xfer(p)
       livstemnstor(ivt(p),c) = livestemn_storage(p)
       livstemntrans(ivt(p),c) = livestemn_xfer(p)
       deadstemnstor(ivt(p),c) = deadstemn_storage(p)
       deadstemntrans(ivt(p),c) = deadstemn_xfer(p)
       livcorsrtnstor(ivt(p),c) = livecrootn_storage(p)
       livcorsrtntrans(ivt(p),c) = livecrootn_xfer(p)
       deadcorsrtnstor(ivt(p),c) = deadcrootn_storage(p)
       deadcorsrtntrans(ivt(p),c) = deadcrootn_xfer(p)
       nretrans(ivt(p),c) = retransn(p)
       photonpool(ivt(p),c) = npool(p)
       ntruncpft(ivt(p),c) = pft_ntrunc(p)
       litrfall(ivt(p),c) = litfall(p) * dt
       photosynth(ivt(p),c) = fpsn(p)
       intco2(ivt(p),c) = cisun(p) * laisun(p) + cisha(p) * laisha(p)
       stomresist(ivt(p),c) = rssun(p) * laisun(p) + rssha(p) * laisha(p)
       abspar(ivt(p),c) = parsun(p) * laisun(p) + parsha(p) * laisha(p)
   end do

end subroutine vic2clmtype
