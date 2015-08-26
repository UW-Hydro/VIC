module lnd_comp_mct

! !USES:

  use mct_mod
  use esmf
  use seq_cdata_mod
  use seq_infodata_mod
  use vic_cesm_interface
  use, intrinsic :: iso_c_binding

!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: lnd_init_mct
  public :: lnd_run_mct
  public :: lnd_final_mct
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  integer :: iulog   ! vic log file unit number
  integer :: shrlogunit   ! generic log unit

  ! lnd -> drv

  integer :: nflds_l2x = 0
  integer :: index_l2x_Sl_t           = 0 ! temperature
  integer :: index_l2x_Sl_tref        = 0 ! 2m reference temperature
  integer :: index_l2x_Sl_qref        = 0 ! 2m reference specific humidity
  integer :: index_l2x_Sl_avsdr       = 0 ! albedo: direct , visible
  integer :: index_l2x_Sl_anidr       = 0 ! albedo: direct , near-ir
  integer :: index_l2x_Sl_avsdf       = 0 ! albedo: diffuse, visible
  integer :: index_l2x_Sl_anidf       = 0 ! albedo: diffuse, near-ir
  integer :: index_l2x_Sl_snowh       = 0 ! snow height
  integer :: index_l2x_Sl_u10         = 0 ! 10m wind
  integer :: index_l2x_Sl_ddvel       = 0 ! dry deposition velocities (optional)
  integer :: index_l2x_Sl_fv          = 0 ! friction velocity  
  integer :: index_l2x_Sl_ram1        = 0 ! aerodynamical resistance
  integer :: index_l2x_Sl_soilw       = 0 ! volumetric soil water
  integer :: index_l2x_Sl_logz0       = 0 ! log z0 
  integer :: index_l2x_Fall_taux      = 0 ! wind stress, zonal
  integer :: index_l2x_Fall_tauy      = 0 ! wind stress, meridional
  integer :: index_l2x_Fall_lat       = 0 ! latent          heat flux
  integer :: index_l2x_Fall_sen       = 0 ! sensible        heat flux
  integer :: index_l2x_Fall_lwup      = 0 ! upward longwave heat flux
  integer :: index_l2x_Fall_evap      = 0 ! evaporation     water flux
  integer :: index_l2x_Fall_swnet     = 0 ! heat flux       shortwave net       
  integer :: index_l2x_Fall_fco2_lnd  = 0 ! co2 flux **For testing set to 0
  integer :: index_l2x_Fall_flxdst1   = 0 ! dust flux size bin 1    
  integer :: index_l2x_Fall_flxdst2   = 0 ! dust flux size bin 2    
  integer :: index_l2x_Fall_flxdst3   = 0 ! dust flux size bin 3    
  integer :: index_l2x_Fall_flxdst4   = 0 ! dust flux size bin 4
  integer :: index_l2x_Fall_flxvoc    = 0 ! MEGAN fluxes
  integer :: index_l2x_Flrl_rofliq    = 0 ! lnd->rtm input fluxes
  integer :: index_l2x_Flrl_rofice    = 0 ! lnd->rtm input fluxes


  ! drv -> lnd

  integer :: nflds_x2l = 0
  integer :: index_x2l_Sa_z           = 0 ! bottom atm level height
  integer :: index_x2l_Sa_u           = 0 ! bottom atm level zon wind
  integer :: index_x2l_Sa_v           = 0 ! bottom atm level mer wind
  integer :: index_x2l_Sa_ptem        = 0 ! bottom atm level pot temp
  integer :: index_x2l_Sa_shum        = 0 ! bottom atm level spec hum
  integer :: index_x2l_Sa_pbot        = 0 ! bottom atm level pressure
  integer :: index_x2l_Sa_tbot        = 0 ! bottom atm level temp
  integer :: index_x2l_Faxa_lwdn      = 0 ! downward lw heat flux
  integer :: index_x2l_Faxa_rainc     = 0 ! prec: liquid "convective"
  integer :: index_x2l_Faxa_rainl     = 0 ! prec: liquid "large scale"
  integer :: index_x2l_Faxa_snowc     = 0 ! prec: frozen "convective"
  integer :: index_x2l_Faxa_snowl     = 0 ! prec: frozen "large scale"
  integer :: index_x2l_Faxa_swndr     = 0 ! sw: nir direct  downward
  integer :: index_x2l_Faxa_swvdr     = 0 ! sw: vis direct  downward
  integer :: index_x2l_Faxa_swndf     = 0 ! sw: nir diffuse downward
  integer :: index_x2l_Faxa_swvdf     = 0 ! sw: vis diffuse downward
  integer :: index_x2l_Sa_co2prog     = 0 ! bottom atm level prognostic co2
  integer :: index_x2l_Sa_co2diag     = 0 ! bottom atm level diagnostic co2
  integer :: index_x2l_Faxa_bcphidry  = 0 ! flux: Black Carbon hydrophilic dry deposition
  integer :: index_x2l_Faxa_bcphodry  = 0 ! flux: Black Carbon hydrophobic dry deposition
  integer :: index_x2l_Faxa_bcphiwet  = 0 ! flux: Black Carbon hydrophilic wet deposition
  integer :: index_x2l_Faxa_ocphidry  = 0 ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2l_Faxa_ocphodry  = 0 ! flux: Organic Carbon hydrophobic dry deposition
  integer :: index_x2l_Faxa_ocphiwet  = 0 ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2l_Faxa_dstwet1   = 0 ! flux: Size 1 dust -- wet deposition
  integer :: index_x2l_Faxa_dstwet2   = 0 ! flux: Size 2 dust -- wet deposition
  integer :: index_x2l_Faxa_dstwet3   = 0 ! flux: Size 3 dust -- wet deposition
  integer :: index_x2l_Faxa_dstwet4   = 0 ! flux: Size 4 dust -- wet deposition
  integer :: index_x2l_Faxa_dstdry1   = 0 ! flux: Size 1 dust -- dry deposition
  integer :: index_x2l_Faxa_dstdry2   = 0 ! flux: Size 2 dust -- dry deposition
  integer :: index_x2l_Faxa_dstdry3   = 0 ! flux: Size 3 dust -- dry deposition
  integer :: index_x2l_Faxa_dstdry4   = 0 ! flux: Size 4 dust -- dry deposition 
  integer :: index_x2l_Flrr_flood     = 0 ! rtm->lnd rof (flood) flux

! tcraig - this can probably be removed
#if (1 == 0)
  ! sno -> drv (only if land-ice model is NOT a stub model)

  integer :: nflds_s2x = 0
  integer :: index_s2x_Ss_tsrf(glc_nec_max)   = 0 ! glc MEC temperature
  integer :: index_s2x_Ss_topo(glc_nec_max)   = 0 ! glc MEC topo height
  integer :: index_s2x_Fgss_qice(glc_nec_max) = 0 ! glc MEC ice flux

  ! drv -> sno (only if land-ice model is NOT a stub model)

  integer :: nflds_x2s = 0
  integer :: index_x2s_Sg_frac(glc_nec_max)   = 0 ! Fraction of glacier in glc MEC class 1
  integer :: index_x2s_Sg_topo(glc_nec_max)   = 0 ! Topo height in glc MEC class 1
  integer :: index_x2s_Fsgg_rofi(glc_nec_max) = 0 ! Ice runoff from glc model
  integer :: index_x2s_Fsgg_rofl(glc_nec_max) = 0 ! Liquid runoff from glc model
  integer :: index_x2s_Fsgg_hflx(glc_nec_max) = 0
#endif

CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_init_mct
!
! !DESCRIPTION:
!     vic lnd model init
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine lnd_init_mct( EClock, cdata, x2l, l2x, cdata_s, x2s, s2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2l, l2x
    type(seq_cdata)             , intent(inout) :: cdata_s
    type(mct_aVect)             , intent(inout) :: x2s, s2x
    character(len=*), optional  , intent(in)    :: NLFilename

    ! Local Variables
    integer(C_INT) :: errno
    type(vic_clock_type) :: vic_clock
    character(len=*, kind=C_CHAR)   :: vic_global_param_file
    integer                         :: LNDID
    integer                         :: mpicom_lnd
    integer                         :: mytask,ierr
    type(mct_gsMap),  pointer       :: GSMap_lnd
    type(mct_gGrid),  pointer       :: dom_lnd
    type(seq_infodata_type), pointer:: infodata
    character(len=CL) :: caseid
    character(len=CL) :: starttype
    character(len=CL) :: calendar
    integer :: dtime
    integer :: start_ymd
    integer :: start_yr
    integer :: start_mon
    integer :: start_day
    integer :: start_tod
    character(len=*), parameter     :: subname = "lnd_init_mct"

!EOP
!-------------------------------------------------------------------------------

   !--- initialize the field index values, subroutine is below in this file
   call cpl_indices_set()

   !--- point/get data from cdata datatype
   call seq_cdata_setptrs(cdata, ID=LNDID, mpicom=mpicom_lnd, &
        gsMap=GSMap_lnd, dom=dom_lnd, infodata=infodata)
   call MPI_COMM_RANK(mpicom_lnd,mytask,ierr)

   !--- copy/hand the mpicom from the driver to the vic model
   call ??(mpicom_lnd)

   !--- get unit number for vic log file
   !--- setup vic log file
   !--- set shr unit number to vic log file
   call shr_file_getLogUnit (shrlogunit)
   if (masterproc) then
      inquire(file='lnd_modelio.nml'//trim(inst_suffix),exist=exists)
      if (exists) then
         iulog = shr_file_getUnit()
         call shr_file_setIO('lnd_modelio.nml'//trim(inst_suffix),iulog)
      end if
      write(iulog,format) "VIC land model initialization"
   else
      iulog = shrlogunit
   end if
   call shr_file_setLogUnit (iulog)

   !--- get the casename, needed for output file names
   call seq_infodata_GetData( infodata, case_name=caseid)

   !--- get the starttype, seq_infodata_start_type_[start,cont,brnch] = clean start, restart, branch
   call seq_infodata_GetData( infodata, start_type=starttype)

   !--- get some clock info
   call seq_timemgr_EClockGetData(EClock,                               &
                                  curr_ymd=start_ymd,                  &
                                  curr_tod=start_tod,
                                  curr_yr=start_yr, curr_mon=start_mon, curr_day=start_day, &
                                  calendar=calendar )
   call seq_timemgr_EClockGetData(EClock, dtime=dtime )

   write(iulog,*)'EClock start_ymd = ',start_ymd,start_tod
   write(iulog,*)'EClock dtime = ',dtime

   !--- VIC global parameter file (namelist)
   if(present(NLFilename)) then
     vic_global_param_file = NLFilename
   else
     vic_global_param_file = "vic.globalconfig.txt"
   endif

   !--- Call the VIC init function
   errno = vic_cesm_init(vic_clock, vic_global_param_file)
   if (errno /= 0) then
     call shr_sys_abort(subname//':: vic_cesm_init returned a errno /= 0' )
   endif

   ! initialize the gsmap, inputs are
   !   gindex = 1d array of global indices on the local mpi task
   !   lsize = the size of gindex (the number of gridcells on this local mpi task)
   !   gsize = global size of grid (nx_global * ny_global)
   !   mpicom_lnd and LNDID from cdata above
   ! outputs are gsmap_lnd

   gindex(:) = from_vic_somewhere(:)
   lsize = size(gindex)
   gsize = global_nx * global_ny
   call mct_gsMap_init( gsMap_lnd, gindex, mpicom_lnd, LNDID, lsize, gsize )

   !--- initialize the dom, data in the dom is just local data of size lsize

   call mct_gGrid_init( GGrid=dom_lnd, CoordChars=trim(seq_flds_dom_coord), &
      OtherChars=trim(seq_flds_dom_other), lsize=lsize )
   call mct_gsMap_orderedPoints(gsMap_lnd, mytask, idata)
   call mct_gGrid_importIAttr(dom_lnd,'GlobGridNum',idata,lsize)
   call mct_gGrid_importRattr(dom_lnd,"lon", lon_data,lsize)
   call mct_gGrid_importRattr(dom_lnd,"lat", lat_data,lsize)
   call mct_gGrid_importRattr(dom_lnd,"lat",area_data,lsize)
   call mct_gGrid_importRattr(dom_lnd,"lat",mask_data,lsize)
   call mct_gGrid_importRattr(dom_lnd,"lat",frac_data,lsize)

   !--- intialize the attribute vectors

   call mct_aVect_init(x2l, rList=seq_flds_x2l_fields, lsize=lsize)
   call mct_aVect_zero(x2l)
   call mct_aVect_init(l2x, rList=seq_flds_l2x_fields, lsize=lsize)
   call mct_aVect_zero(l2x)

   !--- fill the l2x attribute vector with export data

   call lnd_export_mct( l2x )

   !--- fill some scalar export data

   call seq_infodata_PutData(cdata%infodata, &
        lnd_present=.true., lnd_prognostic=.true., &
        sno_present=.false., sno_prognostic=.false.)
   call seq_infodata_PutData( infodata, lnd_nx = nx_global, lnd_ny = ny_global)

   !--- set share log unit back to generic one
   call shr_file_setLogLevel(shrloglev)

end subroutine lnd_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_run_mct
!
! !DESCRIPTION:
!     vic lnd model run
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine lnd_run_mct( EClock, cdata, x2l, l2x, cdata_s, x2s, s2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2l
   type(mct_aVect)             ,intent(inout) :: l2x
   type(seq_cdata)             ,intent(inout) :: cdata_s
   type(mct_aVect)             ,intent(inout) :: x2s
   type(mct_aVect)             ,intent(inout) :: s2x

   ! Local Variables
   integer(C_INT) :: errno
   type(vic_clock_type) :: vic_clock
   character(len=*), parameter     :: subname = "lnd_run_mct"

!EOP
!-------------------------------------------------------------------------------

   !--- set vic log unit
   call shr_file_getLogUnit (shrlogunit)
   call shr_file_setLogUnit (iulog

   !--- point to cdata
   call seq_cdata_setptrs(cdata_l, infodata=infodata)

   !--- get time information
   call seq_timemgr_EClockGetData(EClock, &
        curr_ymd=ymd, curr_tod=tod_sync,  &
        curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)
   write(iulog,*)'EClock run',ymd,tod_sync,yr_sync,mon_sync,day_sync

   !--- get time of next radiation calc for time dependent albedo calc
   call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )

   !--- import data from coupler
   call lnd_import_mct( x2l)

   !--- run vic
   errno = vic_cesm_run(vic_clock)
   if (errno /= 0) then
     call shr_sys_abort(subname//':: vic_cesm_run returned a errno /= 0' )
   endif

   !--- export data to coupler
   call lnd_export_mct( l2x )

   !--- verify driver and vic are in sync, ymd and tod are vic yyyymmdd and sec
   if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
      call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
      write(iulog,*)' vic ymd=',ymd     ,'  clm tod= ',tod
      write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
      call shr_sys_flush(iulog)
      call shr_sys_abort( subname//":: VIC clock not in sync with Master Sync clock" )
   end if

   !--- set share log unit back to generic one
   call shr_file_setLogLevel(shrloglev)

end subroutine lnd_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_final_mct
!
! !DESCRIPTION:
!     vic lnd model finalize
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine lnd_final_mct( EClock, cdata, x2l, l2x, cdata_s, x2s, s2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2l
    type(mct_aVect)             ,intent(inout) :: l2x
    type(seq_cdata)             ,intent(inout) :: cdata_s
    type(mct_aVect)             ,intent(inout) :: x2s
    type(mct_aVect)             ,intent(inout) :: s2x

   ! Local Variables
   integer(C_INT) :: errno
   character(len=*), parameter     :: subname = "lnd_run_mct"
!EOP
!-------------------------------------------------------------------------------

    ! clean up
    errno = vic_cesm_final()

  if (errno /= 0) then
    call shr_sys_abort(subname//':: vic_cesm_final returned a errno /= 0' )
  endif

 end subroutine lnd_final_mct

!====================================================================================

  subroutine lnd_export_mct( l2x )

    !-----------------------------------------------------
    implicit none

    type(mct_aVect)   , intent(inout) :: l2x
    integer :: i,lsize
    character(len=*),parameter :: subname = 'lnd_export_mct'

    !-----------------------------------------------------

    lsize = mct_avect_lsize(l2x)
    l2x%rAttr(:,:) = 0.0_r8

    ! ccsm sign convention is that fluxes are positive downward

!dir$ concurrent
    do i = 1,lsize
!      l2x%rAttr(index_l2x_Sl_landfrac,i) =  adomain%frac(i)
       l2x%rAttr(index_l2x_Sl_t,i)        =  vic_data%t_rad(i)
       l2x%rAttr(index_l2x_Sl_snowh,i)    =  vic_data%h2osno(i)
       l2x%rAttr(index_l2x_Sl_avsdr,i)    =  vic_data%albvsdr(i)
       l2x%rAttr(index_l2x_Sl_anidr,i)    =  vic_data%albnidr(i)
       l2x%rAttr(index_l2x_Sl_avsdf,i)    =  vic_data%albvsdf(i)
       l2x%rAttr(index_l2x_Sl_anidf,i)    =  vic_data%albnidf(i)
       l2x%rAttr(index_l2x_Sl_tref,i)     =  vic_data%t_ref2m(i)
       l2x%rAttr(index_l2x_Sl_qref,i)     =  vic_data%q_ref2m(i)
       l2x%rAttr(index_l2x_Sl_logz0,i)    =  vic_data%logz0(i)
       l2x%rAttr(index_l2x_Flrl_rofice,i) =  vic_data%rofice(i)
       l2x%rAttr(index_l2x_Flrl_rofliq,i) =  vic_data%rofliq(i)
       l2x%rAttr(index_l2x_Fall_taux,i)   = -vic_data%taux(i)
       l2x%rAttr(index_l2x_Fall_tauy,i)   = -vic_data%tauy(i)
       l2x%rAttr(index_l2x_Fall_lat,i)    = -vic_data%eflx_lh_tot(i)
       l2x%rAttr(index_l2x_Fall_sen,i)    = -vic_data%eflx_sh_tot(i)
       l2x%rAttr(index_l2x_Fall_lwup,i)   = -vic_data%eflx_lwrad_out(i)
       l2x%rAttr(index_l2x_Fall_evap,i)   = -vic_data%qflx_evap_tot(i)
       l2x%rAttr(index_l2x_Fall_swnet,i)  =  vic_data%fsa(i)

       !--- optional fields ---
       if (index_l2x_Fall_fco2_lnd /= 0) l2x%rAttr(index_l2x_Fall_fco2_lnd,i) = -vic_data%nee(i)
       if (index_l2x_Sl_fv /= 0 ) l2x%rAttr(index_l2x_Sl_fv,i)= vic_data%fv(i)
       if (index_l2x_Sl_ram1 /= 0 )  l2x%rAttr(index_l2x_Sl_ram1,i) = vic_data%ram1(i)
       if (index_l2x_Sl_fv   /= 0 )  l2x%rAttr(index_l2x_Sl_fv,i)   = vic_data%fv(i)
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x%rAttr(index_l2x_Fall_flxdst1,i)= -vic_data%flxdst(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x%rAttr(index_l2x_Fall_flxdst2,i)= -vic_data%flxdst(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x%rAttr(index_l2x_Fall_flxdst3,i)= -vic_data%flxdst(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x%rAttr(index_l2x_Fall_flxdst4,i)= -vic_data%flxdst(g,4)
    end do

  end subroutine lnd_export_mct

!====================================================================================

  subroutine lnd_import_mct( x2l)

    !-----------------------------------------------------
    implicit none
    !
    ! Arguments
    !
    type(mct_aVect)   , intent(inout) :: x2l
    !
    ! Local Variables
    !
    integer  :: i,lsize
    real(r8) :: e                    !vapor pressure (Pa)
    real(r8) :: qsat                 !saturation specific humidity (kg/kg)
    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: forc_rainl           ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc           ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl           ! snowfxl Atm flux  mm/s
    real(r8) :: esatw                !saturation vapor pressure over water (Pa)
    real(r8) :: esati                !saturation vapor pressure over ice (Pa)
    real(r8) :: a0,a1,a2,a3,a4,a5,a6 !coefficients for esat over water
    real(r8) :: b0,b1,b2,b3,b4,b5,b6 !coefficients for esat over ice
    real(r8) :: tdc, t               !Kelvins to Celcius function and its input
    character(len=32), parameter :: sub = 'lnd_import_mct'

    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
               a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
               a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
               a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
               b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
               b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
               b6=1.838826904e-10_r8)
!
! function declarations
!
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))

    !-----------------------------------------------------

    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.

    lsize = mct_avect_lsize(x2l)

!dir$ concurrent
    do i = 1,lsize
        ! Determine required receive fields

        vic_data%forc_hgt(g)     = x2l%rAttr(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
        vic_data%forc_u(g)       = x2l%rAttr(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
        vic_data%forc_v(g)       = x2l%rAttr(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
        vic_data%forc_th(g)      = x2l%rAttr(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
        vic_data%forc_q(g)       = x2l%rAttr(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
        vic_data%forc_pbot(g)    = x2l%rAttr(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
        vic_data%forc_t(g)       = x2l%rAttr(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
        vic_data%forc_lwrad(g)   = x2l%rAttr(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2
        forc_rainc          = x2l%rAttr(index_x2l_Faxa_rainc,i)   ! mm/s
        forc_rainl          = x2l%rAttr(index_x2l_Faxa_rainl,i)   ! mm/s
        forc_snowc          = x2l%rAttr(index_x2l_Faxa_snowc,i)   ! mm/s
        forc_snowl          = x2l%rAttr(index_x2l_Faxa_snowl,i)   ! mm/s
        vic_data%forc_rain(g)    = forc_rainc + forc_rainl
        vic_data%forc_snow(g)    = forc_snowc + forc_snowl
        vic_data%forc_solndr(g)  = x2l%rAttr(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
        vic_data%forc_solvdr(g)  = x2l%rAttr(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
        vic_data%forc_solndf(g)  = x2l%rAttr(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
        vic_data%forc_solvdf(g)  = x2l%rAttr(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

! tcraig - not sure how much of this is needed for vic

        ! Determine optional receive fields

        if (index_x2l_Sa_co2prog /= 0) then
           co2_ppmv_prog = x2l%rAttr(index_x2l_Sa_co2prog,i)   ! co2 atm state prognostic
        else
           co2_ppmv_prog = co2_ppmv
        end if

        if (index_x2l_Sa_co2diag /= 0) then
           co2_ppmv_diag = x2l%rAttr(index_x2l_Sa_co2diag,i)   ! co2 atm state diagnostic
        else
           co2_ppmv_diag = co2_ppmv
        end if

        ! Determine derived quantities for required fields
        ! First, set forcing height to maximum of atmospheric model forcing height
        ! and a prescribed minimum height

	vic_data%forc_hgt(g)   = max(vic_data%forc_hgt(g), forc_hgt_min)
        vic_data%forc_hgt_u(g) = vic_data%forc_hgt(g)    !observational height of wind [m]
        vic_data%forc_hgt_t(g) = vic_data%forc_hgt(g)    !observational height of temperature [m]
        vic_data%forc_hgt_q(g) = vic_data%forc_hgt(g)    !observational height of humidity [m]
        vic_data%forc_vp(g)    = vic_data%forc_q(g) * vic_data%forc_pbot(g) &
                            / (0.622_r8 + 0.378_r8 * vic_data%forc_q(g))
        vic_data%forc_rho(g)   = (vic_data%forc_pbot(g) - 0.378_r8 * vic_data%forc_vp(g)) &
                            / (rair * vic_data%forc_t(g))
        vic_data%forc_po2(g)   = o2_molar_const * vic_data%forc_pbot(g)
        vic_data%forc_wind(g)  = sqrt(vic_data%forc_u(g)**2 + vic_data%forc_v(g)**2)
        vic_data%forc_solar(g) = vic_data%forc_solad(g,1) + vic_data%forc_solai(g,1) + &
                            vic_data%forc_solad(g,2) + vic_data%forc_solai(g,2)
        vic_data%rainf    (g)  = vic_data%forc_rain(g) + vic_data%forc_snow(g)

        if (vic_data%forc_t(g) > SHR_CONST_TKFRZ) then
           e = esatw(tdc(vic_data%forc_t(g)))
        else
           e = esati(tdc(vic_data%forc_t(g)))
        end if
        qsat           = 0.622_r8*e / (vic_data%forc_pbot(g) - 0.378_r8*e)
        vic_data%forc_rh(g) = 100.0_r8*(vic_data%forc_q(g) / qsat)
        ! Make sure relative humidity is properly bounded
        ! vic_data%forc_rh(g) = min( 100.0_r8, vic_data%forc_rh(g) )
        ! vic_data%forc_rh(g) = max(   0.0_r8, vic_data%forc_rh(g) )

        ! Determine derived quantities for optional fields
        ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
        ! Note that forc_pbot is in Pa

        if (co2_type_idx == 1) then
           co2_ppmv_val = co2_ppmv_prog
        else if (co2_type_idx == 2) then
           co2_ppmv_val = co2_ppmv_diag
        else
           co2_ppmv_val = co2_ppmv
        end if
        vic_data%forc_pco2(g)   = co2_ppmv_val * 1.e-6_r8 * vic_data%forc_pbot(g)
        vic_data%forc_pc13o2(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * vic_data%forc_pbot(g)

     end do

   end subroutine lnd_import_mct

!===============================================================================
!===============================================================================
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: cpl_indices_set
!
! !CALLED FROM: lnd_comp_mct or lnd_comp_esmf
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
! !INTERFACE:
  subroutine cpl_indices_set( )

!
! !DESCRIPTION: 
! Set the coupler indices needed by the land model coupler
! interface.
!
! !USES:
  use seq_flds_mod  , only: seq_flds_x2l_fields, seq_flds_l2x_fields,     &
                            seq_flds_x2s_fields, seq_flds_s2x_fields
  use mct_mod       , only: mct_aVect, mct_aVect_init, mct_avect_indexra, &
                            mct_aVect_clean, mct_avect_nRattr
  use seq_drydep_mod, only: drydep_fields_token, lnd_drydep
  use shr_megan_mod,  only: shr_megan_fields_token, shr_megan_mechcomps_n

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! 01/2011, Erik Kluzek:         Added protex headers
!
! !LOCAL VARIABLES:
    type(mct_aVect)   :: l2x      ! temporary, land to coupler
    type(mct_aVect)   :: x2l      ! temporary, coupler to land
    type(mct_aVect)   :: s2x      ! temporary, glacier to coupler
    type(mct_aVect)   :: x2s      ! temporary, coupler to glacier
    integer           :: num 
    character(len= 2) :: cnum
    character(len=64) :: name
    character(len=32) :: subname = 'cpl_indices_set'  ! subroutine name
!EOP
!
!-----------------------------------------------------------------------

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2l, rList=seq_flds_x2l_fields, lsize=1)
    call mct_aVect_init(l2x, rList=seq_flds_l2x_fields, lsize=1)

    !-------------------------------------------------------------
    ! vic -> drv 
    !-------------------------------------------------------------

    index_l2x_Flrl_rofliq   = mct_avect_indexra(l2x,'Flrl_rofliq')
    index_l2x_Flrl_rofice   = mct_avect_indexra(l2x,'Flrl_rofice')

    index_l2x_Sl_t          = mct_avect_indexra(l2x,'Sl_t')
    index_l2x_Sl_snowh      = mct_avect_indexra(l2x,'Sl_snowh')
    index_l2x_Sl_avsdr      = mct_avect_indexra(l2x,'Sl_avsdr')
    index_l2x_Sl_anidr      = mct_avect_indexra(l2x,'Sl_anidr')
    index_l2x_Sl_avsdf      = mct_avect_indexra(l2x,'Sl_avsdf')
    index_l2x_Sl_anidf      = mct_avect_indexra(l2x,'Sl_anidf')
    index_l2x_Sl_tref       = mct_avect_indexra(l2x,'Sl_tref')
    index_l2x_Sl_qref       = mct_avect_indexra(l2x,'Sl_qref')
    index_l2x_Sl_u10        = mct_avect_indexra(l2x,'Sl_u10')
    index_l2x_Sl_ram1       = mct_avect_indexra(l2x,'Sl_ram1')
    index_l2x_Sl_fv         = mct_avect_indexra(l2x,'Sl_fv')
    index_l2x_Sl_soilw      = mct_avect_indexra(l2x,'Sl_soilw',perrwith='quiet')
    index_l2x_Sl_logz0      = mct_avect_indexra(l2x,'Sl_logz0')
    if ( lnd_drydep )then
       index_l2x_Sl_ddvel = mct_avect_indexra(l2x, trim(drydep_fields_token))
    else
       index_l2x_Sl_ddvel = 0
    end if

    index_l2x_Fall_taux     = mct_avect_indexra(l2x,'Fall_taux')
    index_l2x_Fall_tauy     = mct_avect_indexra(l2x,'Fall_tauy')
    index_l2x_Fall_lat      = mct_avect_indexra(l2x,'Fall_lat')
    index_l2x_Fall_sen      = mct_avect_indexra(l2x,'Fall_sen')
    index_l2x_Fall_lwup     = mct_avect_indexra(l2x,'Fall_lwup')
    index_l2x_Fall_evap     = mct_avect_indexra(l2x,'Fall_evap')
    index_l2x_Fall_swnet    = mct_avect_indexra(l2x,'Fall_swnet')
    index_l2x_Fall_flxdst1  = mct_avect_indexra(l2x,'Fall_flxdst1')
    index_l2x_Fall_flxdst2  = mct_avect_indexra(l2x,'Fall_flxdst2')
    index_l2x_Fall_flxdst3  = mct_avect_indexra(l2x,'Fall_flxdst3')
    index_l2x_Fall_flxdst4  = mct_avect_indexra(l2x,'Fall_flxdst4')

    index_l2x_Fall_fco2_lnd = mct_avect_indexra(l2x,'Fall_fco2_lnd',perrwith='quiet')

    ! MEGAN fluxes
    if (shr_megan_mechcomps_n>0) then
       index_l2x_Fall_flxvoc = mct_avect_indexra(l2x,trim(shr_megan_fields_token))
    else
       index_l2x_Fall_flxvoc = 0
    endif

    nflds_l2x = mct_avect_nRattr(l2x)

    !-------------------------------------------------------------
    ! drv -> vic
    !-------------------------------------------------------------

    index_x2l_Sa_z          = mct_avect_indexra(x2l,'Sa_z')
    index_x2l_Sa_u          = mct_avect_indexra(x2l,'Sa_u')
    index_x2l_Sa_v          = mct_avect_indexra(x2l,'Sa_v')
    index_x2l_Sa_ptem       = mct_avect_indexra(x2l,'Sa_ptem')
    index_x2l_Sa_pbot       = mct_avect_indexra(x2l,'Sa_pbot')
    index_x2l_Sa_tbot       = mct_avect_indexra(x2l,'Sa_tbot')
    index_x2l_Sa_shum       = mct_avect_indexra(x2l,'Sa_shum')
    index_x2l_Sa_co2prog    = mct_avect_indexra(x2l,'Sa_co2prog',perrwith='quiet')
    index_x2l_Sa_co2diag    = mct_avect_indexra(x2l,'Sa_co2diag',perrwith='quiet')

    index_x2l_Faxa_lwdn     = mct_avect_indexra(x2l,'Faxa_lwdn')
    index_x2l_Faxa_rainc    = mct_avect_indexra(x2l,'Faxa_rainc')
    index_x2l_Faxa_rainl    = mct_avect_indexra(x2l,'Faxa_rainl')
    index_x2l_Faxa_snowc    = mct_avect_indexra(x2l,'Faxa_snowc')
    index_x2l_Faxa_snowl    = mct_avect_indexra(x2l,'Faxa_snowl')
    index_x2l_Faxa_swndr    = mct_avect_indexra(x2l,'Faxa_swndr')
    index_x2l_Faxa_swvdr    = mct_avect_indexra(x2l,'Faxa_swvdr')
    index_x2l_Faxa_swndf    = mct_avect_indexra(x2l,'Faxa_swndf')
    index_x2l_Faxa_swvdf    = mct_avect_indexra(x2l,'Faxa_swvdf')
    index_x2l_Faxa_bcphidry = mct_avect_indexra(x2l,'Faxa_bcphidry')
    index_x2l_Faxa_bcphodry = mct_avect_indexra(x2l,'Faxa_bcphodry')
    index_x2l_Faxa_bcphiwet = mct_avect_indexra(x2l,'Faxa_bcphiwet')
    index_x2l_Faxa_ocphidry = mct_avect_indexra(x2l,'Faxa_ocphidry')
    index_x2l_Faxa_ocphodry = mct_avect_indexra(x2l,'Faxa_ocphodry')
    index_x2l_Faxa_ocphiwet = mct_avect_indexra(x2l,'Faxa_ocphiwet')
    index_x2l_Faxa_dstdry1  = mct_avect_indexra(x2l,'Faxa_dstdry1')
    index_x2l_Faxa_dstdry2  = mct_avect_indexra(x2l,'Faxa_dstdry2')
    index_x2l_Faxa_dstdry3  = mct_avect_indexra(x2l,'Faxa_dstdry3')
    index_x2l_Faxa_dstdry4  = mct_avect_indexra(x2l,'Faxa_dstdry4')
    index_x2l_Faxa_dstwet1  = mct_avect_indexra(x2l,'Faxa_dstwet1')
    index_x2l_Faxa_dstwet2  = mct_avect_indexra(x2l,'Faxa_dstwet2')
    index_x2l_Faxa_dstwet3  = mct_avect_indexra(x2l,'Faxa_dstwet3')
    index_x2l_Faxa_dstwet4  = mct_avect_indexra(x2l,'Faxa_dstwet4')

    index_x2l_Flrr_flood    = mct_avect_indexra(x2l,'Flrr_flood')

    nflds_x2l = mct_avect_nRattr(x2l)

    call mct_aVect_clean(x2l)
    call mct_aVect_clean(l2x)

! tcraig - this can probably be removed
#if (1 == 0)
    !-------------------------------------------------------------
    ! drv->sno (for cism coupling)
    !-------------------------------------------------------------

    glc_nec = 0

    if (seq_flds_x2s_fields /= ' ') then
       call mct_aVect_init(x2s, rList=seq_flds_x2s_fields, lsize=1)

       do num = 1,glc_nec_max
          write(cnum,'(i2.2)') num
          name = 'Sg_frac' // cnum
          index_x2s_Sg_frac(num)   = mct_avect_indexra(x2s,trim(name),perrwith='quiet') 
          name = 'Sg_topo' // cnum
          index_x2s_Sg_topo(num)   = mct_avect_indexra(x2s,trim(name),perrwith='quiet')
          name = 'Fsgg_rofi' // cnum
          index_x2s_Fsgg_rofi(num) = mct_avect_indexra(x2s,trim(name),perrwith='quiet')
          name = 'Fsgg_rofl' // cnum
          index_x2s_Fsgg_rofl(num) = mct_avect_indexra(x2s,trim(name),perrwith='quiet')
          name = 'Fsgg_hflx' // cnum
          index_x2s_Fsgg_hflx(num) = mct_avect_indexra(x2s,trim(name),perrwith='quiet')
          if ( index_x2s_Sg_frac(num)   == 0 .and. &
               index_x2s_Sg_topo(num)   == 0 .and. &
               index_x2s_Fsgg_rofi(num) == 0 .and. &
               index_x2s_fsgg_rofl(num) == 0 .and. &
               index_x2s_Fsgg_hflx(num) == 0 ) then
             exit
          end if
          glc_nec = num
       end do
       if (glc_nec == glc_nec_max) then
          call shr_sys_abort (subname // 'error: glc_nec_cpl cannot equal glc_nec_max')
       end if

       nflds_x2s = mct_avect_nRattr(x2s)
       call mct_aVect_clean(x2s)
    end if

    !-------------------------------------------------------------
    ! sno -> drv (for cism coupling)
    !-------------------------------------------------------------

    if (seq_flds_s2x_fields /= ' ') then
       call mct_aVect_init(s2x, rList=seq_flds_s2x_fields, lsize=1)

       do num = 1,glc_nec
          write(cnum,'(i2.2)') num

          name = 'Ss_tsrf' // cnum
          index_s2x_Ss_tsrf(num)   = mct_avect_indexra(s2x,trim(name))
          name = 'Ss_topo' // cnum
          index_s2x_Ss_topo(num)   = mct_avect_indexra(s2x,trim(name))
          name = 'Fgss_qice' // cnum
          index_s2x_Fgss_qice(num) = mct_avect_indexra(s2x,trim(name))
       end do

       nflds_s2x = mct_avect_nRattr(s2x)
       call mct_aVect_clean(s2x)
    end if
#endif

  end subroutine cpl_indices_set

!=======================================================================

end module lnd_comp_mct
