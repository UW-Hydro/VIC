MODULE lnd_comp_mct

  !--------------------------------------------------------------------------
  ! USES:
  !--------------------------------------------------------------------------
  USE mct_mod
  USE esmf
  USE seq_cdata_mod
  USE seq_infodata_mod
  USE shr_kind_mod, ONLY : r8 => shr_kind_r8
  USE vic_cesm_interface
  USE vic_cesm_def_mod
  USE vic_cesm_print_library
  USE, INTRINSIC :: iso_c_binding

  ! PUBLIC TYPES:
  IMPLICIT NONE
  SAVE
  PRIVATE

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------
  PUBLIC :: lnd_init_mct
  PUBLIC :: lnd_run_mct
  PUBLIC :: lnd_final_mct

  !--------------------------------------------------------------------------
  ! Module level variables
  !--------------------------------------------------------------------------
  INTEGER :: iulog   ! vic log file unit number
  INTEGER :: shrlogunit, shrloglev   ! generic log unit
  INTEGER :: year, month, day, seconds
  TYPE(vic_clock), TARGET :: vclock
  TYPE(case_metadata), TARGET :: cmeta

  ! Access global variables in VIC C
  TYPE(domain_struct), bind(C, name='local_domain') :: local_domain
  TYPE(domain_struct), bind(C, name='global_domain') :: global_domain

  TYPE(C_PTR), bind(C, name='x2l_vic') :: x2l_vic
  TYPE(C_PTR), bind(C, name='l2x_vic') :: l2x_vic

  TYPE(x2l_data_struct), DIMENSION(:), POINTER :: x2l_vic_ptr
  TYPE(l2x_data_struct), DIMENSION(:), POINTER :: l2x_vic_ptr

  !--- lnd -> drv
  INTEGER :: nflds_l2x = 0
  INTEGER :: index_l2x_Sl_t           = 0 ! temperature
  INTEGER :: index_l2x_Sl_tref        = 0 ! 2m reference temperature
  INTEGER :: index_l2x_Sl_qref        = 0 ! 2m reference specific humidity
  INTEGER :: index_l2x_Sl_avsdr       = 0 ! albedo: direct , visible
  INTEGER :: index_l2x_Sl_anidr       = 0 ! albedo: direct , near-ir
  INTEGER :: index_l2x_Sl_avsdf       = 0 ! albedo: diffuse, visible
  INTEGER :: index_l2x_Sl_anidf       = 0 ! albedo: diffuse, near-ir
  INTEGER :: index_l2x_Sl_snowh       = 0 ! snow height
  INTEGER :: index_l2x_Sl_u10         = 0 ! 10m wind
  INTEGER :: index_l2x_Sl_ddvel       = 0 ! dry deposition velocities (optional)
  INTEGER :: index_l2x_Sl_fv          = 0 ! friction velocity
  INTEGER :: index_l2x_Sl_ram1        = 0 ! aerodynamical resistance
  INTEGER :: index_l2x_Sl_logz0       = 0 ! log z0
  INTEGER :: index_l2x_Fall_taux      = 0 ! wind stress, zonal
  INTEGER :: index_l2x_Fall_tauy      = 0 ! wind stress, meridional
  INTEGER :: index_l2x_Fall_lat       = 0 ! latent          heat flux
  INTEGER :: index_l2x_Fall_sen       = 0 ! sensible        heat flux
  INTEGER :: index_l2x_Fall_lwup      = 0 ! upward longwave heat flux
  INTEGER :: index_l2x_Fall_evap      = 0 ! evaporation     water flux
  INTEGER :: index_l2x_Fall_swnet     = 0 ! heat flux       shortwave net
  INTEGER :: index_l2x_Fall_fco2_lnd  = 0 ! co2 flux **For testing set to 0
  INTEGER :: index_l2x_Fall_flxdst1   = 0 ! dust flux size bin 1
  INTEGER :: index_l2x_Fall_flxdst2   = 0 ! dust flux size bin 2
  INTEGER :: index_l2x_Fall_flxdst3   = 0 ! dust flux size bin 3
  INTEGER :: index_l2x_Fall_flxdst4   = 0 ! dust flux size bin 4
  INTEGER :: index_l2x_Fall_flxvoc    = 0 ! MEGAN fluxes
  INTEGER :: index_l2x_Flrl_rofliq    = 0 ! lnd->rtm input fluxes
  INTEGER :: index_l2x_Flrl_rofice    = 0 ! lnd->rtm input fluxes

  !--- drv -> lnd
  INTEGER :: nflds_x2l = 0
  INTEGER :: index_x2l_Sa_z           = 0 ! bottom atm level height
  INTEGER :: index_x2l_Sa_u           = 0 ! bottom atm level zon wind
  INTEGER :: index_x2l_Sa_v           = 0 ! bottom atm level mer wind
  INTEGER :: index_x2l_Sa_ptem        = 0 ! bottom atm level pot temp
  INTEGER :: index_x2l_Sa_shum        = 0 ! bottom atm level spec hum
  INTEGER :: index_x2l_Sa_pbot        = 0 ! bottom atm level pressure
  INTEGER :: index_x2l_Sa_tbot        = 0 ! bottom atm level temp
  INTEGER :: index_x2l_Faxa_lwdn      = 0 ! downward lw heat flux
  INTEGER :: index_x2l_Faxa_rainc     = 0 ! prec: liquid 'convective'
  INTEGER :: index_x2l_Faxa_rainl     = 0 ! prec: liquid 'large scale'
  INTEGER :: index_x2l_Faxa_snowc     = 0 ! prec: frozen 'convective'
  INTEGER :: index_x2l_Faxa_snowl     = 0 ! prec: frozen 'large scale'
  INTEGER :: index_x2l_Faxa_swndr     = 0 ! sw: nir direct  downward
  INTEGER :: index_x2l_Faxa_swvdr     = 0 ! sw: vis direct  downward
  INTEGER :: index_x2l_Faxa_swndf     = 0 ! sw: nir diffuse downward
  INTEGER :: index_x2l_Faxa_swvdf     = 0 ! sw: vis diffuse downward
  INTEGER :: index_x2l_Sa_co2prog     = 0 ! bottom atm level prognostic co2
  INTEGER :: index_x2l_Sa_co2diag     = 0 ! bottom atm level diagnostic co2
  INTEGER :: index_x2l_Faxa_bcphidry  = 0 ! flux: Black Carbon hydrophilic dry deposition
  INTEGER :: index_x2l_Faxa_bcphodry  = 0 ! flux: Black Carbon hydrophobic dry deposition
  INTEGER :: index_x2l_Faxa_bcphiwet  = 0 ! flux: Black Carbon hydrophilic wet deposition
  INTEGER :: index_x2l_Faxa_ocphidry  = 0 ! flux: Organic Carbon hydrophilic dry deposition
  INTEGER :: index_x2l_Faxa_ocphodry  = 0 ! flux: Organic Carbon hydrophobic dry deposition
  INTEGER :: index_x2l_Faxa_ocphiwet  = 0 ! flux: Organic Carbon hydrophilic dry deposition
  INTEGER :: index_x2l_Faxa_dstwet1   = 0 ! flux: Size 1 dust -- wet deposition
  INTEGER :: index_x2l_Faxa_dstwet2   = 0 ! flux: Size 2 dust -- wet deposition
  INTEGER :: index_x2l_Faxa_dstwet3   = 0 ! flux: Size 3 dust -- wet deposition
  INTEGER :: index_x2l_Faxa_dstwet4   = 0 ! flux: Size 4 dust -- wet deposition
  INTEGER :: index_x2l_Faxa_dstdry1   = 0 ! flux: Size 1 dust -- dry deposition
  INTEGER :: index_x2l_Faxa_dstdry2   = 0 ! flux: Size 2 dust -- dry deposition
  INTEGER :: index_x2l_Faxa_dstdry3   = 0 ! flux: Size 3 dust -- dry deposition
  INTEGER :: index_x2l_Faxa_dstdry4   = 0 ! flux: Size 4 dust -- dry deposition
  INTEGER :: index_x2l_Flrr_flood     = 0 ! rtm->lnd rof (flood) flux

CONTAINS

  !--------------------------------------------------------------------------
  !> @brief   vic lnd model init
  !--------------------------------------------------------------------------
  SUBROUTINE lnd_init_mct(EClock, cdata, x2l, l2x, cdata_s, x2s, s2x, NLFilename)

    !--- Uses:
    USE shr_file_mod          , ONLY : shr_file_getUnit, shr_file_setLogLevel, &
                                       shr_file_getLogUnit, shr_file_getLogUnit, &
                                       shr_file_setLogUnit, shr_file_setIO
    USE seq_timemgr_mod       , ONLY : seq_timemgr_EClockGetData
    USE mct_mod               , ONLY : mct_gsMap, mct_gGrid, mct_gGrid_importIAttr, &
                                       mct_gGrid_importRAttr, mct_gGrid_init, &
                                       mct_gsMap_orderedPoints
    USE seq_infodata_mod      , ONLY : seq_infodata_PutData, seq_infodata_PutData
    USE seq_flds_mod
    USE vic_cesm_def_mod
    USE vic_cesm_print_library, ONLY : set_print_library_iulog, print_vic_clock

    IMPLICIT NONE

    !--- INPUT/OUTPUT PARAMETERS:
    TYPE(ESMF_Clock)            , INTENT(in)    :: EClock
    TYPE(seq_cdata)             , INTENT(inout) :: cdata
    TYPE(mct_aVect)             , INTENT(inout) :: x2l, l2x
    TYPE(seq_cdata)             , INTENT(inout) :: cdata_s
    TYPE(mct_aVect)             , INTENT(inout) :: x2s, s2x
    CHARACTER(len=*), OPTIONAL  , INTENT(in)    :: NLFilename

    !--- Local Variables
    INTEGER                                     :: LNDID
    INTEGER                                     :: mpicom_lnd
    INTEGER                                     :: mytask, ierr
    INTEGER                                     :: lsize, gsize
    INTEGER                                     :: nx_global, ny_global
    INTEGER                                     :: dtime  ! timestep (seconds)
    INTEGER, POINTER                            :: idata(:)
    INTEGER, DIMENSION(:), ALLOCATABLE          :: gindex  ! Number the local grid points
    REAL(r8), POINTER                           :: lon_data(:)
    REAL(r8), POINTER                           :: lat_data(:)
    REAL(r8), POINTER                           :: area_data(:)
    REAL(r8), POINTER                           :: mask_data(:)
    REAL(r8), POINTER                           :: frac_data(:)
    TYPE(mct_gsMap),  POINTER                   :: GSMap_lnd
    TYPE(mct_gGrid),  POINTER                   :: dom_lnd
    TYPE(seq_infodata_type), POINTER            :: infodata
    CHARACTER(KIND=C_CHAR, len=VICMAXSTRING)    :: caseid
    CHARACTER(KIND=C_CHAR, len=VICMAXSTRING)    :: casedesc
    CHARACTER(KIND=C_CHAR, len=VICMAXSTRING)    :: starttype
    CHARACTER(KIND=C_CHAR, len=VICMAXSTRING)    :: version
    CHARACTER(KIND=C_CHAR, len=VICMAXSTRING)    :: hostname
    CHARACTER(KIND=C_CHAR, len=VICMAXSTRING)    :: username
    CHARACTER(KIND=C_CHAR, len=VICMAXSTRING)    :: calendar
    LOGICAL                                     :: exists  ! true if file exists
    CHARACTER(len=*), PARAMETER                 :: subname = '(lnd_init_mct)'

    !--- Local Variables for C/Fortran Interface
    INTEGER(C_INT)                              :: errno
    CHARACTER(len=VICMAXSTRING, kind=C_CHAR)    :: vic_global_param_file

    !--- initialize the VIC logging to stderr (for now)
    CALL initialize_log()

    !--- initialize the field index values, subroutine is below in this file
    CALL cpl_indices_set()

    !--- point/get data from cdata datatype
    CALL seq_cdata_setptrs(cdata, ID=LNDID, mpicom=mpicom_lnd, &
         gsMap=GSMap_lnd, dom=dom_lnd, infodata=infodata)
    CALL MPI_COMM_RANK(mpicom_lnd, mytask, ierr)

    !--- copy/hand the mpicom from the driver to the vic model
    CALL initialize_vic_cesm_mpi(mpicom_lnd)

    !--- get unit number for vic log file
    !--- setup vic log file
    !--- set shr unit number to vic log file
    CALL shr_file_getLogUnit(shrlogunit)
    IF (mytask == 0) THEN
       iulog = shr_file_getUnit()
       CALL shr_file_setIO('lnd_modelio.nml', iulog)
       WRITE(iulog, *) subname,': VIC land model initialization'
    ELSE
       iulog = shrlogunit
    ENDIF
    CALL set_print_library_iulog(iulog)
    CALL shr_file_setLogUnit(iulog)

    !--- get the casename, needed for output file names
    CALL seq_infodata_GetData(infodata, case_name=caseid)

    !--- get the starttype, seq_infodata_start_type_[start,cont,brnch] = clean start, restart, branch
    CALL seq_infodata_GetData(infodata, start_type=starttype)

    !--- get other metadata
    CALL seq_infodata_GetData(infodata, case_desc=casedesc, model_version=version, &
                              hostname=hostname, username=username)

    cmeta%caseid = Copy_s2a(caseid)
    cmeta%casedesc = Copy_s2a(casedesc)
    cmeta%starttype = Copy_s2a(starttype)
    cmeta%model_version = Copy_s2a(version)
    cmeta%hostname = Copy_s2a(hostname)
    cmeta%username = Copy_s2a(username)

    !--- get some clock info
    CALL seq_timemgr_EClockGetData(EClock, &
         dtime=dtime, &
         curr_yr=year, &
         curr_mon=month, &
         curr_day=day, &
         curr_tod=seconds, &
         calendar=calendar)

    !--- initialize clock and add starting date / time info
    vclock%timestep = dtime
    vclock%current_year = year
    vclock%current_month = month
    vclock%current_day = day
    vclock%current_dayseconds = seconds
    vclock%state_flag = .FALSE.
    vclock%stop_flag = .FALSE.
    vclock%calendar = Copy_s2a(calendar)

    !--- VIC global parameter file (namelist)
    ! commenting out because it was yielding a drv_in <-- need to rethink this
    !IF(PRESENT(NLFilename)) THEN
    !   vic_global_param_file = NLFilename
    !ELSE
    !   vic_global_param_file = 'vic.globalconfig.txt'
    !ENDIF
    vic_global_param_file = 'vic.globalconfig.txt'

    !--- Call the VIC init function
    errno = vic_cesm_init(vclock, cmeta)
    IF (errno /= 0) THEN
       CALL shr_sys_abort(subname//' ERROR: vic_cesm_init returned a errno /= 0')
    ENDIF

    nx_global = global_domain%n_nx
    ny_global = global_domain%n_ny

    ! initialize the gsmap, inputs are
    !   gindex = 1d array of global indices on the local mpi task
    !   lsize = the size of gindex (the number of gridcells on this local mpi task)
    !   gsize = global size of grid (nx_global * ny_global)
    !   mpicom_lnd and LNDID from cdata above
    ! outputs are gsmap_lnd
    lsize = local_domain%ncells_active
    gsize = nx_global * ny_global

    ! allocate domain arrays
    ALLOCATE(gindex(lsize))
    ALLOCATE(lon_data(lsize))
    ALLOCATE(lat_data(lsize))
    ALLOCATE(area_data(lsize))
    ALLOCATE(mask_data(lsize))
    ALLOCATE(frac_data(lsize))

    !--- Unpack the local domain structure
    CALL unpack_vic_domain(gindex, lon_data, lat_data, area_data, mask_data, frac_data)

    CALL mct_gsMap_init(gsMap_lnd, gindex, mpicom_lnd, LNDID, lsize, gsize)

    !-- setup mappings for l2x and x2l structures
    CALL c_f_pointer(l2x_vic, l2x_vic_ptr, [local_domain%ncells_active])
    CALL c_f_pointer(x2l_vic, x2l_vic_ptr, [local_domain%ncells_active])

    !--- initialize the dom, data in the dom is just local data of size lsize
    CALL mct_gGrid_init(GGrid=dom_lnd, CoordChars=TRIM(seq_flds_dom_coord), &
                        OtherChars=TRIM(seq_flds_dom_other), lsize=lsize)
    CALL mct_gsMap_orderedPoints(gsMap_lnd, mytask, idata)
    CALL mct_gGrid_importIAttr(dom_lnd, 'GlobGridNum', idata, lsize)
    CALL mct_gGrid_importRattr(dom_lnd, 'lon', lon_data, lsize)
    CALL mct_gGrid_importRattr(dom_lnd, 'lat', lat_data, lsize)
    CALL mct_gGrid_importRattr(dom_lnd, 'area', area_data, lsize)
    CALL mct_gGrid_importRattr(dom_lnd, 'mask', mask_data, lsize)
    CALL mct_gGrid_importRattr(dom_lnd, 'frac', frac_data, lsize)

    !--- intialize the attribute vectors
    CALL mct_aVect_init(x2l, rList=seq_flds_x2l_fields, lsize=lsize)
    CALL mct_aVect_zero(x2l)
    CALL mct_aVect_init(l2x, rList=seq_flds_l2x_fields, lsize=lsize)
    CALL mct_aVect_zero(l2x)

    !--- fill the l2x attribute vector with export data
    CALL lnd_export_mct(l2x)

    !--- fill some scalar export data
    CALL seq_infodata_PutData(cdata%infodata, &
         lnd_present=.TRUE., lnd_prognostic=.TRUE., &
         sno_present=.FALSE., sno_prognostic=.FALSE.)
    CALL seq_infodata_PutData(infodata, lnd_nx = nx_global, lnd_ny = ny_global)

    !--- deallocate domain arrays
    DEALLOCATE(gindex)
    DEALLOCATE(lon_data)
    DEALLOCATE(lat_data)
    DEALLOCATE(area_data)
    DEALLOCATE(mask_data)
    DEALLOCATE(frac_data)

    !--- set share log unit back to generic one
    CALL shr_file_setLogUnit(shrlogunit)

  END SUBROUTINE lnd_init_mct

  !--------------------------------------------------------------------------
  !> @brief   vic lnd model run
  !--------------------------------------------------------------------------
  SUBROUTINE lnd_run_mct(EClock, cdata, x2l, l2x, cdata_s, x2s, s2x)

    USE shr_file_mod   , ONLY : shr_file_getUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogUnit, &
                                shr_file_setLogUnit
    USE seq_timemgr_mod, ONLY : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn

    IMPLICIT NONE

    !--- INPUT/OUTPUT PARAMETERS:
    TYPE(ESMF_Clock), INTENT(in)     :: EClock
    TYPE(seq_cdata),  INTENT(inout)  :: cdata
    TYPE(mct_aVect),  INTENT(inout)  :: x2l
    TYPE(mct_aVect),  INTENT(inout)  :: l2x
    TYPE(seq_cdata),  INTENT(inout)  :: cdata_s
    TYPE(mct_aVect),  INTENT(inout)  :: x2s
    TYPE(mct_aVect),  INTENT(inout)  :: s2x

    !--- Local Variables
    TYPE(seq_infodata_type), POINTER :: infodata ! CESM information from the driver
    INTEGER(C_INT)                   :: errno
    CHARACTER(len=*), PARAMETER      :: subname = '(lnd_run_mct)'

    !--- set vic log unit
    CALL shr_file_getLogUnit(shrlogunit)
    CALL shr_file_setLogUnit(iulog)

    !--- point to cdata
    CALL seq_cdata_setptrs(cdata, infodata=infodata)

    !--- get time information
    CALL seq_timemgr_EClockGetData(EClock, &
         curr_yr=year, &
         curr_mon=month, &
         curr_day=day, &
         curr_tod=seconds)

    !--- update the clock
    vclock%current_year = year
    vclock%current_month = month
    vclock%current_day = day
    vclock%current_dayseconds = seconds
    vclock%state_flag = seq_timemgr_RestartAlarmIsOn(EClock)
    vclock%stop_flag = seq_timemgr_StopAlarmIsOn(EClock)

    !--- get time of next radiation calc for time dependent albedo calc.
    !--- needed for time (solar angle) dependent albedo calculations.
    !--- this is the julian day of the next radiation calculation
    !--- defined by the atmosphere model.
    !--- fyi, a negative number means don't compute albedos now.
    !--- that will happen when the radiation is computed less frequently
    !--- than the coupling frequency.
    !CALL seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday)

    !--- import data from coupler
    CALL lnd_import_mct(x2l)

    !--- run vic
    errno = vic_cesm_run(vclock)
    IF (errno /= 0) THEN
       CALL shr_sys_abort(subname//' ERROR: vic_cesm_run returned a errno /= 0')
    ENDIF

    !--- export data to coupler
    CALL lnd_export_mct(l2x)

    !--- set share log unit back to generic one
    CALL shr_file_setLogUnit(shrlogunit)

  END SUBROUTINE lnd_run_mct

  !--------------------------------------------------------------------------
  !> @brief   vic lnd model finalize
  !--------------------------------------------------------------------------
  SUBROUTINE lnd_final_mct(EClock, cdata, x2l, l2x, cdata_s, x2s, s2x)

    IMPLICIT NONE

    !--- INPUT/OUTPUT PARAMETERS:
    TYPE(ESMF_Clock)            ,INTENT(inout) :: EClock
    TYPE(seq_cdata)             ,INTENT(inout) :: cdata
    TYPE(mct_aVect)             ,INTENT(inout) :: x2l
    TYPE(mct_aVect)             ,INTENT(inout) :: l2x
    TYPE(seq_cdata)             ,INTENT(inout) :: cdata_s
    TYPE(mct_aVect)             ,INTENT(inout) :: x2s
    TYPE(mct_aVect)             ,INTENT(inout) :: s2x

    !--- Local Variables
    INTEGER(C_INT) :: errno
    CHARACTER(len=*), PARAMETER     :: subname = '(lnd_final_mct)'

    !--- clean up
    errno = vic_cesm_final(vclock)

    IF (errno /= 0) THEN
       CALL shr_sys_abort(subname//' ERROR: vic_cesm_final returned a errno /= 0')
    ENDIF

  END SUBROUTINE lnd_final_mct

  !--------------------------------------------------------------------------
  !> @brief   export fields to coupler
  !--------------------------------------------------------------------------
  SUBROUTINE lnd_export_mct(l2x)

    IMPLICIT NONE

    TYPE(mct_aVect), INTENT(inout) :: l2x
    INTEGER :: i, lsize
    CHARACTER(len=*), PARAMETER :: subname = '(lnd_export_mct)'

    lsize = mct_avect_lsize(l2x)
    l2x%rAttr(:, :) = 0.0_r8

    !--- Copy values to attribute vector
    ! Sign convension and units handeld in VIC driver
    DO i = 1, lsize
       l2x%rAttr(index_l2x_Sl_t, i) = l2x_vic_ptr(i)%l2x_Sl_t
       l2x%rAttr(index_l2x_Sl_tref, i) = l2x_vic_ptr(i)%l2x_Sl_tref
       l2x%rAttr(index_l2x_Sl_qref, i) = l2x_vic_ptr(i)%l2x_Sl_qref
       l2x%rAttr(index_l2x_Sl_avsdr, i) = l2x_vic_ptr(i)%l2x_Sl_avsdr
       l2x%rAttr(index_l2x_Sl_anidr, i) = l2x_vic_ptr(i)%l2x_Sl_anidr
       l2x%rAttr(index_l2x_Sl_avsdf, i) = l2x_vic_ptr(i)%l2x_Sl_avsdf
       l2x%rAttr(index_l2x_Sl_anidf, i) = l2x_vic_ptr(i)%l2x_Sl_anidf
       l2x%rAttr(index_l2x_Sl_snowh, i) = l2x_vic_ptr(i)%l2x_Sl_snowh
       l2x%rAttr(index_l2x_Sl_u10, i) = l2x_vic_ptr(i)%l2x_Sl_u10
       l2x%rAttr(index_l2x_Sl_logz0, i) = l2x_vic_ptr(i)%l2x_Sl_logz0
       l2x%rAttr(index_l2x_Fall_taux, i) = l2x_vic_ptr(i)%l2x_Fall_taux
       l2x%rAttr(index_l2x_Fall_tauy, i) = l2x_vic_ptr(i)%l2x_Fall_tauy
       l2x%rAttr(index_l2x_Fall_lat, i) = l2x_vic_ptr(i)%l2x_Fall_lat
       l2x%rAttr(index_l2x_Fall_sen, i) = l2x_vic_ptr(i)%l2x_Fall_sen
       l2x%rAttr(index_l2x_Fall_lwup, i) = l2x_vic_ptr(i)%l2x_Fall_lwup
       l2x%rAttr(index_l2x_Fall_evap, i) = l2x_vic_ptr(i)%l2x_Fall_evap
       l2x%rAttr(index_l2x_Fall_swnet, i) = l2x_vic_ptr(i)%l2x_Fall_swnet
       l2x%rAttr(index_l2x_Flrl_rofliq, i) = l2x_vic_ptr(i)%l2x_Flrl_rofliq
       l2x%rAttr(index_l2x_Flrl_rofice, i) = l2x_vic_ptr(i)%l2x_Flrl_rofice

       !--- optional fields ---
       IF (index_l2x_Sl_ddvel /= 0) THEN
          l2x%rAttr(index_l2x_Sl_ddvel, i) = l2x_vic_ptr(i)%l2x_Sl_ddvel
       ENDIF
       IF (index_l2x_Fall_fco2_lnd /= 0) THEN
          l2x%rAttr(index_l2x_Fall_fco2_lnd, i) = l2x_vic_ptr(i)%l2x_Fall_fco2_lnd
       ENDIF
       IF (index_l2x_Sl_fv /= 0) THEN
          l2x%rAttr(index_l2x_Sl_fv, i) = l2x_vic_ptr(i)%l2x_Sl_fv
       ENDIF
       IF (index_l2x_Sl_ram1 /= 0) THEN
          l2x%rAttr(index_l2x_Sl_ram1, i) = l2x_vic_ptr(i)%l2x_Sl_ram1
       ENDIF
       IF (index_l2x_Fall_flxdst1 /= 0) THEN
          l2x%rAttr(index_l2x_Fall_flxdst1, i) = l2x_vic_ptr(i)%l2x_Fall_flxdst1
       ENDIF
       IF (index_l2x_Fall_flxdst2 /= 0) THEN
          l2x%rAttr(index_l2x_Fall_flxdst2, i) = l2x_vic_ptr(i)%l2x_Fall_flxdst2
       ENDIF
       IF (index_l2x_Fall_flxdst3 /= 0) THEN
          l2x%rAttr(index_l2x_Fall_flxdst3, i) = l2x_vic_ptr(i)%l2x_Fall_flxdst3
       ENDIF
       IF (index_l2x_Fall_flxdst4 /= 0) THEN
          l2x%rAttr(index_l2x_Fall_flxdst4, i) = l2x_vic_ptr(i)%l2x_Fall_flxdst4
       ENDIF
       IF (index_l2x_Fall_flxvoc /= 0) THEN
          l2x%rAttr(index_l2x_Fall_flxvoc, i) = l2x_vic_ptr(i)%l2x_Fall_flxvoc
       ENDIF

       IF (.NOT. l2x_vic_ptr(i)%l2x_vars_set) THEN
          CALL shr_sys_abort(subname//' ERROR: l2x export vars not set')
       ENDIF
    END DO

  END SUBROUTINE lnd_export_mct

  !--------------------------------------------------------------------------
  !> @brief   import fields from coupler
  !--------------------------------------------------------------------------
  SUBROUTINE lnd_import_mct(x2l)

    IMPLICIT NONE

    !--- Arguments
    TYPE(mct_aVect), INTENT(inout) :: x2l

    !--- Local Variables
    INTEGER  :: i, lsize
    CHARACTER(len=*), PARAMETER :: subname = '(lnd_import_mct)'

    lsize = mct_avect_lsize(x2l)

    DO i = 1, lsize
       x2l_vic_ptr(i)%x2l_Sa_z = x2l%rAttr(index_x2l_Sa_z, i)
       x2l_vic_ptr(i)%x2l_Sa_u = x2l%rAttr(index_x2l_Sa_u, i)
       x2l_vic_ptr(i)%x2l_Sa_v = x2l%rAttr(index_x2l_Sa_v, i)
       x2l_vic_ptr(i)%x2l_Sa_ptem = x2l%rAttr(index_x2l_Sa_ptem, i)
       x2l_vic_ptr(i)%x2l_Sa_shum = x2l%rAttr(index_x2l_Sa_shum, i)
       x2l_vic_ptr(i)%x2l_Sa_pbot = x2l%rAttr(index_x2l_Sa_pbot, i)
       x2l_vic_ptr(i)%x2l_Sa_tbot = x2l%rAttr(index_x2l_Sa_tbot, i)
       x2l_vic_ptr(i)%x2l_Faxa_lwdn = x2l%rAttr(index_x2l_Faxa_lwdn, i)
       x2l_vic_ptr(i)%x2l_Faxa_rainc = x2l%rAttr(index_x2l_Faxa_rainc, i)
       x2l_vic_ptr(i)%x2l_Faxa_rainl = x2l%rAttr(index_x2l_Faxa_rainl, i)
       x2l_vic_ptr(i)%x2l_Faxa_snowc = x2l%rAttr(index_x2l_Faxa_snowc, i)
       x2l_vic_ptr(i)%x2l_Faxa_snowl = x2l%rAttr(index_x2l_Faxa_snowl, i)
       x2l_vic_ptr(i)%x2l_Faxa_swndr = x2l%rAttr(index_x2l_Faxa_swndr, i)
       x2l_vic_ptr(i)%x2l_Faxa_swvdr = x2l%rAttr(index_x2l_Faxa_swvdr, i)
       x2l_vic_ptr(i)%x2l_Faxa_swndf = x2l%rAttr(index_x2l_Faxa_swndf, i)
       x2l_vic_ptr(i)%x2l_Faxa_swvdf = x2l%rAttr(index_x2l_Faxa_swvdf, i)
       x2l_vic_ptr(i)%x2l_Faxa_bcphidry = x2l%rAttr(index_x2l_Faxa_bcphidry, i)
       x2l_vic_ptr(i)%x2l_Faxa_bcphodry = x2l%rAttr(index_x2l_Faxa_bcphodry, i)
       x2l_vic_ptr(i)%x2l_Faxa_bcphiwet = x2l%rAttr(index_x2l_Faxa_bcphiwet, i)
       x2l_vic_ptr(i)%x2l_Faxa_ocphidry = x2l%rAttr(index_x2l_Faxa_ocphidry, i)
       x2l_vic_ptr(i)%x2l_Faxa_ocphodry = x2l%rAttr(index_x2l_Faxa_ocphodry, i)
       x2l_vic_ptr(i)%x2l_Faxa_ocphiwet = x2l%rAttr(index_x2l_Faxa_ocphiwet, i)
       x2l_vic_ptr(i)%x2l_Faxa_dstwet1 = x2l%rAttr(index_x2l_Faxa_dstwet1, i)
       x2l_vic_ptr(i)%x2l_Faxa_dstwet2 = x2l%rAttr(index_x2l_Faxa_dstwet2, i)
       x2l_vic_ptr(i)%x2l_Faxa_dstwet3 = x2l%rAttr(index_x2l_Faxa_dstwet3, i)
       x2l_vic_ptr(i)%x2l_Faxa_dstwet4 = x2l%rAttr(index_x2l_Faxa_dstwet4, i)
       x2l_vic_ptr(i)%x2l_Faxa_dstdry1 = x2l%rAttr(index_x2l_Faxa_dstdry1, i)
       x2l_vic_ptr(i)%x2l_Faxa_dstdry2 = x2l%rAttr(index_x2l_Faxa_dstdry2, i)
       x2l_vic_ptr(i)%x2l_Faxa_dstdry3 = x2l%rAttr(index_x2l_Faxa_dstdry3, i)
       x2l_vic_ptr(i)%x2l_Faxa_dstdry4 = x2l%rAttr(index_x2l_Faxa_dstdry4, i)
       x2l_vic_ptr(i)%x2l_Flrr_flood = x2l%rAttr(index_x2l_Flrr_flood, i)

       !--- Determine optional receive fields
       IF (index_x2l_Sa_co2prog /= 0) THEN
          x2l_vic_ptr(i)%x2l_Sa_co2prog = x2l%rAttr(index_x2l_Sa_co2prog, i)
       ENDIF
       IF (index_x2l_Sa_co2diag /= 0) THEN
          x2l_vic_ptr(i)%x2l_Sa_co2diag = x2l%rAttr(index_x2l_Sa_co2diag, i)
       ENDIF

       x2l_vic_ptr(i)%x2l_vars_set =  .TRUE.

    END DO

  END SUBROUTINE lnd_import_mct

  !--------------------------------------------------------------------------
  !> @brief   Set the coupler indices needed by the land model coupler
  !--------------------------------------------------------------------------
  SUBROUTINE cpl_indices_set()

    !--- USES:
    USE seq_flds_mod  , ONLY: seq_flds_x2l_fields, seq_flds_l2x_fields,     &
                              seq_flds_x2s_fields, seq_flds_s2x_fields
    USE mct_mod       , ONLY: mct_aVect, mct_aVect_init, mct_avect_indexra, &
                              mct_aVect_clean, mct_avect_nRattr
    USE seq_drydep_mod, ONLY: drydep_fields_token, lnd_drydep
    USE shr_megan_mod , ONLY: shr_megan_fields_token, shr_megan_mechcomps_n

    !--- ARGUMENTS:
    IMPLICIT NONE

    ! LOCAL VARIABLES:
    TYPE(mct_aVect)   :: l2x      ! temporary, land to coupler
    TYPE(mct_aVect)   :: x2l      ! temporary, coupler to land
    TYPE(mct_aVect)   :: s2x      ! temporary, glacier to coupler
    TYPE(mct_aVect)   :: x2s      ! temporary, coupler to glacier
    INTEGER           :: num
    CHARACTER(len= 2) :: cnum
    CHARACTER(len=64) :: name
    CHARACTER(len=32) :: subname = '(cpl_indices_set)'  ! subroutine name

    !--- Determine attribute vector indices

    !--- create temporary attribute vectors
    CALL mct_aVect_init(x2l, rList=seq_flds_x2l_fields, lsize=1)
    CALL mct_aVect_init(l2x, rList=seq_flds_l2x_fields, lsize=1)

    !-------------------------------------------------------------
    ! vic -> drv
    !-------------------------------------------------------------
    index_l2x_Flrl_rofliq   = mct_avect_indexra(l2x, 'Flrl_rofliq')
    index_l2x_Flrl_rofice   = mct_avect_indexra(l2x, 'Flrl_rofice')

    index_l2x_Sl_t          = mct_avect_indexra(l2x, 'Sl_t')
    index_l2x_Sl_snowh      = mct_avect_indexra(l2x, 'Sl_snowh')
    index_l2x_Sl_avsdr      = mct_avect_indexra(l2x, 'Sl_avsdr')
    index_l2x_Sl_anidr      = mct_avect_indexra(l2x, 'Sl_anidr')
    index_l2x_Sl_avsdf      = mct_avect_indexra(l2x, 'Sl_avsdf')
    index_l2x_Sl_anidf      = mct_avect_indexra(l2x, 'Sl_anidf')
    index_l2x_Sl_tref       = mct_avect_indexra(l2x, 'Sl_tref')
    index_l2x_Sl_qref       = mct_avect_indexra(l2x, 'Sl_qref')
    index_l2x_Sl_u10        = mct_avect_indexra(l2x, 'Sl_u10')
    index_l2x_Sl_ram1       = mct_avect_indexra(l2x, 'Sl_ram1')
    index_l2x_Sl_fv         = mct_avect_indexra(l2x, 'Sl_fv')
    index_l2x_Sl_logz0      = mct_avect_indexra(l2x, 'Sl_logz0')
    IF (lnd_drydep) THEN
       index_l2x_Sl_ddvel = mct_avect_indexra(l2x, TRIM(drydep_fields_token))
    ELSE
       index_l2x_Sl_ddvel = 0
    ENDIF

    index_l2x_Fall_taux     = mct_avect_indexra(l2x, 'Fall_taux')
    index_l2x_Fall_tauy     = mct_avect_indexra(l2x, 'Fall_tauy')
    index_l2x_Fall_lat      = mct_avect_indexra(l2x, 'Fall_lat')
    index_l2x_Fall_sen      = mct_avect_indexra(l2x, 'Fall_sen')
    index_l2x_Fall_lwup     = mct_avect_indexra(l2x, 'Fall_lwup')
    index_l2x_Fall_evap     = mct_avect_indexra(l2x, 'Fall_evap')
    index_l2x_Fall_swnet    = mct_avect_indexra(l2x, 'Fall_swnet')
    index_l2x_Fall_flxdst1  = mct_avect_indexra(l2x, 'Fall_flxdst1')
    index_l2x_Fall_flxdst2  = mct_avect_indexra(l2x, 'Fall_flxdst2')
    index_l2x_Fall_flxdst3  = mct_avect_indexra(l2x, 'Fall_flxdst3')
    index_l2x_Fall_flxdst4  = mct_avect_indexra(l2x, 'Fall_flxdst4')
    ! co2 var, perrwith='quiet' supresses errors
    index_l2x_Fall_fco2_lnd = mct_avect_indexra(l2x, 'Fall_fco2_lnd', perrwith='quiet')

    !--- MEGAN fluxes
    ! MEGAN = Model of Emissions of Gases and Aerosols from Nature
    ! Emissions from land surface
    IF (shr_megan_mechcomps_n>0) THEN
       index_l2x_Fall_flxvoc = mct_avect_indexra(l2x, TRIM(shr_megan_fields_token))
    ELSE
       index_l2x_Fall_flxvoc = 0
    ENDIF

    nflds_l2x = mct_avect_nRattr(l2x)

    !-------------------------------------------------------------
    ! drv -> vic
    !-------------------------------------------------------------
    index_x2l_Sa_z          = mct_avect_indexra(x2l, 'Sa_z')
    index_x2l_Sa_u          = mct_avect_indexra(x2l, 'Sa_u')
    index_x2l_Sa_v          = mct_avect_indexra(x2l, 'Sa_v')
    index_x2l_Sa_ptem       = mct_avect_indexra(x2l, 'Sa_ptem')
    index_x2l_Sa_pbot       = mct_avect_indexra(x2l, 'Sa_pbot')
    index_x2l_Sa_tbot       = mct_avect_indexra(x2l, 'Sa_tbot')
    index_x2l_Sa_shum       = mct_avect_indexra(x2l, 'Sa_shum')
    ! co2 vars, perrwith='quiet' supresses errors
    index_x2l_Sa_co2prog    = mct_avect_indexra(x2l, 'Sa_co2prog', perrwith='quiet')
    index_x2l_Sa_co2diag    = mct_avect_indexra(x2l, 'Sa_co2diag', perrwith='quiet')

    index_x2l_Faxa_lwdn     = mct_avect_indexra(x2l, 'Faxa_lwdn')
    index_x2l_Faxa_rainc    = mct_avect_indexra(x2l, 'Faxa_rainc')
    index_x2l_Faxa_rainl    = mct_avect_indexra(x2l, 'Faxa_rainl')
    index_x2l_Faxa_snowc    = mct_avect_indexra(x2l, 'Faxa_snowc')
    index_x2l_Faxa_snowl    = mct_avect_indexra(x2l, 'Faxa_snowl')
    index_x2l_Faxa_swndr    = mct_avect_indexra(x2l, 'Faxa_swndr')
    index_x2l_Faxa_swvdr    = mct_avect_indexra(x2l, 'Faxa_swvdr')
    index_x2l_Faxa_swndf    = mct_avect_indexra(x2l, 'Faxa_swndf')
    index_x2l_Faxa_swvdf    = mct_avect_indexra(x2l, 'Faxa_swvdf')
    index_x2l_Faxa_bcphidry = mct_avect_indexra(x2l, 'Faxa_bcphidry')
    index_x2l_Faxa_bcphodry = mct_avect_indexra(x2l, 'Faxa_bcphodry')
    index_x2l_Faxa_bcphiwet = mct_avect_indexra(x2l, 'Faxa_bcphiwet')
    index_x2l_Faxa_ocphidry = mct_avect_indexra(x2l, 'Faxa_ocphidry')
    index_x2l_Faxa_ocphodry = mct_avect_indexra(x2l, 'Faxa_ocphodry')
    index_x2l_Faxa_ocphiwet = mct_avect_indexra(x2l, 'Faxa_ocphiwet')
    index_x2l_Faxa_dstdry1  = mct_avect_indexra(x2l, 'Faxa_dstdry1')
    index_x2l_Faxa_dstdry2  = mct_avect_indexra(x2l, 'Faxa_dstdry2')
    index_x2l_Faxa_dstdry3  = mct_avect_indexra(x2l, 'Faxa_dstdry3')
    index_x2l_Faxa_dstdry4  = mct_avect_indexra(x2l, 'Faxa_dstdry4')
    index_x2l_Faxa_dstwet1  = mct_avect_indexra(x2l, 'Faxa_dstwet1')
    index_x2l_Faxa_dstwet2  = mct_avect_indexra(x2l, 'Faxa_dstwet2')
    index_x2l_Faxa_dstwet3  = mct_avect_indexra(x2l, 'Faxa_dstwet3')
    index_x2l_Faxa_dstwet4  = mct_avect_indexra(x2l, 'Faxa_dstwet4')

    index_x2l_Flrr_flood    = mct_avect_indexra(x2l, 'Flrr_flood')

    nflds_x2l = mct_avect_nRattr(x2l)

    CALL mct_aVect_clean(x2l)
    CALL mct_aVect_clean(l2x)

  END SUBROUTINE cpl_indices_set


  !--------------------------------------------------------------------------
  !> @brief   Unpack the VIC domain structure
  !--------------------------------------------------------------------------
  SUBROUTINE unpack_vic_domain(gindex, lon_data, lat_data, area_data, &
                               mask_data, frac_data)

    IMPLICIT NONE

    !--- INPUT/OUTPUT PARAMETERS:
    INTEGER,  DIMENSION(:), INTENT(inout) :: gindex  ! Number the local grid points
    REAL(r8), DIMENSION(:), INTENT(inout) :: lon_data
    REAL(r8), DIMENSION(:), INTENT(inout) :: lat_data
    REAL(r8), DIMENSION(:), INTENT(inout) :: area_data
    REAL(R8), DIMENSION(:), INTENT(inout) :: mask_data
    REAL(r8), DIMENSION(:), INTENT(inout) :: frac_data

    !--- LOCAL VARIABLES:
    INTEGER :: i
    TYPE(location_struct), DIMENSION(:), POINTER :: locations
    CHARACTER(len=*), PARAMETER :: subname = '(unpack_vic_domain)'

    !--- Associate locations pointer in local_domain structure
    CALL c_f_pointer(local_domain%locations, locations, [local_domain%ncells_active])

    !--- Get domain information for each cell in local domain
    DO i = 1, local_domain%ncells_active
       gindex(i) = locations(i)%io_idx + 1
       lon_data(i) = locations(i)%longitude
       lat_data(i) = locations(i)%latitude
       area_data(i) = locations(i)%area
       frac_data(i) = locations(i)%frac
       mask_data(i) = 1  ! 1 for all active VIC grid cells
    END DO

  END SUBROUTINE unpack_vic_domain

  !--------------------------------------------------------------------------
  !> @brief   Copy fortran string to character array
  !--------------------------------------------------------------------------
  PURE FUNCTION Copy_s2a(s)  RESULT (a)   ! copy s(1:Clen(s)) to char array
    CHARACTER(*),INTENT(IN) :: s
    CHARACTER(kind=C_CHAR)  :: a(LEN(s))
    INTEGER                 :: i

    DO i = 1,LEN(s) - 1 ! last char is null
       a(i) = s(i:i)
    END DO
    a(LEN(s)) = C_NULL_CHAR
  END FUNCTION Copy_s2a

END MODULE lnd_comp_mct
