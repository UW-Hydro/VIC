!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_initialize2
!
! !INTERFACE:
  subroutine clm_initialize2(dt, nlevgrnd, begg, endg, begc, endc, begp, &
	endp, nveg, vegtype, vegfrac, npfts)
!
! !DESCRIPTION:
! Land model initialization.
! o Initializes run control variables via the [clm_inparm] namelist.
! o Reads surface data on model grid.
! o Defines the multiple plant types and fraction areas for each surface type.
! o Builds the appropriate subgrid <-> grid mapping indices and weights.
! o Set up parallel processing.
! o Initializes time constant variables.
! o Reads restart data for a restart or branch run.
! o Reads initial data and initializes the time variant variables for an initial run.
! o Initializes history file output.
! o Initializes river routing model.
! o Initializes accumulation variables.
!
! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
!    use clm_varctl      , only : finidat, fpftdyn
    use clm_varpar      , only : mxveg, mxpft
    use clmtypeInitMod  , only : initClmtype
!    use mkarbinitMod    , only : mkarbinit
    use iniTimeConstMod , only : iniTimeConst
    use pftvarcon,         only : pftconrd
    use ndepStreamMod    , only : ndep_init
    use CNiniTimeVarMod  , only : CNiniTimeVar
    use CNEcosystemDynMod, only : CNEcosystemDynInit
    use CNBalanceCheckMod, only : BeginCBalance, BeginNBalance
    use CNSetValueMod,     only : CNZeroFluxes_dwt
!#if (defined CNDV)
!    use pftdynMod             , only : pftwt_init
!    use CNDVEcosystemDyniniMod, only : CNDVEcosystemDynini
!#endif
!    use STATICEcosysDynMod , only : EcosystemDynini, readAnnualVegetation
!    use STATICEcosysDynMod , only : interpMonthlyVeg
!    use clm_time_manager, only : get_curr_date, get_nstep, advance_timestep, &
!                                 timemgr_init, timemgr_restart_io, timemgr_restart
!    use clm_time_manager, only : get_step_size, get_curr_calday
!    use shr_orb_mod        , only : shr_orb_decl
!    use clm_varorb         , only : eccen, mvelpp, lambm0, obliqr

! !Arguments    
    implicit none

    real(r8) :: dt                    ! timestep
    integer  :: nlevgrnd              ! # soil layers
    integer  :: begp, endp            ! clump beg and ending pft indices
    integer  :: begc, endc            ! clump beg and ending column indices
    integer  :: begg, endg            ! clump beg and ending gridcell indices
    integer  :: vegtype(1:mxveg+1)       ! VIC vegetation types
    real(r8) :: vegfrac(1:mxveg+1)       ! VIC vegetation fractions
!    real(r8) :: pftfrac(1:mxveg+1,0:mxpft)
    integer  :: nveg                  ! # VIC vegetation types
    integer  :: npfts                 ! # CLM PFTs

!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
! 12/10/12: Adapted for use in VIC by Michael Brunke
! 05/29/14: Added pftfrac for ingesting PFT fractions, Michael Brunke
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: nl,nlg,na,nag         ! indices
    integer  :: i,j,k,n1,n2           ! indices
    integer  :: yr                    ! current year (0, ...)
    integer  :: mon                   ! current month (1 -> 12)
    integer  :: day                   ! current day (1 -> 31)
    integer  :: ncsec                 ! current time of day [seconds]
    integer  :: nc                    ! clump index
    integer  :: nclumps               ! number of clumps on this processor
    integer  :: begg_atm, endg_atm    ! proc beg and ending gridcell indices
    character(len=256) :: fnamer      ! name of netcdf restart file 
    character(len=256) :: pnamer      ! full pathname of netcdf restart file
    real(r8) :: dtime                 ! time step increment (sec)
    integer  :: nstep                 ! model time step
    real(r8) :: calday                ! calendar day for nstep
    real(r8) :: caldaym1              ! calendar day for nstep-1
    real(r8) :: declin                ! solar declination angle in radians for nstep
    real(r8) :: declinm1              ! solar declination angle in radians for nstep-1
    real(r8) :: eccf                  ! earth orbit eccentricity factor
    character(len=32) :: subname = 'initialize2' ! subroutine name
!----------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Allocate memory and initialize values of clmtype data structures
    ! ------------------------------------------------------------------------

    call initClmtype(begg, endg, begc, endc, begp, endp, nveg, vegtype, &
	vegfrac, npfts, nlevgrnd)

    call pftconrd

    ! ------------------------------------------------------------------------
    ! Initialize time constant variables 
    ! ------------------------------------------------------------------------

    ! Initialize Ecosystem Dynamics 

!#if (defined CNDV)
!    call CNDVEcosystemDynini()
!#elif (!defined CN)
!    call EcosystemDynini()
!#endif
!#if (defined CN) || (defined CNDV)

    ! --------------------------------------------------------------
    ! Initialize CLMSP ecosystem dynamics when drydeposition is used
    ! so that estimates of monthly differences in LAI can be computed
    ! --------------------------------------------------------------

!    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
!       call EcosystemDynini()
!    end if
!#endif

    call iniTimeConst(nlevgrnd, begg, endg, begc, endc, begp, endp)

    ! ------------------------------------------------------------------------
    ! Initialize CN Ecosystem Dynamics (must be after time-manager initialization)
    ! ------------------------------------------------------------------------
    call CNEcosystemDynInit( begc, endc, begp, endp, dt )

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! Initialize accumulator fields to be time accumulated for various purposes.
    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed.

!    call t_startf('init_accflds')
!    call initAccFlds()
!    call t_stopf('init_accflds')

    ! ------------------------------------------------------------------------
    ! Set arbitrary initial conditions for time varying fields 
    ! used in coupled carbon-nitrogen code
    ! ------------------------------------------------------------------------

!    if (nsrest == nsrStartup) then
       call CNiniTimeVar(begg, endg, begc, endc, begp, endp)
!    end if

    ! ------------------------------------------------------------------------
    ! Initialization of dynamic pft weights
    ! ------------------------------------------------------------------------

    ! Determine correct pft weights (interpolate pftdyn dataset if initial run)
    ! Otherwise these are read in for a restart run

!#if (defined CNDV)
!    call pftwt_init()
!#else
!    if (fpftdyn /= ' ') then
!       call pftdyn_init()
!       call pftdyn_interp( )
!    end if
!#endif


    ! ------------------------------------------------------------------------
    ! Initialize nitrogen deposition
    ! ------------------------------------------------------------------------

!#ifdef CN
    call ndep_init(begg, endg)
!?    call ndep_interp()
!#endif
    
    ! ------------------------------------------------------------------------
    ! Initialize accumator buffers
    ! ------------------------------------------------------------------------

    ! Initialize clmtype variables that are obtained from accumulated fields.
    ! This routine is called in an initial run at nstep=0
    ! This routine is also always called for a restart run and must 
    ! therefore be called after the restart file is read in

!    call initAccClmtype()

    ! Initialize C and N balance, added MAB 10/25/13

     call BeginCBalance(endc)
     call BeginNBalance(endc)

    ! Set conversion and product pool fluxes to 0, added MAB 10/25/13

     call CNZeroFluxes_dwt(begc, endc, begp, endp)

    ! --------------------------------------------------------------
    ! Note - everything below this point needs updated weights
    ! --------------------------------------------------------------

    !
    ! Even if CN is on, and dry-deposition is active, read CLMSP annual vegetation to get estimates of monthly LAI
    !
!    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
!       call readAnnualVegetation()
!    end if

    ! End initialization

!    call t_startf('init_wlog')
!    if (masterproc) then
!       write(iulog,*) 'Successfully initialized the land model'
!       if (nsrest == nsrStartup) then
!          write(iulog,*) 'begin initial run at: '
!       else
!          write(iulog,*) 'begin continuation run at:'
!       end if
!       call get_curr_date(yr, mon, day, ncsec)
!       write(iulog,*) '   nstep= ',get_nstep(), ' year= ',yr,' month= ',mon,&
!            ' day= ',day,' seconds= ',ncsec
!       write(iulog,*)
!       write(iulog,'(72a1)') ("*",i=1,60)
!       write(iulog,*)
!    endif
!    call t_stopf('init_wlog')

!    if (get_nstep() == 0 .or. nsrest == nsrStartup) then
       ! Initialize albedos (correct pft filters are needed)

!       if (finidat == ' ' .or. do_initsurfalb) then
!          call t_startf('init_orb')
!          calday = get_curr_calday()
!          call t_startf('init_orbd1')
!          call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf )
!          call t_stopf('init_orbd1')
          
!          dtime = get_step_size()
!          caldaym1 = get_curr_calday(offset=-int(dtime))
!          call t_startf('init_orbd2')
!          call shr_orb_decl( caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf )
!          call t_stopf('init_orbd2')
          
!          call t_startf('init_orbSA')
!          call initSurfAlb( calday, declin, declinm1 )
!          call t_stopf('init_orbSA')
!          call t_stopf('init_orb')
!       else if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
          ! Call interpMonthlyVeg for dry-deposition so that mlaidiff will be calculated
          ! This needs to be done even if CN or CNDV is on!
!          call interpMonthlyVeg()
!       end if

       ! Determine gridcell averaged properties to send to atm

!    end if

  end subroutine clm_initialize2
