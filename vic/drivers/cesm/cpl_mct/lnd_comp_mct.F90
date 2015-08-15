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

  subroutine lnd_init_mct( EClock, cdata, x2d, d2x, cdata_s, x2s, s2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    type(seq_cdata)             , intent(inout) :: cdata_s
    type(mct_aVect)             , intent(inout) :: x2s, s2x
    character(len=*), optional  , intent(in)    :: NLFilename

    ! Local Variables
    integer(C_INT) :: errno
    type(vic_clock_type) :: vic_clock
    character(len=*), parameter     :: subname = "lnd_init_mct"

!EOP
!-------------------------------------------------------------------------------

   call seq_infodata_PutData(cdata%infodata, &
        lnd_present=.true., lnd_prognostic=.true., &
        sno_present=.false., sno_prognostic=.false.)

  errno = vic_cesm_init(vic_clock)

  if (errno /= 0) then
    call shr_sys_abort(subname//':: vic_cesm_init returned a errno /= 0' )
  endif


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

subroutine lnd_run_mct( EClock, cdata, x2d, d2x, cdata_s, x2s, s2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d
   type(mct_aVect)             ,intent(inout) :: d2x
   type(seq_cdata)             ,intent(inout) :: cdata_s
   type(mct_aVect)             ,intent(inout) :: x2s
   type(mct_aVect)             ,intent(inout) :: s2x

   ! Local Variables
   integer(C_INT) :: errno
   type(vic_clock_type) :: vic_clock
   character(len=*), parameter     :: subname = "lnd_run_mct"

!EOP
!-------------------------------------------------------------------------------

  errno = vic_cesm_run(vic_clock)

  if (errno /= 0) then
    call shr_sys_abort(subname//':: vic_cesm_run returned a errno /= 0' )
  endif

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
subroutine lnd_final_mct( EClock, cdata, x2d, d2x, cdata_s, x2s, s2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d
    type(mct_aVect)             ,intent(inout) :: d2x
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

!===============================================================================

end module lnd_comp_mct
