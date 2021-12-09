!>
!! @section DESCRIPTION
!!
!! Fortran interface for CESM driver. Uses the Fortran/C ISO_C_BINDING.

MODULE vic_cesm_interface

  ! !PUBLIC TYPES:
  IMPLICIT NONE
  SAVE
  PRIVATE ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------
  PUBLIC :: initialize_log
  PUBLIC :: initialize_vic_cesm_mpi
  PUBLIC :: vic_cesm_init
  PUBLIC :: vic_cesm_run
  PUBLIC :: vic_cesm_final

  !--------------------------------------------------------------------------
  !> @brief   Initialize logging to stderr until after mpi is initialized
  !--------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE initialize_log() BIND(C, name='initialize_log')
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT NONE
     END SUBROUTINE initialize_log
  END INTERFACE

  !--------------------------------------------------------------------------
  !> @brief   Init MPI Interface
  !--------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE initialize_vic_cesm_mpi(MPI_COMM_VIC_F) BIND(C, name='initialize_vic_cesm_mpi')
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT NONE
       INTEGER, INTENT(inout) :: MPI_COMM_VIC_F
     END SUBROUTINE initialize_vic_cesm_mpi
  END INTERFACE

  !--------------------------------------------------------------------------
  !> @brief   Init Interface
  !--------------------------------------------------------------------------
  INTERFACE
     INTEGER(C_INT) FUNCTION vic_cesm_init(vclock, cmeta) BIND(C, name='vic_cesm_init')
       USE, INTRINSIC :: ISO_C_BINDING
       USE vic_cesm_def_mod
       IMPLICIT NONE
       TYPE(vic_clock), INTENT(in) :: vclock
       TYPE(case_metadata), INTENT(in) :: cmeta
     END FUNCTION vic_cesm_init
  END INTERFACE

  !--------------------------------------------------------------------------
  !> @brief   Run Interface
  !--------------------------------------------------------------------------
  INTERFACE
     INTEGER(C_INT) FUNCTION vic_cesm_run(vclock) BIND(C, name='vic_cesm_run')
       USE, INTRINSIC :: ISO_C_BINDING
       USE vic_cesm_def_mod
       IMPLICIT NONE
       TYPE(vic_clock), INTENT(in) :: vclock
     END FUNCTION vic_cesm_run
  END INTERFACE

  !--------------------------------------------------------------------------
  !> @brief   Finalize Interface
  !--------------------------------------------------------------------------
  INTERFACE
     INTEGER(C_INT) FUNCTION vic_cesm_final(vclock) BIND(C, name='vic_cesm_final')
       USE, INTRINSIC :: ISO_C_BINDING
       USE vic_cesm_def_mod
       IMPLICIT NONE
       TYPE(vic_clock), INTENT(in) :: vclock
     END FUNCTION vic_cesm_final
  END INTERFACE

END MODULE vic_cesm_interface
