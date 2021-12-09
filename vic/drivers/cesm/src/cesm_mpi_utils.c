/******************************************************************************
 * @section DESCRIPTION
 *
 * MPI support routines for VIC's CESM driver
 *****************************************************************************/

#include <vic_driver_cesm.h>

/******************************************************************************
 * @brief   Wrapper around VIC's initialize mpi function, first translating the
            Fortran MPI Communicator to C.
 *****************************************************************************/
void
initialize_vic_cesm_mpi(MPI_Fint *MPI_COMM_VIC_F)
{
    extern MPI_Comm MPI_COMM_VIC;

    MPI_COMM_VIC = MPI_Comm_f2c(*MPI_COMM_VIC_F);

    initialize_mpi();
}
