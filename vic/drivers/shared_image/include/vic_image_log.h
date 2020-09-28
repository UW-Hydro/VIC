/******************************************************************************
 *
 * @section DESCRIPTION
 *
 * Logging macros for netcdf routines
 *
 * @section LICENSE
 *
 * Copyright (c) 2010, Zed A. Shaw and Mongrel2 Project Contributors.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 *     * Neither the name of the Mongrel2 Project, Zed A. Shaw, nor the names
 *       of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *****************************************************************************/

#ifndef __vic_image_log_h__
#define __vic_image_log_h__

#include <vic_def.h>
#include <vic_mpi.h>

// Macros for logging
#define clean_ncerrno(e) (e == NC_NOERR ? "None" : nc_strerror(e))

// Error Level is always active
#ifdef NO_LINENOS
#define log_ncerr(e, M, ...) print_trace(); fprintf(LOG_DEST, \
                                                    "[ERROR] errno: %s: " M "\n", \
                                                    clean_ncerrno(e), \
                                                    ## __VA_ARGS__); \
    exit(EXIT_FAILURE);
#else
#define log_ncerr(e, M, ...) print_trace(); fprintf(LOG_DEST, \
                                                    "[ERROR] %s:%d: errno: %s: " M " \n", \
                                                    __FILE__, __LINE__, \
                                                    clean_ncerrno(e), \
                                                    ## __VA_ARGS__); \
    exit(EXIT_FAILURE);
#endif

#define check_nc_status(A, M, ...) if (A != NC_NOERR) {log_ncerr(A, M, \
                                                                 ## __VA_ARGS__); \
                                                       errno = 0; exit( \
                                                           EXIT_FAILURE);}

#define log_mpi_err(e, M, ...) print_trace(); \
    print_mpi_error_str(e); fprintf(LOG_DEST, \
                                    "[ERROR] %s:%d: errno: %d: " M " \n", \
                                    __FILE__, __LINE__, e, \
                                    ## __VA_ARGS__); \
    MPI_Abort(MPI_COMM_VIC, e);

#define check_mpi_status(A, M, ...) if (A != MPI_SUCCESS) {log_mpi_err(A, M, \
                                                                       ## __VA_ARGS__); \
                                                           errno = 0; MPI_Abort( \
                                                               MPI_COMM_VIC, A); \
}

#endif
