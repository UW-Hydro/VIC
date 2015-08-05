/******************************************************************************
 *
 * @section DESCRIPTION
 *
 * Logging macros
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

#ifndef __vic_log_h__
#define __vic_log_h__

#include <stdio.h>
#include <errno.h>
#include <string.h>

// Set the log level
// To turn off warning statments, set LOG_LVL >= 30
// | Level     | Numeric value    |
// |---------  |---------------   |
// | ERROR     | Always Active    |
// | WARNING   | < 30             |
// | INFO      | < 20             |
// | DEBUG     | < 10             |
#ifndef LOG_LVL
#define LOG_LVL 25
#endif

FILE *LOG_DEST;

void finalize_logging(void);
void get_current_datetime(char *cdt);
void get_logname(const char *path, int id, char *filename);
void initialize_log(void);
void setup_logging(int id);

// Macros for logging
#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define clean_ncerrno(e) (nc_strerror(e))

// Debug Level
#if LOG_LVL < 10
#define debug(M, ...) fprintf(LOG_DEST, "[DEBUG] %s:%d: " M "\n", __FILE__, \
                              __LINE__, ## __VA_ARGS__)
#else
#define debug(M, ...)

#endif

// Info Level
#if LOG_LVL < 20
#ifdef NO_LINENOS
#define log_info(M, ...) fprintf(LOG_DEST, "[INFO] " M "\n", ## __VA_ARGS__)
#else
#define log_info(M, ...) fprintf(LOG_DEST, "[INFO] %s:%d: " M "\n", __FILE__, \
                                 __LINE__, ## __VA_ARGS__)
#endif
#else
#define log_info(M, ...)
#endif

// Warn Level
#if LOG_LVL < 30
#ifdef NO_LINENOS
#define log_warn(M, ...) fprintf(LOG_DEST, "[WARN] errno: %s: " M "\n", \
                                 clean_errno(), ## __VA_ARGS__); errno = 0
#else
#define log_warn(M, ...) fprintf(LOG_DEST, "[WARN] %s:%d: errno: %s: " M "\n", \
                                 __FILE__, __LINE__, \
                                 clean_errno(), ## __VA_ARGS__); errno = 0
#endif
#else
#define log_warn(M, ...)

#endif

// Error Level is always active
#ifdef NO_LINENOS
#define log_err(M, ...) fprintf(LOG_DEST, "[ERROR] errno: %s: " M "\n", \
                                clean_errno(), ## __VA_ARGS__); exit(1);
#define log_ncerr(e) fprintf(LOG_DEST, "[ERROR] errno: %s \n", \
                             clean_ncerrno(e)); exit(1);
#else
#define log_err(M, ...) fprintf(LOG_DEST, "[ERROR] %s:%d: errno: %s: " M "\n", \
                                __FILE__, __LINE__, \
                                clean_errno(), ## __VA_ARGS__); exit(1);
#define log_ncerr(e) fprintf(LOG_DEST, "[ERROR] %s:%d: errno: %s \n", __FILE__, \
                             __LINE__, clean_ncerrno(e)); exit(1);
#endif

// These depend on previously defined macros

// do not try to be smart and make this go away on LOG_LVL, the _debug
// here means that it just doesn't print a message, it still does the
// check.  MKAY?
#define check_debug(A, M, ...) if (!(A)) {debug(M, ## __VA_ARGS__); errno = 0; \
                                          exit(1); }

#define check(A, M, ...) if (!(A)) {log_err(M, ## __VA_ARGS__); errno = 0; exit( \
                                        1); }

#define sentinel(M, ...)  {log_err(M, ## __VA_ARGS__); errno = 0; exit(1); }

#define check_mem(A) check((A), "Out of memory.")

#define TRACE(C, E) debug("--> %s(%s:%d) %s:%d ", "" # C, State_event_name( \
                              E), E, __FUNCTION__, __LINE__)

#define error_response(F, C, M, ...)  {Response_send_status(F, &HTTP_ ## C); \
                                       sentinel(M, ## __VA_ARGS__); }

#define error_unless(T, F, C, M, ...) if (!(T)) \
        error_response(F, C, M, ## __VA_ARGS__)

#endif
