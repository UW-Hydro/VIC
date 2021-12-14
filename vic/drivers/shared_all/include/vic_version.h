/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for build time metadata
 *****************************************************************************/

#ifndef VIC_VERSION_H
#define VIC_VERSION_H

#define STR_HELPER(x) # x
#define STR(x) STR_HELPER(x)

#ifndef VERSION
#define VERSION "5.0.1 February 1, 2017"
#endif

#ifndef SHORT_VERSION
#define SHORT_VERSION "5.0.1"
#endif

#ifndef GIT_VERSION
#define GIT_VERSION "unset"
#endif

#ifndef USERNAME
#define USERNAME "unset"
#endif

#ifndef HOSTNAME
#define HOSTNAME "unset"
#endif

#define BUILD_DATE __DATE__
#define BUILD_TIME __TIME__

/* Get compiler metadata */
#if defined(__clang__)
/* Clang/LLVM. ---------------------------------------------- */
# define COMPILER "clang"
# define COMPILER_VERSION __clang_version__

#elif defined(__ICC) || defined(__INTEL_COMPILER)
/* Intel ICC/ICPC. ------------------------------------------ */
# define COMPILER "icc"
# define COMPILER_VERSION __VERSION__

#elif defined(__GNUC__) || defined(__GNUG__)
/* GNU GCC/G++. --------------------------------------------- */
# define COMPILER "gcc"
# define COMPILER_VERSION STR(__GNUC__) "." STR(__GNUC_MINOR__) "." STR( \
        __GNUC_PATCHLEVEL__)

#elif defined(__PGI)
/* Portland Group PGCC/PGCPP. ------------------------------- */
# define COMPILER "pgcc"
# define COMPILER_VERSION __PGIC__ __PGIC_MINOR __PGIC_PATCHLEVEL__

#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)
/* Oracle Solaris Studio. ----------------------------------- */
# define COMPILER "pgcc"
# define COMPILER_VERSION __PGIC__ __PGIC_MINOR __PGIC_PATCHLEVEL__
#endif
#ifndef COMPILER
# define COMPILER "unknown"
# define COMPILER_VERSION "unknown"
#endif

/* C Standard */
#ifdef __STDC_VERSION__
# define CSTANDARD __STDC_VERSION__
#endif
#ifndef CSTANDARD
# define CSTANDARD "unknown"
#endif

/* Platform */

#ifdef __APPLE__
# define PLATFORM "APPLE"
#elif __linux__
# define PLATFORM "LINUX"
#elif __unix__ // all unices, not all compilers
# define PLATFORM "UNIX"
#endif
#ifndef PLATFORM
# define PLATFORM "unknown"
#endif

#endif
