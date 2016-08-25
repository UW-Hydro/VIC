cmake_minimum_required(VERSION 3.0)
project(vic)

# ============================================================================
# Settings.
# ----------------------------------------------------------------------------

set(OTHER_INCLUDE_PATH  your_include_path)
set(OTHER_LIBRARY_PATH  your_library_path)
set(MPI_COMPILER        gcc)

# ============================================================================

# Cmake configurations options.
option(CLASSIC      "Compile VIC Classic Driver"    ON)
option(IMAGE        "Compile VIC Image Driver"      ON)
option(SHARED_LIB   "With Vic shared library"       ON)
option(INSTALL      "Install VIC to system path"    ON)

set(CMAKE_INSTALL_PREFIX /usr/local/vic)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
set(CMAKE_C_FLAGS_DEBUG  ${CMAKE_C_FLAGS_DEBUG} )

include_directories(./vic_run/include ./drivers/shared_all/include
        ./drivers/image/include ./drivers/shared_image/include
        ./drivers/classic/include ${OTHER_INCLUDE_PATH})
link_directories(${OTHER_LIBRARY_PATH})

file(GLOB VIC_RUN_SRC ./vic_run/src/*.c)
file(GLOB VIC_CLASSIC_SRC ./drivers/classic/src/*.c)
file(GLOB VIC_IMAGE_SRC ./drivers/image/src/*.c)

file(GLOB VIC_RUN_HDR ./vic_run/include/*.h)
file(GLOB VIC_CLASSIC_HDR ./drivers/classic/include/*.h)
file(GLOB VIC_IMAGE_HDR ./drivers/image/include/*.h)

file(GLOB VIC_SHARE_SRC ./drivers/shared_all/src/*.c)
file(GLOB VIC_SHARE_HDR ./drivers/shared_all/hdr/*.h)

file(GLOB VIC_SHARE_IMG_SRC ./drivers/shared_image/src/*.c)
file(GLOB VIC_SHARE_IMG_HDR ./drivers/shared_image/hdr/*.h)

# ----------------------------------------------------------------------------
# VIC Classic Driver
# ----------------------------------------------------------------------------
if(CLASSIC)
    set(CLASSIC_SRC ${VIC_CLASSIC_SRC} ${VIC_SHARE_SRC} ${VIC_RUN_SRC}
                    ${VIC_CLASSIC_HDR} ${VIC_SHARE_HDR} ${VIC_RUN_HDR})

    add_executable(vic_classic ${CLASSIC_SRC})
    target_link_libraries(vic_classic m)
endif()

# ----------------------------------------------------------------------------
# VIC Image Driver
# ----------------------------------------------------------------------------

#           Find out the include and library directory of MPI.
if(IMAGE)
    set(MPI_INCLUDE_SEARCH_PATHS
            /usr/include
            /usr/include/openmpi
            /usr/local/include
            /usr/local/include/openmpi
            /usr/local/openmpi/include)
    set(MPI_LIB_SEARCH_PATHS
            /lib/
            /lib/openmpi
            /lib64/
            /usr/lib
            /usr/lib/openmpi
            /usr/lib/x86_64-linux-gnu
            /usr/lib64
            /usr/local/lib
            /usr/local/lib64
            /usr/local/openmpi/lib)
        
    find_path(MPI_INCLUDE_DIR NAMES mpi.h PATHS ${MPI_INCLUDE_SEARCH_PATHS})
    find_path(MPI_LIB_DIR NAMES libmpi.so PATHS ${MPI_LIB_SEARCH_PATHS})

    if(NOT MPI_INCLUDE_DIR)
        message(STATUS "Could not find MPI include. VIC Image Driver would not be built.")
        set(BUILD_IMAGE OFF)
    else()
        message(STATUS "MPI include at ${MPI_INCLUDE_DIR}")
        set(BUILD_IMAGE ON)
        include_directories(${MPI_INCLUDE_DIR})
    endif()

    if(NOT MPI_LIB_DIR)
        message(STATUS "Could not find MPI library. VIC Image Driver would not be built.")
        set(BUILD_IMAGE OFF)
    else()
        message(STATUS "MPI library at ${MPI_LIB_DIR}")
        set(BUILD_IMAGE ON)
        link_directories(${MPI_LIB_DIR})
    endif()
endif()

#           Find out the include and library directory of netCDF.
if(IMAGE)
    set(NETCDF_INCLUDE_SEARCH_PATHS
            /usr/include
            /usr/include/netcdf
            /usr/local/include
            /usr/local/include/netcdf
            /usr/local/netcdf/include)
    set(NETCDF_LIB_SEARCH_PATHS
            /lib/
            /lib64/
            /usr/lib
            /usr/lib/netcdf
            /usr/lib/x86_64-linux-gnu
            /usr/lib64
            /usr/local/lib
            /usr/local/lib64
            /usr/local/netcdf/lib)

    find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h PATHS ${NETCDF_INCLUDE_SEARCH_PATHS})
    find_path(NETCDF_LIB_DIR NAMES libnetcdf.so PATHS ${NETCDF_LIB_SEARCH_PATHS})

    if(NOT NETCDF_INCLUDE_DIR)
        message(STATUS "Could not find netCDF include. VIC Image Driver would not be built.")
        set(BUILD_IMAGE OFF)
    else()
        message(STATUS "NetCDF include at ${NETCDF_LIB_DIR}")
        set(BUILD_IMAGE ON)
        include_directories(${NETCDF_INCLUDE_DIR})
    endif()

    if(NOT NETCDF_LIB_DIR)
        message(STATUS "Could not find netCDF library. VIC Image Driver would not be built.")
        set(BUILD_IMAGE OFF)
    else()
        message(STATUS "NetCDF library at ${NETCDF_LIB_DIR}")
        set(BUILD_IMAGE ON)
        link_directories(${NETCDF_LIB_DIR})
    endif()
endif()

#           Build Image Driver if MPI and netCDF were found.
if(BUILD_IMAGE)
    set(IMAGE_SRC ${VIC_IMAGE_SRC} ${VIC_SHARE_SRC} ${VIC_SHARE_IMG_SRC} ${VIC_RUN_SRC}
                  ${VIC_IMAGE_HDR} ${VIC_SHARE_HDR} ${VIC_SHARE_IMG_HDR} ${VIC_RUN_HDR}
                  ${MPI_INCLUDE_DIR}/mpi.h ${NETCDF_INCLUDE_DIR}/netcdf.h)
    set(CMAKE_C_COMPILER ${MPI_COMPILER})

    add_executable(vic_image ${IMAGE_SRC})
    target_link_libraries(vic_image m netcdf mpi)
endif()

# ----------------------------------------------------------------------------
# VIC Shared library and it's install.
# ----------------------------------------------------------------------------
if(SHARED_LIB)
    add_library(vic SHARED ${VIC_RUN_SRC} ${VIC_RUN_HDR})
    target_link_libraries(vic m)
endif()

if(INSTALL)
    if(SHARED_LIB)
        file(GLOB_RECURSE SHLIB libvic*.so)
        install(FILES ${SHLIB} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/)
        install(FILES ${VIC_RUN_HDR} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/)
    endif()
endif()
