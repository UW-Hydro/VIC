cmake_minimum_required(VERSION 3.0)
project(vic)

# ============================================================================
# Settings. Add your include and library path of dependences here.
# ----------------------------------------------------------------------------
set(INCLUDE_PATH    /usr/local/openmpi/include)
set(LIB_PATH        /usr/lib/x86_64-linux-gnu/hdf5/serial /usr/local/openmpi/lib)

set(MPI_COMPILER mpicc)
set(LOG_LVL 5)
# ============================================================================

# Cmake configurations options.
option(CLASSIC      "Compile VIC Classic Driver"    ON)
option(IMAGE        "Compile VIC Image Driver"      ON)
option(SHARED_LIB   "With Vic shared library"       ON)
option(INSTALL      "Install VIC to system path"    ON)

set(CMAKE_INSTALL_PREFIX /usr/local/vic)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} -DLOG_LVL=${LOG_LVL} -O3)
set(CMAKE_C_FLAGS_DEBUG  ${CMAKE_C_FLAGS_DEBUG})

include_directories(./vic_run/include ./drivers/shared_all/include
        ./drivers/shared_image/include ./drivers/image/include
        ./drivers/classic/include ${INCLUDE_PATH})

link_directories(${LIB_PATH})

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

if(CLASSIC)
    set(CLASSIC_SRC ${VIC_CLASSIC_SRC} ${VIC_SHARE_SRC} ${VIC_RUN_SRC}
                    ${VIC_CLASSIC_HDR} ${VIC_SHARE_HDR} ${VIC_RUN_HDR})

    add_executable(vic_classic ${CLASSIC_SRC})
    target_link_libraries(vic_classic m)
endif()

if(IMAGE)
    set(IMAGE_SRC ${VIC_IMAGE_SRC} ${VIC_SHARE_SRC} ${VIC_SHARE_IMG_SRC} ${VIC_RUN_SRC}
                  ${VIC_IMAGE_HDR} ${VIC_SHARE_HDR} ${VIC_SHARE_IMG_HDR} ${VIC_RUN_HDR})
    set(CMAKE_C_COMPILER ${MPI_COMPILER})

    add_executable(vic_image ${IMAGE_SRC})
    target_link_libraries(vic_image m netcdf mpi)
endif()

if(SHARED_LIB)
    add_library(vic SHARED ${VIC_RUN_SRC} ${VIC_RUN_HDR})
    target_link_libraries(vic_classic m)
endif()

if(INSTALL)
    if(SHARED_LIB)
        file(GLOB_RECURSE SHLIB libvic*.so)
        install(FILES ${SHLIB} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/)
        install(FILES ${VIC_RUN_HDR} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/)
    endif()
endif()
