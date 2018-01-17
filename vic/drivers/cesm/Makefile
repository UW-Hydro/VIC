 ##############################################################################
 # @section DESCRIPTION
 #
 # VIC CESM Driver Makefile
 # FOR BUILD TESTING ONLY
 #
 # @section LICENSE
 #
 # The Variable Infiltration Capacity (VIC) macroscale hydrological model
 # Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
 # and Environmental Engineering, University of Washington.
 #
 # The VIC model is free software; you can redistribute it and/or
 # modify it under the terms of the GNU General Public License
 # as published by the Free Software Foundation; either version 2
 # of the License, or (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License along with
 # this program; if not, write to the Free Software Foundation, Inc.,
 # 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 ##############################################################################

# -----------------------------------------------------------------------
# SET ENVIRONMENT-SPECIFIC OPTIONS HERE
# -----------------------------------------------------------------------

# VIC DRIVER PATH (include and src are subdirs of this)
DRIVERPATH = .

# VIC RUN PATH
SHAREDPATH = ../shared_all
SHAREDIMAGEPATH = ../shared_image

# VIC RUN PATH
VICPATH = ../../vic_run

# NETCDF include .. change the path depending on travis or other location
ifeq (true, ${TRAVIS})
NETCDFPATH = ${TRAVIS_NETCDFPATH}
else ifeq (hydro-c1, $(shell hostname))
NETCDFPATH = /usr/local/netcdf-4.3.0+gcc-4.7.2
else
NETCDFPATH = /opt/local
endif

# MPI include .. change the path depending on travis or other location
ifeq (true, ${TRAVIS})
MPIPATH = ${TRAVIS_MPIPATH}
else ifeq (hydro-c1, $(shell hostname))
MPIPATH = /usr/local/netcdf-4.3.0+gcc-4.7.2
else
MPIPATH = /opt/local
endif

# VIC EXTENTIONS PATH
EXTPATH = ../../extensions

# Set the routing model extension
include ${EXTPATH}/rout_stub/rout.mk

# Set SHELL = your shell here
SHELL = /bin/bash

# Set CC = your compiler here
ifeq (true, ${TRAVIS})
CC = ${MPIPATH}/bin/mpicc
else
CC = mpicc
endif

ifeq (true, ${TRAVIS})
FC = ${MPIPATH}/bin/mpif90
else
FC = mpifort
endif

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
HOSTNAME := $(shell uname -n)

# set includes
INCLUDES = -I ${DRIVERPATH}/include \
		   -I ${VICPATH}/include \
		   -I ${SHAREDPATH}/include \
		   -I ${SHAREDIMAGEPATH}/include \
		   -I ${NETCDFPATH}/include \
		   -I ${EXTPATH}/rout_stub/include

# Set libraries
LIBRARY = -lm -L${NETCDFPATH}/lib -lnetcdf

# Set compiler flags
CFLAGS  =  ${INCLUDES} -ggdb -O0 -Wall -Wextra -fPIC \
					 -fopenmp \
					 -DLOG_LVL=$(LOG_LVL) \
					 -DGIT_VERSION=\"$(GIT_VERSION)\" \
					 -DUSERNAME=\"$(USER)\" \
					 -DHOSTNAME=\"$(HOSTNAME)\"
FFLAGS  = ${INCLUDES} -fPIC

# Set the log level
# To turn off warning statements, set LOG_LVL >= 30
# | Level     | Numeric value    |
# |---------  |---------------   |
# | ERROR     | Always Active    |
# | WARNING   | < 30             |
# | INFO      | < 20             |
# | DEBUG     | < 10             |
LOG_LVL = 5

COMPLIB = lndlib
EXT = .a

LDFLAGGS = -shared


# -----------------------------------------------------------------------
# MOST USERS DO NOT NEED TO MODIFY BELOW THIS LINE
# -----------------------------------------------------------------------

HDRS = \
	$(wildcard ${VICPATH}/include/*.h) \
	$(wildcard ${DRIVERPATH}/include/*.h) \
	$(wildcard ${SHAREDPATH}/include/*.h) \
	$(wildcard ${SHAREDIMAGEPATH}/include/*.h) \
	$(wildcard ${EXTPATH}/rout_stub/include/*.h)

CSRCS = \
	$(wildcard ${VICPATH}/src/*.c) \
	$(wildcard ${DRIVERPATH}/src/*.c) \
	$(wildcard ${SHAREDPATH}/src/*.c) \
	$(wildcard ${SHAREDIMAGEPATH}/src/*.c) \
	$(wildcard ${EXTPATH}/rout_stub/src/*.c)

FSRCS = \
	$(wildcard ${DRIVERPATH}/src/*.F90)

COBJS = $(CSRCS:%.o=%.c)
FOBJS = $(FSRCS:%.o=%.F90)

all:
	make depend
	make model

default:
	make depend
	make model

full:
	make clean
	make depend
	make tags
	make model

clean::
	\rm -f .depend

model: $(OBJS)
	$(FC) $(LDFLAGGS) $(COBJS) $(FOBJS) $(FFLAGS) $(LIBRARY) -o ${COMPLIB}${EXT}

# -------------------------------------------------------------
# tags
# so we can find our way around
# -------------------------------------------------------------
tags:	TAGS
TAGS:	$(SRCS) $(HDRS)
	etags $(SRCS) $(HDRS)
clean::
	rm -f TAGS


# -------------------------------------------------------------
# depend
# -------------------------------------------------------------
depend: .depend
.depend: $(CSRCS) $(FSRCS) $(HDRS)
	$(CC) $(CFLAGS) $(INCLUDES) -M $(CSRCS) -DLOG_LVL=$(LOG_LVL) > $@
	$(FC) $(FFLAGS) $(INCLUDES) -M $(FSRCS) > $@

clean::
	rm -f .depend
