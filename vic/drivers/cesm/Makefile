 ##############################################################################
 # @section DESCRIPTION
 #
 # VIC CESM Driver Makefile
 # FOR BUILD TESTING ONLY
 #
 # @section LICENSE
 #
 # The Variable Infiltration Capacity (VIC) macroscale hydrological model
 # Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
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
SHAREDPATH = ../shared

# VIC RUN PATH
VICPATH = ../../vic_run

# NETCDF include .. change the path depending on travis or other location
ifeq (true, ${TRAVIS})
NETCDFPATH = /usr
else ifeq (hydro-c1, $(shell hostname))
NETCDFPATH = /usr/local/netcdf-4.3.0+gcc-4.7.2
else
NETCDFPATH = /opt/local
endif

# MPI include .. change the path depending on travis or other location
ifeq (true, ${TRAVIS})
MPIPATH = /usr
else ifeq (hydro-c1, $(shell hostname))
MPIPATH = /usr/local/netcdf-4.3.0+gcc-4.7.2
else
MPIPATH = /opt/local
endif

# Set SHELL = your shell here
SHELL = /bin/bash

# Set CC = your compiler here
ifeq (true, ${TRAVIS})
CC = mpicc
else
CC = mpicc-mpich-mp
endif

# set includes
INCLUDES = -I ${DRIVERPATH}/include -I ${VICPATH}/include -I ${SHAREDPATH}/include -I ${NETCDFPATH}/include -I ${MPIPATH}/include

# Uncomment for normal optimized code flags (fastest run option)
#CFLAGS  = -O3 -Wall -Wno-unused
LIBRARY = -lm -L${NETCDFPATH}/lib -lnetcdf -L${MPIPATH}/lib -lmpi

# Uncomment to include debugging information
CFLAGS  =  ${INCLUDES} -ggdb -O0 -Wall -Wextra -fPIC
#LIBRARY = -lm

# Uncomment to include execution profiling information
#CFLAGS  = ${INCLUDES} -O3 -pg -Wall -Wno-unused
#LIBRARY = -lm

# Uncomment to debug memory problems using electric fence (man efence)
#CFLAGS  = ${INCLUDES} -g -Wall -Wno-unused
#LIBRARY = -lm -lefence -L/usr/local/lib

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
	$(wildcard ${SHAREDPATH}/include/*.h)

SRCS = \
	$(wildcard ${VICPATH}/src/*.c) \
	$(wildcard ${DRIVERPATH}/src/*.c) \
	$(wildcard ${SHAREDPATH}/src/*.c)

OBJS = $(SRCS:%.o=%.c)

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
	rm -f ${DRIVERPATH}/src/*.o
	rm -f ${SHAREDPATH}/src/*.o
	rm -f ${VICPATH}/src/*.o
	rm -f core log
	rm -f ${DRIVERPATH}/src/*~ ${SHAREDPATH}/src/*.~ ${VICPATH}/src/*~
	rm -rf ${COMPLIB} ${COMPLIB}.dSYM

model: $(OBJS)
	$(CC) $(LDFLAGGS) $(OBJS) $(CFLAGS) $(LIBRARY) -o ${COMPLIB}${EXT}

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
.depend: $(SRCS) $(HDRS)
	$(CC) $(CFLAGS) -M $(SRCS) -DLOG_LVL=$(LOG_LVL) > $@

clean::
	rm -f .depend
