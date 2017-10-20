 ##############################################################################
 # @section DESCRIPTION
 #
 # VIC Image Driver Makefile
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

# VIC EXTENTIONS PATH
EXTPATH = ../../extensions

# Set the routing model extension
# | Level      | Value                              |
# |----------  |---------------------------------   |
# | rout_stub  | Stub routing model (no routing)    |
# | rout_rvic  | Use RVIC routing scheme            |
ROUT=rout_stub

include ${EXTPATH}/${ROUT}/rout.mk

ifndef NC_LIBS
NC_LIBS = $(shell nc-config --libs)
endif
$(info $$NC_LIBS is [${NC_LIBS}])

ifndef NC_CFLAGS
NC_CFLAGS = $(shell nc-config --cflags)
endif
$(info $$NC_CFLAGS is [${NC_CFLAGS}])

ifndef MPICC
MPICC = mpicc
endif
$(info $$MPICC is [${MPICC}])

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
UNAME_S := $(shell uname -s)
HOSTNAME := $(shell uname -n)

# Set the log level
# To turn off warning statements, set LOG_LVL >= 30
# | Level     | Numeric value    |
# |---------  |---------------   |
# | ERROR     | Always Active    |
# | WARNING   | < 30             |
# | INFO      | < 20             |
# | DEBUG     | < 10             |
LOG_LVL = 5

# set includes
INCLUDES = -I ${DRIVERPATH}/include \
		   -I ${VICPATH}/include \
		   -I ${SHAREDPATH}/include \
		   -I ${SHAREDIMAGEPATH}/include \
		   -I ${EXTPATH}/${ROUT}/include

# Uncomment to include debugging information
CFLAGS  =  ${INCLUDES} ${NC_CFLAGS} -ggdb -O0 -Wall -Wextra -std=c99 \
					 -fopenmp \
					 -DLOG_LVL=$(LOG_LVL) \
					 -DGIT_VERSION=\"$(GIT_VERSION)\" \
					 -DUSERNAME=\"$(USER)\" \
					 -DHOSTNAME=\"$(HOSTNAME)\"

ifeq (true, ${TRAVIS})
# Add extra debugging for builds on travis
CFLAGS += -rdynamic -Wl,-export-dynamic
endif

LIBRARY = -lm ${NC_LIBS}

COMPEXE = vic_image
EXT = .exe

# -----------------------------------------------------------------------
# MOST USERS DO NOT NEED TO MODIFY BELOW THIS LINE
# -----------------------------------------------------------------------

HDRS = \
	$(wildcard ${VICPATH}/include/*.h) \
	$(wildcard ${DRIVERPATH}/include/*.h) \
	$(wildcard ${SHAREDPATH}/include/*.h) \
	$(wildcard ${SHAREDIMAGEPATH}/include/*.h) \
	$(wildcard ${EXTPATH}/${ROUT}/include/*.h)

SRCS = \
	$(wildcard ${VICPATH}/src/*.c) \
	$(wildcard ${DRIVERPATH}/src/*.c) \
	$(wildcard ${SHAREDPATH}/src/*.c) \
	$(wildcard ${SHAREDIMAGEPATH}/src/*.c) \
	$(wildcard ${EXTPATH}/${ROUT}/src/*.c)

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
	\rm -f core log
	\rm -rf ${COMPEXE}${EXT} ${COMPEXE}${EXT}.dSYM

model: $(OBJS)
	$(MPICC) -o ${COMPEXE}${EXT} $(OBJS) $(CFLAGS) $(LIBRARY)

# -------------------------------------------------------------
# tags
# so we can find our way around
# -------------------------------------------------------------
tags:	TAGS
TAGS:	$(SRCS) $(HDRS)
	etags $(SRCS) $(HDRS)
clean::
	\rm -f TAGS


# -------------------------------------------------------------
# depend
# -------------------------------------------------------------
depend: .depend
.depend: $(SRCS) $(HDRS)
	$(MPICC) $(CFLAGS) -M $(SRCS) > $@

clean::
	\rm -f .depend
