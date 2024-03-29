 ##############################################################################
 # @section DESCRIPTION
 #
 # VIC Classic Driver Makefile
 ##############################################################################

# -----------------------------------------------------------------------
# SET ENVIRONMENT-SPECIFIC OPTIONS HERE
# -----------------------------------------------------------------------

# VIC DRIVER PATH (include and src are subdirs of this)
DRIVERPATH = .

# VIC RUN PATH
SHAREDPATH = ../shared_all

# VIC RUN PATH
VICPATH = ../../vic_run

# Set SHELL = your shell here
SHELL = /bin/bash

# Set CC = your compiler here
ifndef CC
CC=gcc
endif

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
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

# set include
INCLUDES = -I ${DRIVERPATH}/include -I $(SHAREDPATH)/include -I ${VICPATH}/include

# Uncomment for normal optimized code flags (fastest run option)
#CFLAGS  = -O3 -Wall -Wno-unused
# LIBRARY = -lm

# Uncomment to include debugging information
CFLAGS  =  ${INCLUDES} -g -Wall -Wextra -std=c99 \
					 -DLOG_LVL=$(LOG_LVL) \
					 -DGIT_VERSION=\"$(GIT_VERSION)\" \
					 -DUSERNAME=\"$(USER)\" \
					 -DHOSTNAME=\"$(HOSTNAME)\"
LIBRARY = -lm

# Uncomment to include execution profiling information
#CFLAGS  = ${INCLUDES} -O3 -pg -Wall -Wno-unused -DLOG_LVL=$(LOG_LVL)
#LIBRARY = -lm

# Uncomment to debug memory problems using electric fence (man efence)
#CFLAGS  = ${INCLUDES} -g -Wall -Wno-unused -DLOG_LVL=$(LOG_LVL)
#LIBRARY = -lm -lefence -L/usr/local/lib

COMPEXE = vic_classic
EXT = .exe

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
	\rm -f core log
	\rm -rf ${COMPEXE}${EXT} ${COMPEXE}${EXT}.dSYM

model: $(OBJS)
	$(CC) -o ${COMPEXE}${EXT} $(OBJS) $(CFLAGS) $(LIBRARY)

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
	$(CC) $(CFLAGS) -M $(SRCS) > $@

clean::
	\rm -f .depend
