 ##############################################################################
 # @section DESCRIPTION
 #
 # VIC Routing Extension Makefile Include
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

# make a list of all *.h files in routing folder
INCL_ROUT := $(wildcard ${EXTPATH}/${ROUT}/include/*.h)

# make a list of all *.c files in routing folder
SRCS_ROUT := $(wildcard ${EXTPATH}/${ROUT}/src/*.c)

# convert the list of all *.c to a list of *.o files (object files)
OBJS_ROUT = $(SRCS_ROUT:%.o=%.c)
