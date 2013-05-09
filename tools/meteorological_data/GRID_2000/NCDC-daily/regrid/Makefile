# Makefile for regrid code
#
# NOTES
# If the code is to be run on another operating system "make clean" before
# recompiling.
# It is recommended that this code is NOT optimised above O2, due to the
# dynamic memory allocation.
# The -C flag will not work on this code due to the dynamic memory allocation.

F77 = f77 

FFLAGS  = -O
 
TARGET = ../regrd
 
OBJECTS = \
	indexx.o \
	isaleap.o \
	iso_weight.o \
	latlong.o \
	near_neigh.o \
	symap.o\
	chgtemp.o\
	grd.o 

SOURCE = \
	indexx.f \
	isaleap.f \
	iso_weight.f \
	latlong.f \
	near_neigh.f \
	symap.f\
	chgtemp.f\
	grd.f

LIBS = 
 
.f.o:
	$(F77) $(FFLAGS) -c $< -o $@
 
$(TARGET): $(OBJECTS)
	$(F77) $(FFLAGS) $(OBJECTS) -o $(TARGET) $(LIBS)

clean:
	/bin/rm -f *.o $(TARGET)
