# make a list of all *.h files in routing folder
INCL_ROUT := $(wildcard ${EXTPATH}/${ROUT}/include/*.h)

# make a list of all *.c files in routing folder
SRCS_ROUT := $(wildcard ${EXTPATH}/${ROUT}/src/*.c)

# convert the list of all *.c to a list of *.o files (object files)
OBJS_ROUT = $(SRCS_ROUT:%.o=%.c)
