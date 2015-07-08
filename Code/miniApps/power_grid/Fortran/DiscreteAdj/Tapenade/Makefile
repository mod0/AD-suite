#
# OpenAD executable
#
OAD = openad

#
# Fortran compiler and flags
#
FC = gfortran
FFLAGS = -g

#
# C compiler and flags
#
CC = gcc
CFLAGS = -g
DEFINE_FLAG = -DUSE_TAPENADE
TARGET = powergrid
DEPENDENCIES = blas.o linpack.o timer.o \
	       lbfgsb.o print_active.o \
	       gnufor2.o numerics_tapenade.o \
	       adStack.o adBuffer.o main.o

all: ${TARGET}
	@echo Compiling powergrid
	${FC} ${FFLAGS} ${DEPENDENCIES} -o powergrid

powergrid: ${DEPENDENCIES}


%.o : %.f90
	${FC} ${FFLAGS} -c $< -o $@

%.o : %.f
	${FC} ${FFLAGS} -c $< -o $@

main.o:main.f90
	${FC} -cpp ${FFLAGS} ${DEFINE_FLAG} -c $< -o $@

clean:
	rm -f *.o *.mod powergrid main.s
