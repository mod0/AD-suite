# setup the compiler and flags
FORT=gfortran
FFLAGS=-g -O0

TARGET = pm

OBJS90=fluid.o grid.o mathutil.o \
     matrix.o mgmres.o linsolve.o \
     fvm.o gnufor2.o runspe10.o

OBJS77=print_active.o

%.o : %.f
	${FORT} -cpp ${FFLAGS} -c $< -o $@

%.o : %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o $@


all: ${TARGET}

pm: ${OBJS77} ${OBJS90}
	${FORT} ${FFLAGS} $^ -o porous_media

clean:
	rm *.o *.mod porous_media

