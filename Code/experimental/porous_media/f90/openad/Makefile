# setup the compiler and flags
FORT=gfortran
FFLAGS=-g -O0

TARGET = pm

OBJS90=fluid.o grid.o mathutil.o \
     matrix.o mgmres.o linsolve.o \
     fvm.o gnufor2.o runspe10.o

%.o : %.f
	${FORT} -cpp ${FFLAGS} -c -g $< -o $@

%.o : %.f90
	${FORT} -cpp ${FFLAGS} -c -g $< -o $@


all: ${TARGET}

pm: ${OBJS90}
	${FORT} ${FFLAGS} $^ -g -o porous_media

clean:
	rm *.o *.mod porous_media
