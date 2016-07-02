# setup the compiler and flags
FORT=gfortran
FFLAGS=-g -O0

OBJS90=gnufor2.o utils.o mgmres.o numerics.o
MAINS=main.o main_oil_vs_mu_sigma.o

all:  ${OBJS90} ${MAINS}
	${FORT} ${FFLAGS} ${OBJS90} main.o -o porous_media
	${FORT} ${FFLAGS} ${OBJS90} main_oil_vs_mu_sigma.o -o porous_oil_mu_vs_sigma

%.o : %.f
	${FORT} -cpp ${FFLAGS} -c $< -o $@

%.o : %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o $@

clean:
	rm *.o *.mod porous_media porous_oil_mu_vs_sigma

