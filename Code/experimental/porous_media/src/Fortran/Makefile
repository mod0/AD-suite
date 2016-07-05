# setup the compiler and flags
FORT=gfortran
FFLAGS=-g -O0
CC=gcc
CFLAGS=-g -O0
OBJS_SUPP = gnufor2.o utils.o mgmres.o
OBJS_SIM = numerics.o
OBJS= ${OBJS_SUPP} ${OBJS_SIM}
MAINS=driver.o 

all:  numerics.f90 ${OBJS} ${MAINS}
	${FORT} ${FFLAGS} ${OBJS} driver.o -o porous_media

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o $@

%.o: %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o $@

numerics.f90: parameters.f90 mathutil.f90 matrix.f90 linsolve.f90 finitevolume.f90 simulation.f90
	paste --delimiter=\\n --serial $^ > $@
	sed -i -ne '/!%>LINSOLVE/ {p; r sparse_pmgmres_template.f90' -e ':a; n; /!%<LINSOLVE/ {p; b}; ba}; p' numerics.f90
	sed -i 's/call sparse_dummy_method/call sparse_pmgmres_method/' numerics.f90


clean:
	rm -f *.o *.mod porous_media 
	rm -i *.png
	rm -i *.txt

