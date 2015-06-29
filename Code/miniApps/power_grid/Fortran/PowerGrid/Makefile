FC = gfortran
FFLAGS = -cpp

TARGET = powergrid
DEPENDENCIES = blas.o linpack.o timer.o \
	       lbfgsb.o constants.o print_active.o \
	       gnufor2.o main.o

%.o : %.f
	${FC} ${FFLAGS} -c $< -o $@

%.o : %.f90
	${FC} ${FFLAGS} -c $< -o $@

all: ${TARGET}
	@echo Compiling powergrid
	${FC} ${FFLAGS} ${DEPENDENCIES} -o powergrid

powergrid: ${DEPENDENCIES}
	

clean:
	rm *.o *.mod powergrid
