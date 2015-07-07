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

# At the moment, removing this flag will 
# not work anyway. Use -f Makefile_continuous_adjoint
# to build the continuous adjoint vertion
DEFINE_FLAG=-DUSE_OPENAD
TARGET = powergrid
RTSUPP = w2f__types OAD_active OAD_cp OAD_tape OAD_rev
DEPENDENCIES = blas.o linpack.o timer.o \
	       lbfgsb.o print_active.o \
	       gnufor2.o numerics_openad.pre.xb.x2w.w2f.post.o main.o

all: ${TARGET}
	@echo Compiling powergrid
	${FC} ${FFLAGS} $(addsuffix .o, $(RTSUPP) iaddr) ${DEPENDENCIES} -o powergrid

powergrid: $(addsuffix .o, $(RTSUPP) iaddr) ad_inline.o  ${DEPENDENCIES}


numerics_openad.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90, $(RTSUPP)) iaddr.c: toolChain

toolChain: numerics_openad.f90
	openad -c -m rj $<

numerics_openad.f90: numerics.f90
	 gcc -cpp -E -C -P ${DEFINE_FLAG} numerics.f90 > numerics_openad.f90

ad_inline.f: toolChain

%.o : %.f90
	${FC} ${FFLAGS} -c $< -o $@

%.o : %.f
	${FC} ${FFLAGS} -c $< -o $@

main.o:main.f90
	${FC} -cpp ${FFLAGS} ${DEFINE_FLAG} -c $< -o $@

clean:
	rm -f *.o *.mod powergrid openad.*.log~ numerics_openad.* *.mod-whirl $(addsuffix .f90, $(RTSUPP)) ad_inline.f ad_template.f iaddr.c numerics.s
