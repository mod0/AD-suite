include ../../../make.config

OBJ_DIR = ../obj
MOD_DIR = ../modules
OBJS_SUPP = supporting_routines.o
OBJS_AD_SUPP = adBuffer.o adStack.o
OBJS_SIM = numerics_d.o numerics_b.o
OBJS_EXT = sparse_pmgmres_template_tangent_linear.o sparse_pmgmres_template_adjoint.o 
OBJS= ${OBJS_SUPP} ${OBJS_AD_SUPP} ${OBJS_SIM} ${OBJS_EXT}

.PHONY: all numerics.f90 numerics_d.f90 numerics_b.f90 supporting_routines.f90 clean

all:   numerics.f90 numerics_d.f90 numerics_b.f90 supporting_routines.f90 ${OBJS} ${MAINS}

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.c
	${CC} ${CFLAGS} -c $< -o ${OBJ_DIR}/$@ 

numerics.f90: parameters.f90 mathutil.f90 matrix.f90 linsolve.f90 finitevolume.f90 simulation.f90
	@echo "Creating a single numerics file from all the other application files"
	paste --delimiter=\\n --serial $^ > $@

numerics_d.f90: numerics.f90
	${TAPENADE} -ext sparse_pmgmres_template_lib -head "wrapper(mu, sigma)\(oil)" -forward numerics.f90 -o numerics -noisize 

numerics_b.f90: numerics.f90
	${TAPENADE} -ext sparse_pmgmres_template_lib -head "wrapper(mu, sigma)\(oil)" -reverse numerics.f90 -o numerics -noisize 

supporting_routines.f90: gnufor2.f90 utils.f90 mgmres.f90 
	@echo "Creating a single supporting routines file from supporting files"
	paste --delimiter=\\n --serial $^ > $@

clean:
	rm -f ${OBJ_DIR}/*
	rm -f ${MOD_DIR}/*  
