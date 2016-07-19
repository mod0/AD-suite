include ../../../make.config

OBJ_DIR = ../obj
MOD_DIR = ../modules
OBJS_SUPP = supporting_routines.o 
OBJS_SIM = numerics.o
OBJS_EXT = sparse_pmgmres_template.o
OBJS= ${OBJS_SUPP} ${OBJS_SIM} ${OBJS_EXT}

.PHONY: all numerics.f90 supporting_routines.f90 clean

all:  numerics.f90 supporting_routines.f90 ${OBJS} ${MAINS}

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

numerics.f90: parameters.f90 mathutil.f90 matrix.f90 linsolve.f90 finitevolume.f90 simulation.f90
	@echo "Creating a single numerics file from all the other application files"
	paste --delimiter=\\n --serial $^ > $@

supporting_routines.f90: gnufor2.f90 utils.f90 mgmres.f90
	@echo "Creating a single supporting routines file from supporting files"
	paste --delimiter=\\n --serial $^ > $@

clean:
	rm -f ${OBJ_DIR}/*
	rm -f ${MOD_DIR}/*  

