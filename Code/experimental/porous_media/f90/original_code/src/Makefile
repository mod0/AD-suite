include ../../../make.config

OBJ_DIR = ../obj
MOD_DIR = ../modules
PARAMS_FILE = ../../../data/data_${SCENARIO}/parameters.f90
OBJS_SUPP = supporting_routines.o
OBJS_SIM = numerics.o
OBJS= ${OBJS_SUPP} ${OBJS_SIM}

.PHONY: all numerics.f90 supporting_routines.f90 clean

all:  numerics.f90 supporting_routines.f90 ${OBJS} ${MAINS}

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

numerics.f90: ${PARAMS_FILE} mathutil.f90 matrix.f90 linsolve.f90 finitevolume.f90 simulation.f90
	@echo "Creating a single numerics file from all the other application files"
	paste --delimiter=\\n --serial $^ > $@
	@echo "Replacing the dummy linear solver inside the numerics file with that from the template"
	sed -i -ne '/!%>LINSOLVE/ {p; r sparse_pmgmres_template.f90' -e ':a; n; /!%<LINSOLVE/ {p; b}; ba}; p' numerics.f90
	@echo "Replacing the call to the dummy linear solver with the call to sparse pmgmres"
	sed -i 's/call sparse_dummy_method/call sparse_pmgmres_method/' numerics.f90

supporting_routines.f90: gnufor2.f90 utils.f90 mgmres.f90
	@echo "Creating a single supporting routines file from supporting files"
	paste --delimiter=\\n --serial $^ > $@

clean:
	rm -f ${OBJ_DIR}/*
	rm -f ${MOD_DIR}/*  

