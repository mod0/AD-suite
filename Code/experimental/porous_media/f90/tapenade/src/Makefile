include ../../../make.config

OBJ_DIR = ../obj
MOD_DIR = ../modules
OBJS_SUPP = supporting_routines.o
OBJS_SIM = numerics_d.o numerics_b.o
OBJS= ${OBJS_SUPP} ${OBJS_SIM}

PARAMS_FILE = ../../../data/data_${SCENE}/parameters.f90

.PHONY: all numerics.f90 numerics_d.f90 numerics_b.f90 supporting_routines.f90 clean

all:  numerics.f90 numerics_d.f90 numerics_b.f90 supporting_routines.f90 ${OBJS} ${MAINS}

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

numerics.f90: ${PARAMS_FILE} mathutil.f90 matrix.f90 linsolve.f90 finitevolume.f90 simulation.f90
	@echo "Creating a single numerics file from all the other application files"
	paste --delimiter=\\n --serial $^ > $@

numerics_d.f90: numerics.f90
	${TAPENADE} -head "wrapper(mu, sigma)\(oil)" -forward numerics.f90 -o numerics -noisize 
	@echo "Replacing the differentiated dummy linear solver inside the numerics file with that from the template"
	sed -i -ne '/subroutine sparse_dummy_method_d/I{r sparse_pmgmres_template_tangentlin.f90' -e ':a; n; /end subroutine sparse_dummy_method_d/I{d; b}; ba}; p' numerics_d.f90
	@echo "Replacing the call to the differentiated dummy linear solver with the call to tangent linear version of sparse pmgmres"
	sed -i 's/call sparse_dummy_method_d/call sparse_pmgmres_method_d/I' numerics_d.f90
	@echo "Replacing the dummy linear solver inside the numerics file with that from the template"
	sed -i -ne '/subroutine sparse_dummy_method\b/I{r sparse_pmgmres_template.f90' -e ':a; n; /end subroutine sparse_dummy_method\b/I{d; b}; ba}; p' numerics_d.f90
	@echo "Replacing the call to the dummy linear solver with the call to sparse pmgmres"
	sed -i 's/call sparse_dummy_method\b/call sparse_pmgmres_method/I' numerics_d.f90
	sed -i -e 's/use mathutil\b/use mathutil_d/I' numerics_d.f90

numerics_b.f90: numerics.f90
	${TAPENADE} -head "wrapper(mu, sigma)\(oil)" -reverse numerics.f90 -o numerics -noisize 
	@echo "Replacing the differentiated dummy linear solver inside the numerics file with that from the template"
	sed -i -ne '/subroutine sparse_dummy_method_b/I{r sparse_pmgmres_template_adjoint.f90' -e ':a; n; /end subroutine sparse_dummy_method_b/I{d; b}; ba}; p' numerics_b.f90                                       
	@echo "Replacing the call to the differentiated dummy linear solver with the call to adjoint version of sparse pmgmres"
	sed -i 's/call sparse_dummy_method_b/call sparse_pmgmres_method_b/I' numerics_b.f90
	@echo "Replacing the dummy linear solver inside the numerics file with that from the template"
	sed -i -ne '/subroutine sparse_dummy_method\b/I{r sparse_pmgmres_template.f90' -e ':a; n; /end subroutine sparse_dummy_method\b/I{d; b}; ba}; p' numerics_b.f90                                       
	@echo "Replacing the call to the dummy linear solver with the call to sparse pmgmres"
	sed -i 's/call sparse_dummy_method/call sparse_pmgmres_method/I' numerics_b.f90
	sed -i -e 's/use mathutil\b/use mathutil_b/I' numerics_b.f90

supporting_routines.f90: gnufor2.f90 utils.f90 mgmres.f90
	@echo "Creating a single supporting routines file from supporting files"
	paste --delimiter=\\n --serial $^ > $@

clean:
	rm -f ${OBJ_DIR}/*
	rm -f ${MOD_DIR}/*  

