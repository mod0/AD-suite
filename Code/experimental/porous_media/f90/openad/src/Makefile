include ../../../make.config

OBJ_DIR = ../obj
MOD_DIR = ../modules
OBJS_SUPP = supporting_routines.o
OBJS_AD_SUPP = w2f__types.o OAD_active.o OAD_tape.o OAD_rev.o OAD_cp.o iaddr.o
OBJS_SIM = numerics_f.o numerics_rj.o numerics_rs.o
OBJS= ${OBJS_SUPP} ${OBJS_AD_SUPP} ${OBJS_SIM}

.PHONY: all numerics.f90 numerics_f.f90 numerics_rj.f90 numerics_rs.f90 supporting_routines.f90 clean

all:   numerics.f90 numerics_f.f90 numerics_rj.f90 numerics_rs.f90 supporting_routines.f90 ${OBJS} 

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.c
	${CC} ${CFLAGS} -c $< -o ${OBJ_DIR}/$@ 

numerics.f90: parameters.f90 mathutil.f90 matrix.f90 linsolve.f90 finitevolume.f90 simulation.f90
	@echo "Creating a single numerics file from all the other application files"
	paste --delimiter=\\n --serial $^ > $@

numerics_f.f90: numerics.f90
	${OPENAD} -c -m f numerics.f90 
	cp numerics.pre.xb.x2w.w2f.post.f90 numerics_f.f90
	@echo "Replace sparse_dummy_method by sparse_pmgmres_method"
	sed -i -ne '/subroutine sparse_dummy_method/I{r sparse_pmgmres_template_tangentlin.f90' -e ':a; n; /end subroutine/I{d; b}; ba}; p' numerics_f.f90
	@echo "Replacing the call to the dummy linear solver with the call to tangent linear version of sparse pmgmres"
	sed -i 's/call sparse_dummy_method/call sparse_pmgmres_method/I' numerics_f.f90
	rm -f numerics.pre.*
	rm -f *.log~
	rm -f *.mod-whirl

numerics_rj.f90: numerics.f90
	${OPENAD} -c -m rj numerics.f90 
	cp numerics.pre.xb.x2w.w2f.post.f90 numerics_rj.f90
	rm -f numerics.pre.*
	rm -f *.log~
	rm -f *.mod-whirl

numerics_rs.f90: numerics.f90
	${OPENAD} -c -m rs numerics.f90 
	cp numerics.pre.xb.x2w.w2f.post.f90 numerics_rs.f90
	rm -f numerics.pre.*
	rm -f *.log~
	rm -f *.mod-whirl

supporting_routines.f90: gnufor2.f90 utils.f90 mgmres.f90 
	@echo "Creating a single supporting routines file from supporting files"
	paste --delimiter=\\n --serial $^ > $@

clean:
	rm -f ${OBJ_DIR}/*
	rm -f ${MOD_DIR}/*  
