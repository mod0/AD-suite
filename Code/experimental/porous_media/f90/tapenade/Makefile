include ../../make.config

OBJ_DIR = obj
MOD_DIR = modules
APP_SRC_DIR = src
OBJS_SUPP = supporting_routines.o
OBJS_SIM = numerics_d.o numerics_b.o
OBJS = ${OBJS_SUPP} ${OBJS_SIM}
OBJS_MAIN = eval_deriv_tapenade_1_forward.o eval_deriv_tapenade_1_reverse.o  

# Determine if this directory needs to be built
# Refer: http://stackoverflow.com/a/28790801/1220493
IS_FORTRAN=$(shell echo ${LANG}|grep -i -c F90)
IS_TAPENADE=$(shell echo ${DERIV}|grep -i -c TAPENADE)
IS_ALL_LANG=$(shell echo ${LANG}|grep -i -c ALL)
IS_ALL_DERIV=$(shell echo ${DERIV}|grep -i -c ALL)
FLAG=$(shell echo $$(((${IS_FORTRAN} + ${IS_ALL_LANG}) * (${IS_TAPENADE} * ${IS_ALL_DERIV}))))

# Add targets if directory needs to be built
CONDITIONAL_BUILD =
ifneq	(${FLAG}, 0)
	CONDITIONAL_BUILD += cdir cpmod src ${OBJS_MAIN} link
endif

.PHONY: cdir cpmod src build link clean

build: ${CONDITIONAL_BUILD}

link:
	${FORT} ${FFLAGS} $(addprefix ${OBJ_DIR}/, ${OBJS}) ${OBJ_DIR}/eval_deriv_tapenade_1_forward.o -o ../eval_deriv_tapenade_1_forward -J${MOD_DIR} $(addprefix -, ${F90_LIBS})
	${FORT} ${FFLAGS} $(addprefix ${OBJ_DIR}/, ${OBJS}) ${OBJ_DIR}/eval_deriv_tapenade_1_reverse.o -o ../eval_deriv_tapenade_1_forward -J${MOD_DIR} $(addprefix -, ${F90_LIBS})

cdir: 
	mkdir -p ${OBJ_DIR}
	mkdir -p ${MOD_DIR}

cpmod:
	cp ${NETCDF_F90_MOD} ${MOD_DIR}

src:
	make -C ${APP_SRC_DIR}

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

clean:
	make -C ${APP_SRC_DIR} clean
	rm -f ../eval_deriv_tapenade_1_forward.o
	rm -f ../eval_deriv_tapenade_1_reverse.o
