include ../../make.config

OBJ_DIR = obj
MOD_DIR = modules
APP_SRC_DIR = src
OBJS_SUPP = supporting_routines.o 
OBJS_AD_SUPP = w2f__types.o OAD_active.o OAD_tape.o OAD_rev.o OAD_cp.o iaddr.o 
OBJS_SIM_TL = numerics_f.o

OBJS_EXT_TL = 

OBJS_MAIN_TL = eval_deriv_openad_1_forward.o

OBJS_TL = ${OBJS_SUPP} ${OBJS_AD_SUPP} ${OBJS_SIM_TL} ${OBJS_EXT_TL}
OBJS_ADJ = ${OBJS_SUPP} ${OBJS_AD_SUPP} ${OBJS_SIM_ADJ} ${OBJS_EXT_ADJ}

# Determine if this directory needs to be built
# Refer: http://stackoverflow.com/a/28790801/1220493
IS_FORTRAN=$(shell echo ${LANG}|grep -i -c F90)
IS_TAPENADE=$(shell echo ${DERIV}|grep -i -c OPENAD)
IS_ALL_LANG=$(shell echo ${LANG}|grep -i -c ALL)
IS_ALL_DERIV=$(shell echo ${DERIV}|grep -i -c ALL)
FLAG=$(shell echo $$(((${IS_FORTRAN} + ${IS_ALL_LANG}) * (${IS_TAPENADE} + ${IS_ALL_DERIV}))))

# Add targets if directory needs to be built
CONDITIONAL_BUILD =
ifneq	(${FLAG}, 0)
	CONDITIONAL_BUILD += cdir cpmod src ${OBJS_MAIN_TL} ${OBJS_MAIN_ADJ} link
endif

.PHONY: cdir cpmod src build link clean

build: ${CONDITIONAL_BUILD}


link:
		${FORT} ${FFLAGS} $(addprefix ${OBJ_DIR}/, ${OBJS_TL}) ${OBJ_DIR}/${OBJS_MAIN_TL} -o ../eval_deriv_openad_1_forward -J${MOD_DIR} $(addprefix -, ${F90_LIBS})

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
	rm -f ../eval_deriv_openad_1_forward

