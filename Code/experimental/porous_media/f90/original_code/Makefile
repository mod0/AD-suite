include ../../make.config

OBJ_DIR = obj
MOD_DIR = modules
APP_SRC_DIR = src
OBJS_SUPP = supporting_routines.o
OBJS_SIM = numerics.o
OBJS = ${OBJS_SUPP} ${OBJS_SIM}
OBJS_MAIN = eval_original_code.o  

.PHONY: cdir build src clean

build: cdir src ${OBJS_MAIN}
	${FORT} ${FFLAGS} $(addprefix ${OBJ_DIR}/, ${OBJS}) ${OBJ_DIR}/eval_original_code.o -o ../eval_original_code -J${MOD_DIR}

cdir: 
	mkdir -p ${OBJ_DIR}
	mkdir -p ${MOD_DIR}

src:
	make -C ${APP_SRC_DIR}

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

clean:
	make -C ${APP_SRC_DIR} clean
	rm -f ../eval_original_code


