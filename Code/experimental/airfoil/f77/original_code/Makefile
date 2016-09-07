include ../../make.config

OBJ_DIR = obj
MOD_DIR = modules
APP_SRC_DIR = src
OBJS_SUPP = io.o 
OBJS_SIM = routines.o
OBJS_EXT = 
OBJS = ${OBJS_SUPP} ${OBJS_SIM} ${OBJS_EXT}
OBJS_MAIN = eval_original_code.o  

.PHONY: cdir src build link clean

build: cdir src ${OBJS_MAIN} link

link:
	${FORT} ${FFLAGS} $(addprefix ${OBJ_DIR}/, ${OBJS}) ${OBJ_DIR}/eval_original_code.o -o ../eval_original_code -J${MOD_DIR} $(addprefix -, ${F77_LIBS})

cdir: 
	mkdir -p ${OBJ_DIR}
	mkdir -p ${MOD_DIR}

src:
	make -C ${APP_SRC_DIR}

%.o: %.F
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

clean:
	make -C ${APP_SRC_DIR} clean
	rm -f ../eval_original_code
