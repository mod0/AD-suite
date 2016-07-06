OBJ_DIR = obj
MOD_DIR = modules
APP_SRC_DIR = src
OBJS_SUPP = supporting_routines.o
OBJS_SIM = numerics.o
OBJS= ${OBJS_SUPP} ${OBJS_SIM}
MAINS=driver1.o 

.PHONY: all src clean

all: src ${MAINS}
	${FORT} ${FFLAGS} $(addprefix ${OBJ_DIR}/, ${OBJS}) ${OBJ_DIR}/driver1.o -o driver1 -J${MOD_DIR}

src:
	make -C ${APP_SRC_DIR}

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.f90
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

clean:
	make -C ${APP_SRC_DIR} clean
	rm -f driver1
	rm -if *.txt
	rm -if *.png

