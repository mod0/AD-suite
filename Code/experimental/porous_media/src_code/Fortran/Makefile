OBJ_DIR = app_obj
MOD_DIR = app_mod
APP_SRC_DIR = app_src
OBJS_SUPP = supporting_routines.o
OBJS_SIM = numerics.o
OBJS= ${OBJS_SUPP} ${OBJS_SIM}
MAINS=driver1.o 

.PHONY: all app_src clean

all: app_src ${MAINS}
	${FORT} ${FFLAGS} $(addprefix ${OBJ_DIR}/, ${OBJS}) driver1.o -o driver1 -J${MOD_DIR}

app_src:
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

