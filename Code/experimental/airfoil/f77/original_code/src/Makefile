include ../../../make.config


OBJ_DIR = ../obj
MOD_DIR = ../modules
OBJS_SUPP = io.o 
OBJS_SIM = routines.o
OBJS_EXT = 
OBJS= ${OBJS_SUPP} ${OBJS_SIM} ${OBJS_EXT}

.PHONY: all routines.F io.F clean

all:  routines.F io.F ${OBJS} ${MAINS}

%.o: %.F
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

%.o: %.f
	${FORT} -cpp ${FFLAGS} -c $< -o ${OBJ_DIR}/$@ -J${MOD_DIR}

routines.F: 
	@echo "Including routines.F"

input.F:
	@echo "Including input.F"

clean:
	rm -f ${OBJ_DIR}/*
	rm -f ${MOD_DIR}/*  
