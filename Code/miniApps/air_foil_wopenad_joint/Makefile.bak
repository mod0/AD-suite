
#
# TAPENADE executable
#

#TPN  = ~giles/bin/tapenade

TPN  = ${TAPENADE_HOME}/bin/tapenade

#
# Intel 8.0 flags
#
 
FC = ifort 
CC = icc

FFLAGS   = -O3 -ipo -static
FFLAGS2  = -DCOMPLEX -O3 -ipo -static
CFLAGS   = -O3 -ipo -static
CPPFLAGS = 
LFLAGS   = -Vaxlib -lm -static-libcxa

#
# gfortran flags
#
 
FC = gfortran
CC = gcc

FFLAGS   = 
FFLAGS2  = -DCOMPLEX
CFLAGS   =
CPPFLAGS = 
LFLAGS   = -lm -static-libgcc

TARGETS =  airfoil air_lin air_adj testlinadj

#
# transformation rules
#

.F.o:     ; /bin/rm -f $*.o; ${FC} ${FFLAGS} -c $*.F ; /bin/rm -f $*.cpp
.c.o:     ; /bin/rm -f $*.o; ${CC} ${CFLAGS} -c $*.c ; /bin/rm -f $*.cpp

#
# program dependencies
#

AIRFOIL = airfoil.o input.o routines.o

AIR_LIN = air_lin.o input.o \
          flux_face_d.o  flux_wall_d.o  time_cell_d.o  lift_wall_d.o \
          flux_face_dx.o flux_wall_dx.o time_cell_dx.o lift_wall_dx.o

AIR_ADJ = air_adj.o input.o routines.o adStack.o \
          flux_face_b.o  flux_wall_b.o  time_cell_b.o  lift_wall_b.o \
          flux_face_bx.o flux_wall_bx.o time_cell_bx.o lift_wall_bx.o

TESTLINADJ = testlinadj.o routines.o Croutines.o adStack.o \
          flux_face_d.o  flux_wall_d.o  time_cell_d.o  lift_wall_d.o \
          flux_face_b.o  flux_wall_b.o  time_cell_b.o  lift_wall_b.o \
          flux_face_dx.o flux_wall_dx.o time_cell_dx.o lift_wall_dx.o \
          flux_face_bx.o flux_wall_bx.o time_cell_bx.o lift_wall_bx.o

#
# programs
#

ALL:		$(TARGETS)

airfoil:	                           ${AIRFOIL}
		${FC} ${FFLAGS} -o airfoil ${AIRFOIL} ${LFLAGS}

air_lin:	                           ${AIR_LIN}
		${FC} ${FFLAGS} -o air_lin ${AIR_LIN} ${LFLAGS}

air_adj:	                           ${AIR_ADJ}
		${FC} ${FFLAGS} -o air_adj ${AIR_ADJ} ${LFLAGS}

testlinadj:	                              ${TESTLINADJ}
		${FC} ${FFLAGS} -o testlinadj ${TESTLINADJ} ${LFLAGS}

#
# Complex version of nonlinear routines
#

Croutines.o:	                                   routines.F
		${FC} ${FFLAGS2} -c -o Croutines.o routines.F

#
# Linear versions of nonlinear routines
#

time_cell_d.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -forward              \
		       -head       time_cell \
		       -output     time_cell \
		       -vars       "q"       \
		       -outvars    "q adt"   \
		        routines.f;
		${FC} ${FFLAGS} -c time_cell_d.f;
	/bin/rm routines.f time_cell_d.f *.msg

time_cell_dx.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -forward                          \
		       -head         time_cell           \
		       -output       time_cell           \
		       -vars         "x1 x2 x3 x4 q"     \
		       -outvars      "x1 x2 x3 x4 q adt" \
		       -difffuncname "_dx"               \
		        routines.f;
		${FC} ${FFLAGS} -c time_cell_dx.f;
		/bin/rm routines.f time_cell_dx.f *.msg

flux_face_d.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -forward                                \
		       -head       flux_face                   \
		       -output     flux_face                   \
		       -vars       "q1 q2 adt1 adt2 res1 res2" \
		       -outvars    "q1 q2 adt1 adt2 res1 res2" \
		        routines.f;
		${FC} ${FFLAGS} -c flux_face_d.f;
		/bin/rm routines.f flux_face_d.f *.msg

flux_face_dx.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -forward                                        \
		       -head         flux_face                         \
		       -output       flux_face                         \
		       -vars         "x1 x2 q1 q2 adt1 adt2 res1 res2" \
		       -outvars      "x1 x2 q1 q2 adt1 adt2 res1 res2" \
		       -difffuncname "_dx"                             \
		        routines.f;
		${FC} ${FFLAGS} -c flux_face_dx.f;
		/bin/rm routines.f flux_face_dx.f *.msg

flux_wall_d.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -forward              \
		       -head       flux_wall \
		       -output     flux_wall \
		       -vars       "q res"   \
		       -outvars    "q res"   \
		        routines.f;
		${FC} ${FFLAGS} -c flux_wall_d.f;
		/bin/rm routines.f flux_wall_d.f *.msg

flux_wall_dx.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -forward                    \
		       -vars         "x1 x2 q res" \
		       -outvars      "x1 x2 q res" \
		       -head         flux_wall     \
		       -output       flux_wall     \
		       -difffuncname "_dx"         \
		        routines.f;
		${FC} ${FFLAGS} -c flux_wall_dx.f;
		/bin/rm routines.f flux_wall_dx.f *.msg

lift_wall_d.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -forward              \
		       -head       lift_wall \
		       -output     lift_wall \
		       -vars       "q lift"  \
		       -outvars    "q lift"  \
		        routines.f;
		${FC} ${FFLAGS} -c lift_wall_d.f;
		/bin/rm routines.f lift_wall_d.f *.msg

lift_wall_dx.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -forward                     \
		       -head         lift_wall      \
		       -output       lift_wall      \
		       -vars         "x1 x2 q lift" \
		       -outvars      "x1 x2 q lift" \
		       -difffuncname "_dx"          \
		        routines.f;
		${FC} ${FFLAGS} -c lift_wall_dx.f;
		/bin/rm routines.f lift_wall_dx.f *.msg

#
# Adjoint versions of nonlinear routines
#

time_cell_b.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -backward             \
		       -head       time_cell \
		       -output     time_cell \
		       -vars       "q"       \
		       -outvars    "q adt"   \
		        routines.f;
		${FC} ${FFLAGS} -c time_cell_b.f;
		/bin/rm routines.f time_cell_b.f *.msg

time_cell_bx.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -backward                         \
		       -head         time_cell           \
		       -output       time_cell           \
		       -vars         "x1 x2 x3 x4 q"     \
		       -outvars      "x1 x2 x3 x4 q adt" \
		       -difffuncname "_bx"               \
		        routines.f;
		${FC} ${FFLAGS} -c time_cell_bx.f;
		/bin/rm routines.f time_cell_bx.f *.msg

flux_face_b.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -backward                               \
		       -head       flux_face                   \
		       -output     flux_face                   \
		       -vars       "q1 q2 adt1 adt2 res1 res2" \
		       -outvars    "q1 q2 adt1 adt2 res1 res2" \
		        routines.f;
		${FC} ${FFLAGS} -c flux_face_b.f;
		/bin/rm routines.f flux_face_b.f *.msg

flux_face_bx.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -backward                                       \
		       -head         flux_face                         \
		       -output       flux_face                         \
		       -vars         "x1 x2 q1 q2 adt1 adt2 res1 res2" \
		       -outvars      "x1 x2 q1 q2 adt1 adt2 res1 res2" \
		       -difffuncname "_bx"                             \
		        routines.f;
		${FC} ${FFLAGS} -c flux_face_bx.f;
		/bin/rm routines.f flux_face_bx.f *.msg

flux_wall_b.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -backward             \
		       -head       flux_wall \
		       -output     flux_wall \
		       -vars       "q res"   \
		       -outvars    "q res"   \
		        routines.f;
		${FC} ${FFLAGS} -c flux_wall_b.f;
		/bin/rm routines.f flux_wall_b.f *.msg

flux_wall_bx.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -backward                   \
		       -head         flux_wall     \
		       -output       flux_wall     \
		       -vars         "x1 x2 q res" \
		       -outvars      "x1 x2 q res" \
		       -difffuncname "_bx"         \
		        routines.f;
		${FC} ${FFLAGS} -c flux_wall_bx.f;
		/bin/rm routines.f flux_wall_bx.f *.msg

lift_wall_b.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -backward             \
		       -head       lift_wall \
		       -output     lift_wall \
		       -vars       "q lift"  \
		       -outvars    "q lift"  \
		        routines.f;
		${FC} ${FFLAGS} -c lift_wall_b.f;
		/bin/rm routines.f lift_wall_b.f *.msg

lift_wall_bx.o:	routines.F
		${CC} -E -C -P routines.F > routines.f;
		${TPN} -backward                    \
		       -head         lift_wall      \
		       -output       lift_wall      \
		       -vars         "x1 x2 q lift" \
		       -outvars      "x1 x2 q lift" \
		       -difffuncname "_bx"          \
		        routines.f;
		${FC} ${FFLAGS} -c lift_wall_bx.f;
		/bin/rm routines.f lift_wall_bx.f *.msg

#
# clean things up
#

clean:	
	/bin/rm -f *.o *.cpp *.f *.msg

realclean:	
	/bin/rm -f *.o *.cpp *.f *.msg $(TARGETS)

#
# create tar file
#

tar:	
	tar cvf airfoil.tar README *.F *.c *.inc Makefile *.dat

