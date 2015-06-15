#
# OpenAD executable
#
OAD = openad 

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

FFLAGS   = -O0 -g
F90FLAGS =  -O0 -g
FFLAGS2  = -DCOMPLEX
CFLAGS   = -g 
CPPFLAGS = 
LFLAGS   = -lm -static-libgcc

RTSUPP  =	w2f__types OAD_active

RTSUPP_ADJ = OAD_cp OAD_tape OAD_rev 

TARGETS =  airfoil air_lin air_adj #testlinadj

AIRFOIL = airfoil.o input.o routines.o flux_face.o flux_wall.o time_cell.o lift_wall.o 

AIR_LIN = 	time_cell_d.o  time_cell_dx.o\
			flux_face_d.o  flux_face_dx.o \
			flux_wall_d.o flux_wall_dx.o \
			lift_wall_d.o  lift_wall_dx.o\
			input.o print_active.o\
			$(addsuffix .o, $(RTSUPP))\
			air_lin.o

AIR_ADJ = 	time_cell.o time_cell_b.o  time_cell_bx.o \
			flux_face_b.o flux_face_bx.o \
			flux_wall_b.o flux_wall_bx.o \
			lift_wall_b.o lift_wall_bx.o \
			input.o print_active.o\
			$(addsuffix .o, $(RTSUPP))\
			$(addsuffix .o, $(RTSUPP_ADJ))\
			air_adj.o

all: $(TARGETS)

airfoil:	                           ${AIRFOIL}
		@echo Compiling Airfoil;
		${FC} ${FFLAGS} -o airfoil ${AIRFOIL} ${LFLAGS}

air_lin:	                           ${AIR_LIN}
		@echo Compiling Airfoil Linear;
		${FC} ${FFLAGS} -o air_lin ${AIR_LIN} ${LFLAGS}

air_adj:	                           ${AIR_ADJ}
		@echo Compiling Airfoil Adjoint;
		${FC} ${FFLAGS} -o air_adj ${AIR_ADJ} ${LFLAGS}

%.o : %.F
	${FC} ${FFLAGS} -o $@ -c $<

%.o : %.f90
	${FC} ${F90FLAGS} -o $@ -c $< 

%.o : %.f
	${FC} ${F90FLAGS} -o $@ -c $< 

%.o : %.c
	${CC} -o $@ -c $< 

#---------------------------------------------------------------------------
# Linear versions of the routines
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------	
# Time Cell
#---------------------------------------------------------------------------
time_cell_d.o:	time_cell.F
	${CC} -E -C -P -DOPENAD_D time_cell.F > time_cell.f;
	${OAD} -c -m f time_cell.f
	@echo Rename function name in place
	@sed -i'' "s/time_cell/time_cell_d/" time_cell.pre.xb.x2w.w2f.post.f  
	@echo Rename file to time_cell_d.f
	@mv time_cell.pre.xb.x2w.w2f.post.f time_cell_d.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} ${FFLAGS} -c time_cell_d.f

time_cell_dx.o:	time_cell.F
	${CC} -E -C -P -DOPENAD_DX time_cell.F > time_cell.f;
	${OAD} -c -m f time_cell.f
	@echo Rename function name in place
	@sed -i'' "s/time_cell/time_cell_dx/" time_cell.pre.xb.x2w.w2f.post.f  
	@echo Rename file to time_cell_dx.f
	@mv time_cell.pre.xb.x2w.w2f.post.f time_cell_dx.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} ${FFLAGS} -c time_cell_dx.f

#---------------------------------------------------------------------------	
# Flux Face
#---------------------------------------------------------------------------

flux_face_d.o:	flux_face.F
	${CC} -E -C -P -DOPENAD_D flux_face.F > flux_face.f;
	${OAD} -c -m f flux_face.f
	@echo Rename function name in place
	@sed -i'' "s/flux_face/flux_face_d/" flux_face.pre.xb.x2w.w2f.post.f  
	@echo Rename file to flux_face_d.f
	@mv flux_face.pre.xb.x2w.w2f.post.f flux_face_d.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} ${FFLAGS} -c flux_face_d.f

flux_face_dx.o:	flux_face.F
	@${CC} -E -C -P -DOPENAD_DX flux_face.F > flux_face.f;
	@${OAD} -c -m f flux_face.f
	@echo Rename function name in place
	@sed -i'' "s/flux_face/flux_face_dx/" flux_face.pre.xb.x2w.w2f.post.f  
	@echo Rename file to flux_face_dx.f
	@mv flux_face.pre.xb.x2w.w2f.post.f flux_face_dx.f
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} ${FFLAGS} -c flux_face_dx.f

#---------------------------------------------------------------------------	
# Flux Wall
#---------------------------------------------------------------------------
flux_wall_d.o:	flux_wall.F
	${CC} -E -C -P -DOPENAD_D flux_wall.F > flux_wall.f;
	${OAD} -c -m f flux_wall.f
	@echo Rename function name in place
	@sed -i'' "s/flux_wall/flux_wall_d/" flux_wall.pre.xb.x2w.w2f.post.f  
	@echo Rename file to flux_wall_d.f
	@mv flux_wall.pre.xb.x2w.w2f.post.f flux_wall_d.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} ${FFLAGS} -c flux_wall_d.f

flux_wall_dx.o:	flux_wall.F
	@${CC} -E -C -P -DOPENAD_DX flux_wall.F > flux_wall.f;
	@${OAD} -c -m f flux_wall.f
	@echo Rename function name in place
	@sed -i'' "s/flux_wall/flux_wall_dx/" flux_wall.pre.xb.x2w.w2f.post.f  
	@echo Rename file to flux_wall_dx.f
	@mv flux_wall.pre.xb.x2w.w2f.post.f flux_wall_dx.f
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} ${FFLAGS} -c flux_wall_dx.f

#---------------------------------------------------------------------------	
# Lift Wall
#---------------------------------------------------------------------------
lift_wall_d.o:	lift_wall.F
	${CC} -E -C -P -DOPENAD_D lift_wall.F > lift_wall.f;
	${OAD} -c -m f lift_wall.f
	@echo Rename function name in place
	@sed -i'' "s/lift_wall/lift_wall_d/" lift_wall.pre.xb.x2w.w2f.post.f  
	@echo Rename file to lift_wall_d.f
	@mv lift_wall.pre.xb.x2w.w2f.post.f lift_wall_d.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} ${FFLAGS} -c lift_wall_d.f

lift_wall_dx.o:	lift_wall.F
	@${CC} -E -C -P -DOPENAD_DX lift_wall.F > lift_wall.f;
	@${OAD} -c -m f lift_wall.f
	@echo Rename function name in place
	@sed -i'' "s/lift_wall/lift_wall_dx/" lift_wall.pre.xb.x2w.w2f.post.f  
	@echo Rename file to lift_wall_dx.f
	@mv lift_wall.pre.xb.x2w.w2f.post.f lift_wall_dx.f
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} ${FFLAGS} -c lift_wall_dx.f

#---------------------------------------------------------------------------
# Adjoint versions of the routines
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------	
# Time Cell
#---------------------------------------------------------------------------
time_cell_b.o:	time_cell.F
	${CC} -E -C -P -DOPENAD_D time_cell.F > time_cell.f;
	${OAD} -c -m rs time_cell.f
	@echo Rename function name in place
	@sed -i'' "s/time_cell/time_cell_b/" time_cell.pre.xb.x2w.w2f.post.f  
	@echo Rename file to time_cell_b.f
	@mv time_cell.pre.xb.x2w.w2f.post.f time_cell_b.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} -c OAD_cp.f90  
	${FC} -c OAD_tape.f90
	${FC} -c OAD_rev.f90
	${FC} -c ad_inline.f
	${FC} ${FFLAGS} -c time_cell_b.f

time_cell_bx.o:	time_cell.F
	${CC} -E -C -P -DOPENAD_DX time_cell.F > time_cell.f;
	${OAD} -c -m rs time_cell.f
	@echo Rename function name in place
	@sed -i'' "s/time_cell/time_cell_bx/" time_cell.pre.xb.x2w.w2f.post.f  
	@echo Rename file to time_cell_bx.f
	@mv time_cell.pre.xb.x2w.w2f.post.f time_cell_bx.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} -c OAD_cp.f90  
	${FC} -c OAD_tape.f90
	${FC} -c OAD_rev.f90
	${FC} -c ad_inline.f
	${FC} ${FFLAGS} -c time_cell_bx.f

#---------------------------------------------------------------------------	
# Flux Face
#---------------------------------------------------------------------------
flux_face_b.o:	flux_face.F
	${CC} -E -C -P -DOPENAD_D flux_face.F > flux_face.f;
	${OAD} -c -m rs flux_face.f
	@echo Rename function name in place
	@sed -i'' "s/flux_face/flux_face_b/" flux_face.pre.xb.x2w.w2f.post.f  
	@echo Rename file to flux_face_b.f
	@mv flux_face.pre.xb.x2w.w2f.post.f flux_face_b.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} ${FFLAGS} -c flux_face_b.f

flux_face_bx.o:	flux_face.F
	${CC} -E -C -P -DOPENAD_DX flux_face.F > flux_face.f;
	${OAD} -c -m rs flux_face.f
	@echo Rename function name in place
	@sed -i'' "s/flux_face/flux_face_bx/" flux_face.pre.xb.x2w.w2f.post.f  
	@echo Rename file to flux_face_bx.f
	@mv flux_face.pre.xb.x2w.w2f.post.f flux_face_bx.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} ${FFLAGS} -c flux_face_bx.f

#---------------------------------------------------------------------------	
# Flux Wall
#---------------------------------------------------------------------------
flux_wall_b.o:	flux_wall.F
	${CC} -E -C -P -DOPENAD_D flux_wall.F > flux_wall.f;
	${OAD} -c -m rs flux_wall.f
	@echo Rename function name in place
	@sed -i'' "s/flux_wall/flux_wall_b/" flux_wall.pre.xb.x2w.w2f.post.f  
	@echo Rename file to flux_wall_b.f
	@mv flux_wall.pre.xb.x2w.w2f.post.f flux_wall_b.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} -c OAD_cp.f90  
	${FC} -c OAD_tape.f90
	${FC} -c OAD_rev.f90
	${FC} -c ad_inline.f
	${FC} ${FFLAGS} -c flux_wall_b.f

flux_wall_bx.o:	flux_wall.F
	${CC} -E -C -P -DOPENAD_DX flux_wall.F > flux_wall.f;
	${OAD} -c -m rs flux_wall.f
	@echo Rename function name in place
	@sed -i'' "s/flux_wall/flux_wall_bx/" flux_wall.pre.xb.x2w.w2f.post.f  
	@echo Rename file to flux_wall_bx.f
	@mv flux_wall.pre.xb.x2w.w2f.post.f flux_wall_bx.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} -c OAD_cp.f90  
	${FC} -c OAD_tape.f90
	${FC} -c OAD_rev.f90
	${FC} -c ad_inline.f
	${FC} ${FFLAGS} -c flux_wall_bx.f

#---------------------------------------------------------------------------	
# Lift Wall
#---------------------------------------------------------------------------
lift_wall_b.o:	lift_wall.F
	${CC} -E -C -P -DOPENAD_D lift_wall.F > lift_wall.f;
	${OAD} -c -m rs lift_wall.f
	@echo Rename function name in place
	@sed -i'' "s/lift_wall/lift_wall_b/" lift_wall.pre.xb.x2w.w2f.post.f  
	@echo Rename file to lift_wall_b.f
	@mv lift_wall.pre.xb.x2w.w2f.post.f lift_wall_b.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} -c OAD_cp.f90  
	${FC} -c OAD_tape.f90
	${FC} -c OAD_rev.f90
	${FC} -c ad_inline.f
	${FC} ${FFLAGS} -c lift_wall_b.f

lift_wall_bx.o:	lift_wall.F
	${CC} -E -C -P -DOPENAD_DX lift_wall.F > lift_wall.f;
	${OAD} -c -m rs lift_wall.f
	@echo Rename function name in place
	@sed -i'' "s/lift_wall/lift_wall_bx/" lift_wall.pre.xb.x2w.w2f.post.f  
	@echo Rename file to lift_wall_bx.f
	@mv lift_wall.pre.xb.x2w.w2f.post.f lift_wall_bx.f 
	${CC} -c iaddr.c
	${FC} -c w2f__types.f90
	${FC} -c OAD_active.f90
	${FC} -c OAD_cp.f90  
	${FC} -c OAD_tape.f90
	${FC} -c OAD_rev.f90
	${FC} -c ad_inline.f
	${FC} ${FFLAGS} -c lift_wall_bx.f

clean:
	rm -f *.o *.f *.B *.xaif *.s *.log~ *.f90 *.mod-whirl *.mod airfoil air_lin air_adj



