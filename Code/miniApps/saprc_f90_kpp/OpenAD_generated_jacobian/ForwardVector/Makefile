#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User: Set here the F90 compiler and options
#       Pedefined compilers: INTEL, PGF, HPUX, LAHEY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#OpenAD executable
OAD = openad

#COMPILER = G95
#COMPILER = LAHEY
#COMPILER = INTEL
#COMPILER = PGF
#COMPILER = HPUX
COMPILER = GFORTRAN

FC_G95     = g95
FOPT_G95   = -cpp -O -pg -fbounds-check -fimplicit-none  -Wall -ftrace=full

CC_G95 = gcc
COPT_G95 = -O0 -g

FC_LAHEY   = lf95
# More aggressive for production runs:
#FOPT_LAHEY = -Cpp --pca -O
# More checking for debugging:
FOPT_LAHEY = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap --trace --chkglobal

FC_INTEL   = ifort
# More aggressive for production runs:
#FOPT_INTEL = -cpp -O -fp-model precise -pc80 -prec_div
# More checking for debugging:
FOPT_INTEL = -cpp -O0 -fp-model strict -implicitnone -ftrapuv \
              -debug all -warn all # -check all

CC_INTEL = icc
COPT_INTEL = -O0 -g

FC_PGF     = pgf90
# More aggressive for production runs:
FOPT_PGF   = -Mpreprocess -O -fast -pc 80 -Kieee
# More checking for debugging:
#FOPT_PGF   = -Mpreprocess -O0 -Mbounds -Mchkfpstk -Mchkptr -Mchkstk \
#             -Ktrap=fp -pc 80 -Kieee

FC_HPUX    = f90
FOPT_HPUX  = -O -u +Oall +check=on

FC_GFORTRAN     = gfortran
FOPT_GFORTRAN   = -cpp -O

CC_GFORTRAN = gcc
COPT_GFORTRAN = -O0 -g

# define FULL_ALGEBRA for non-sparse integration
FC   = $(FC_$(COMPILER))
FOPT = $(FOPT_$(COMPILER)) # -DFULL_ALGEBRA

CC = $(CC_$(COMPILER))
COPT = $(COPT_$(COMPILER))

LIBS =
#LIBS = -llapack -lblas

# Command to create Matlab mex gateway routines
# Note: use $(FC) as the mex Fortran compiler
MEX  = mex

GENSRC =   saprc2006_Global.o #\
	 saprc2006_Parameters.o       \
	 saprc2006_Precision.o

GENOBJ =  saprc2006_Global.o #\
	 saprc2006_Parameters.o       \
	 saprc2006_Precision.o

FUNSRC = saprc2006_Function.f90
FUNOBJ = saprc2006_Function.o

JACSRC = saprc2006_JacobianSP.f90  # saprc2006_Jacobian.f90
JACOBJ = saprc2006_JacobianSP.o    # saprc2006_Jacobian.o

HESSRC = saprc2006_HessianSP.f90   saprc2006_Hessian.f90
HESOBJ = saprc2006_HessianSP.o     saprc2006_Hessian.o

STMSRC = saprc2006_StoichiomSP.f90 saprc2006_Stoichiom.f90
STMOBJ = saprc2006_StoichiomSP.o   saprc2006_Stoichiom.o

UTLSRC = saprc2006_Rates.f90 saprc2006_Util.f90 saprc2006_Monitor.f90
UTLOBJ = saprc2006_Rates.o   saprc2006_Util.o   saprc2006_Monitor.o

LASRC  = saprc2006_LinearAlgebra.f90
LAOBJ  = saprc2006_LinearAlgebra.o

STOCHSRC = saprc2006_Stochastic.f90
STOCHOBJ = saprc2006_Stochastic.o

OADSRC = w2f__types.f90 OAD_active.f90 iaddr.c saprc2006_Func_OpenAD.pre.xb.x2w.w2f.post.f90
OADOBJ = w2f__types.o OAD_active.o iaddr.o saprc2006_Func_OpenAD.pre.xb.x2w.w2f.post.o

MAINSRC = saprc2006_Main.f90   saprc2006_Initialize.f90   saprc2006_Integrator.f90 saprc2006_Model.f90
MAINOBJ = saprc2006_Main.o     saprc2006_Initialize.o     saprc2006_Integrator.o   saprc2006_Model.o

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User: modify the line below to include only the
#       objects needed by your application
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ALLOBJ = $(OADOBJ) $(GENOBJ) $(JACOBJ) $(HESOBJ) $(STMOBJ) \
	 $(UTLOBJ) $(LAOBJ)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User: modify the line below to include only the
#       executables needed by your application
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all:    exe

$(OADSRC): saprc2006_Precision.f90 saprc2006_Parameters.f90 saprc2006_Function.f90
	cat saprc2006_Precision.f90 saprc2006_Parameters.f90 saprc2006_Function.f90 > saprc2006_Func_OpenAD.f90
	sed -i'' 's/SUBROUTINE Fun/SUBROUTINE Fun_OAD/g' saprc2006_Func_OpenAD.f90
	$(OAD) -c -m fv saprc2006_Func_OpenAD.f90

exe:	$(ALLOBJ) $(MAINOBJ)
	$(FC) $(FOPT) $(ALLOBJ) $(MAINOBJ) $(LIBS) -o saprc2006.exe

stochastic:$(ALLOBJ) $(STOCHOBJ) $(MAINOBJ)
	$(FC) $(FOPT) $(ALLOBJ) $(STOCHOBJ) $(MAINOBJ) $(LIBS) \
	-o saprc2006_stochastic.exe

mex:    $(ALLOBJ)
	$(MEX) FC#$(FC) -fortran -O saprc2006_mex_Fun.f90     $(ALLOBJ)
	$(MEX) FC#$(FC) -fortran -O saprc2006_mex_Jac_SP.f90  $(ALLOBJ)
	$(MEX) FC#$(FC) -fortran -O saprc2006_mex_Hessian.f90 $(ALLOBJ)

clean:
	rm -f saprc2006*.o saprc2006*.mod \
	saprc2006*.dat saprc2006.exe saprc2006*.mexglx \
	saprc2006.map OAD_active.f90 w2f__types.f90 iaddr.c \
	saprc2006_Func_OpenAD.pre.* oad_active.mod \
	oad_intrinsics.mod w2f__types.mod \
	iaddr.o OAD_active.o w2f__types.o \
	*.mod-whirl openad.*.log~ \
	saprc2006_Func_OpenAD.f90

distclean:
#	rm -f saprc2006*.o saprc2006*.mod \
	saprc2006*.dat saprc2006.exe saprc2006.map \
	saprc2006*.f90 saprc2006_*.mexglx

saprc2006_Func_OpenAD.pre.xb.x2w.w2f.post.o: saprc2006_Func_OpenAD.pre.xb.x2w.w2f.post.f90
	$(FC) $(FOPT) -c $<

OAD_active.o: OAD_active.f90
	$(FC) $(FOPT) -c $<

w2f__types.o: w2f__types.f90
	$(FC) $(FOPT) -c $<

iaddr.o: iaddr.c
	$(CC) $(COPT) -c $<

saprc2006_Precision.o: saprc2006_Precision.f90
	$(FC) $(FOPT) -c $<

saprc2006_Parameters.o: saprc2006_Parameters.f90 \
	            saprc2006_Precision.o
	$(FC) $(FOPT) -c $<

saprc2006_Monitor.o: saprc2006_Monitor.f90 \
	             saprc2006_Precision.o
	$(FC) $(FOPT) -c $<

saprc2006_Global.o: saprc2006_Global.f90 \
	            saprc2006_Parameters.o saprc2006_Precision.o
	$(FC) $(FOPT) -c $<

saprc2006_Initialize.o: saprc2006_Initialize.f90  $(GENOBJ)
	$(FC) $(FOPT) -c $<

saprc2006_Function.o: saprc2006_Function.f90  $(GENOBJ)
	$(FC) $(FOPT) -c $<

saprc2006_Stochastic.o: saprc2006_Stochastic.f90  $(GENOBJ)
	$(FC) $(FOPT) -c $<

saprc2006_JacobianSP.o: saprc2006_JacobianSP.f90 $(GENOBJ)
	$(FC) $(FOPT) -c $<

saprc2006_Jacobian.o: saprc2006_Jacobian.f90  $(GENOBJ) saprc2006_JacobianSP.o
	$(FC) $(FOPT) -c $<

saprc2006_LinearAlgebra.o: saprc2006_LinearAlgebra.f90 $(GENOBJ) saprc2006_JacobianSP.o
	$(FC) $(FOPT) -c $<

saprc2006_Rates.o: saprc2006_Rates.f90  $(GENOBJ)
	$(FC) $(FOPT) -c $<

saprc2006_HessianSP.o: saprc2006_HessianSP.f90  $(GENOBJ)
	$(FC) $(FOPT) -c $<

saprc2006_Hessian.o:  saprc2006_Hessian.f90 $(GENOBJ) saprc2006_HessianSP.o
	$(FC) $(FOPT) -c $<

saprc2006_StoichiomSP.o: saprc2006_StoichiomSP.f90 $(GENOBJ)
	$(FC) $(FOPT) -c $<

saprc2006_Stoichiom.o: saprc2006_Stoichiom.f90  $(GENOBJ) saprc2006_StoichiomSP.o
	$(FC) $(FOPT) -c $<

saprc2006_Util.o: saprc2006_Util.f90  $(GENOBJ) saprc2006_Monitor.o
	$(FC) $(FOPT) -c $<

saprc2006_Main.o: saprc2006_Main.f90  $(ALLOBJ) saprc2006_Initialize.o saprc2006_Model.o saprc2006_Integrator.o
	$(FC) $(FOPT) -c $<

saprc2006_Model.o: saprc2006_Model.f90  $(ALLOBJ) saprc2006_Integrator.o
	$(FC) $(FOPT) -c $<

saprc2006_Integrator.o: saprc2006_Integrator.f90  $(ALLOBJ)
	$(FC) $(FOPT) -c $<
