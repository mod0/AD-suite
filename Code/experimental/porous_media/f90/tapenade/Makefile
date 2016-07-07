OBJS= gnufor2.o utils.o mgmres.o adBuffer.o adStack.o
DIFFS= numerics_d.o numerics_b.o
MAINS= main_solver_tuning_tanglin.o main_deriv_check_tanglin.o \
       main_tanglin.o main_solver_tuning_adj.o main_deriv_check_adj.o\
       main_oil_vs_mu_sigma.o
LIBS= -lgfortran -lm

all:$(OBJS) $(DIFFS) $(MAINS)
	${FORT} ${FFLAGS} ${OBJS} numerics_d.o main_tanglin.o  ${LIBS} -o porous_tanglin
	${FORT} ${FFLAGS} ${OBJS} numerics_d.o main_solver_tuning_tanglin.o ${LIBS} -o porous_solver_tuning_tanglin
	${FORT} ${FFLAGS} ${OBJS} numerics_d.o main_deriv_check_tanglin.o ${LIBS} -o porous_deriv_check_tanglin
	${FORT} ${FFLAGS} ${OBJS} numerics_b.o main_solver_tuning_adj.o ${LIBS} -o porous_solver_tuning_adj
	${FORT} ${FFLAGS} ${OBJS} numerics_b.o main_deriv_check_adj.o ${LIBS} -o porous_deriv_check_adj
	${FORT} ${FFLAGS} ${OBJS} numerics_d.o main_oil_vs_mu_sigma.o ${LIBS} -o porous_mu_sigma_oil

%.o : %.f90
	${FORT} ${FFLAGS} -c $< -o $@

%.o : %.f
	${FORT} ${FFLAGS} -c $< -o $@

clean:
	rm -f *.o
	rm -f *.mod-whirl
	rm -f *.mod
	rm -f porous_tanglin
	rm -f porous_solver_tuning_tanglin
	rm -f porous_deriv_check_tanglin
	rm -f porous_solver_tuning_adj
	rm -f porous_deriv_check_adj
	rm -f porous_mu_sigma_oil
	rm -i *~
	rm -i *.msg
	rm -i *test.txt
	rm -i test*.txt
