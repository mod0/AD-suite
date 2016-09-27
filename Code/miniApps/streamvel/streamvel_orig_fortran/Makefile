ifndef F90C
F90C=gfortran
endif
ifndef CC
CC=gcc
endif
driver: stream_vel_variables.o conj_grad.o stream_vel_timedep.o driver.o 
	${F90C} -g -o $@ $^

%.o : %.f90
	${F90C} ${F90FLAGS} -g  -O -o $@ -c $< 

%.o : %.f
	${F90C} ${F90FLAGS} -g -O -o $@ -c $< 

%.o : %.c
	${CC} -g -O -o $@ -c $< 

clean: 
	rm -f  *.o *.mod* driver *~ 

.PHONY: clean 

