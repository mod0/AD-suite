ifndef F90C
F90C=gfortran
endif
ifndef CC
CC=gcc
endif
RTSUPP=w2f__types OAD_active OAD_cp OAD_tape OAD_rev ad_inline
driverADMjoint:   $(addsuffix .o, $(RTSUPP) iaddr) ad_inline.o stream_vel_variables_passive.o conj_grad.o conj_grad_ad.o numCore.pre.xb.x2w.w2f.post.o driver.o 
	${F90C} -g -o $@ $^
numCore.pre.xb.x2w.w2f.post.f90 $(addsuffix .f90, $(RTSUPP)) iaddr.c : toolChain 
toolChain : numCore.f90
#	echo "Faking openadcall"
	openad -c -m rj $<
#	svn revert OAD_rev.f90

        
numCore.f90: stream_vel_variables.f90 conj_gradStubs.f90 stream_vel_timedep.f90
	cat $^ | sed 's/use conj_grad/use conj_gradStub/' > $@

stream_vel_variables_passive.f90: stream_vel_variables.f90
	cat $< | sed 's/stream_vel_variables/stream_vel_variables_passive/' > $@

ad_inline.f:toolChain

%.o : %.f90
	${F90C} ${F90FLAGS} -g  -O -o $@ -c $< 

%.o : %.f
	${F90C} ${F90FLAGS} -g -O -o $@ -c $< 

%.o : %.c
	${CC} -g -O0 -o $@ -c $< 

#head.f90 lu.f90: 
#	ln -s ../$@ ./

clean: 
	rm -f  numCore.* *.o *.mod* driverADMjoint *~ 
	rm -f ad_template* ad_inline.f OAD_* w2f__*  iaddr* head.f90 lu.f90
#	rm -f  *.o *.mod* driverADMjoint *~ 

.PHONY: clean toolChain
