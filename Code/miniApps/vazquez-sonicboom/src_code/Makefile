SRCDIR=sources
OBJDIR=generated
EXEDIR=data

ADSOURCENAMES=conddirflux.f  fluroe.f  gradnod.f \
psiroe.f  transpiration.f  vcurvm.f gradfb.f normcoq.f calcnormpeau.f
INCLUDENAMES=Param3D.h  Paramopt3D.h param.h
INCLDIFFNAMES=DiffCall_flow_d.h DiffCall_flow_b.h
FSOURCENAMES=OptDes3D.f testAD.f second.f don3d.f calprestuyere.f mail3d.f \
 cmvvno.f cmvfac.f calcfac.f recherche.f mail3d2d.f mail2d.f \
 jacobi.f jacobiroe.f excsub.f inversion.f \
 normcoq.f transpiration.f vcurvm.f gradfb.f \
 calcnormpeau.f conddirflux.f fluroe.f gradnod.f rept3d.f \
 tilt.f seg3d.f calfro.f caldeg.f initcontrol.f inim3d.f \
 research.f cdmat.f diagajout.f trins.f imptranspiration.f \
 condbords.f matvl.f resu3d.f aerof1.f clippin.f improe.f \
 dfsomme.f dfdecentre.f dfexact.f rotation.f trirap.f \
 insert.f partit.f ecriture.f cdlns.f ecritvigie.f \
 psiroe.f fluxdt.f flux.f zzoutsur.f calcaires.f suctrl.f \
 rhoopt.f smlesgr.f lissage.f fullvcycle.f formmagique.f \
 gradient.f etatadjointO2.f testadjoint.f testdercout.f \
 portance.f adjointe.f fnellecout_gp.f fcnelj.f vanleer.f \
 gtolocal.f dcoutdgamm.f dcoutdw_gp.f testetat.f airpeauinit.f \
 smlesdi.f calcvno.f vnonorm.f physcons.f interpol.f fdecfex.f \
 verifcontr.f zonag.f airagnormag.f showdi.f clcible.f

FSOURCES=$(FSOURCENAMES:%=$(SRCDIR)/%)
FOBJECTS=$(FSOURCENAMES:%.f=$(OBJDIR)/%.o)
ADSOURCES=$(ADSOURCENAMES:%=$(SRCDIR)/%)
INCLUDES=$(INCLUDENAMES:%=$(SRCDIR)/%)
INCLDIFFS=$(INCLDIFFNAMES:%=$(SRCDIR)/%)

AD=$(DIFFCOMMAND)
#FC=pgf77 -g # -Mprof=func # -pg
FC=f77 -O
#CC=pgcc -D_BYTESWAP_H # -Mprof=func # -pg
CC=cc


../generated :
	mkdir -p ../generated

all : $(OBJDIR)/Paramopt3D.h $(EXEDIR)/executable
	cd $(EXEDIR); ./executable > ../Results ; cd .. ; cat Results

$(EXEDIR)/executable : $(FOBJECTS) $(OBJDIR)/adBuffer.o $(OBJDIR)/adStack.o $(OBJDIR)/allDiff_d.o $(OBJDIR)/allDiff_b.o
	$(FC) $^ -o $@

$(FOBJECTS):$(OBJDIR)/%.o:$(SRCDIR)/%.f $(INCLUDES) $(INCLDIFFS)
	$(FC) -c $< -o $@

$(OBJDIR)/adBuffer.o : ../../ADFirstAidKit/adBuffer.f
	$(FC) -c $< -o $@

$(OBJDIR)/adStack.o : ../../ADFirstAidKit/adStack.c
	$(CC) -c $< -o $@

$(OBJDIR)/allDiff_d.o : $(OBJDIR)/allDiff_d.f
	$(FC) -c $< -o $@

$(OBJDIR)/allDiff_b.o : $(OBJDIR)/allDiff_b.f
	$(FC) -c $< -o $@

$(OBJDIR)/allDiff_d.f : $(ADSOURCES) $(INCLUDES)
	time $(AD) -d -msglevel 30 -debug 30 -head psiroe -vars "ua" -outvars "ce" -O generated/ -o allDiff $(ADSOURCES) 2>&1 > $(OBJDIR)/traceADtgt

$(OBJDIR)/allDiff_b.f : $(ADSOURCES) $(INCLUDES)
	time $(AD) -b -msglevel 30 -debug 30 -head psiroe -vars "ua" -outvars "ce" -O generated/ -o allDiff $(ADSOURCES) 2>&1 > $(OBJDIR)/traceADadj

$(OBJDIR)/Paramopt3D.h : $(INCLUDES)
	cp $(SRCDIR)/*.h $(OBJDIR)

diffresults:
	diff Results refs; echo ""

diff:
	diff $(OBJDIR)/allDiff_d.msg refs ; diff $(OBJDIR)/allDiff_b.msg refs ; \
	diff $(OBJDIR)/allDiff_d.f refs ; diff $(OBJDIR)/allDiff_b.f refs ; diff Results refs ; echo ""

accept:
	cp Results refs
	cp $(OBJDIR)/allDiff_*.f refs
	cp $(OBJDIR)/allDiff_*.msg refs

clean:
	/bin/rm -f *~ .*~ */*~ */.*~ $(OBJDIR)/* $(EXEDIR)/executable Results
	/bin/rm -f $(EXEDIR)/fort.59 $(EXEDIR)/fort.21 $(EXEDIR)/volu.res
	/bin/rm -f $(EXEDIR)/surf.res $(EXEDIR)/conv.data
	/bin/rm -f $(EXEDIR)/fort.88 $(EXEDIR)/fort.20

# PROFILE
profile: $(EXEDIR)/executable-profile
	cd $(EXEDIR); ./executable-profile > ../Results; cd ..; cat Results

$(EXEDIR)/executable-profile : $(OBJDIR)/allDiff_profile_p.o $(OBJDIR)/profile.o $(OBJDIR)/adBuffer.o $(OBJDIR)/adStack.o # $(OBJDIR)/dpStack.o  $(OBJDIR)/dpTest.o
	$(FC) $(FFLAGS) $^ -o $@

$(OBJDIR)/profile.o : ../../ADFirstAidKit/profile.c
	gcc -c $< -o $@

$(OBJDIR)/allDiff_profile_p.o : $(OBJDIR)/allDiff_profile_p.f
	$(FC) -c $< -o $@

$(OBJDIR)/allDiff_profile_b.f: $(ADSOURCES) $(INCLUDES)
	time $(AD) -profile -b -msglevel 30 -head psiroe -vars "ua" -outvars "ce" -O generated/ -o allDiff_profile $(ADSOURCES) 2>&1 > $(OBJDIR)/traceADadjprofile

$(OBJDIR)/allDiff_profile_d.f : $(ADSOURCES) $(INCLUDES)
	time $(AD) -profile -d -msglevel 30 -head psiroe -vars "ua" -outvars "ce" -O generated/ -o allDiff_profile $(ADSOURCES) 2>&1 > $(OBJDIR)/traceADtgtprofile

$(OBJDIR)/allDiff_profile_p.f: $(FSOURCES) $(OBJDIR)/allDiff_profile_b.f $(OBJDIR)/allDiff_profile_d.f 
	time $(AD) -p -profile -O $(OBJDIR) -o allDiff_profile $^ 2>&1 > $(OBJDIR)/traceADparseprofile

