OBJDIR=../generated
ADDIR=../src_deriv
EXEDIR=../data

ADSOURCENAMES=conddirflux.f  fluroe.f  gradnod.f \
psiroe.f  transpiration.f  vcurvm.f gradfb.f normcoq.f calcnormpeau.f
INCLUDENAMES=Param3D.h Paramopt3D.h param.h
FSOURCENAMES=OptDes3D.f second.f don3d.f calprestuyere.f mail3d.f \
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
 gradient.f testadjoint.f testdercout.f \
 portance.f adjointe.f fnellecout_gp.f fcnelj.f vanleer.f \
 gtolocal.f dcoutdgamm.f dcoutdw_gp.f testetat.f airpeauinit.f \
 smlesdi.f calcvno.f vnonorm.f physcons.f interpol.f fdecfex.f \
 verifcontr.f zonag.f airagnormag.f showdi.f

FSOURCES=$(FSOURCENAMES:%=%)
FOBJECTS=$(FSOURCENAMES:%.f=$(OBJDIR)/%.o)
ADSOURCES=$(ADSOURCENAMES:%=%)
INCLUDES=$(INCLUDENAMES:%=%)

AD=tapenade -msglevel 5 -debug 30
FC=f77 -O
CC=cc

all: runprimal runtapenade


runtapenade : $(EXEDIR)/testADtapenade.exe $(OBJDIR)
	cd $(EXEDIR) ; ./testADtapenade.exe > $(OBJDIR)/Results ; cd $(OBJDIR) ; cat Results

$(EXEDIR)/testADtapenade.exe : $(OBJDIR) $(FOBJECTS) $(OBJDIR)/testAD_tapenade.o $(OBJDIR)/allDiff_d.o $(OBJDIR)/allDiff_b.o $(OBJDIR)/adBuffer.o $(OBJDIR)/adStack.o $(INCLUDES)
	$(FC) $(FOBJECTS) $(OBJDIR)/testAD_tapenade.o $(OBJDIR)/allDiff_d.o $(OBJDIR)/allDiff_b.o $(OBJDIR)/adBuffer.o $(OBJDIR)/adStack.o -o $@

$(OBJDIR)/testAD_tapenade.o : testAD_tapenade.f $(INCLUDES)
	$(FC) -c $< -o $@

$(OBJDIR)/allDiff_d.o : $(ADDIR)/allDiff_d.f $(ADDIR)/Paramopt3D.h
	$(FC) -c $< -o $@

$(ADDIR)/allDiff_d.f : $(ADSOURCES) $(INCLUDES)
	time $(AD) -d -head psiroe -vars "ua" -outvars "ce" -O $(ADDIR)/ -o allDiff $(ADSOURCES) 2>&1 > $(ADDIR)/traceADtgt

$(OBJDIR)/allDiff_b.o : $(ADDIR)/allDiff_b.f $(ADDIR)/Paramopt3D.h
	$(FC) -c $< -o $@

$(ADDIR)/allDiff_b.f : $(ADSOURCES) $(INCLUDES)
	time $(AD) -b -head psiroe -vars "ua" -outvars "ce" -O $(ADDIR)/ -o allDiff $(ADSOURCES) 2>&1 > $(ADDIR)/traceADadj

$(OBJDIR)/adBuffer.o : ../src_deriv/adBuffer.f
	$(FC) -c $< -o $@

$(OBJDIR)/adStack.o : ../src_deriv/adStack.c
	$(CC) -c $< -o $@

$(ADDIR)/Paramopt3D.h : $(INCLUDES)
	cp *.h $(ADDIR)



runprimal : $(EXEDIR)/primal.exe $(OBJDIR)
	cd $(EXEDIR) ; ./primal.exe > $(OBJDIR)/Results ; cd $(OBJDIR) ; cat Results

$(EXEDIR)/primal.exe : $(OBJDIR) $(FOBJECTS) $(OBJDIR)/testAD_primal.o $(INCLUDES)
	$(FC) $(FOBJECTS) $(OBJDIR)/testAD_primal.o -o $@

$(OBJDIR)/testAD_primal.o : testAD_primal.f $(INCLUDES)
	$(FC) -c $< -o $@



$(FOBJECTS):$(OBJDIR)/%.o: %.f $(INCLUDES)
	$(FC) -c $< -o $@



$(OBJDIR) :
	mkdir -p $(OBJDIR)

clean :
	rm -rf $(OBJDIR)/* $(EXEDIR)/*.exe $(ADDIR)/Param3D.h $(ADDIR)/Paramopt3D.h $(ADDIR)/param.h





