      PROGRAM OPTDES3D
c---------------------------------------------------------------------   
c 
c       Programme d'optimisation 3D. On optimise une aile d'avion ou bien 
c       un demi-falcon.
c       Methode de transpiration.
c
c       Algorithme :
c       -------------
c
c           1- Initialisation :
c                 Lecture des donnees
c                 Lecture du maillage 3D
c                 Construction du maillage de peau
c                 Construction de la normale VNO a la peau, des aires
c                 Normalisation de la normale VNO (stockee dans VNON)
c                 Initialisation du control CTRL
c                 Construction des zones grossieres 
c                 Calcul de la pression cible ou bien de la trainee cible 
c           2- Optimisation :
c                 Resolution des equations d'Euler :
c                              initialisation
c                              calcul du flux explicite
c                              conditions aux bords avec transpiration
c                              construction de la matrice implicite
c                              conditions aux bords avec transpiration
c                              resolution par Jacobi
c                              obtention d'un etat W.
c                 Calcul de la fonctionnelle cout J(W,Gam)
c                 Derivees de J : DJ/DW et DJ/DGAM
c                 Calcul de l'etat-adjoint Pi
c                 Calcul du gradient Grad = DJ/DGAM - <Pi,DPSI/DGAM>
c                 Methode multi-niveau
c                 Recherche du pas optimal rhoopt
c                 CTRL = CTRL - rhoopt LPP*L*Grad
c                 Visualisation des nouvelles coordonnees de la peau :
c                         coorpbar = coorp + ctrl*VNON
c                 On remonte en 2.
c 
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c---------------------------------------------------------------------

      REAL*8 CTRL(nnsp), CTRL1(nnsp), COUTINIT, GMAX, Pasopt
      REAL*8 CD, CL, Pi, fadd,volcontrol
      INTEGER isp,i,is,k,kniv,iniv,kvar,jtp,nnsgp,nspp,ifac
      INTEGER nunit,numfic, fic,LISGRAD,itopt1,nsp1, numvigie
      INTEGER is1, is2, is3,izzout,izztask
      integer ktopt

C [llh Nan]      call second(ct0)

C [llh Nan]      cttot=0.0
C [llh Nan]      ctini=0.0
C [llh Nan]      ctoto=0.0
C [llh Nan]      cteto=0.0
C [llh Nan]      ctato=0.0
C [llh Nan]      ctjac=0.0
C [llh Nan]      ctjad=0.0
C [llh Nan]      ctjet=0.0      
c
c Lecture des donnees
c
      
      CALL DON3D

C [llh]      write (6,*)
C [llh]     .  'Old resu.data, conv.data and surf.data files erased..'
      call system('rm -f volu.res ; touch volu.res;')
      call system('rm -f conv.data ; touch conv.data;')
      call system('rm -f surf.res ; touch surf.res;')
      call system('rm -f smm*.res')
      call system('rm -f vmm*.res')
      do i=1,30
        reazz(i)= 0.0d0
      end do
      
      
c      
      IF ((COEFM1.EQ.10).AND.(DEFORTUY.EQ.1).AND.(NCONTOPT.EQ.0)) THEN
c
c  Si DEFORTUY = 1, on est alors dans le cas d'un calcul d'optimisation sur
c  une tuyere, et on calcule la pression desiree sur un maillage deforme.
C  Calcul de la pression desiree sur un maillage 3D de tuyere deforme.
C  On fait le calcul sans transpiration (ITRANS=0 et CTRL=0) et en partant 
C  d'une solution uniforme. A la fin, on remet la variable ITRANS a 1.
C
         CALL CALPRESTUYERE(CTRL)
C
      ENDIF
c
c Lecture du maillage 3D 
c
      CALL MAIL3D(faczz(7),faczz(8))
c
c Conditions de Dirichlet
c

c      If (diric.eq.1) then
c         CALL CONDDIR
c      endif

c
c Construction des segments, des normales.
c
      CALL CMVVNO
c 
      CALL CMVFAC
c
c Verfification du maillage :
c ---------------------------
c Verification de la somme des normales exterieures
c (Verification globale)
c         
      CALL CALCFAC
c
c Verification de l'orientation des normales
c
      CALL RECHERCHE
c
      CALL CMVFAC

c
c Verification de la somme des normales autour d'une cellule
c (Verification locale)
c
c      CALL RESEARCH
c
C***********************************************************************
c Reprise d'un calcul d'optimisation : ncontopt = 1, sinon, ncontopt = 0
c***********************************************************************
c
              IF (NCONTOPT.EQ.0) THEN
c             =======================
c
c Construction du maillage de peau
c
         numfic = 20
c
         CALL MAIL3D2D(numfic)
c
         CALL MAIL2D
c
c
c Construction des normales a la peau initiale.
c
         CALL CALCVNO
c
c Normalisation de VNO 
c
         CALL VNONORM
c
c Calcul des aires des cellules sur la peau.
c

         CALL CALCAIRES
c
c Initialisation aux aires initiales.
c
         CALL AIRPEAUINIT

c Smoothing the normals
c
C         call smooles     !this subroutine should be tested...

c
         CALL INITCONTROL(CTRL)

c

            

         
c
c  Si ITRANS = 1 alors :
c                  si COEFM1 = 10 alors : (cas du travail sur la tuyere)
c                              contr=1 : ctrl = forme desiree par transpiration
c                              contr=2 : initialisation a 0
c                  si COEFM1<>10 alors : (cas de l'aile ou du falcon)
c                              contr<>1 et 2 : initialisation a 0
c  Si ITRANS<>1 alors : contr peut avoir n'importe quelle valeur : 
c                                                         initialisation a 0
c
c     Fichiers ou les differentes geometries de la coque s'inscrivent
c     ainsi que l'ecoulement.
c     nunit : pour zvis a chaque iteration d'optimisation.
c     -----
c     fic : pour zvis dans la recherche du pas optimal.
c     ---
c     numvigie : pour vigie : visualisation de l'ecoulement a chaque
c     --------                iteration d'optimisation.
c     
c
         nunit = 200
C [llh]         PRINT*,'OPTDES3D.f:  NUNIT=',NUNIT
         fic = 300
         numvigie = 22
c
c*** Initialisation de la variable itopt pour la boucle d'optimisation
c
         ITOPT = 0
C
                  ELSE  ! Le cas Ncontopt.ne.0
c                ======
c
c*** Lecture du maillage de peau :
c     nombre de noeuds : nsp
c     nombre de triangles : ntp
c     coordonnees de la nouvelle geometrie : coorp
c     tableau des sommets : nup
c     aires de la nouvelle geometrie : airesp
c     aires de la geometrie initiale : airesp0
C     tableau de correspondance entre les noeuds du maillages 3D et les
c     noeuds de la coque : node3d2d
C     pression desiree : pdesp
c     trainee de choc initiale, portance : CLTARGET, CD
c

c
c*** These are the files needed for a restart:
c     
c        54: geometry, cl/d, cl/d targets
c        55: normals, deformation (i.e. ctrl)                    
c        35: nunit, itopt, ncf, COUTP, Roinit,faczz(29)            
c        59: kt0, t0  (read in etat.f) 
c        21: last solution (read in etat.f)
c                    
c                    
c     
c
c                    
c                    
         READ(54,145) NSP,NTP,NSEGP,nseg
C [llh]         write(6,*) NSP,NTP,NSEGP,nseg
         READ(54,146) ((coorp(i,isp),i=1,3),isp=1,nsp)
c         write(6,*) ((coorp(i,isp),i=1,3),isp=1,nsp)
         READ(54,148) ((nup(k,jtp),k=1,3),jtp=1,ntp)
c         write(6,*) ((nup(k,jtp),k=1,3),jtp=1,ntp)
         READ(54,148) ((nubop(k,nnsgp),k=1,2),nnsgp=1,nsegp)
c         write(6,*) ((nubop(k,nnsgp),k=1,2),nnsgp=1,nsegp)
         READ(54,146) (AIRESP(ISP),ISP=1,NSP)
c         write(6,*) (AIRESP(ISP),ISP=1,NSP)
         READ(54,146) (AIRESP0(ISP),ISP=1,NSP)
c         write(6,*) (AIRESP0(ISP),ISP=1,NSP)
         READ(54,148) (node2d3d(isp),isp=1,nsp)
c         write(6,*)  (node2d3d(isp),isp=1,nsp)
         READ(54,150) CLTARGET,CDTARGET,CL,CD
c         write(6,*) CLTARGET,CDTARGET,CL,CD
c         READ(54,149) (Pdesp(is),is=1,ns)
C
c*** Lecture des normales a la geometrie initiale : vn0
c    Lecture des normales de la nouvelle geometrie : vncoq
c    Lecture du controle : ctrl
c    
         READ(55,146) ((VNO(i,isp),i=1,3),isp=1,nsp)
         READ(55,146) ((VNON(i,isp),i=1,3),isp=1,nsp)
         READ(55,146) (CTRL(ISP),ISP=1,NSP)
C
c*** Lecture du numero de fichier courant dans lequel s'inscrira la
c            nouvelle geometrie : nunit
c    Lecture du numero de l'iteration d'optimisation courante : itopt
c    Lecture du nombre de fonctionnelles cout calculees : ncf
c    Lecture de la fonctionnelle cout a l'iteration precedente : coutp
c    Lecture du pas optimal a l'iteration precedente : Roinit
c    faczz(29) : normalization residual
c    
c         READ(35,147) nunit, fic, itopt, ncf, COUTP, Roinit
         READ(35,147) nunit, itopt, ncf, COUTP, Roinit,faczz(29)
         coutp=1.0e5
c
C [llh]         PRINT*,'OPTDES3D.f: LECTURE SUR FORT.35: NUNIT=',NUNIT
C
                     ENDIF
c                   =======
c
c
C [llh]         PRINT*,'OPTDES3D.f(3): NUNIT=',NUNIT
C
C*** Construction des zones par agglomeration
C
      CALL ZONAG
c
      CALL AIRAGNORMAG
C
c*** Initialisation de la solution soit par un ecoulement uniforme (ncont=0),
c    soit par un ecoulement non uniforme (ncont=1).
c    Si ncontopt=1 alors mettre ncont=1 ==> calcul a partir de la solution
c    precedente.
c
      CALL INIM3D

C [llh]      call showdi
      

c
                    IF (NCONTOPT.EQ.0) THEN
c                   =======================
C
C  Si on travaille sur la tuyere, et que l'on ne veut pas calculer la 
c  pression desiree sur le maillage deforme, alors, la variable DEFORTUY=0
c  et on calcule la pression desiree par transpiration simulant la deformation
c  (COEFM1=10, ITRANS=1, contr=1)
C
      IF ((COEFM1.EQ.10).AND.(DEFORTUY.EQ.0)) THEN
C
C     Ecriture dans le listing du controle initial
c
C [llh]         WRITE(6, *) ' '
C [llh]         WRITE(6, *) 'Le controle cible a pour valeurs :'
C [llh]         WRITE(6, *) ' '
c
         DO ISP = 1,NSP
            IS = NODE2D3D(ISP)
            PRINT 1515,ISP,CTRL(ISP),COOR(1,IS)
         END DO
c
 1515       FORMAT('   ISP= ',I3,'  CTRL= ', F10.4,'   X= ',F10.4)
c
C [llh]         WRITE(6, *) ' '
C [llh]         WRITE(6, *) ' '
C
C [llh]         WRITE(6, *) ' Resolution des equations d"Euler afin de'
C [llh]         WRITE(6, *) ' calculer la pression cible'
c         CALL ETAT(CTRL,KTMAX)
c         CALL PRESDESTUY(CTRL)
         if (igicc.eq.0) CALL ETAT(CTRL,ktmax)
c
      ENDIF
c
C  Si on travaille sur l'aile d'avion ou bien le falcon, on calcule la pression
c  desiree ; soit par lissage (npres=1 et ITRANS=0), soit avec une incidence 
c  simulee par transpiration (npres=0 et ITRANS=1 et contr=0).
c
C      IF ((COEFM1.NE.10).AND.(NPRES.NE.2)) THEN
C         CALL PRESCIBLE(ctrl)
C      ENDIF
c
c*** Modification le 24.05.95 ***
c    ------------------------
c    On calcule une autre fonctionnelle cout : elle va minimiser la trainee
c    de choc. J(W,gamma) = CD + omega*(CL - CLTARGET)^2
c    On doit donc calculer CLTARGET : il suffit de calculer un ecoulement et
c    de calculer CL stocke dans CLTARGET.
c    NPRES=2, CONTR=3 et ITRANS=0
c
      IF (COEFM1.NE.10) THEN

         if (igicc.eq.0) CALL ETAT(CTRL,ktmax)

      ENDIF
      STOP
c
c*** Verification qu'avant l'entree dans la boucle d'optimisation le contole
c    CTRl est bien egal a 0.
c
      CALL VERIFCONTR(CTRL)
c
c**** Results output for GiD, initial state
c
      izztask= 2
c      call zzoutsur(ctrl,izztask)
c
ccccccc

      
C [llh Nan]      call second(ct1)
C [llh Nan]      ctini= ct1 - ct0
      

      if (intzz(1).gt.0) then

C [llh]        WRITE(6, *) '************************************************'
C [llh]        WRITE(6, *) ' '
C [llh]        WRITE(6, *) '        NO OPTIMIZATION IS DEMANDED '
C [llh]        WRITE(6, *) '        ****************************'
C [llh]        WRITE(6, *) ' '
C [llh]        WRITE(6, *) ' '
C [llh]        WRITE(6, *) '************************************************'
        
C [llh]        PRINT *,'********************************'
C [llh]        PRINT *,'    F I N   N O R M A L E'
C [llh]        PRINT *,'********************************'

        stop

      end if
      
c
                                ENDIF
c                               =====
c
c***  Methode multi-niveau (Travail sur une seule grille : vcyc=0 ou
c                           Travail avec V-cycles en dents de scie : vcyc=1)
c
      KNIV = NIVG
c
c--------------------------------------------------------------
c***                    Boucle d'optimisation               ***
c--------------------------------------------------------------
c

c         WRITE(54,145) NSP,NTP,NSEGP,nseg
c         WRITE(54,146) ((coorp(i,isp),i=1,3),isp=1,nsp)
c         WRITE(54,148) ((nup(k,jtp),k=1,3),jtp=1,ntp)
c         WRITE(54,148) ((nubop(k,nnsgp),k=1,2),nnsgp=1,nsegp)
c         WRITE(54,146) (AIRESP(ISP),ISP=1,NSP)
c         WRITE(54,146) (AIRESP0(ISP),ISP=1,NSP)
c         WRITE(54,148) (node2d3d(isp),isp=1,nsp)
c         cl=0.0d0
c         cd=0.0d0
c         WRITE(54,150) CLTARGET,CDINIT,CL,CD
c
c         WRITE(55,146) ((VNO(i,isp),i=1,3),isp=1,nsp)
c         WRITE(55,146) ((VNON(i,isp),i=1,3),isp=1,nsp)
c         WRITE(55,146) (CTRL(ISP),ISP=1,NSP)
c
c*** Ecriture du numero de fichier courant dans lequel s'inscrira la
c    nouvelle geometrie : nunit
c    Ecriture du numero de l'iteration d'optimisation courante : itopt
c    Ecriture du nombre de fonctionnelles cout calculees : ncf
c    Ecriture de la fonctionnelle cout a l'iteration courante : coutp
c    Ecriture du pas optimal a l'iteration courante : Roinit
c    
c         WRITE(35,147) nunit, itopt, ncf, COUTP, Roinit
c
c         call flunow(54)
c         call flunow(55)
c         call flunow(35)
c         stop
         

      kmov   = 0
      ktopt  = 0

      
      WRITE(6, *) '************************************************'
      WRITE(6, *) ' '
      WRITE(6, *) '        6- ENTREE DANS LA BOUCLE D"OPTIMISATION '
      WRITE(6, *) '        *************************************** '
      WRITE(6, *) ' '
      WRITE(6, *) ' '
      WRITE(6, *) '************************************************'
c
1000  CONTINUE
c
      ITOPT = ITOPT + 1
      ktopt = ktopt + 1
      imvng = 0                                             !bruno wont move the mesh (default)
c
c*** Resolution des Equations d'Euler. Sortie : W
c
      
      WRITE(6, *) ' '
      WRITE(6, *) ' '
      WRITE(6, *) ' '
      WRITE(6, *) ' '
      WRITE(6, *) '************************************************'
      WRITE(6, *) ' '
      WRITE(6, *) '          ITERATION D"OPTIMISATION No',itopt
      WRITE(6, *) '          **********************************'
      WRITE(6, *) ' '
      WRITE(6, *) '************************************************'

      if (imvmsh.gt.0   .and. 
     .  mod(ktopt,kmvcou).eq.0 .and. itopt.gt.1) then        
c
c*** Move mesh
        
c... 1. Smooth surface displacement vector
c... add surface displacement to coorp and erase ctrl

        kmov=kmov+1                                         !move-mesh counter
        
        call smlesdi(ctrl)
                
c... 2. Move the mesh by calling movemesh library        

c [llh] BECAUSE OF BAD DECLARATION OF NSMAX
c        call Movelib_mov(
c     .    coor,coorp,node2d3d,nu,logfac,nubo,nsfac,ns,nt,nfac,nseg,nsp)
        
        imvng=1                                             !the mesh has moved now
        
c... 3. Recalculate the new skin stuff and other animals

        call cmvvno
        call cmvfac
        call calcfac
        call recherche
        call cmvfac
        call calcvno
        call vnonorm
        call calcaires
        call airpeauinit
        
        WRITE(6, *) '    --->> Mesh updated, keep on working...'
        WRITE(6, *) '************************************************'
        WRITE(6, *) '************************************************'
      
      end if

      
      WRITE(6, *) ' '
      WRITE(6, *) 'Entree dans la resolution de l"etat '
      WRITE(6, *) ' '
      WRITE(6, *) ' '
c
c*** Verification de la derivee de W par rapport a gamma
c      
      if (testderW.eq.1) then
         CALL TESTETAT(CTRL)
         stop
      endif

c
c  On entre dans l'optimisation (a itopt=1) avec ctrl=0
c
      
      IF (ONESHOT.EQ.1) THEN
         IF (ITOPT.EQ.1) THEN
            CALL ETAT(CTRL,KTMAX)
         ELSE
            CALL ETAT(CTRL,KTMAXOS)
         ENDIF
      ELSE
         CALL ETAT(CTRL,KTMAX)
      ENDIF

c
c*** Calcul de la fonctionnelle cout. Sortie : Cout
c
      WRITE(6, *) 'Calcul de la fonctionnelle cout '
      WRITE(6, *) '------------------------------- '
      WRITE(6, *) ' '
c
c  Calcul de la fonctionnelle cout avec uniquement la pression cible (cas
c  de la tuyere).
c
c         CALL FONCCOUT
c
c
c  Calcul d'une autre fonctionnelle cout : elle minimise la trainee de choc.
c  Procedure FnelleCout.f. Cette fois, elle depend du controle ctrl.
c  Calcul de CD et CL ==> portance.f
c  Calcul de J ==> FnelleCout.f
c
         CALL PORTANCE(CTRL,CD,CL)

         call flunow(6)
         
         if (itopt.eq.1) then
           pmold=pmnew
           dixol=0.0d0
         end if
         
         CALL FNELLECOUT(CD,CL)

c
         WRITE(6, *) ' '
         WRITE(6,*) 'CD = ',CD,' CL = ',CL
         WRITE(6,*) 'CDTARGET = ',CDTARGET,' CLTARGET = ',CLTARGET
         WRITE(6,*) 'Iteration optimisation = ',ITOPT
         WRITE(6, *) ' '
c
         ncf = ncf + 1
c
         OPEN(7,ACCESS='APPEND')
         WRITE(7,*)ITOPT,COUT
         CLOSE(7)
c        
      WRITE(6, *) ' '
      WRITE(6, *) '*** LA FONCTIONNELLE COUT VAUT : ',Cout
      WRITE(6, *) ' '

      reazz(2)= cout
      reazz(3)= cl
      reazz(4)= cd
      reazz(5)= cltarget
      reazz(6)= cdtarget

c
c*** Calcul de la derivee de J par rapport a W. Sortie : DJDW
      LISGRAD=0
      IF(ITOPT.eq.1)LISGRAD=0
      IF(LISGRAD.eq.0)   THEN
c                        ****
        CALL DCOUTDW(CTRL,CL,CD)
c
c*** Si dJ/dGamma est non nul alors,
c*** calcul de la derivee de J par rapport a Gamma. Sortie : DJDGAM
c

        CALL DCOUTDGAMM(CTRL)

            if (testder.eq.1) then 
               CALL TESTDERCOUT(CTRL)
               stop
            endif
c
         if (testadj.eq.1) then
            CALL TESTADJOINT
            stop
         endif

C         call flunow(6)

C         CALL ETATADJOINT(ctrl)
c
c***  Calcul du gradient. Sortie : Grad
c
      WRITE(6, *) ' '
      WRITE(6, *) 'Calcul du gradient '
      WRITE(6, *) '------------------ '
c
C         CALL GRADIENT(CTRL)
c
c
                         ELSE
c                        ****
c On lit le gradient sur fort.701
c
      read (701) itopt1, nsp1, (grad(isp),isp=1,nsp1)
      if(nsp1.ne.nsp)then
                     print *,'incompatibilite fort.701,nsp1=',nsp1
                     stop
                     endif
      print*,'gradient lu sur fort.701:'
c
                        ENDIF
c                       *****
c
c    Test de verification du gradient, par differences finies
c
C [llh Nan]         call second(ctf)
C [llh Nan]         cttot= ctf-ct0

         if (testgrad.eq.1) then
            COUTINIT = COUT
            CALL FORMMAGIQUE(CTRL,COUTINIT)
                 print *,'fin normale apres verification du gradient'
                 stop
         endif
C
C*** Passage au multi-niveau : calcul de LPP*L*Grad (si multiniv = 1)
C    Calcul du nouveau parametre de controle de forme CTRL
C    Visualisation de la nouvelle coque dans Ecriture.f
C    Calcul des nouvelles aires aux cellules qui serviront pour le multiniveau
c
         IF (COEFM1.EQ.10) THEN
C
C     Les points du bord ne doivent pas bouger
c
            DO ISP = 1,NSP
               IS = NODE2D3D(ISP)
               IF(COOR(1,IS).LT.0.2D0.OR.COOR(1,IS).GT.1.8D0)
     $              GRAD(ISP)=0.d0
            END DO
C
      ENDIF
c
C     Si methode Full-V-Cycles alors :
C     (elle a ete programmee pour 4 niveaux au maximum)
C
      if ((multiniv.eq.1).and.(fvc.eq.1).and.(nivg.le.4)) then
        if (iliss.eq.0) CALL FULLVCYCLE(KNIV)
      endif
      
      WRITE(6,*)' '
      if (iliss.eq.0) then
        WRITE(6,*) 'Travail sur le niveau ',kniv
      else if (iliss.eq.1) then
        WRITE(6,*) 'Multiniveau aditive'
      end if
      WRITE(6,*)' '
        
c      izztask= 1
c      itopt=2
c      call zzoutsur(ctrl,izztask)
c      write(6,*) itopt
      
      if (multiniv.eq.1) then
        
        if (iliss.eq.0) then
          

          CALL LISSAGE(KNIV)
          DO 10 ISP = 1,NSP
            GRAD(ISP) = DELTA(ISP)
   10     CONTINUE

c          call restrgr
          call smlesgr                                      !then, smooth it
          
        else if (iliss.eq.1) then
          
c...   MULTINIVEAU ADDITIVE:
          
          CALL LISSAGE(nivg)
          do isp=1,nsp
            deltaz(isp)=delta(isp)
          end do

          fadd= 1.0d0/4.0d0

          do iniv=nivg-1,1,-1
            do isp=1,nsp
              deltaz(isp) = deltaz(isp) - fadd*delta(isp)
            end do             
            CALL LISSAGE(iniv)
            do isp=1,nsp
              deltaz(isp) = deltaz(isp) + fadd*delta(isp)
            end do
            fadd= fadd*fadd
          end do

          do isp = 1,nsp
            grad(isp) = deltaz(isp)
          end do

        end if

c        call restrgr                                        !first, restrict gradient
        call smlesgr                                        !then, smooth it
        
      endif

c
c      izztask= 1
c      itopt=3
c      call zzoutsur(ctrl,izztask)
c      write(6,*) itopt

      IF (COEFM1.EQ.10) THEN
C
c     Les points du bord ne doivent pas bouger
c
            DO ISP = 1,NSP
               IS = NODE2D3D(ISP)
               IF(COOR(1,IS).LT.0.2D0.OR.COOR(1,IS).GT.1.8D0)
     $              GRAD(ISP)=0.d0
            END DO
C
         ENDIF
c
c*** Calcul du pas optimal rhoopt.
c

      IF (ONESHOT.EQ.0) THEN
C
C Pour faire des economies de cout, on calcule le pas optimal a la 1ere
C iteration d'optimisation. Ensuite, on prend des pas 2 fois plus petits pour
C aller sur les niveaux fins. Si on ne travaille que sur un niveau, ne pas
C oublier de mettre NIVG=NIVEAU.
C
         IF (NIVG.EQ.NIVEAU .or. iliss.eq.1) THEN
            CALL RHOOPT(CTRL,fic)
            PASOPT = ROOPT
         ELSE
            IF (ITOPT.LT.100) THEN 
               CALL RHOOPT(CTRL,fic)
               PASOPT = ROOPT
            ELSE
              IF (KNIV.EQ.NIVG) PASOPT = ROOPT
              IF ((KNIV.GE.NIVEAU).AND.(KNIV.LT.NIVG)) PASOPT=PASOPT/2.
            ENDIF
         ENDIF
C
      ELSE
c
c*** Pas optimaux sur chaque grille dans le cas de la tuyere
c

c         IF (KNIV.EQ.NIVG) PASOPT = 12655          !  8361.
c         IF (KNIV.EQ.2) PASOPT = 6327.5            !  4181.
c         IF (KNIV.EQ.NIVEAU) PASOPT = 3163.75      !  2090.

c
c         IF (ITOPT.EQ.1) THEN
c            GMAX = 0.
c            DO ISP = 1,NSP
c               GMAX = MAX(GMAX,GRAD(ISP))
c            END DO
c            PASOPT = ROOS/GMAX
c         ENDIF
C
      ENDIF
C
c            CALL RHOOPT(CTRL,fic)
c
c            pasopt = roopt
c         pasopt = roinit

      WRITE(6, *) ' '
      WRITE(6, *) 'Apres recherche du pas, le pas optimal vaut : '
      WRITE(6, *)  pasopt,' a l"iteration d"optimisation No ',itopt
      WRITE(6,*) 'sur le niveau',kniv
      WRITE(6, *) ' '
C

      call suctrl(volcontrol)

      write(6,*) 'volcontrol= ',volcontrol
      
      do isp = 1,nsp
        ctrl1(isp)= ctrl( isp) - pasopt * (grad(isp) - volcontrol)
        ctrl( isp)= ctrl1(isp)
      end do
      
c
c**** Results output for GiD
c
      izzout = 1
      izztask= 0
      
      if (mod(itopt,izzout).eq.0) izztask=1
C [llh]      call zzoutsur(ctrl,izztask)

      call flunow(6)

c
ccccccc

c
c**** Calcul des nouvelles coordonnees de la coque
c
      
      IF (itmax.ge.100) then
c
        if ((itopt.eq.1).or.(mod(itopt,5).eq.0)) then
C [llh]          PRINT*,'OPTDES3D.f(4): NUNIT=',NUNIT
          CALL ECRITURE(ctrl,nunit)
        endif
c
      ELSE
c
C [llh]        PRINT*,'OPTDES3D.f(5): NUNIT=',NUNIT
        CALL ECRITURE(ctrl,nunit)
c
      ENDIF
c
      nunit = nunit + 1
C
      CALL CALCAIRES
C

      if ((multiniv.eq.1).and.(vcyc.eq.1)) then
c         KNIV = KNIV - 1
c         if (KNIV.EQ.0) KNIV = NIVG
      endif

c
      CALL ECRITVIGIE(numvigie)
c
      numvigie = numvigie + 1
c
c      dixol= xman
c      pmold= pmnew

      
c
c
c
c             SORTIES
c             *******
c
c
c
C
c*** Ecriture au format vigie de l'ecoulement autour de l'obstacle.
c
c*** Sauvegarde de plusieurs variables afin de pouvoir continuer un calcul
c    d'optimisation.
c
c*** Ecriture du maillage de peau :
c     nombre de noeuds : nsp
c     nombre de triangles : ntp
c     coordonnees de la nouvelle geometrie : coorp
c     tableau des sommets : nup
c     aires de la nouvelle geometrie : airesp
c     aires de la geometrie initiale : airesp0
C     tableau de correspondance entre les noeuds du maillages 3D et les
c     noeuds de la coque : node3d2d
C     pression desiree : pdesp
c     coefficients de portance et de trainee : CLTARGET,CD,CDINIT,CL
C
         if (itopt.ge.1) then

           rewind(54)
           rewind(55)
           rewind(35)

           
           WRITE(54,145) NSP,NTP,NSEGP,nseg
           WRITE(54,146) ((coorp(i,isp),i=1,3),isp=1,nsp)
           WRITE(54,148) ((nup(k,jtp),k=1,3),jtp=1,ntp)
           WRITE(54,148) ((nubop(k,nnsgp),k=1,2),nnsgp=1,nsegp)
           WRITE(54,146) (AIRESP(ISP),ISP=1,NSP)
           WRITE(54,146) (AIRESP0(ISP),ISP=1,NSP)
           WRITE(54,148) (node2d3d(isp),isp=1,nsp)
           WRITE(54,150) CLTARGET,CDINIT,CL,CD
c         WRITE(54,149) (Pdesp(is),is=1,ns)
C
c*** Ecriture des normales a la geometrie initiale : vn0
c    Ecriture des normales de la nouvelle geometrie : vncoq
c    Ecriture du controle : ctrl
c    
           WRITE(55,146) ((VNO(i,isp),i=1,3),isp=1,nsp)
           WRITE(55,146) ((VNON(i,isp),i=1,3),isp=1,nsp)
           WRITE(55,146) (CTRL(ISP),ISP=1,NSP)
c
c*** Ecriture du numero de fichier courant dans lequel s'inscrira la
c    nouvelle geometrie : nunit
c    Ecriture du numero de l'iteration d'optimisation courante : itopt
c    Ecriture du nombre de fonctionnelles cout calculees : ncf
c    Ecriture de la fonctionnelle cout a l'iteration courante : coutp
c    Ecriture du pas optimal a l'iteration courante : Roinit
c    faczz(29) : normalization residual
c    
           WRITE(35,147) nunit, itopt, ncf, COUTP, Roinit,faczz(29)

         end if
           
         IF (ITOPT.LT.ITMAX) GOTO 1000
C
 146     FORMAT(8e15.8)
 149     FORMAT(8e15.8)
 150     FORMAT(4e15.8)
 145     FORMAT(3i8)
 147     FORMAT(3i8,2e15.8)
 148     FORMAT(8i8)

c 145     FORMAT(3i5)
c 147     FORMAT(3i5,2e15.8)
c 148     FORMAT(8i5)
C

C [llh Nan]         call second(ctf)
C [llh Nan]         cttot= ctf-ct0

C [llh Nan]         call zzcputev
         
         
         
      PRINT *,'********************************'
      PRINT *,'    F I N   N O R M A L E'
      PRINT *,'********************************'

      
      STOP
      END
