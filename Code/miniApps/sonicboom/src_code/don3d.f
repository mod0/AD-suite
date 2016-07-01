      SUBROUTINE DON3D
c---------------------------------------------------------------------   
c Reads the application parameters on file fort.10
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c---------------------------------------------------------------------
c     Global variables definition
      INTEGER itrs  , itrt  
      COMMON/comtri/itrs, itrt
c
c     Local variables definition  
      REAL*8    uin, tinf  , sound 
      REAL*8    xmcin  , xmcout
      REAL*8    tetadeg
      REAL*8    pi
c
      INTEGER i, nalpha
c
      pi                           = 4.0*ATAN(1.0)
c
      gam                          = 1.4
      gam1                         = gam - 1.0
c
c     Initialisation du nombre de fonctionnelles cout calculees
c
      ncf = 0
c
c     Prandtl number
c
      pr                           = 0.72
c 

      
      REWIND(10)
c
c     Test case identificator
c
      READ(10, *) coefm1
c
c     Initial deformation, like that produced by aero in the aeroelastic case
c     To be read from file disp.aero
c
      READ(10, *) intzz(3)
c
c     Flag for viscous/inviscid simulation
c
      READ(10, *) ivis
c
c     Reynolds number and characteristic length for the 
c     computation of the viscosity
c
      READ(10, *) rey
c
c     Free stream Mach number
c
      READ(10, *) xmach
c
c     Angles of attack on planes xy, xz (attack and yaw angles)
c
      READ(10, *) tetaxy, tetaxz
c
c     Determine from the attack and yaw angles the span-wise direction
c
      faczz(24)=   -1.0                                     !default direction z
      if (tetaxz.lt.1.e-10 .and. tetaxz.gt.-1.e-10)
     .  faczz(24)=  1.0                                     !change to direction y
c
c     Flag for fixed/moving mesh simulation
c
      READ(10, *) idefor
c
      IF ((ivis .EQ. 1) .AND. (idefor .EQ. 1)) THEN 
         WRITE(6, *) 'Viscous simulation on a deforming mesh'
         WRITE(6, *) 'This option is not yet possible'
         CALL TILT
      ENDIF
c
c     Reinitialisation strategy
c
      READ(10, *) ncont
c
c     Maximum number of non-linear iterations (time steps) for the 
c     simulation
c
      READ(10, *) ktmax
c
c     Maximal time for the simulation
c     Steady state tolerance
c
      READ(10, *) tmax, resf
c
c     Explicit/implicit flag
c
      READ(10, *) nexp
c
c     Storage strategy for the implicit Jacobian matrix
c     Linear solver identifier
c
      READ(10, *) irlax
c
c     Maximum number of relaxations
c     Tolerance for the linear system resolution
c
      READ(10, *) nbrel, err
c
c     CFL law coefficients (cflmax is used for explicit computation)
c
      READ(10,*) xcfl, ycfl, zcfl, cflmax
c
c     Convective flux solver identifier
c
      READ(10, *) iflux

c
c     Spatial approximation order
c     Predictor/corrector flag
c
      READ(10, *) nordre, ipred
c
      IF ((nordre .NE. 1) .AND. (nexp .EQ. 1)) ipred = 1
      IF (nexp .EQ. 0) ipred = 0
c
c     Local time step strategy
c
      READ(10, *) iloc
c
      READ(10, *) itrs, itrt
c
      READ(10, *) ient, epsiim, epsiex
c
c     Characteristic quantities of the flow
c
      READ(10, *) rhoref, lref, vref
c
c     Wing surface (zz) if negative, it is used the calculated skin surface
c     that is stored in faczz(25) and calculated in calcaires
c
      READ(10, *) faczz(30)

c
c     Conditions aux bords: force all st-w to be dirichlets (zz)
c
c      READ(10,*) diric, bound

      READ(10,*) diric
c
c     Verification de dpsi/dw
c
      READ(10,*) validpsi
c
c     Verification de la methode de Newton
c
      READ(10,*) Newton
c
c     Verification de la derivee de W par rapport au controle
c
      READ(10,*) testderW
c
c      Verification de la derivee de la fonctionnelle cout
c
      READ(10,*) testder
c
c     Verification de l'adjoint
c
      READ(10,*) testadj
c
c     Verification du gradient
c
      READ(10,*) testgrad
c      
c     Conditions de transpiration
c
      READ(10,*) itrans
c
c     Choix du calcul de la pression desiree dans le cas ou l'on travaille
c     sur la tuyere
c
      READ(10,*) defortuy
c
c     Choix de la construction de la pression cible si travail sur l'aile
c     ou le Falcon
c
      READ(10,*) npres
c
c     Choix sur le calcul d'initialisation du controle
c
      READ(10,*) contr
c
c     Redemarrage d'un calcul d'optimisation ou pas
c
      READ(10,*) ncontopt
c
c     Nombre maximal d'iterations d'optimisation
c
      READ(10,*) itmax
c
c     Residu a atteindre pour calculer l'etat-adjoint
c
      READ(10,*) errjacPi
c
c     Nombre de relaxations de Jacobi
c
      READ(10,*) nbrelPi
c
c     Methode One-Shot
c
      READ(10,*) oneshot
c
c     Parametre afin de calculer le pas constant pour le one-shot
c
      READ(10,*) roos
c
c     Nombre d'iterations en temps pour resoudre l'etat en One-Shot
c
      READ(10,*) ktmaxos
c
c     Nombre de relaxations dans Jacobi pour l'etat et l'etat-adjoint en
c     One-Shot
c
      READ(10,*) nbrelos
c
c     Pas de descente initial pour la methode de gradient
c
      READ(10,*) roinit
      faczz(9)=roinit
c
c     Methode multiniveau
c
      READ(10,*) multiniv
c
c     V-cycles en dents de scie
c
      READ(10,*) vcyc
c
c     Methode Full-V-Cycle
c
      READ(10,*) fvc
c
c     Nombre de fois ou un V-Cycle est repete dans la methode Full V_Cycle
c
      READ(10,*) nbvc
c
c     Nombre de V_cycles comprenant tous les niveaux (de nivg a niveau) dans 
c     la methode Full_V_Cycle, afin d'arreter le calcul sur le niveau fin.
c
      READ(10,*) nbvcfin
c
c     Niveau de grille le plus fin
c
      READ(10,*) niveau
c
c     Niveau de grille le plus grossier
c
      READ(10,*) nivg
c
c     Lissage type (zz)
c
      READ(10,*) iliss
c
c     Parametre de lissage
c
      READ(10,*) theta
c
c     VPGP: constant volume optimization by projection (zz)
c
      READ(10,*) isuct
c
c     VPGP: spanwise sections, streamwise sections and leading edge toothpaste tube effect (zz)
c
      READ(10,*) intzz(6),intzz(7),intzz(8)
c
c     Optimiztion parameters (zz) edrag,elift,epres,frali
c
      READ(10,*) faczz(1),faczz(2),faczz(3),faczz(4)
c
c     Spatial position (zz) for the "below plane" pressure gradient calculation in j
c
      READ(10,*) faczz(6), faczz(5)
c
c     Spatial position (zz) for extracting the optimization surface (typical of symetric problems)
c
      READ(10,*) faczz(8), faczz(7)
c
c     Calcul de la trainee de choc
c
      READ(10,*) coeftrainee
c
c     Lecture du coefficient de deviation de P-Pcible
c
      READ(10,*) coefpres
c
c     Angle d'incidence pour calculer les coefficients CD et CL
c
      READ(10,*) tetacdcl
c
c     Cltarget, Cdtarget and flag for use them or not
c
      READ(10,*) igicc,cltarget,cdtarget
c
c     mesh movement parameters (DIFFERENT THAN IDEFOR)
c
      imvmsh = 0                                            !default value
      kmvcou = 1000
      
      READ(10,*) imvmsh,kmvcou
c
c     calculate flow, do not optimize and maybe read init.res file
c
c
      intzz(1)=0                                            !default value
      READ(10,*) intzz(1)

c
c     correct normals in accute edges
c
      intzz(2)=0                                            !default value
      READ(10,*) intzz(2)

c
c     update transpiration condition reference normals (the "guide normals")
c
      intzz(4)=0                                            !default value
      READ(10,*) intzz(4)

c
c     keep constant incidence angle (0) or  update it to keep constant lift
c
      intzz(5)=0                                            !default value
      READ(10,*) intzz(5)

cccccccccc
cccccccccc
cccc  End reading, close fort.10
cccccccccc
cccccccccc

      CLOSE(10)

cccccccccc

      
      tetacdcl = tetacdcl*Pi/180.
      faczz(19)=-10.0d0                                     !normalization residual

      pref                         = rhoref*vref*vref

      xmu                          = 1.0

      IF (coefm1 .EQ. 1) THEN
c
c        Shock tube problem
c
         idefor                    = 0
         iloc                      = 0
         roin                      = 1.0
         uxin                      = 0.0
         uyin                      = 0.0
         uzin                      = 0.0
         pin                       = 1.0  
         roout                     = 0.125
         uxout                     = 0.0
         uyout                     = 0.0
         uzout                     = 0.0
         pout                      = 0.1
c 
         uin                       = SQRT(uxin*uxin + uyin*uyin + 
     &                                    uzin*uzin)
c
      ENDIF
c
      IF ((coefm1 .EQ. 2) .OR. (coefm1 .EQ. 5)) THEN
c
c        External flow 
c
         teta                      = tetadeg*pi/180.0
         tetaxy                    = tetaxy*pi/180.0
         tetaxz                    = tetaxz*pi/180.0

         tetaT= tetaxz
         if (faczz(24).gt.0.d0) tetaT=tetaxy

         tinf                      = 1.0/(gam*gam1*xmach*xmach)
         tbrd                      = tinf*(1.0 + (0.5*gam1*xmach*xmach))

         roin                      = 1.0

         uxin                      = cos(tetaxz) * cos(tetaxy) 
         uyin                      = sin(tetaxy) 
         uzin                      = sin(tetaxz) * cos(tetaxy)

         
         pin                       = (roin*(uxin*uxin + uyin*uyin +
     &                                      uzin*uzin))/
     &                               (gam*xmach*xmach)
c
      ENDIF
c
      IF (COEFM1.EQ.10) THEN
c
c     External flow without incidence
C
         tinf                      = 1.0/(gam*gam1*xmach*xmach)
         tbrd                      = tinf*(1.0 + (0.5*gam1*xmach*xmach))
c
         roin                      = 1.0
         uxin                      = 1.0
         uyin                      = 0.0
         uzin                      = 0.0
         pin                       = (roin*(uxin*uxin + uyin*uyin +
     &                                uzin*uzin))/(gam*xmach*xmach)
c
      ENDIF
c
c
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) '         ---->>>> A L Y A <<<<----'
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) ' PROGRAMME D"OPTIMISATION TRIDIMENSIONNELLE'
C [llh]      WRITE(6,*) ' ------------------------------------------'
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) ' INRIA-SOPHIA-ANTIPOLIS'
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) '               1. LES DONNEES DU PROBLEME'
C [llh]      WRITE(6,*) '               ************************** '
C [llh]      WRITE(6,*) ' '
c
C [llh]      IF (coefm1 .EQ. 1) WRITE(6, *) 'Travail sur le tube a choc'
C [llh]      IF ((coefm1 .EQ. 2).or.(coefm1.eq.5).or.(coefm1.eq.10))
C [llh]     $     WRITE(6, *) 'Ecoulement autour d"un obstacle'
c
      IF (iviS.EQ. 0) THEN 
C [llh]         WRITE(6,*) 'Ecoulement Eulerien'
C [llh]         WRITE(6,*) ' '
      ENDIF
c
C [llh]      IF (ivis   .EQ. 1) WRITE(6, *) 'Navier-Stokes computation'
c
C [llh]      IF (idefor .EQ. 1) WRITE(6, *) 'Deforming mesh computation'
c
C [llh]      WRITE(6, *) ' '
c
C [llh]      IF ((coefm1 .EQ. 2) .OR. (coefm1 .EQ. 5))
C [llh]     &   WRITE(6, 1000) ' Free stream Mach number : ', xmach
c
C [llh]      WRITE(6, *) ' '
c
C [llh]      IF (nexp   .EQ. 0) WRITE(6, *) 'Implicit time integration'
C [llh]      IF (nexp   .EQ. 1) WRITE(6, *) 'Explicit time integration'
c
C [llh]      WRITE(6,*) ' '
c
C [llh]      IF (iflux  .EQ. 1) WRITE(6, *) 'van Leer flux vector splitting'
C [llh]      IF (iflux  .EQ. 2) WRITE(6, *) 'Roe approximate Riemann solver'
c
C [llh]      WRITE(6,*) ' '
c
C [llh]      IF (nordre .EQ. 1) 
C [llh]     &   WRITE(6, *) 'First order spatial approximation' 
C [llh]      IF (nordre .GE. 2) 
C [llh]     &   WRITE(6, *) 'Second order spatial approximation' 
c
C [llh]      WRITE(6,*) ' '
c
      IF (diric.eq.1) then
C [llh]        WRITE(6, *) 'Forcing all Steger-Warming'
C [llh]        WRITE(6, *) 'conditions to be  Dirichlet'
      end if
c
C [llh]      WRITE(6,*) ' '
c
c      IF (bound.eq.1)
c     &   WRITE(6, *) 'Conditions de Steger-Warming'
c
C [llh]      WRITE(6,*) ' '
c
C [llh]      IF (itrans.eq.1) 
C [llh]     $     write(6,*) 'Inclusion des conditions de transpiration'
c
C [llh]      WRITE(6,*) ' '
C [llh]      IF (COEFM1.NE.10) THEN
c        WRITE(6,*) 'Travail avec la fonctionnelle minimisant'
c        WRITE(6,*) 'la trainee de choc. La variable ITRANS vaut'
c        WRITE(6,*) 'alors 1 et le controle est initialise a 0.'
c        WRITE(6,*) 'La portance cible est calculee par la '
c        WRITE(6,*) 'simulation de l"incidence en transpiration.'
c        WRITE(6,*) 'Le coefficient coeftrainee vaut :',coeftrainee
c        WRITE(6,*) 'Le coefficient coefpres vaut :',coefpres
c        WRITE(6,*) 'L"angle d"incidence vaut :',tetacdcl,'en radians.'
c        
C [llh]        write(6,*) 
C [llh]        write(6,*) 'Values at infinite:'
C [llh]        write(6,*) '-------------------'
C [llh]        write(6,*) 'X-velocity           =   ',uxin
C [llh]        write(6,*) 'Y-velocity           =   ',uyin
C [llh]        write(6,*) 'Z-velocity           =   ',uzin
C [llh]        write(6,*) 'Density              =   ',roin
C [llh]        write(6,*) 'Pressure             =   ',pin
C [llh]        write(6,*) 'Dynamic Pressure     =   ',
C [llh]     .    (roin*(uxin*uxin + uyin*uyin + uzin*uzin))/2.d0
C [llh]        write(6,*) 'Temperature          =   ',tinf
C [llh]        write(6,*) 'Mach number          =   ',xmach
C [llh]        write(6,*) 
C [llh]        write(6,*) 'Incidence angle: '
C [llh]        write(6,*)
C [llh]     .    '   Theta X-Y         =   ',tetaxy*180.d0/pi,'  degrees'
C [llh]        write(6,*)
C [llh]     .    '   Theta X-Z         =   ',tetaxz*180.d0/pi,'  degrees'
C [llh]        if(intzz(5).eq.1) then
C [llh]          write(6,*)
C [llh]     .      '   Condition         :     Variable (lift dependent) '
C [llh]        else if (intzz(5).eq.0) then
C [llh]          write(6,*)
C [llh]     .      '   Condition         :     Fixed'
C [llh]        end if
C [llh]      ENDIF

c     
C*****                    Fin des verifications              *****C
C          
C [llh]      IF ((coefm1.eq.10).and.(defortuy.eq.1).and.(itrans.eq.0)) THEN
C [llh]         write(6,*) 'Calcul de la pression desiree sur le maillage'
C [llh]         write(6,*) 'deforme geometriquement de la tuyere.'
C [llh]         WRITE(6,*) ' '
C [llh]      ENDIF
c
C [llh]      IF ((coefm1.eq.10).and.(defortuy.eq.0).and.(itrans.eq.1)) THEN
C [llh]         write(6,*) 'Calcul de la pression desiree de la tuyere par'
C [llh]         write(6,*) 'transpiration et avec un controle initialise'
C [llh]         write(6,*) 'a la tuyere cible virtuelle'
C [llh]         WRITE(6,*) ' '
C [llh]      ENDIF
c
      IF (coefm1.ne.10) then
         IF ((npres.eq.1).and.(itrans.eq.0)) then
C [llh]            write(6,*) 'Calcul de la pression desiree par lissage'
C [llh]            WRITE(6,*) ' '
         endif
c
         if ((npres.eq.0).and.(itrans.eq.1)) then
c            write(6,*) 'Pression cible calculee avec une petite '
c            write(6,*) 'incidence de ',tetacdcl,' radians simulee '
c            write(6,*) 'par transpiration.'
            WRITE(6,*) ' '
         endif
      ENDIF
c
      if ((coefm1.eq.10).and.(contr.eq.1)) then
C [llh]         write(6,*) 'Initialisation du controle a la forme desiree'
C [llh]         write(6,*) 'virtuelle'
C [llh]         WRITE(6,*) ' '
      endif
c
      if ((coefm1.eq.10).and.(contr.eq.2)) then
C [llh]         write(6,*) 'Initialisation du controle sur la tuyere a 0'
      endif
c
C [llh]      WRITE(6,*) ' '
c
c
C [llh]      WRITE(6,*) ' '
c
      if (oneshot.eq.0) then
C [llh]         WRITE(6,*) 'Optimisation par une methode de gradient :'
c         WRITE(6,*) 'recherche du pas optimal par une minimisation 1D.'
C [llh]         write(6,*) 'Le pas initial de descente rhoopt vaut',roinit
C [llh]         WRITE(6,*) ' '
      endif
c
      if (ncontopt.eq.1) then
C [llh]         write(6,*) 'Redemarrage de l"optimisation'
      endif
c
C [llh]      WRITE(6,*) ' '
c
      IF (oneshot.eq.1) THEN
C [llh]         write(6,*) 'Methode One-Shot'
C [llh]         write(6,*) '----------------'
C [llh]         write(6,*) 'Le nombre d"iterations en temps pour le'
C [llh]         write(6,*) 'calcul de l"etat vaut',ktmaxos
C [llh]         write(6,*) 'Le nombre de relaxations de Jacobi pour '
C [llh]         write(6,*) 'etat et etat-adjoint vaut',nbrelos
      ENDIF
c
      IF (multiniv.eq.1) THEN
C [llh]      write(6,*) 'Methode multiniveau'
C [llh]      WRITE(6,*) '-------------------'
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) 'Parametre de lissage theta =',theta
C [llh]      WRITE(6,*) ' '
c
      if ((fvc.eq.1).and.(vcyc.eq.1)) then
C [llh]         write(6,*)' Attention il y a incompatibilite dans le choix de
C [llh]     $        la strategie V_cycles.'
         stop
      endif
c

      if (iliss.eq.0) then
      
      if (vcyc.eq.1) then
C [llh]         write(6,*) 'V-cycles en dents de scie'
C [llh]         write(6,*) 'Le niveau le plus fin est egal a',niveau
C [llh]         write(6,*) 'Le niveau le plus grossier est egal a',nivg
      endif
c
      if ((fvc.eq.1).and.(nivg.gt.4)) then
C [llh]         write(6,*) 'Attention on ne peut pas faire du Full_V_Cycle'
C [llh]         write(6,*) 'On fait alors du V_cycles en dents de scie.'
         fvc=0
         vcyc=1
      endif
c
      if (fvc.eq.1) then
C [llh]         write(6,*) 'Methode Full V-Cycle'
C [llh]         write(6,*) ' '
C [llh]         write(6,*) 'chaque V_cycle sera repete ',nbvc,' fois.'
C [llh]         write(6,*) nbvcfin,' cycles seront effectues en alternant sur'
C [llh]         write(6,*) 'tous les niveaux.'
C [llh]         write(6,*) ' '
c
c On va alors calculer le nombre d'iterations maximales en fonction de nbvc 
c et nbvcfin afin que le calcul s'arrete sur le niveau fin.
c
         nalpha = 0
         do i = 1,nivg-1
            nalpha = nalpha + nbvc*i
         end do
         nalpha = nalpha+1
         itmax = nalpha + nbvcfin*nivg
C [llh]         write(6,*)'Le nombre maximal d"iterations est egal a :',itmax
      endif
c
      if ((fvc.eq.0).and.(vcyc.eq.0)) then
C [llh]         write(6,*) 'Travail sur le niveau',nivg
      endif
c
      else if (iliss.eq.1) then

C [llh]        write (6,*) 'Additive Multilevel'

      else if (iliss.eq.2) then

C [llh]        write (6,*)
C [llh]     .    'Additive Multilevel + Least Square Gradient Smoothing'

      end if
      ENDIF
c
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) 'Fin de la declaration des donnees du probleme '
C [llh]      WRITE(6,*) '********************************************* '


1000  FORMAT(a27,e9.3)
c

      RETURN
      END
