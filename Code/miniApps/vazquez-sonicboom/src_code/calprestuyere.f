


      SUBROUTINE CALPRESTUYERE(CTRL)
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
c*** Cette procedure calcule la pression desiree sur un maillage de tuyere
c    deforme. (cas ou DEFORTUY=1 et COEFM1=10)
c
      INTEGER IS, IFAC, ISEG, is1, is2, is3, NUMFIC, NUNIT
      INTEGER k, jt, isp, c, d, e,f,i,g
      REAL*8 PI, ymincoq, CTRL(NNSP)
      REAL*8 XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
C
C   Lecture du maillage initial de la tuyere : c'est la geometrie d'un tube 
C   a choc :
c
c     y
c     ^
c     |   /-------------------------/
c     | / |                       / |    z
c  0.5|---|---------------------|   |   /
c     |   |                     |   | /
c     |   |---------------------|---|   0.5
c     | /                       | /
c     |-------------------------|---------->
c     -2                        4          x
c
c
         READ(60,*) ns, nt, nfac
         READ(60, *) ((coor(i,is)   , i=1,3), is=1,ns)
         READ(60, *) ((nu(k,jt)     , k=1,4), jt=1,nt)
         READ(60, *) (logfac(ifac)  , ifac=1,nfac)
         READ(60, *) ((nsfac(k,ifac), k=1,3), ifac=1,nfac)
c
c   A partir du maillage initial, on construit le maillage desire :
c   on change le plan y=5 en y=f(x)
c
         PI = acos(-1.0)
c
         ymincoq = 1.0e+16
         f=0
         g=0
c
         DO IS = 1,NS
c
            IF ((((COOR(1,IS).GE.0.).AND.(COOR(1,IS).LE.2.)).AND.
     $           (COOR(2,IS).EQ.5.))) THEN
               f=f+1
               COOR(2,IS) = 4.75 + 0.25 *SIN(PI*(COOR(1,IS)+0.5)) 
               ymincoq = MIN(ymincoq, coor(2,is))
            ENDIF
c
            IF ((((COOR(1,IS).GE.0.).AND.(COOR(1,IS).LE.2.)).AND.
     $           (COOR(2,IS).EQ.4.))) THEN
               g=g+1
               COOR(2,IS) = 3.75 + 0.25 *SIN(PI*(COOR(1,IS)+0.5)) 
            ENDIF
c
         END DO
c
         print*,'Le nombre de noeuds ayant pour ordonnee y=f(x) est :'
         print*,f,g
c
c   Redefinition des logfac : 
c   sur les plans x=-2 et x=4, conditions d'entree-sortie : logfac=4
c   partout ailleurs, conditions de glissement, mais on definit a part, 
c   la partie ou y = f(x) (logfac=7)
c
      xmin                         = 1.0e+16
      ymin                         = 1.0e+16
      zmin                         = 1.0e+16
      xmax                         =-1.0e+16
      ymax                         =-1.0e+16
      zmax                         =-1.0e+16
c
      DO 100 is=1,ns
         xmin                      = MIN(xmin, coor(1,is)) 
         ymin                      = MIN(ymin, coor(2,is)) 
         zmin                      = MIN(zmin, coor(3,is)) 
         xmax                      = MAX(xmax, coor(1,is)) 
         ymax                      = MAX(ymax, coor(2,is)) 
         zmax                      = MAX(zmax, coor(3,is)) 
100   CONTINUE 
c
         c=0
         d=0
         e=0
c
         DO IFAC = 1,NFAC
            IF (((COOR(1,NSFAC(1,IFAC)).EQ.XMIN).AND.
     $           (COOR(1,NSFAC(2,IFAC)).EQ.XMIN).AND.
     $           (COOR(1,NSFAC(3,IFAC)).EQ.XMIN)).OR.
     $           (COOR(1,NSFAC(1,IFAC)).EQ.XMAX).AND.
     $           (COOR(1,NSFAC(2,IFAC)).EQ.XMAX).AND.
     $           (COOR(1,NSFAC(3,IFAC)).EQ.XMAX)) THEN
               LOGFAC(IFAC)=4
               c = c+1
            ELSE
               IF (((COOR(2,NSFAC(1,IFAC)).ge.ymincoq).AND.
     $              (COOR(2,NSFAC(2,IFAC)).ge.ymincoq).AND.
     $              (COOR(2,NSFAC(3,IFAC)).ge.ymincoq)).and.
     $         (((COOR(1,NSFAC(1,IFAC)).GE.0.).AND.
     $              (COOR(1,NSFAC(1,IFAC)).LE.2.))
     $              .AND.
     $          ((COOR(1,NSFAC(2,IFAC)).GE.0.).AND.
     $              (COOR(1,NSFAC(2,IFAC)).LE.2.))
     $              .AND.
     $          ((COOR(1,NSFAC(3,IFAC)).GE.0.).AND.
     $              (COOR(1,NSFAC(3,IFAC)).LE.2.)))) THEN
                  LOGFAC(IFAC) = 7
                  d = d+1
               ELSE
                  LOGFAC(IFAC) = 2
                  e=e+1
               ENDIF
            ENDIF
         END DO
c
         print*,'Nombre de faces non glissantes = ',c
         print*,'Nombre de faces glissantes sur la coque = ',d
         print*,'Nombre de faces glissantes pas sur la coque = ',e
c
      CALL CALFRO
c
      CALL SEG3D
c
      DO 105 is=1,ns
         ndeg(is)                  = 0
105   CONTINUE  
c
c     Forms the list of edges attached to each vertex 
c
      DO 5 iseg=1,nseg
c
         is1                       = nubo(1,iseg)
         is2                       = nubo(2,iseg)
         ndeg(is1)                 = ndeg(is1) + 1
         IF (ndeg(is1) .GT. ndmax) GOTO 15
         ndeg(is2)                 = ndeg(is2) + 1
         IF (ndeg(is2) .GT. ndmax) GOTO 15
         jaret(is1,ndeg(is1))      = iseg
         jaret(is2,ndeg(is2))      = iseg
c
5     CONTINUE
c
      GOTO 1010
c
15    CONTINUE
c
      WRITE(6, *) 'Error'
      WRITE(6, *) 'Increase NDMAX : ', ndmax
c
1010  CONTINUE
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
C
         DO 85 is=1,ns
            logfr(is)              = 0            
85       CONTINUE
c
         do ifac = 1,nfac
            if ((logfac(ifac).eq.2).or.(logfac(ifac).eq.7)) then
               logfr(nsfac(1,ifac))   = 2
               logfr(nsfac(2,ifac))   = 2
               logfr(nsfac(3,ifac))   = 2
            endif
         end do
c
         do ifac = 1,nfac
            if (logfac(ifac).eq.4) then
               logfr(nsfac(1,ifac))   = 4
               logfr(nsfac(2,ifac))   = 4
               logfr(nsfac(3,ifac))   = 4
            endif
         end do
c
c Verification de l'orientation des normales
c
      CALL RECHERCHE
C
      CALL CMVFAC
c
      CALL RESEARCH
c
c Construction du maillage de peau
c
         numfic = 199
c
         CALL MAIL3D2D(numfic)
c
         CALL MAIL2D
c
         CALL INIM3D
c
        print*,'ns = ',ns,' nt = ',nt,' nseg = ',nseg,' nfac = ',nfac
C
c   Calcul de la pression desiree, sans transpiration : on aura mis dans don3d
c   la variable itrans a 0. On part d'une solution uniforme (on a mis ncont
c   a 0).
c
        CALL INITCONTROL(CTRL)
c
         CALL ETAT(CTRL,KTMAX)
C
         DO IS = 1,NS
            PDESP(IS) = gam1*(ua(5,is) - 0.5*(ua(2,is)**2 +
     $              ua(3,is)**2 + ua(4,is)**2)/ua(1,is))
            write(24,*) PDESP(IS)
         END DO
C
c   Ecriture dans le fichier fort.21, de la pression desiree au format Vigie
c
c         NUNIT = 21
c         CALL RESPRES(NUNIT)
c
c         CALL ECRITVIGIE
C 
c   Apres avoir calcule la pression desiree sur la geometrie finale, on remet
c   la variable ITRANS a 1. De plus, pour commencer le calcul d'optimisation,
c   on met la variable CONTR a 2 (cas ou COEFM1=10, ITRANS=1 et CTRL=0), de 
c   facon a commencer avec un controle CTRL=0.
c
         ITRANS = 1
         CONTR = 2
C
         RETURN
         END
