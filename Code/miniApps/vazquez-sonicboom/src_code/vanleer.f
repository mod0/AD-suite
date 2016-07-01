
C     -----------------------------------------------
      SUBROUTINE VANLEER( NEQUATION, NBLOC, ORDRE, COORLOC, VNLOC, RHO, 
     $                    RHOU, RHOV, RHOW, ENERGIE, GRADLOC,
     $     GAMMA, FLULOC )
C     -----------------------------------------------
C
C
C     Cette procedure calcule le flux de Van-Leer 3-D
C     sur un segment.
C     
C     -----------------------------------------------
C     Description de la liste d'appel
C     -----------------------------------------------
C
C     ENTREE :
C
C     NEQUATION : le nombre d'equations resolues (5
C         en 3-D)
C     ORDRE : Odre de precision de la methode.
C     COORLOC : Tableau de reels contenant les  coor-
C         -donnees des points de  la molecule de cal-
C         -cul.
C     RHO : La masse volumique sur les deux points du
C         segment.
C     RHOU : La quantite de mouvement horizontale sur
C         les deux points du segment.
C     RHOV : La quantite de  mouvement verticale  sur
C         les deux points du segment.
C     RHOW : La quantite de  mouvement transverse sur
C         les deux points du segment.
C     ENERGIE : L'energie  totale sur les deux points 
C         du segment.
C     GAMMA : Rapport  des  chaleurs   specifiques  a
C          pression et volume constants.
C
C     SORTIE :
C
C     FLULOC : Les flux convectifs explicites sur le 
C          segment courant .
C
C     -----------------------------------------------
C     Debut des declarations
C     -----------------------------------------------
C
      IMPLICIT NONE
C
C     inclusion du header
C
C     variables d'appel
      INTEGER NEQUATION, ORDRE, NBLOC
      REAL*8 COORLOC(NBLOC,3), VNLOC(3)
      REAL*8 GAMMA
      REAL*8 RHO(NBLOC), RHOU(NBLOC), RHOV(NBLOC) 
      REAL*8 RHOW(NBLOC), ENERGIE(NBLOC)
      REAL*8 FLULOC(NEQUATION)
      REAL*8 GRADLOC(NBLOC,NEQUATION,3)
C
C     Variables locales
C
C     Tableaux de travail
      REAL*8 VPHYS(5,2), VINTA(5), VINTB(5)
      REAL*8 MV(2), U(2), V(2), W(2), PRESSION(2)
      REAL*8 FGM(5), FGP(5), FGMP(5)
C
C     Indices de boucle
C
C     Divers
      REAL*8 GAM4
      REAL*8 UN, VN, WN
      REAL*8 XNN, YNN, ZNN, NORME, RAYON, CMULT
      INTEGER ISRAYNUL
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
      GAM4= 1./ (GAMMA - 1.)
C
C     Passage variables conservatives -> variables physiques
C
      CALL PHYSCONS( NEQUATION, GAMMA, 2, RHO, RHOU, RHOV, RHOW, 
     $               ENERGIE, VPHYS )
C
C     Calcul des variables physiques interpolees de part et
C     d'autre du milieu du segment.
C
      CALL INTERPOL( NEQUATION, ORDRE, COORLOC(1,1), COORLOC(1,2),
     $               COORLOC(1,3), VPHYS(1,1), VPHYS(1,2), GRADLOC,
     $               VINTA, VINTB )
C
CJMM$$      CALL INTERPOL( NEQUATION, ORDRE, COORLOC(1,1), COORLOC(1,2),
CJMM$$     $               VPHYS(1,1), VPHYS(1,2), VPHYS(1,3), VPHYS(1,4),
CJMM$$     $               VPHYS(1,5), VPHYS(1,6), VINTA, VINTB )
C
      MV(1) = VINTA(1)
      U(1)  = VINTA(2)
      V(1)  = VINTA(3)
      W(1)  = VINTA(4)
      PRESSION(1) = VINTA(5)
C
      MV(2) = VINTB(1)
      U(2)  = VINTB(2)
      V(2)  = VINTB(3)
      W(2)  = VINTB(4)
      PRESSION(2) = VINTB(5)
c
c  debug
c      print*,'vanleer:',VINTA(1), VINTA(2),VINTA(3),VINTA(4), VINTA(5)
c      print*,'vanleer:',VINTb(1), VINTb(2),VINTb(3),VINTb(4), VINTb(5)
C
C     Calcul de l'energie interpolee.
C
      ENERGIE(1) = GAM4* PRESSION(1) + 0.5 * MV(1) * 
     $                    ( U(1)**2 + V(1)**2 + W(1)**2 )
      ENERGIE(2) = GAM4* PRESSION(2) + 0.5 * MV(2) * 
     $                    ( U(2)**2 + V(2)**2 + W(2)**2 )
C
C     Calcul de la normale au bisegment
C
      NORME = SQRT(VNLOC(1)**2 + VNLOC(2)**2 + VNLOC(3)**2)
C
      XNN = VNLOC(1)/ NORME
      YNN = VNLOC(2)/ NORME
      ZNN = VNLOC(3)/ NORME
C
C     Flux Van Leer : F = F+(NUBO1) + F-(NUBO2)
C                     G = G+(NUBO1) + G-(NUBO2)
C
C     1/ Contribution en NUBO1
C     -------------------------
C
C     --  Rotation (XN,YN) --
C
      RAYON = SQRT(XNN**2 + YNN**2)
C
C     Si RAYON <> 0, les vecteurs (- YNN,XNN,0) et 
C        (-XNN*ZNN,-YNN*ZNN,RAYON) forment une base du plan
C        orthogonal au vecteur VNOCL.
C     Sinon, ces deux vecteurs sont egaux au vecteur nul, et
C        l'on doit utiliser les vecteurs (0,1,0) et (-1,0,0)
C        pour avoir une base acceptable.
C
C     ISRAYNUL = 1 si RAYON = 0, 0 sinon.
C
      ISRAYNUL = (1.+ SIGN(1.,-RAYON))/2.
      UN = XNN* U(1) + YNN* V(1) + ZNN* W(1)
C
C     Donc, si RAYON == 0, on a VN = U et WN = V. 
C
      VN = (- YNN* U(1) + XNN* V(1) ) / (RAYON + ISRAYNUL) +
     $     ISRAYNUL* V(1)
      WN = -( XNN * U(1) + YNN * V(1) )* ZNN/ (RAYON + ISRAYNUL) +
     $     RAYON * W(1) - ISRAYNUL* U(1)
C
C     -- Decentrage
C
      cmult=-1.d0  
      CALL FDECFEX( cmult, GAMMA, MV(1), UN, VN, WN, PRESSION(1),
     $               ENERGIE(1), FGM )
C      CALL FDECFEX( -1., GAMMA, MV(1), UN, VN, WN, PRESSION(1),
C     $               ENERGIE(1), FGM )
C     
C     2/ Contribution en NUBO2
C     -------------------------
C     
C     -- Valeurs centrees + rotation (XN,YN)
C
      UN = XNN* U(2) + YNN* V(2) + ZNN* W(2)
C
C     Donc, si RAYON == 0, on a VN = U et WN = V. 
C
      VN = (- YNN* U(2) + XNN* V(2) ) / (RAYON + ISRAYNUL) +
     $     ISRAYNUL* V(2)
      WN = -( XNN * U(2) + YNN * V(2) )* ZNN/ (RAYON + ISRAYNUL) +
     $     RAYON * W(2) - ISRAYNUL* U(2)
c debug
c      print*,'vanleer3:',mv(2), un, vn, wn,gamma, PRESSION(2),
c     $ ENERGIE(2)
C     
C     -- Decentrage
C 
      cmult=1.d0    
      CALL FDECFEX( cmult, GAMMA, MV(2), UN, VN, WN, PRESSION(2),
     $               ENERGIE(2), FGP )
C      CALL FDECFEX( 1., GAMMA, MV(2), UN, VN, WN, PRESSION(2),
C     $               ENERGIE(2), FGP )
c debug
c       print*,'vanleer4:',FGP(1)    
C     3/ Report sur les noeuds des segments et rotation inverse.
C     ======================================
C
c  debug
c      print*,'vanleer5:',FGM(1), FGM(2),FGM(3),FGM(4), FGM(5)
c      print*,'vanleer5:',Fgp(1), Fgp(2),Fgp(3),Fgp(4), Fgp(5)
C
      FGMP(2) = FGM(2) + FGP(2)
      FGMP(3) = FGM(3) + FGP(3)
      FGMP(4) = FGM(4) + FGP(4)
C
      FLULOC(1) = ( FGM(1) + FGP(1) )* NORME
      FLULOC(2) = (XNN* FGMP(2) -
     $             (YNN* FGMP(3) + XNN* ZNN* FGMP(4))/
     $             (RAYON + ISRAYNUL) + FGMP(4)* ISRAYNUL)* NORME
      FLULOC(3) = (YNN* FGMP(2) +
     $             (XNN* FGMP(3) - YNN* ZNN* FGMP(4))/
     $             (RAYON + ISRAYNUL) + FGMP(3)* ISRAYNUL)* NORME
      FLULOC(4) = (ZNN* FGMP(2) + RAYON* FGMP(4) -
     $            FGMP(2)* ISRAYNUL)* NORME
      FLULOC(5) = ( FGM(5) + FGP(5) )* NORME
C
      RETURN
      END
