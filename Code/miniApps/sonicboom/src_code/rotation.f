C     -----------------------------------------------
      SUBROUTINE ROTATION(NEQUATION, RHO, RHOU, RHOV, RHOW, ENERGIE,
     $                     GAMMA, VNLOC, UN, VN, WN, C, NORME, ROT )
C     -----------------------------------------------
C     CETTE PROCEDURE CALCULE LES TERMES DE LA ROTATION 3D :
C
C      _                                                                  _
C     |  1              0                   0                 0         0 |
C     |  0      Cos(theta)Cos(phi)   Sin(theta)Cos(phi)    sin(phi)     0 |
C  R= |  0        -sin(theta)          cos(theta)             0         0 |
C     |  0     -cos(theta)sin(phi)   -sin(theta)sin(phi)   Cos(phi)     0 |
C     |  0              0                   0                 0         1 |
C     |_                                                                 _|
C
C     Cos(theta)Cos(phi) = XNN ( = eta  / {N} )
C                                     x
C     Sin(theta)Cos(phi) = YNN ( = eta  / {N} )
C                                     y
C     sin(phi)           = ZNN ( = eta  / {N} )
C                                     z
C                2     2     2  1/2   
C     {N} = ( eta + eta + eta  )
C                x     y     z
C
C                 2      2  1/2
C     RAYON = (XNN  + YNN  )
C
C     -----------------------------------------------
C     Debut des declarations
C     -----------------------------------------------
C
      IMPLICIT NONE
C
C     variables d'appel
C
      INTEGER NEQUATION
      REAL*8 RHO, RHOU, RHOV, RHOW, ENERGIE
      REAL*8 GAMMA
      REAL*8 VNLOC(3)
      REAL*8 UN, VN, WN, C, ROT(NEQUATION,NEQUATION)
C
C     Divers
      REAL*8 U, V, W, P, NORME, XNN, YNN, ZNN, RAYON, ISRAYNUL
      REAL*8 RUN,RVN,RWN
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
      U  = RHOU/ RHO
      V  = RHOV/ RHO
      W  = RHOW/ RHO
      P  = (GAMMA-1.)* (ENERGIE - 0.5* RHO* (U*U + V*V + W*W))
C
C     Vitesse du son
C
      C  = SQRT( GAMMA* P/ RHO )
C
C     Calcul de la normale au segment courant
C
C     Norme du vecteur normal a la face de la cellule.
C
      NORME = SQRT( VNLOC(1)**2 + VNLOC(2)**2 + VNLOC(3)**2 )
C
C     XNN, YNN, ZNN : composantes du vecteur normal normalisees
C
      XNN = VNLOC(1)/ NORME
      YNN = VNLOC(2)/ NORME
      ZNN = VNLOC(3)/ NORME
C       
C     Rotation (U,V,W) --> (UN,VN,WN) et variables auxilliaires
C
      RAYON = SQRT(XNN*XNN + YNN*YNN)
C
      UN =   XNN* U + YNN* V + ZNN* W
C
C     * Si RAYON <> 0, les vecteurs (- YNN,XNN,0) et 
C        (-XNN*ZNN,-YNN*ZNN,RAYON) forment une base du plan
C        orthogonal au vecteur VNOCL.
C     * Sinon, ces deux vecteurs sont egaux au vecteur nul, et
C        l'on doit utiliser les vecteurs (0,1,0) et (-1,0,0)
C        pour avoir une base acceptable. C'est-a-dire que :
C        U se transforme en W, V reste V et W se transforme en -U,
C        d'ou la matrice de rotation :
C                             1  0  0  0  0 
C                             0  0  0  1  0
C                             0  0  1  0  0
C                             0 -1  0  0  0
C                             0  0  0  0  1
C
C      *  ISRAYNUL = 1 si RAYON = 0, 0 sinon.
C
      IF ( RAYON .LT. 1.E-8 ) THEN
         ISRAYNUL = 1.
      ELSE
         ISRAYNUL = 0.
      ENDIF
C
C     Donc, si RAYON == 0, on a VN = V et WN = -U. 
C
      VN = (- YNN* U + XNN* V ) / (RAYON + ISRAYNUL) + ISRAYNUL* V
      WN = -( XNN * U + YNN * V )* ZNN/ (RAYON + ISRAYNUL) +
     $     RAYON * W - ISRAYNUL* U
c
C     On calcule les coefficients de la matrice de rotation et
C     on les stocke dans le tableau ROT.
C
      ROT(1,1) = 1.
      ROT(1,2) = 0.
      ROT(1,3) = 0.
      ROT(1,4) = 0.
      ROT(1,5) = 0.
C
      ROT(2,1) = 0.
      ROT(2,2) = XNN
      ROT(2,3) = YNN
      ROT(2,4) = ZNN 
      ROT(2,5) = 0.
C
      ROT(3,1) = 0.
      ROT(3,2) = - YNN/ (RAYON + ISRAYNUL) 
      ROT(3,3) = XNN/ (RAYON + ISRAYNUL) + ISRAYNUL
      ROT(3,4) = 0.
      ROT(3,5) = 0.
C
      ROT(4,1) = 0.
      ROT(4,2) = - XNN* ZNN/ (RAYON + ISRAYNUL) - ISRAYNUL
      ROT(4,3) = - YNN* ZNN/ (RAYON + ISRAYNUL) 
      ROT(4,4) = RAYON
      ROT(4,5) = 0.
C
      ROT(5,1) = 0.
      ROT(5,2) = 0.
      ROT(5,3) = 0.
      ROT(5,4) = 0.
      ROT(5,5) = 1.
c                            t    _
c*** Test : On veut voir si   R * W = W
C
      RUN = ROT(2,2)*UN + ROT(3,2)*VN + ROT(4,2)*WN
      RVN = ROT(2,3)*UN + ROT(3,3)*VN + ROT(4,3)*WN
      RWN = ROT(2,4)*UN + ROT(3,4)*VN + ROT(4,4)*WN
C
      RUN = ABS(RUN-U)
      RVN = ABS(RVN-V)
      RWN = ABS(RWN-W)
c
c      if (RUN.gt.1d-10) print*,'diff1 = ',RUN
c      if (RVN.gt.1d-10) print*,'diff2 = ',RVN
c      if (RWN.gt.1d-10) print*,'diff3 = ',RWN
c
      RETURN
      END
