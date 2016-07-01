C     -----------------------------------------------
      SUBROUTINE INTERPOL( NEQUATION, ORDRE, X, Y, Z, WI, WJ,
     $                     GRADLOC, WIJ, WJI )
C     -----------------------------------------------
C    
C     Cette procedure calcule les valeurs des incon-
C     -nues physiques au milieu du segment
C
C     Ordre 1 et Ordre 2 sont programmes.
C
C     -----------------------------------------------
C     Description de la liste d'appel
C     -----------------------------------------------
C
C     ENTREE : X, Y, Z coordonnees locales
C              GRADLOC(2,5,3) pentes locales (dx,dy,dz)
C     SORTIE : WIJ et WJI variables d'etat locales
C
C     -----------------------------------------------
C     Debut des declarations
C     -----------------------------------------------
C
      IMPLICIT NONE
C
C     variables d'appel
C
      INTEGER NEQUATION, ORDRE
      REAL*8 X(2), Y(2), Z(2)
      REAL*8 WI(NEQUATION), WJ(NEQUATION)
      REAL*8 GRADLOC(2,NEQUATION,3)
      REAL*8 WIJ(NEQUATION), WJI(NEQUATION)
C
C     Variables locales
C
      REAL*8 flur(5), fltr(5), WMIN,WMAX
      REAL*8 beta, beta2, beta3, aix, aiy, aiz
      INTEGER K
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
      beta                         = 0.5
      beta2                        = beta
      beta3                        = 0.5*(1.0 - 2.0*beta)
c
         flur(1)                     = 0.0
         flur(2)                     = 0.0
         flur(3)                     = 0.0
         flur(4)                     = 0.0
         flur(5)                     = 0.0
c         
         fltr(1)                     = 0.0
         fltr(2)                     = 0.0
         fltr(3)                     = 0.0
         fltr(4)                     = 0.0
         fltr(5)                     = 0.0
c
         IF (ordre .EQ. 2) then
c
         aix                       = x(2) - x(1)
         aiy                       = y(2) - y(1)
         aiz                       = z(2) - z(1)
c
         flur(1)                     = beta2*
     &         (aix*gradloc(1,1,1) + aiy*gradloc(1,1,2) + 
     &          aiz*gradloc(1,1,3)) + beta3*(wj(1) - wi(1))
c
         flur(2)                     = beta2*
     &         (aix*gradloc(1,2,1) + aiy*gradloc(1,2,2) + 
     &          aiz*gradloc(1,2,3)) + beta3*(wj(2) - wi(2))
c
         flur(3)                     = beta2*
     &         (aix*gradloc(1,3,1) + aiy*gradloc(1,3,2) + 
     &          aiz*gradloc(1,3,3)) + beta3*(wj(3) - wi(3))
c
         flur(4)                     = beta2*
     &         (aix*gradloc(1,4,1) + aiy*gradloc(1,4,2) + 
     &          aiz*gradloc(1,4,3)) + beta3*(wj(4) - wi(4))
c
         flur(5)                     = beta2*
     &         (aix*gradloc(1,5,1) + aiy*gradloc(1,5,2) + 
     &          aiz*gradloc(1,5,3)) + beta3*(wj(5) - wi(5))
c
         fltr(1)                     = beta2*
     &         (aix*gradloc(2,1,1) + aiy*gradloc(2,1,2) + 
     &          aiz*gradloc(2,1,3)) + beta3*(wj(1) - wi(1))
c
         fltr(2)                     = beta2*
     &         (aix*gradloc(2,2,1) + aiy*gradloc(2,2,2) + 
     &          aiz*gradloc(2,2,3)) + beta3*(wj(2) - wi(2))
c
         fltr(3)                     = beta2*
     &         (aix*gradloc(2,3,1) + aiy*gradloc(2,3,2) + 
     &          aiz*gradloc(2,3,3)) + beta3*(wj(3) - wi(3))
c
         fltr(4)                     = beta2*
     &         (aix*gradloc(2,4,1) + aiy*gradloc(2,4,2) + 
     &          aiz*gradloc(2,4,3)) + beta3*(wj(4) - wi(4))
c
         fltr(5)                     = beta2*
     &         (aix*gradloc(2,5,1) + aiy*gradloc(2,5,2) + 
     &          aiz*gradloc(2,5,3)) + beta3*(wj(5) - wi(5))
c
      ENDIF
c
      DO 100 K = 1, NEQUATION, 1
C
         WIJ(K) = WI(K) + FLUR(K)
         WJI(K) = WJ(K) - FLTR(K)
c
c limitation des valeurs interpolees a l'intervalle
c limite par les valeurs nodales
c
         WMAX=MAX(WI(K),WJ(K))
         WMIN=MIN(WI(K),WJ(K))
c
         WIJ(K)=MAX(WIJ(K),WMIN)
         WIJ(K)=MIN(WIJ(K),WMAX)
c
         WJI(K)=MAX(WJI(K),WMIN)
         WJI(K)=MIN(WJI(K),WMAX)
C
  100 CONTINUE
C
      RETURN
      END
