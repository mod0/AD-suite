
C     -----------------------------------------------
      SUBROUTINE PHYSCONS( NEQUATION, GAMMA, SIZE, RHO, RHOU, RHOV,
     $                     RHOW, ENERGIE, VPHYS )
C     -----------------------------------------------
C
C     Passage variables conservatives --> variables physiques
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
C
      INTEGER NEQUATION, SIZE
      REAL*8 GAMMA
      REAL*8 RHO(SIZE), RHOU(SIZE), RHOV(SIZE), RHOW(SIZE)
      REAL*8 ENERGIE(SIZE)
      REAL*8 VPHYS(NEQUATION,SIZE)
C
C     Variables locales
C
C     Tableaux de travail
C
C     Indices de boucle
      INTEGER K
C
C     Divers
C
C     Procedures et Fonctions
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
      DO 100 K = 1, SIZE, 1
         VPHYS(1,K) = ABS(RHO(K))
         VPHYS(2,K) = RHOU(K)/ RHO(K)
         VPHYS(3,K) = RHOV(K)/ RHO(K)
         VPHYS(4,K) = RHOW(K)/ RHO(K)
         VPHYS(5,K) = ABS((GAMMA - 1.)*( ENERGIE(K) -
     $     0.5*(RHOU(K)**2 + RHOV(K)**2 + RHOW(K)**2)/RHO(K) ))
c
  100 CONTINUE
C
      RETURN
      END
