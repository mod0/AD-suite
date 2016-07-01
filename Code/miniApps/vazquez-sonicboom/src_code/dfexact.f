
C     -----------------------------------------------
      SUBROUTINE DFEXACT( NEQUATION, RHO, UN, VN, WN, E, GAMMA, GAM1,
     $                    GAM3, DFX )
C     -----------------------------------------------
C
C     Calcul de DFX jacobien exact de F en "le point courant".
C
C              _             _
C             |              |
C             |   Rho*UN     |
C             |         2    |
C             |   Rho*UN + P |
C             |              |
C     F(W)=   |   Rho*UN*VN  |
C             |              |
C             |   Rho*UN*WN  |
C             |              |
C             |   UN*(E+P)   |
C             |_            _|
C                                              2     2     2
C      P(W) = (gamma - 1) * (E - 0.5 * Rho *(UN  + VN  + WN )
C                                                2          2          2
C                                      ( (Rho*UN) + (Rho*VN) + (Rho*WN) ) 
C           = (gamma - 1) * (E - 0.5 * ( -------------------------------) )
C                                      (            Rho                 )
C
C     gam3 = gamma - 3
C     gam1 = gamma - 1
C
C     Derivation par rapport aux variables conservatives :
C                 (Rho, Rho*UN, Rho*VN, Rho*WN, E)
C
C     -----------------------------------------------
C     Debut des declarations
C     -----------------------------------------------
C
      IMPLICIT NONE
C
C     variables d'appel
      INTEGER NEQUATION
      REAL*8 UN, VN, WN, E, RHO, GAMMA, GAM1, GAM3
      REAL*8 DFX(NEQUATION,NEQUATION)
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
      DFX(1,1)= 0.
      DFX(1,2)= 1.
      DFX(1,3)= 0.
      DFX(1,4)= 0.
      DFX(1,5)= 0.
C
      DFX(2,1)= 0.5*( GAM3*UN**2 + GAM1*VN**2 + GAM1*WN**2)
      DFX(2,2)= - GAM3*UN
      DFX(2,3)= - GAM1*VN
      DFX(2,4)= - GAM1*WN
      DFX(2,5)= GAM1
C
      DFX(3,1)= - UN*VN
      DFX(3,2)= VN
      DFX(3,3)= UN
      DFX(3,4)= 0.
      DFX(3,5)= 0.
C
      DFX(4,1)= - UN*WN
      DFX(4,2)= WN
      DFX(4,3)= 0.
      DFX(4,4)= UN
      DFX(4,5)= 0.
C
      DFX(5,1)= - GAMMA*E*UN/RHO + GAM1*UN*( UN**2 + VN**2 + WN**2 ) 
      DFX(5,2)= GAMMA*E/RHO - 0.5*GAM1*( 3.*UN**2 + VN**2 + WN**2 )
      DFX(5,3)= - GAM1*UN*VN
      DFX(5,4)= - GAM1*UN*WN
      DFX(5,5)= GAMMA*UN
C
      RETURN
      END
