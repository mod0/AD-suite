      SUBROUTINE VERIFCONTR(CTRL)
C
C**** Verification qu'avant l'entree dans la boucle d'optimisation le contole
c    CTRl est bien egal a 0.
c
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      REAL*8 CTRL(NNSP)
      INTEGER ISP
C
      DO ISP = 1,NSP
            CTRL(ISP) = 0.
      END DO
C
      RETURN
      END    
