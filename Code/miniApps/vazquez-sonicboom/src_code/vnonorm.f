      SUBROUTINE VNONORM
C
C*** Normalisation de la normale VNO a la coque : stockee dans VNON
c    VNON est utilise dans le calcul des nouvelles coordonnees de la coque a
c    chaque iteration d'optimisation.
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER ISP, I
      REAL*8 NORMALIS
C
      DO ISP = 1,NSP
      NORMALIS = 0.
         DO I = 1,3
            NORMALIS = NORMALIS + VNO(I,ISP)*VNO(I,ISP)
         END DO
         NORMALIS = SQRT(NORMALIS)
         DO I = 1,3
            VNON(I,ISP) = VNO(I,ISP)/NORMALIS
         END DO
      END DO
C
      RETURN
      END
