      SUBROUTINE CALCCONTROL(CTRL,FIC)
C
c*** Remise a jour du control CTRL apres une iteration d'optimisation
c    Calcul des nouvelles coordonnees de la coque (COORP)
c
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER IS, ISP, KVAR, FIC
      REAL*8 CTRL(NNSP)
C
         DO 30 ISP = 1,NSP
            CTRL(ISP) = CTRL(ISP) - ROOPT*GRAD(ISP)
 30      CONTINUE
C
         CALL PLOTCONTROL(CTRL,FIC)
C
         DO ISP = 1,NSP
            DO KVAR = 1,3
               COORP(KVAR,ISP) = COORP(KVAR,ISP) + CTRL(ISP)
     $                           *VNO(KVAR,ISP)
            END DO
         END DO
C
         RETURN
         END
