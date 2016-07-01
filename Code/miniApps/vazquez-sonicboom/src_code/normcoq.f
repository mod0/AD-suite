      SUBROUTINE NORMCOQ(CTRL,VNCOQ)
C
c     Cette procedure calcule les normales VNCOQ a la coque.
c     Ces normales sont utilisees pour la transpiration explicite
c     (procedure Transpiration.f) et la transpiration implicite
c     (procedure ImpTranspiration.f)
c
      INCLUDE 'Paramopt3D.h'
C
      REAL*8 CTRL(NNSP), VNCOQ(3,NNSP), COORPBAR(3,NNSP)
      INTEGER ISP, K
C
C$AD II-LOOP
      DO ISP = 1,NSP
         DO K = 1,3
            COORPBAR(K,ISP) = COORP(K,ISP) + CTRL(ISP)*VNON(K,ISP)
         END DO
      END DO
C        
      CALL CALCNORMPEAU(VNCOQ,COORPBAR)    
C
      RETURN
      END
