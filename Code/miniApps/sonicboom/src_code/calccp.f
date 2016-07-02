      SUBROUTINE CALCCP
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER IS
      REAL*8 RHO, U, V, W, E, VCARRE, PRESSION, DENOM, CP(NSMAX)
C
      DO IS = 1,NS
         RHO = UA(1,IS)
         U = UA(2,IS)/RHO
         V = UA(3,IS)/RHO
         W = UA(4,IS)/RHO
         E = UA(5,IS)
         VCARRE = U**2 + V**2 + W**2
         PRESSION = GAM1*(E - 0.5*RHO*VCARRE)
         DENOM = ROIN*0.5*VCARRE
         CP(IS) = (PRESSION - PIN)/DENOM
      END DO
C
      RETURN
      END
