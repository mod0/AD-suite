      SUBROUTINE CALCAIRES
C
C**** Cette procedure calcule les aires des triangles sur la coque, 
C     AIRTP(1:jtp) et les aires des cellules duales de la coque, AIRESP(1:nsp)
c     utilisees dans la procedure Lissage.f.
C
      INCLUDE 'Paramopt3D.h'
C
      INTEGER JTP, ISP, I
C
      DO 10 JTP = 1,NTP
         AIRTP(JTP) = SQRT(VNP(1,JTP)**2+VNP(2,JTP)**2+VNP(3,JTP)**2)
10    CONTINUE
C
      DO 50 ISP = 1,NSP
         AIRESP(ISP) = 0.
 50   CONTINUE
C
      DO 20 ISP = 1,NSP
         DO 30 JTP = 1,NTP
            DO 40 I = 1,3
               IF (NUP(I,JTP).EQ.ISP) THEN
                  AIRESP(ISP) = AIRESP(ISP) + AIRTP(JTP)/3.
               ENDIF
 40         CONTINUE
 30      CONTINUE
 20   CONTINUE
C
      RETURN
      END
