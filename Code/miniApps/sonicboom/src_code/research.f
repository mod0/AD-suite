      SUBROUTINE RESEARCH
C
C     Cette procedure verifie que la somme des vnocl et des vnfac pour
c     un noeud donne est bien egale a 0.
c
      Include 'Param3D.h'
C
      INTEGER IFAC,K,ISEG,I,L,IS,II,NBFAC
      REAL*8 SOMME1(3),SOMME2(3),EPSIL,SUM(3)
      REAL*8 SUM1(3),SUM2(3),S(3),TOTO(3)
C
      OPEN(15)
C
      DO 2000 IS = 1,NS
         DO 2 I = 1,3
            SOMME1(I)=0.
            SOMME2(I)=0.
            SUM1(I)=0.
            SUM2(I)=0.
 2       CONTINUE
C
         DO 3 ISEG = 1,NSEG
            IF ((NUBO(1,ISEG).EQ.IS).OR.(NUBO(2,ISEG).EQ.IS)) THEN
               IF (IS.EQ.NUBO(1,ISEG)) EPSIL = -1.
               IF (IS.EQ.NUBO(2,ISEG)) EPSIL = 1.
               DO K = 1,3
                  SOMME1(K) = SOMME1(K) + VNOCL(K,ISEG)*EPSIL
                  SUM1(K) = SUM1(K) + ABS(VNOCL(K,ISEG))
               END DO
            ENDIF
 3       CONTINUE
C
         NBFAC = 0
         DO 4 IFAC = 1,NFAC
            IF ((NSFAC(1,IFAC).EQ.IS).OR.(NSFAC(2,IFAC).EQ.IS).OR.
     $           (NSFAC(3,IFAC).EQ.IS)) THEN
               NBFAC = NBFAC + 1
               DO L = 1,3
                  SOMME2(L) = SOMME2(L) + VNFAC(L,IFAC)
                  SUM2(L) = SUM2(L) + ABS(VNFAC(L,IFAC))
               END DO
            ENDIF
 4       CONTINUE
c
         DO 5 II = 1,3
            SUM(ii) = somme1(ii) + somme2(ii)
            S(II) = SUM2(II) + SUM1(II)
            TOTO(II) = ABS(SUM(II))/S(II)
            IF (TOTO(II).GT.1.D-10) WRITE(15,*)IS,LOGFR(IS),TOTO(II)
 5       CONTINUE
C
 2000 CONTINUE
c
         RETURN
         END
