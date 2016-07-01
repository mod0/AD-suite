      SUBROUTINE CDMAT(NEQUATION,
     $     NNdim,NNSGdim,NUBO,LOGFR,ZM,DIAG,nn,nnsg)
C
C**** Insertion de conditions de Dirichlet aux bords : 
c     si logfac(ifac)=4, alors logfr(nsfac(i=1,3,ifac)) = 1
C    
C     Les termes extra-diagonaux de la matrice implicite sont mis a 0
C     et le bloc diagonal est reduit a l'identite.
C
      IMPLICIT NONE

      INTEGER NEQUATION, NN, NNSG,NNdim,NNSGdim
      INTEGER ISEG, IA, NUBO1, NUBO2, IS, IB
      REAL*8 ZM(NEQUATION,NEQUATION,2,NNSGdim),
     .  DIAG(NNdim,NEQUATION,NEQUATION)
      INTEGER LOGFR(NNdim), NUBO(2,NNSGdim)
C
      DO ISEG = 1,NNSG
         NUBO1 = NUBO(1,ISEG)
         NUBO2 = NUBO(2,ISEG)
         IF (LOGFR(NUBO1).EQ.1 .or. LOGFR(NUBO1).EQ.5) THEN 
               DO IA = 1,5
                  DO IB = 1,5

                     ZM(IA,IB,2,ISEG) = 0.
                     ZM(IA,IB,2,ISEG) = 0.
                     ZM(IA,IB,2,ISEG) = 0.
                     ZM(IA,IB,2,ISEG) = 0.
                     ZM(IA,IB,2,ISEG) = 0.

                  END DO
               END DO
         ENDIF
         IF (LOGFR(NUBO2).EQ.1 .or. LOGFR(NUBO2).EQ.5) THEN
               DO IA = 1,5
                  DO IB = 1,5
                     ZM(IA,IB,1,ISEG) = 0.
                     ZM(IA,IB,1,ISEG) = 0.
                     ZM(IA,IB,1,ISEG) = 0.
                     ZM(IA,IB,1,ISEG) = 0.
                     ZM(IA,IB,1,ISEG) = 0.
                  END DO
               END DO
         ENDIF
      END DO
C
      DO IS = 1,NN
         IF (LOGFR(IS).EQ.1 .or. LOGFR(IS).EQ.5) THEN
            DIAG(IS,1,1) = 1.
            DIAG(IS,1,2) = 0.
            DIAG(IS,1,3) = 0.
            DIAG(IS,1,4) = 0.
            DIAG(IS,1,5) = 0.
C
            DIAG(IS,2,1) = 0.
            DIAG(IS,2,2) = 1.
            DIAG(IS,2,3) = 0.
            DIAG(IS,2,4) = 0.
            DIAG(IS,2,5) = 0.
C
            DIAG(IS,3,1) = 0.
            DIAG(IS,3,2) = 0.
            DIAG(IS,3,3) = 1.
            DIAG(IS,3,4) = 0.
            DIAG(IS,3,5) = 0.
C
            DIAG(IS,4,1) = 0.
            DIAG(IS,4,2) = 0.
            DIAG(IS,4,3) = 0.
            DIAG(IS,4,4) = 1.
            DIAG(IS,4,5) = 0.
C
            DIAG(IS,5,1) = 0.
            DIAG(IS,5,2) = 0.
            DIAG(IS,5,3) = 0.
            DIAG(IS,5,4) = 0.
            DIAG(IS,5,5) = 1.
C
         ENDIF
      END DO
C
      RETURN
      END
