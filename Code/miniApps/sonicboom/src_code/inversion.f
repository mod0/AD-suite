         
      SUBROUTINE INVERSION
C
c     Inversion of the diagonal terms
c
      INCLUDE 'Param3D.h'
C
C     Variables locales
C
      REAL*8 RMAT(5,5), var
      INTEGER IS, IA, K, I, IB
C
C     Fin des declarations
C

      
C$ DOACROSS LOCAL(IS,IA,IB,K,VAR,I,RMAT)
      DO 800 is=1,ns

C
         DO 805 ia=1,5
C
            Do 100 ib = 1,5

c
            rmat(ia,ib) = diag(is,ia,ib)

c
100         Continue
C            
805     CONTINUE
C         
         DO 815 k=1,5
C            
            var = 1.0/rmat(k,k)
            rmat(k,k) = 1.0
C    
            Do 101 ia = 1,5
c
            rmat(k,ia) = rmat(k,ia)*var
C
101         Continue
C            
            DO 850 i=1,5
C               
               IF (i .EQ. k) GOTO 850
C               
               var = rmat(i,k)
               rmat(i,k) = 0.0
C               
               Do 102 ia = 1,5
c
               rmat(i,ia) = rmat(i,ia) - var*rmat(k,ia)
C
102            Continue
C               
850        CONTINUE
C            
815     CONTINUE
C         
         DO 860 ia=1,5
C          
            Do 103 ib = 1,5
c
c         if (is.eq.3603) write(6,*) 'nuuu',ib,ia,rmat(ia,ib),var
            diag(is,ia,ib) = rmat(ia,ib)
c
103         Continue
C            
860     CONTINUE
C         
800   CONTINUE
C


      RETURN
      END
