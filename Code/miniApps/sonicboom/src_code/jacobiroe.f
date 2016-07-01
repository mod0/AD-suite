      SUBROUTINE JACOBIROE
c---------------------------------------------------------------------   
c     Resout le systeme lineaire : A X = f (1)
c              Vol   dPSI  n         n+1   n               n
c     avec A = --- - ----(W ) , X = W   - W   et f = -PSI(W )
c              dt     dW          
c
c     On pose A = D - L - U
c                     n+1            n
c     Donc, (1) => D X    = (L + U) X  + f
c                     n+1    -1           n
c               =>   X    = D  ( (L + U) X  + f )
c
c                                                       n
c     On va calculer en un premier temps le terme (L+U)X  + f
c     et ensuite, on le multiplie par l'inverse de D.
c
c     Les termes extra-diagonaux sont stockes dans STMAT(nsegs) et
c     les termes diagonaux inverses sont stockes dans DIAG(ns,5,5).
c
c     Pour un segment donne iseg, il y a donc 2 extremites nubo(1,iseg)
c     et nubo(2,iseg). La matrice elementaire STMAT a donc 2 x (5x5) 
c     elements.
c
c     Construction des termes de la matrice implicite :
c
c     En se placant sur la ligne i, le flux a ete calcule de la facon
c     suivante :
c           ___     L  R                L  R
c     PSI = \  Phi(W ,W ,eta) - \  Phi(W ,W ,eta)
c        i  /       i  j        /       j  i
c           ---                 ---
c          iseg                iseg
c         i=nubo1             i=nubo2
c         j=nubo2             j=nubo1
c
c            ___                     ___
c            \  dPHI      dPHI       \  dPHI       dPHI
c     DPSI = /  ---- dW + ---- dW  - /  ---- dW  + ---- dW
c            --- dW    i    dW   j   --- dW    j    dW    i
c         i=nubo1  L          R    i=nubo2 L          R
c
c                           dPHI                  dPHI
c   => Termes diagonaux :   ----(W ,W )    et   - ----(W ,W )
c                            dW   i  j             dW   j  i
c                              L                     R
c             => DIAG
c
c                               dPHI               dPHI
c   => Termes extradiagonaux :  ----(W ,W )  et  - ----(W ,W )
c                                dW   i  j          dW   j  i
c                                  R                  L
c                                    ||                 ||
c                                    \/                 \/
c                               STMAT(1..25)      STMAT(26..51)
c

c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
c     Local variables definition
      INTEGER iseg  , is, ia, i, j, id1, id2
      REAL*8    res   , res0
      REAL*8    x1, x2, x3, x4, x5
c
c     ***************** INITIALISATION ********************
c                    n+1   -1 
c     On initialise X   a D  f . 
c                               n
c     Dans f, est stocke - PSI(W ), c'est-a-dire ce(5,ns)
c
c      DO 5 ia=1,5
c
c         DO 10 is=1,ns
c
c            dx(ia,is)              = diag(is,ia,1)*ce(1,is) + 
c     &                               diag(is,ia,2)*ce(2,is) + 
c     &                               diag(is,ia,3)*ce(3,is) + 
c     &                               diag(is,ia,4)*ce(4,is) +
c     &                               diag(is,ia,5)*ce(5,is)
c
c10       CONTINUE
c
c5     CONTINUE
c     
      IF (nbrel .GT. 0) THEN
c
         DO 300 nit=1,nbrel
c
c********************** MISE A JOUR *********************
c                         n
c     Dans dy est stocke X
c     Dans dx, on met f pour commencer a calculer les termes
c     extradiagonaux.
c

C$DOACROSS LOCAL(is,ia)
           do is=1,ns
             do ia=1,5
               dy(ia,is)        = dx(ia,is)
               dx(ia,is)        = ce(ia,is)
             end do
             dx(5,is)            = adh1(is)*dx(5,is)
           end do

c
c            DO 30 is=1,ns
c               dx(5,is)            = adh1(is)*dx(5,is)
c30          CONTINUE  

c
c            IF (istok .EQ. 1) THEN
c

C$DOACROSS LOCAL(iseg,i,j,id1,id2)
               DO 1000 iseg=1,nseg
c
                  i                = nubo(2,iseg)
                  j                = nubo(1,iseg)
c
                  id1              = 50*(iseg - 1)
                  id2              = 50*(iseg - 1) + 25
c
                  dx(1,i)          = dx(1,i) + stmat(id1+1)*dy(1,j)    +
     &                                         stmat(id1+2)*dy(2,j)    +
     &                                         stmat(id1+3)*dy(3,j)    +
     &                                         stmat(id1+4)*dy(4,j)    +
     &                                         stmat(id1+5)*dy(5,j)
c
                  dx(2,i)          = dx(2,i) + adh1(i)*(
     &                                         stmat(id1+5+1)*dy(1,j)  +
     &                                         stmat(id1+5+2)*dy(2,j)  +
     &                                         stmat(id1+5+3)*dy(3,j)  +
     &                                         stmat(id1+5+4)*dy(4,j)  +
     &                                         stmat(id1+5+5)*dy(5,j))
c     
                  dx(3,i)          = dx(3,i) + adh1(i)*(
     &                                         stmat(id1+10+1)*dy(1,j) +
     &                                         stmat(id1+10+2)*dy(2,j) +
     &                                         stmat(id1+10+3)*dy(3,j) +
     &                                         stmat(id1+10+4)*dy(4,j) +
     &                                         stmat(id1+10+5)*dy(5,j)) 
c
                  dx(4,i)          = dx(4,i) + adh1(i)*(
     &                                         stmat(id1+15+1)*dy(1,j) +
     &                                         stmat(id1+15+2)*dy(2,j) +
     &                                         stmat(id1+15+3)*dy(3,j) +
     &                                         stmat(id1+15+4)*dy(4,j) +
     &                                         stmat(id1+15+5)*dy(5,j))
c
                  dx(5,i)          = dx(5,i) + adh1(i)*(
     &                                         stmat(id1+20+1)*dy(1,j) + 
     &                                         stmat(id1+20+2)*dy(2,j) +
     &                                         stmat(id1+20+3)*dy(3,j) +
     &                                         stmat(id1+20+4)*dy(4,j) +
     &                                         stmat(id1+20+5)*dy(5,j))
c
                  dx(1,j)          = dx(1,j) + stmat(id2+1)*dy(1,i)    +
     &                                         stmat(id2+2)*dy(2,i)    + 
     &                                         stmat(id2+3)*dy(3,i)    + 
     &                                         stmat(id2+4)*dy(4,i)    +
     &                                         stmat(id2+5)*dy(5,i)
c
                  dx(2,j)          = dx(2,j) + adh1(j)*(
     &                                         stmat(id2+5+1)*dy(1,i)  + 
     &                                         stmat(id2+5+2)*dy(2,i)  +
     &                                         stmat(id2+5+3)*dy(3,i)  +
     &                                         stmat(id2+5+4)*dy(4,i)  +
     &                                         stmat(id2+5+5)*dy(5,i))
c     
                  dx(3,j)          = dx(3,j) + adh1(j)*(
     &                                         stmat(id2+10+1)*dy(1,i) + 
     &                                         stmat(id2+10+2)*dy(2,i) +
     &                                         stmat(id2+10+3)*dy(3,i) +
     &                                         stmat(id2+10+4)*dy(4,i) +
     &                                         stmat(id2+10+5)*dy(5,i))
c
                  dx(4,j)          = dx(4,j) + adh1(j)*(
     &                                         stmat(id2+15+1)*dy(1,i) +
     &                                         stmat(id2+15+2)*dy(2,i) +
     &                                         stmat(id2+15+3)*dy(3,i) +
     &                                         stmat(id2+15+4)*dy(4,i) +
     &                                         stmat(id2+15+5)*dy(5,i))
c
                  dx(5,j)          = dx(5,j) + adh1(j)*(
     &                                         stmat(id2+20+1)*dy(1,i) +
     &                                         stmat(id2+20+2)*dy(2,i) +
     &                                         stmat(id2+20+3)*dy(3,i) +
     &                                         stmat(id2+20+4)*dy(4,i) +
     &                                         stmat(id2+20+5)*dy(5,i))
c
1000           CONTINUE
c     
c            ENDIF
c
C$DOACROSS LOCAL(is,x1,x2,x3,x4,x5)
            DO 35 is=1,ns
c
               x1                  = diag(is,1,1)*dx(1,is) + 
     &       diag(is,1,2)*dx(2,is) + diag(is,1,3)*dx(3,is) + 
     &       diag(is,1,4)*dx(4,is) + diag(is,1,5)*dx(5,is)
               x2                  = diag(is,2,1)*dx(1,is) +
     &       diag(is,2,2)*dx(2,is) + diag(is,2,3)*dx(3,is) +
     &       diag(is,2,4)*dx(4,is) + diag(is,2,5)*dx(5,is)
               x3                  = diag(is,3,1)*dx(1,is) +
     &       diag(is,3,2)*dx(2,is) + diag(is,3,3)*dx(3,is) +
     &       diag(is,3,4)*dx(4,is) + diag(is,3,5)*dx(5,is)
               x4                  = diag(is,4,1)*dx(1,is) + 
     &       diag(is,4,2)*dx(2,is) + diag(is,4,3)*dx(3,is) + 
     &       diag(is,4,4)*dx(4,is) + diag(is,4,5)*dx(5,is)
               x5                  = diag(is,5,1)*dx(1,is) +
     &       diag(is,5,2)*dx(2,is) + diag(is,5,3)*dx(3,is) +
     &       diag(is,5,4)*dx(4,is) + diag(is,5,5)*dx(5,is)
c
               dx(1,is)            = x1
               dx(2,is)            = x2
               dx(3,is)            = x3
               dx(4,is)            = x4
               dx(5,is)            = x5
c
35          CONTINUE
c
c           Computing the residual 
c
            res                    = 0.0
c
C$DOACROSS LOCAL(is), REDUCTION(res)
            DO 40 is=1,ns
c 
               res                 = res + 
     &                               (dx(1,is) - dy(1,is))**2 + 
     &                               (dx(2,is) - dy(2,is))**2 +
     &                               (dx(3,is) - dy(3,is))**2 + 
     &                               (dx(4,is) - dy(4,is))**2 +
     &                               (dx(5,is) - dy(5,is))**2
c
40          CONTINUE
c
            IF (nit .EQ. 1) res0   = res
c
            res                    = SQRT(res/res0)
c
            IF ((nit .EQ. nbrel) .OR. (res .LT. err)) GOTO 50
c
300      CONTINUE
c
      ENDIF
c
50    CONTINUE
c
c      WRITE(17, 1200) cfl, nit, res
1200  FORMAT(f12.2,2x,i6,2x,e12.5)
c
c      call flunow(17)
      
      RETURN
      END
