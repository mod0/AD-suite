C     -----------------------------------------------
      SUBROUTINE JACOBI ( NEQUATION, NNdim, NNSGdim, NUBO, RELAX, ZM, 
     $                    DIAG, BVEC, ERRX, X, RES, NITER, IETAT,ctjac,
     $                    nn,nnsg)
C     -----------------------------------------------
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
c                               ZM(1:5,1:5,1,NSG)   ZM(1:5,1:5,2,NSG)
c
C     -----------------------------------------------
C     Debut des declarations
C     -----------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER NN, NNSG, NEQUATION
      integer nndim,nnsgdim
C
C     variables d'appel
C
      INTEGER RELAX, IETAT
      INTEGER NUBO(2,NNSGdim)
      REAL*8 RES, ERRX
      REAL*8 ZM(NEQUATION,NEQUATION,2,NNSGdim)
      REAL*8 DIAG(NNdim,NEQUATION,NEQUATION)
      REAL*8 X(NEQUATION,NNdim), BVEC(NEQUATION,NNdim)
      REAL*4 ctjac,ct1,ct2
C
C     Variables locales
C
C     Tableaux de travail
C
C*** ATTENTION!!!!!!!!!!!!, comme ce sont des tableaux internes a la procedure,
C    il faut les dimensionner au nombre de noeuds exact du maillage
C    sur lequel on travaille.
C
      REAL*8 Y(5,20600), XTMP(5,20600)                      !booemi
c      REAL*8 Y(5,20300), XTMP(5,20300)                      !delwin
c      REAL*8 Y(5,15640), XTMP(5,15640)
c      REAL*8 Y(5,30514), XTMP(5,30514)
c      REAL*8 Y(5,10000), XTMP(5,10000)
C
C     Indices de boucle
      INTEGER IS, K, J, NITER, NSG
C
C     Divers
      INTEGER IA, IB
C
C     Variables memorisees
      REAL*8 RES0
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
c      if (NN.ne.15640)then
c          print*,'Jacobi: NN different de 15640: NN=',NN
c          if(nn.gt.15640)stop
c      endif
c     ***************** INITIALISATION ********************
c                    n+1   -1 
c     On initialise X   a D  f . 
c                               n
c     Dans f, est stocke - PSI(W ), c'est-a-dire ce(5,ns)
c
c      DO 103 K = 1, NEQUATION, 1
c         DO 102 J = 1, NEQUATION, 1
c            DO 101 IS = 1, NN, 1
c               X(K,IS)= X(K,IS) + DIAG(IS,K,J) * BVEC(J,IS)
c 101        CONTINUE
c 102     CONTINUE
c 103  CONTINUE
C
C     Iterations de Jacobi
C
C [llh Nan]      call second(ct1)
      
      DO 121 NITER = 1, RELAX, 1
C
c********************** MISE A JOUR *********************
c                         n
c     Dans dy est stocke X
c     Dans dx, on met f pour commencer a calculer les termes
c     extradiagonaux.
c
cpp         DO 105 K = 1, NEQUATION, 1
cpp            DO 104 IS = 1, NN, 1

         DO IS = 1, NN, 1
            DO K = 1, NEQUATION, 1
               Y(K,IS)= X(K,IS)
               X(K,IS)= BVEC(K,IS)
            END DO
         END DO
C
C        Boucle sur les segments
C


C     
C              Jacobi sur les blocs extra-diagonaux
C

C$DOACROSS LOCAL(NSG,K,J,IA,IB)
         DO NSG = 1, NNSG, 1
            DO J = 1, NEQUATION, 1
               DO K = 1, NEQUATION, 1
                  IA= NUBO(1,NSG)
                  IB= NUBO(2,NSG)
C     
                  X(K,IA)= X(K,IA) - ZM(K,J,2,NSG)* Y(J,IB)
                  X(K,IB)= X(K,IB) - ZM(K,J,1,NSG)* Y(J,IA)
               END DO
            END DO
         END DO
C
C        Boucle sur les noeuds : elements diagonaux
C

C$DOACROSS LOCAL(IS,K)
         DO IS = 1, NN, 1
            DO K = 1, NEQUATION, 1
               XTMP(K,IS) = 0.
            END DO
         END DO
C

C$DOACROSS LOCAL(IS,J,K)
         DO IS = 1 , NN, 1
            DO J = 1, NEQUATION, 1
               DO K = 1, NEQUATION, 1
                  XTMP(K,IS) = XTMP(K,IS) + DIAG(IS,K,J) * X(J,IS)
C
               END DO
            END DO
         END DO
C
C$DOACROSS LOCAL(IS,K)
         DO IS = 1, NN, 1
            DO K = 1, NEQUATION, 1
C     
c        Sous-relaxation de Jacobi de 0.8
c
               X(K,IS) = 0.8*XTMP(K,IS) + 0.2*Y(K,IS)
C
            END DO
         END DO
C
C        Calcul du residu iteratif
C
         RES= 0.
C$DOACROSS LOCAL(IS,K), REDUCTION(RES)
         DO IS = 1, NN, 1
            DO K = 1, NEQUATION, 1
               RES= RES + (X(K,IS) - Y(K,IS))**2
c      if(is.gt.2201)then
c        print*,'Jacobi:        IS=',IS,'  res=',res
c        print*,'Jacobi: diag=',
c     $DIAG(IS,1,1),DIAG(IS,2,2),DIAG(IS,3,3),DIAG(IS,4,4)
c      endif
            END DO
         END DO
C
         RES= SQRT( RES )
c
c     Pour l'etat: normalisation residu
c
         IF(IETAT.EQ.1)THEN
            IF (NITER .EQ. 1) THEN
               RES0= RES
            ENDIF
C
C        Normalisation de l'erreur lineaire
C
            RES = RES/ RES0
         ENDIF
C
c
      print*,'Jacobi: niter=',niter,'  res=',res
c
C        Condition de sortie
C
c         write (6,*) res

         IF (RES .LT. ERRX) THEN
            GOTO 200
         ENDIF
c
c
C
  121 CONTINUE
C
  200 CONTINUE
      print*,'Jacobi: niter=',niter,'  res=',res

C [llh Nan]      call second(ct2)
C [llh Nan]      ctjac=ctjac+ct2-ct1
      
      RETURN
      END
