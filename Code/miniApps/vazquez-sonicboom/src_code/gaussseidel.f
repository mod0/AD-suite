
C     -----------------------------------------------
      SUBROUTINE GAUSSSEIDEL ( NEQU_GS, NN_GS, NNSG_GS, NUBO_GS, 
     $                    RELAX_GS, ZM_GS, 
     $                    DIAG_GS, BVEC_GS, ERRX_GS, X_GS, 
     $                    RES_GS, NITER_GS, IETAT_GS)
C     -----------------------------------------------
c---------------------------------------------------------------------   
c     Resout le systeme lineaire : A X = f (1)
c              Vol   dPSI  n         n+1   n               n
c     avec A = --- - ----(W ) , X = W   - W   et f = -PSI(W )
c              dt     dW          
c
c     On pose A = D - L - U
c                         n+1       n
c     Donc, (1) => (D - L) X    =  U X  + f
c                     n+1    -1     n+1   n
c               =>   X    = D  ( L X + U X  + f )
c
c
c     Les termes extra-diagonaux sont stockes dans STMAT(nsegs) 
c     et les termes diagonaux inverses sont stockes dans 
c     DIAG_GS(ns,5,5).
c
c     Pour un segment donne iseg, il y a donc 2 extremites 
c     NUBO_GS(1,<iseg) et NUBO_GS(2,iseg). La matrice elementaire 
c     STMAT a donc 2 x (5x5) elements.
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
c             => DIAG_GS
c
c   => Termes extradiagonaux : 
c
c                dPHI               dPHI
c                ----(W ,W )  et  - ----(W ,W )
c                dW   i  j          dW   j  i
c                  R                  L
c                     ||                 ||
c                     \/                 \/
c            ZM_GS(1:5,1:5,1,NSG)   ZM_GS(1:5,1:5,2,NSG)
c
C     -----------------------------------------------
C     Debut des declarations
C     -----------------------------------------------
C
c---------------------------------------------------------------------
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------
c
      INTEGER NN_GS, NNSG_GS, NEQU_GS
C
C     variables d'appel
C
      INTEGER RELAX_GS, IETAT_GS
      INTEGER NUBO_GS(2,NNSG_GS)
      REAL*8 RES_GS, ERRX_GS
      REAL*8 ZM_GS(NEQU_GS,NEQU_GS,2,NNSG_GS)
      REAL*8 DIAG_GS(NN_GS,NEQU_GS,NEQU_GS)
      REAL*8 X_GS(NEQU_GS,NN_GS), BVEC_GS(NEQU_GS,NN_GS)
C
C     Variables locales
C
C     Tableaux de travail
C
C*** ATTENTION!!!!!!!!!!, comme ce sont des tableaux internes a la 
C    procedure, il faut les dimensionner au nombre de noeuds 
C    exact du maillage sur lequel on travaille.
C
c      REAL*8 Y_NS(5,15640), XTMP(5,15640)
c      REAL*8 Y_NS(5,30514), XTMP(5,30514)
      REAL*8 Y_NS(5,4151), XTMP(5,4151)

C
C     Indices de boucle
      INTEGER IS, K, J, NITER_GS, NSG, KS
C
C     Divers
      INTEGER isgs, ia, ib, icc, idd, nsg2, incgs, iseg
C
C     Variables memorisees
      REAL*8 RES0
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
      if (NN_GS.ne.15640)print*,'GaussSeidel: ns pas bon'
c     ***************** INITIALISATION ********************
c                    n+1   -1 
c     On initialise X   a D  f . 
c                               n
c     Dans f, est stocke - PSI(W ), c'est-a-dire ce(5,ns)
c
c      DO 103 K = 1, NEQU_GS, 1
c         DO 102 J = 1, NEQU_GS, 1
c            DO 101 IS = 1, NN_GS, 1
c               X_GS(K,IS)= 
c     $   X_GS(K,IS) + DIAG_GS(IS,K,J) * BVEC_GS(J,IS)
c 101        CONTINUE
c 102     CONTINUE
c 103  CONTINUE
C
C     Iterations de Gauss-Seidel
C
      DO 121 NITER_GS = 1, RELAX_GS, 1
C
c********************** MISE A JOUR *********************
c                         n
c     Dans dy est stocke X
c     Dans dx, on met f pour commencer a calculer les termes
c     extradiagonaux.
c
         DO 105 K = 1, NEQU_GS, 1
            DO 104 IS = 1, NN_GS, 1
               Y_NS(K,IS)= X_GS(K,IS)
               X_GS(K,IS)= BVEC_GS(K,IS)
 104        CONTINUE
 105     CONTINUE
c
c   Boucle sur les noeuds
c
            IF (MOD(nit, 2) .EQ. 0) THEN
               isgs                = NN_GS + 1
               incgs               =-1
            ELSE
               isgs                = 0
               incgs               = 1
            ENDIF
c
            DO 200 KS=1,NN_GS
c           ==============
c
               isgs                = isgs + incgs
               nsg2                = ndeg(isgs)
C
C
C        Boucle sur les segments voisins de isgs
C
               DO 106 NSG = 1, NSG2
C
c
                  iseg             = jaret(isgs,nsg)
c
                  icc              = 0
c
                  IF (isgs .EQ. NUBO_GS(1,iseg)) THEN
                     icc           = 2
                     idd           = 2
                  ELSE
                     icc           = 1
                     idd           = 1
                  ENDIF
c
                  ia                = nubo(3-icc,iseg)
                  ib                = nubo(icc,iseg)
c
                  DO K = 1, NEQU_GS, 1
                  DO  J = 1, NEQU_GS, 1
                      X_GS(K,IA)= 
     $                X_GS(K,IA) - 
     $                ZM_GS(K,J,idd,ISEG)* Y_NS(J,IB)
c
                  ENDDO
                  ENDDO
c
  106          CONTINUE
C
C     Multiplication par inverse de D;
C
C
         IS=ISGS
c
         DO K = 1, NEQU_GS, 1
              XTMP(K,IS) = 0.
              DO J = 1, NEQU_GS, 1
C
c
                  XTMP(K,IS) = 
     $            XTMP(K,IS) + 
     $             DIAG_GS(IS,K,J) * X_GS(J,IS)
C
              ENDDO
C
c        Sous-relaxation de 0.99
c
              X_GS(K,IS) = 0.99*XTMP(K,IS) + 0.01*Y_NS(K,IS)
         ENDDO
C
C  FIN BOUCLE SOMMET
C
  200    CONTINUE
c ===============
C
C        Calcul du residu iteratif
C
         RES_GS= 0.
         DO 120 K = 1, NEQU_GS, 1
            DO 119 IS = 1, NN_GS, 1
               RES_GS = RES_GS + 
     $                 (X_GS(K,IS) - Y_NS(K,IS))**2
  119       CONTINUE
  120    CONTINUE
C
         RES_GS= SQRT( RES_GS )
c
c     Pour l'etat: normalisation residu
c
         IF(IETAT_GS.EQ.1)THEN
            IF (NITER_GS .EQ. 1) THEN
               RES0= RES_GS
            ENDIF
C
C        Normalisation de l'erreur lineaire
C
            RES_GS = RES_GS/ RES0
         ENDIF
         write(6,*)'GaussSeidel: residu=',RES_GS
C
C        Condition de sortie
C
         IF (RES_GS .LT. ERRX_GS) THEN
            GOTO 300
         ENDIF
C
  121 CONTINUE
C
  300 CONTINUE
      RETURN
      END
