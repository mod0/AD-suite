

      SUBROUTINE TESTETAT(CTRL)
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER IS,ISP,ISP0,ISPK,K,is1
      REAL*8 EPS, UAGAM(5,NSMAX), DPSIDGAM(5,NSMAX), WPRIM(5,NSMAX)
      REAL*8 CEBAR(5,NSMAX), CTRL(NNSP)
      REAL*8 DIFF, DIFF1, res,err1, dmax
      INTEGER iter, ktmax0, ia, ib, IETAT
      REAL*4 ctjac

C**** 1ere etape : Resolution de l'etat
c                  On obtient W(Gamma)
c                  On part d'une solution deja convergee et on ne fait
c                  qu'une seule iteration pour calculer la matrice implicite
c                  stockee dans zm et diag.
c
      ktmax0 = ktmax
      ktmax = kt0+5
      dmax = 0.
c
      CALL ETAT(CTRL,ktmax)
c
      ktmax = ktmax0
c
      DO IS = 1,NS
         DO K = 1,5
            UAGAM(K,IS) = UA(K,IS)
         END DO
      END DO
C
C**** 2eme etape : Calcul de la derivee de Psi par rapport au controle 
C     gamma (informatique : ctrl(1:nnsp) ) par differences divisees.
C
C     dPsi    Psi(gamma + epsilon) - Psi(gamma)
C    ------ = ---------------------------------
C    dgamma               epsilon
C
C    Attention : dans la variable CE est stockee -Psi.
c
      EPS = 1.D-07
c

C
C
c*** Calcul de Psi(Gamma)
c
      DO 20 IS = 1,NS
         DO 30 K = 1,5
            CE(K,IS) = 0.
 30      CONTINUE
 20   CONTINUE

      CALL TRANSPIRATION(CE,CTRL)
C
      ISP = 1
C
         CTRL(ISP) = CTRL(ISP) + EPS
C
c*** Calcul de Psi(Gamma + deltaGamma)
c
         DO IS = 1,NS
            DO K = 1,5
               CEBAR(K,IS) = 0.
            END DO
         END DO
C
         CALL TRANSPIRATION(CEBAR,CTRL)
C
C*** Calcul de la difference (Psi(gamma+detaGamma)-Psi(Gamma))/epsilon
c
         DO IS = 1,NS
            DO K = 1,5
               DPSIDGAM(K,IS) = 0.
               WPRIM(K,IS) = 0.
            END DO
         END DO
C
         DO ISP0 = 1,NSP
             IS = NODE2D3D(ISP0)
            DO K = 1,5
               DPSIDGAM(K,IS) = (CE(K,IS)-CEBAR(K,IS))/EPS
            END DO
         END DO
C      
c*** 3eme etape : Calcul de (dPsi/dW)   W'  = dpsi/dgamma
c                 avec Jacobi
c
      CALL INVERSION
C
C**** On enleve maintenant le terme vols/DeltaT a la diagonale 
C
      Do 1 ia = 1,5
         Do 2 is = 1,ns
            diag(is,ia,ia) = diag(is,ia,ia)-vols(is)/dtl(is)
2        Continue
1     Continue
c
      CALL INVERSION
c
      IETAT = 1
c
         CALL JACOBI(NEQUATION,nsmax,nsgmax,nubo,nbrelPi,zm,diag,
     $     DPSIDGAM, errjacPi,WPRIM,res,iter,IETAT,ctjac,ns,nseg)
c
             open(19,access='append')
              write(19,180)res,iter
              close(19)
180   FORMAT(e14.7,1x,i5)
C
C*** 4eme etape : Resolution de l'etat avec W(gamma+epsilon)
c
      CALL ETAT(CTRL,ktmax)
C
      DO ISPK = 1,NSP
         IS = NODE2D3D(ISPK)
         DO K = 1,5
            DIFF = -(UA(K,IS) - UAGAM(K,IS))/EPS
            DIFF1 = ABS(DIFF - WPRIM(K,IS))
            write(71,*) DIFF1, DIFF, WPRIM(K,IS)
            dmax = max(dmax,diff1)
         END DO
      END DO
c
      print*,'max = ',dmax
c
         CTRL(ISP) = CTRL(ISP) - EPS
C
      RETURN
      END
