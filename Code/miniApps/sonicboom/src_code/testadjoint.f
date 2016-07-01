




      SUBROUTINE TESTADJOINT
C
      INCLUDE 'Param3D.h'
C
      INTEGER IA,IB,IS,KVAR,JVAR,NSG,ISEG,ITER
      REAL*8 aux, S(5,NSMAX)
      REAL*8 V(5,NSMAX),DJDWTEST(5,NSMAX),PRODSCAL
      REAL*8 ADJ(5,NSMAX),RES
      INTEGER IETAT
      REAL*4 ctjac
C
C**** En ayant resolu l'equation d'etat Psi(W)=0 en implicite, nous avons 
C     utilise une methode de Jacobi ou il a fallu inverser la diagonale
C     Diag.
C     On inverse donc a nouveau les termes diagonaux : 1/DIAG --> DIAG
C     Il faut remarquer que dans Diag, les termes de bord sont inclus.
C
      CALL INVERSION
C
C**** On enleve maintenant le terme vols/DeltaT a la diagonale 
C
      Do 1 ia = 1,5
         Do 2 is = 1,ns
            diag(is,ia,ia) = diag(is,ia,ia)-vols(is)/dtl(is)
2        Continue
1     Continue
C
C*** INITIALISATION
C
      DO IS = 1,NS
         DO KVAR = 1,5
            DJDWTEST(KVAR,IS) = 1.
            S(KVAR,IS) = 1.
            V(KVAR,IS) = 0.
         END DO
      END DO
C
C-------------------------------------------------------
C                           DPSI
C               CALCUL DE   ---- . S = V
C                            DW
C                        
C-------------------------------------------------------
C
C*** TERMES EXTRA-DIAGONAUX
C
      DO KVAR = 1,5
         DO JVAR = 1,5
            DO NSG = 1,NSEG
               IA = NUBO(1,NSG)
               IB = NUBO(2,NSG)
               V(KVAR,IA) = V(KVAR,IA) + ZM(KVAR,JVAR,2,NSG)*S(JVAR,IB)
               V(KVAR,IB) = V(KVAR,IB) + ZM(KVAR,JVAR,1,NSG)*S(JVAR,IA)
            END DO
         END DO
      END DO
C
C*** TERMES DIAGONAUX
C
      DO IS = 1,NS
         DO KVAR = 1,5
            DO JVAR = 1,5
               V(KVAR,IS) = V(KVAR,IS) + DIAG(IS,KVAR,JVAR)*S(JVAR,IS)
            END DO
         END DO
      END DO
C
C-------------------------------------------------------
C                            DPSI*
C                CALCUL DE  ------ . ADJ  = DJDWTEST
C                             DW
C-------------------------------------------------------
C
C**** On transpose la matrice dPSI/dW
C
C Termes diagonaux :
C ------------------
C
      Do 3 is = 1,ns
         Do 4 ia = 1,5
            Do 5 ib = ia+1,5
               aux = diag(is,ia,ib)
               diag(is,ia,ib) = diag(is,ib,ia)
               diag(is,ib,ia) = aux
5           Continue
4        Continue
3     Continue
C
C Termes extra-diagonaux :
C ------------------------
C
      Do 6 iseg = 1,nseg
         Do 7 ia = 1,5
            Do 8 ib = 1,5
               aux = zm(ia,ib,1,iseg)
               zm(ia,ib,1,iseg) = zm(ib,ia,2,iseg)
               zm(ib,ia,2,iseg) = aux
8           Continue
7        Continue
6     Continue
C
C**** On inverse a nouveau les termes diagonaux pour la resolution par Jacobi
C     de ADJ 
C
      CALL INVERSION
c
      IETAT = 0
C
      CALL JACOBI(nequation,nsmax,nsgmax,nubo,nbrelPi,zm,diag,djdwtest,
     $     errjacPi,Adj,res,iter,IETAT,ctjac,ns,nseg)
C
              open(19,access='append')
              write(19,180)res,iter
              close(19)
c
180   FORMAT(e14.7,1x,i5)
C
      PRODSCAL = 0.
      DO IS = 1,NS
         DO KVAR = 1,5
            PRODSCAL = PRODSCAL + ADJ(KVAR,IS)*V(KVAR,IS)
         END DO
      END DO
C
      PRINT*,'LE PRODUIT SCALAIRE EST EGAL A :',PRODSCAL
C
      RETURN
      END
