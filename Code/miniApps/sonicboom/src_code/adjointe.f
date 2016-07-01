      SUBROUTINE ADJOINTE
C
C----------------------------------------------------------------------
C            Calcul de la transposee de dPSI/dW
C----------------------------------------------------------------------
C
      INCLUDE 'Param3D.h'
C
      INTEGER is, ia, ib, iseg
      REAL*8 aux
C
C**** En ayant resolu l'equation d'etat Psi(W)=0 en implicite, nous avons 
C     utilise une methode de Jacobi ou il a fallu inverser la diagonale
C     Diag.
C     On inverse donc a nouveau les termes diagonaux : 1/DIAG --> DIAG
C     Il faut remarquer que dans Diag, les termes de bord sont inclus.
C
c

      write(6,*) 'inversion1',diag(1330,1,1)
      CALL INVERSION
      write(6,*) 'inversion1',diag(1330,1,1)
C
C**** On enleve maintenant le terme vols/DeltaT a la diagonale 
C AD 30 jan
      print *,'Maintien de Masse dans adjoint a .001'


      Do 1 ia = 1,5
         Do 2 is = 1,ns
c            diag(is,ia,ia) = diag(is,ia,ia)-         vols(is)/dtl(is)

c AD 30 jan 
            diag(is,ia,ia) = diag(is,ia,ia)- 0.999 * vols(is)/dtl(is)
c                                            ******
2        Continue
1     Continue
C
C**** On transpose la matrice dPSI/dW
C
C Termes diagonaux :
C ------------------
C
c         write(6,*) 'queee',diag(3603,1,1), vols(3603)/dtl(3603)

      Do is = 1,ns
         Do ia = 1,5
            Do ib = ia+1,5
               aux = diag(is,ia,ib)
               diag(is,ia,ib) = diag(is,ib,ia)
               diag(is,ib,ia) = aux
            end do
         end do
      end do
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
C     de l'etat-adjoint PI : 
C                              dPSI*       dJ
C                              ----- PI = ----
C                                dW        dW
C



      write(6,*) 'inversion2',diag(1330,1,1)
      CALL INVERSION
      write(6,*) 'inversion2',diag(1330,1,1)
C


      RETURN
      END
