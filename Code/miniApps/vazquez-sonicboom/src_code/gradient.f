      SUBROUTINE GRADIENT(CTRL)
C
C**** 1ere etape : Calcul de la derivee de Psi par rapport au controle 
C     gamma (informatique : ctrl(1:nnsp) ) par differences divisees.
C
C     dPsi    Psi(gamma + epsilon) - Psi(gamma)
C    ------ = ---------------------------------
C    dgamma               epsilon
C
C    Attention : dans la variable CE est stockee -Psi.
C
C    2eme etape : on en deduit le Gradient stocke dans Grad(1:nnsp) :
c                 Grad = dJ/dgamma - < Pi,dPsi/dgamma >
c
C*** Entree : CTRL (Gamma)
c
c    Calculs internes : ce --> Psi(gamma)
c                       cebar --> Psi(gamma+eps)
c
C    Sortie : Grad (Gradient)
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER IS, ISP, K, ISPK, I, ISP0, NUGRAD
      REAL*8 EPS, CEBAR(5,NSMAX), NORM, CTRL(NNSP)
C [llh Nan]      REAL*4 t00, t1
C
      EPS = 1.
c
C [llh Nan]      CALL SECOND(t00)
C
c *** Calcul de Psi(gamma)
c
         DO ISP0 = 1,NSP
             IS = NODE2D3D(ISP0)
            DO K = 1,5
               CE(K,IS) = 0.
            END DO
         END DO
C
      CALL TRANSPIRATION(CE,CTRL)
C
      DO 100 ISP = 1,NSP

         GRAD(ISP) = 0.
C
c*** Calcul de Psi(Gamma + eps)
c
         CTRL(ISP) = CTRL(ISP) + EPS
C
         DO ISP0 = 1,NSP
             IS = NODE2D3D(ISP0)
            DO K = 1,5
               CEBAR(K,IS) = 0.
            END DO
         END DO
C
         CALL TRANSPIRATION(CEBAR,CTRL)
C
c*** Calcul du gradient :
c
c               isp
c               --- ___            dPsi
c   Grad(isp) = \   || (1:5,is) * ------(1:5,is)
c               /                 dgamma
c               ---
c             ispk=1
c        is=node2d3d(ispk)
c
c           dPsi            ce(1:5,is) - cebar(1:5,is)
c   avec : ------ (1:5,is)= --------------------------
c          dgamma                    epsilon


         DO 60 ISPK = 1,NSP
            IS = NODE2D3D(ISPK)
            DO 70 K = 1,5
               GRAD(ISP) = GRAD(ISP) - PIADJ(K,IS)*
     $                            (CE(K,IS)-CEBAR(K,IS))/EPS
 70         CONTINUE
 60      CONTINUE
C
c  Cas ou DJ/DGamma est non nul :
c
         GRAD(ISP) = GRAD(ISP) + DJDGAM(ISP)
c
         CTRL(ISP) = CTRL(ISP) - EPS
c
 100  CONTINUE
C
      NORM = 0.
      DO 200 ISP = 1,NSP
         NORM = NORM + GRAD(ISP)**2
 200  CONTINUE
C
      OPEN(8,ACCESS='APPEND')
      WRITE(8,300)ITOPT,SQRT(NORM)
 300  FORMAT(I5,E14.7,I5)
      CLOSE(8)
c
c*** Ecriture du gradient dans les fichiers fort.700.....
c
      nugrad=100+itopt
      nugrad=36
      print*,'Debut ecriture du gradient sur unite=',nugrad
      OPEN(unit=36,access='sequential',form='unformatted')
      write (36) itopt, nsp, (grad(isp),isp=1,nsp)

c      OPEN(unit=36,status='unknown')
c      do isp=1,nsp
c         write (36,*) isp,node2d3d(isp),grad(isp)
c      end do
c      close(nugrad)

      print*,'gradient ecrit sur fort:',nugrad
c
c
      
      reazz(1)=SQRT(NORM) 
      
      WRITE(6, *) ' '
      print *,'  Norme l2 du gradient: ',SQRT(NORM)
      WRITE(6, *) ' '
C
C [llh Nan]      CALL SECOND(t1)
C
C [llh Nan]      WRITE(6,*) 'Temps effectue pour calculer le gradient :',t1-t00
c
      RETURN
      END
