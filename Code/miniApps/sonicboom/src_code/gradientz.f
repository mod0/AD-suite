

      SUBROUTINE GRADIENTZ(KNIV,CTRL)
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
c                 Grad = - < Pi,dPsi/dgamma >
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
      REAL*8 t0, t1
      INTEGER IZONE, KNIV
C
      EPS = 1.D-06
C
      CALL SECOND(t0)
c
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
      DO 100 IZONE = 1,NBZ(KNIV)
C
         GRADZ(IZONE) = 0.
C
c*** Calcul de Psi(Gamma + eps)
c
         DELTAZ(IZONE) = DELTAZ(IZONE) + EPS
C
         CALL LP(KNIV,DELTAZ,DELTA)
C
         DO ISP = 1,NSP
C
            CTRL(ISP) = CTRL(ISP) + DELTA(ISP)
C
         END DO
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
               GRADZ(IZONE) = GRADZ(IZONE) - PIADJ(K,IS)*
     $                            (CE(K,IS)-CEBAR(K,IS))/EPS
 70         CONTINUE
 60      CONTINUE
C
         GRADZ(IZONE) = GRADZ(IZONE) + DJDGAMZ(IZONE)
c
         DELTAZ(IZONE) = DELTAZ(IZONE) - EPS
c
         DO ISP = 1,NSP
C
            CTRL(ISP) = CTRL(ISP) - DELTA(ISP)    
C
         END DO 
C
 100  CONTINUE
C
      CALL LP(KNIV,GRADZ,GRAD)
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
      nugrad=400+itopt
      write (nugrad), itopt, nsp, (grad(isp),isp=1,nsp)
      print*,'gradient ecrit sur fort.',nugrad
c
      WRITE(6, *) ' '
      print *,'  Norme l2 du gradient: ',SQRT(NORM)
      WRITE(6, *) ' '
C
      CALL SECOND(t1)
C
      WRITE(6,*) 'Temps effectue pour calculer le gradient :',t1-t0
c
      RETURN
      END
