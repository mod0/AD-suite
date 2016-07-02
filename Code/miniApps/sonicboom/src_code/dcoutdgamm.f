      SUBROUTINE DCOUTDGAMM(CTRL)
C
C**** Calcul de la derivee de J par rapport au controle 
C     gamma (informatique : ctrl(1:nnsp) ) par differences divisees.
C
C      dJ      J(gamma + epsilon) - J(gamma)
C    ------ = ---------------------------------
C    dgamma               epsilon
C
C*** Entree : CTRL (Gamma)
c
C    Sortie : DJDGAM 
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER IS, ISP, K, ISPK, I, ISP0
      REAL*8 EPS, CTRL(NNSP), COUTINIT
      REAL*8 CDP, CLP
C [llh Nan]      REAL*4 ti, t1
C
      EPS = 1.D-06
c
      write(6,*) 'Entree dans DCoutdGamm.'
      call flunow(6)
c
c
C [llh Nan]      CALL SECOND(ti)
C
c *** Calcul de J(gamma)
c
      CALL PORTANCE(CTRL,CDP,CLP)
c
      CALL FNELLECOUT(CDP,CLP)
C
      COUTINIT = COUT
C
      DO 100 ISP = 1,NSP
C
c*** Calcul de J(Gamma + eps)
c
         CTRL(ISP) = CTRL(ISP) + EPS
C
         CALL PORTANCE(CTRL,CDP,CLP)
c
         CALL FNELLECOUT(CDP,CLP)
C
         DJDGAM(ISP) = (COUT-COUTINIT)/EPS
C
         CTRL(ISP) = CTRL(ISP) - EPS         
C
 100  CONTINUE
C
C [llh Nan]      CALL SECOND(t1)
C
C [llh Nan]      WRITE(6,*) 'Temps effectue pour calculer dJ/dGam :'
C [llh Nan]      WRITE(6,*) t1-ti
C
      WRITE(6,*) 'Sortie de DCoutDGamm.'
c
      RETURN
      END
