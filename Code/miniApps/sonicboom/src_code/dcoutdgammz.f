

      SUBROUTINE DCOUTDGAMMZ(KNIV,CTRL)
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
      INTEGER KNIV, IZONE
      REAL*8 EPS, CTRL(NNSP), COUTINIT
      REAL*8 CDP, CLP, t00, t1
C
      WRITE(6,*) 'Entree dans dCoutdGamm.'
c
      EPS = 1.D-06
C
c *** Calcul de J(gamma)
c
      CALL SECOND(t00)
C
      CALL PORTANCE(CTRL,CDP,CLP)
c
      CALL FNELLECOUT(CDP,CLP)
C
      COUTINIT = COUT
C
      DO IZONE = 1,NBZ(KNIV)
         DELTAZ(IZONE) = 0.
      END DO
c
      DO 100 IZONE = 1,NBZ(KNIV)
C
c*** Calcul de J(Gamma + eps)
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
         CALL PORTANCE(CTRL,CDP,CLP)
c
         CALL FNELLECOUT(CDP,CLP)
C
         DJDGAMZ(IZONE) = (COUT-COUTINIT)/EPS
C
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
      CALL SECOND(t1)
C
      WRITE(6,*) 'Temps effectue pour calculer dJ/dGam :',t1-t00
C
      WRITE(6,*) 'Sortie de DCoutDGamm.'
      RETURN
      END
