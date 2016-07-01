C     -----------------------------------------------
      SUBROUTINE DFDECENTRE( NEQUATION, RHO, UN, VN, WN, E, GAMMA, GAM1, 
     &                       GAM4, GAM5, GAM6, XMM, XMP, C, DFD )
C     -----------------------------------------------
C	CETTE PROCEDURE CALCULE LA DERIVEE DU FLUX DECENTRE
C	 + _   
C       F (W) PAR RAPPORT A L'ETAT W.
C                + _      _     - _
C     A NOTER : F (W) = F(W) - F (W)
C 
C       _                             2   +                            _
C      |  Rho*C/4 * (rho*UN/rho*C + 1) = f   = F1                      |
C      |   +                              1                            |
C      |  f  / gamma * ( (gamma-1)*rho*UN/rho + 2C ) = F2              | 
C      |   1         +                                                 |
C  +_  |            f  * rho*VN.rho = F3                               |
C F(W)=|             1                                                 |
C      |             +                                                 |
C      |            f  * rho*WN/rho = F4                               |
C      |   *         1                        2              2         |
C      |  f  * [ ( (gamma-1)*rho*UN/rho + 2C ) / 2*(gamma -1)          |
C      |   1                          2          2       2             |
C      |                   + ((rho*VN) + (rho*WN) )/2*rho ] = F5       |
C      |                                                               |
C      |_                                                             _|
C
C       DANS LE CAS OU  |UN/C| < 1
C
C                             1/2
C       C = ( gamma * P / Rho )  
C                               1            2          2          2     
C       P = (gamma - 1) * ( E - -* ( (rho*UN) + (rho*VN) + (rho*WN) )/rho )
C                               2 
C
C     XMM = 0.5 * (UN/C - 1)
C     XMP = 0.5 * (UN/C + 1)
C                                  +
C     Derivation du flux decentre F par rapport aux variables conservatives :
C                       (Rho, Rho*UN, Rho*VN, Rho*WN, E)
C
C     -----------------------------------------------
C     Debut des declarations
C     -----------------------------------------------
C
      IMPLICIT NONE
C
C     variables d'appel
      INTEGER NEQUATION
      REAL*8 UN, VN, WN, RHO, C, E, XMM, XMP
      REAL*8 GAM1, GAM4, GAMMA, GAM5, GAM6
      REAL*8 DFD(NEQUATION,NEQUATION)
C
C     Divers
      REAL*8 RHOC, AX, VNWN2, F1, F5
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
C     CONSTANTES DU CALCUL
C
C     GAMMA = 1.4 , GAM1 = GAMMA - 1 , GAM5 = 1/GAMMA
C     GAM4 = GAMMA*GAM1, GAM6 = 1./( GAMMA**2 - 1. )
C
      RHOC   = RHO* C
C     
      AX    =  GAM1* UN + 2.*C
      VNWN2 =  VN**2 + WN**2
C
      F1 = RHO*C/4. * (UN/C + 1.)**2
C
      F5 = (AX**2 * GAM6 + VNWN2)/2.
C
C     Calcul des derivees de  F+
C
c**** Derivees de F1 par rapport a (rho, rhoU, rhoV, rhoW, E)
c
      DFD(1,1)= - GAM4* E * XMM * XMP/ RHOC/ 2.
c
      DFD(1,2)= XMP/4. * ( GAM4*UN*2*XMM/C + 4. )
c
      DFD(1,3)= GAM4/2. * XMP * VN/C * XMM
c
      DFD(1,4)= GAM4/2. * WN/C * XMM * XMP
c
      DFD(1,5)= - GAM4/2./C * XMM * XMP
C    
c**** Derivees de F2 par rapport a (rho, rhoU, rhoV, rhoW, E)
c
      DFD(2,1)= GAM5*AX*DFD(1,1) + F1*GAM5/RHO * (GAM4*E/RHOC - AX)
c
      DFD(2,2)= GAM5*AX*DFD(1,2) + GAM1*F1*GAM5/RHO * (1.-GAMMA*UN/C)
c
      DFD(2,3)= GAM5*AX*DFD(1,3) - GAM1*VN*F1/RHOC
c
      DFD(2,4)= GAM5*AX*DFD(1,4) - GAM1*WN*F1/RHOC
c
      DFD(2,5)= GAM5*AX*DFD(1,5) + GAM1*F1/RHOC
C   
c**** Derivees de F3 par rapport a (rho, rhoU, rhoV, rhoW, E)
c  
      DFD(3,1)= VN*(DFD(1,1) - F1/RHO)
c
      DFD(3,2)= VN*DFD(1,2)
c
      DFD(3,3)= F1/RHO + VN*DFD(1,3)
c
      DFD(3,4)= VN*DFD(1,4)
c
      DFD(3,5)= VN*DFD(1,5)
C     
c**** Derivees de F4 par rapport a (rho, rhoU, rhoV, rhoW, E)
c
      DFD(4,1)= WN*(DFD(1,1) - F1/RHO)
c
      DFD(4,2)= WN*DFD(1,2)
c
      DFD(4,3)= WN*DFD(1,3)
c
      DFD(4,4)= F1/RHO + WN*DFD(1,4)
c
      DFD(4,5)= WN*DFD(1,5)
C     
C**** Derivees de F5 par rapport a (rho, rhoU, rhoV, rhoW, E)
c
      DFD(5,1)= F5*DFD(1,1) + F1/RHO *
     &                    ((AX*(GAM4*E/RHOC - AX)*GAM6)-VNWN2)
c
      DFD(5,2)= F5*DFD(1,2) + GAM1*GAM6*F1/RHO * AX * (1.- GAMMA*UN/C)
c
      DFD(5,3)= F5*DFD(1,3) + F1*VN/RHO * (1. - AX*GAM4*GAM6/C)
c
      DFD(5,4)= F5*DFD(1,4) + F1*WN/RHO * (1. - AX*GAM4*GAM6/C)
c
      DFD(5,5)= F5*DFD(1,5) + GAM4*GAM6*F1*AX/RHOC
C
      RETURN
      END
