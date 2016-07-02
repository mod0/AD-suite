      SUBROUTINE FCNELJ(CTRL,RHO,FJ,fic)
C
c    Cette procedure est appelee par la procedure Rhoopt.f dans le but
c    de trouver le pas optimal de descente a chaque iteration d'optimisation
c
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER ISP,fic
      REAL*8 CTRL1(NNSP), CTRL(NNSP)
      REAL*8 FJ,RHO,CCD,CCL
C
      DO ISP = 1,NSP
        CTRL1(ISP) = CTRL(ISP) - RHO * GRAD(ISP)
      END DO
c
c*** Ecriture dans les fichiers fort.300....de la coque a chaque recherche
c    d'un pas de descente.
c

      
      PRINT*,'FcnelJ.f: NUNIT=',fic
      CALL ECRITURE(CTRL1,fic)
c
      fic = fic + 1
C
      CALL ETAT(CTRL1,ktmax)
c
      PRINT*,'Longueur du pas = ',rho
C
c   Modification le 29.05.95
c   ------------------------
      CALL PORTANCE(CTRL1,CCD,CCL)
c
      CALL FNELLECOUT(CCD,CCL)
c
c      CALL FONCCOUT
c
      ncf = ncf + 1
c
      PRINT*,'Fonctionnelle cout = ',cout,' CD=',CCD,' CL =',CCL
      print*,'Iteration optimisation no : ',itopt
c
      FJ = COUT
c
      RETURN
      END
