      SUBROUTINE FORMMAGIQUE(CTRL,COUTINIT)
C
c     Test de verification du calcul du gradient GRAD(isp) par la methode
c     des differences divisees.
c
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      REAL*8 EPS,COUTINIT,DIFF
      REAL*8 DIFF1, CTRL(NNSP), CCD, CCL
      INTEGER ISP,KVAR,ispk,iteps
C
      WRITE(6, *) '************************************************* '
      WRITE(6, *) '*  VERIFICATION DU GRADIENT PAR LES DIFFERENCES * '
      WRITE(6, *) '*  DIVISEES                                     * '
      WRITE(6, *) '************************************************* '
      WRITE(6, *) ' '
c
      EPS = 1.E-04
c
      do iteps = 1,3

C
c      DO ISP = 1,6
c
      ISP = 1
C
         CTRL(ISP) = CTRL(ISP) + EPS
C
         CALL ETAT(CTRL,ktmax)
c
c         CALL FONCCOUT
c
         CALL PORTANCE(CTRL,CCD,CCL)
C
         CALL FNELLECOUT(CCD,CCL)
C
         DIFF = (COUT - COUTINIT)/EPS
C
         WRITE(6,*) 'epsilon = ',eps
         WRITE(6, *) 'Le gradient vaut grad(',isp,') = ',grad(isp)
         WRITE(6, *) 'Les differences divisees valent :',DIFF
c
         DIFF1 = DIFF-GRAD(ISP)
c
         WRITE(6, *) 'La difference vaut: ',DIFF1
         WRITE(6, *) ' '
C
         CTRL(ISP) = CTRL(ISP) - EPS
C
c      END DO
c
         eps = eps/10.
c
         end do
C
      COUT = COUTINIT
c
      RETURN
      END
