      SUBROUTINE AIRPEAUINIT
C
      Include 'Param3D.h'
      Include 'Paramopt3D.h'
C
C*** Calcul des aires de la coque initiale : stockees dans airespo
c    et utilisees dans le calcul de la fonctionnelle cout.
c
      Integer isp
c
         DO ISP = 1,NSP
            AIRESP0(isp) = AIRESP(ISP)
         END DO
c
         return
         end
