          

      SUBROUTINE PLOTCONTROL(CTRL,fich)
C
      INCLUDE 'Paramopt3D.h'
C
      INTEGER ISP,fich
      REAL*8 CTRL(nnsp)
C
         DO ISP = 1,NSP
            write(fich,110) coorp(1,isp),ctrl(isp)
         END DO
C
 110     format(2e15.6)
c
         RETURN
         END
