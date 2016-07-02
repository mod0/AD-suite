      SUBROUTINE CDLNS
c     -----------------------------------------------------------------
c     No-slip boundary condition in the case of viscous flows
c     -----------------------------------------------------------------
      INCLUDE 'Param3D.h'
c     -----------------------------------------------------------------
c     Local variables definition
      INTEGER is1, is2, is3, if11, ifac
c
      DO 100 if11=1,nf11
c
         ifac                      = noe1(if11)
c
         is1                       = nsfac(1,ifac)
         is2                       = nsfac(2,ifac)
         is3                       = nsfac(3,ifac)
c
         ce(2,is1)                 = 0.0
         ce(3,is1)                 = 0.0
         ce(4,is1)                 = 0.0
         ce(5,is1)                 = ce(1,is1)*tbrd
c
         ce(2,is2)                 = 0.0
         ce(3,is2)                 = 0.0
         ce(4,is2)                 = 0.0
         ce(5,is2)                 = ce(1,is2)*tbrd
c
         ce(2,is3)                 = 0.0
         ce(3,is3)                 = 0.0
         ce(4,is3)                 = 0.0
         ce(5,is3)                 = ce(1,is3)*tbrd
c
100   CONTINUE
c
      RETURN
      END
