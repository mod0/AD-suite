C        Generated by TAPENADE     (INRIA, Tropics team)
C  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
C
C  Differentiation of flux_wall in reverse (adjoint) mode:
C   gradient     of useful results: res q x1 x2
C   with respect to varying inputs: res q x1 x2
C   RW status of diff variables: res:in-out q:incr x1:incr x2:incr
C
C
C
C
C
      SUBROUTINE FLUX_WALL_BX(x1, x1b, x2, x2b, q, qb, res, resb)
      IMPLICIT NONE
      REAL*8 gam, gm1, cfl, eps, mach, alpha
C
      COMMON /constants/ gam, gm1, cfl, eps, mach, alpha
C
C
      INTEGER n
C
      REAL*8 x1(2), x2(2), q(4), res(4), dx, dy, ri, u, v, p
      REAL*8 x1b(2), x2b(2), qb(4), resb(4), dxb, dyb, rib, pb
      REAL*8 tempb
C
C
      dx = x1(1) - x2(1)
      dy = x1(2) - x2(2)
C
      ri = 1.d0/q(1)
      p = gm1*(q(4)-0.5d0*ri*(q(2)**2+q(3)**2))
C
C
      pb = dx*resb(3) - dy*resb(2)
      dxb = p*resb(3)
      dyb = -(p*resb(2))
      tempb = -(gm1*0.5d0*pb)
      qb(4) = qb(4) + gm1*pb
      rib = (q(2)**2+q(3)**2)*tempb
      qb(2) = qb(2) + ri*2*q(2)*tempb
      qb(3) = qb(3) + ri*2*q(3)*tempb
      qb(1) = qb(1) - rib/q(1)**2
      x1b(2) = x1b(2) + dyb
      x2b(2) = x2b(2) - dyb
      x1b(1) = x1b(1) + dxb
      x2b(1) = x2b(1) - dxb
      END
