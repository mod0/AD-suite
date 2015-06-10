C        Generated by TAPENADE     (INRIA, Tropics team)
C  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
C
C  Differentiation of lift_wall in reverse (adjoint) mode:
C   gradient     of useful results: lift q
C   with respect to varying inputs: lift q
C   RW status of diff variables: lift:in-out q:incr
C
C
C
C
C
      SUBROUTINE LIFT_WALL_B(x1, x2, q, qb, lift, liftb)
      IMPLICIT NONE
      REAL*8 gam, gm1, cfl, eps, mach, alpha
C
      COMMON /constants/ gam, gm1, cfl, eps, mach, alpha
C
C
C
      REAL*8 x1(2), x2(2), q(4), lift, dx, ri, u, v, p
      REAL*8 qb(4), liftb, rib, pb
      REAL*8 tempb
C
C
      dx = x1(1) - x2(1)
C
      ri = 1.d0/q(1)
C
C
      pb = dx*liftb
      tempb = -(gm1*0.5d0*pb)
      qb(4) = qb(4) + gm1*pb
      rib = (q(2)**2+q(3)**2)*tempb
      qb(2) = qb(2) + ri*2*q(2)*tempb
      qb(3) = qb(3) + ri*2*q(3)*tempb
      qb(1) = qb(1) - rib/q(1)**2
      END
