C        Generated by TAPENADE     (INRIA, Tropics team)
C  Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
C
C  Differentiation of time_cell in forward (tangent) mode:
C   variations   of useful results: adt
C   with respect to varying inputs: q
C   RW status of diff variables: q:in adt:out
C
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                                                                   c
C      Nonlinear routines used in airfoil calculations              c
C                                                                   c
C      linear/adjoint versions are created on-the-fly by Tapenade   c
C      complex version is generated using compiler flag -DCOMPLEX   c
C                                                                   c
C      Copyright Devendra Ghate and Mike Giles, 2005                c
C      but can be freely used with due acknowledgement              c
C                                                                   c
C      Tapenade developed by Laurent Hascoet and others at INRIA    c
C                                                                   c
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C
C
      SUBROUTINE TIME_CELL_D(x1, x2, x3, x4, q, qd, adt, adtd)
      IMPLICIT NONE
      REAL*8 gam, gm1, cfl, eps, mach, alpha
C
      COMMON /constants/ gam, gm1, cfl, eps, mach, alpha
C
C
      INTEGER i
C
      REAL*8 x1(2), x2(2), x3(2), x4(2), q(4), adt, dx(4), dy(4), ri, u
     +       , v, p, c
      REAL*8 qd(4), adtd, rid, ud, vd, pd, cd
      INTRINSIC SQRT
      INTRINSIC ABS
      REAL*8 arg1
      REAL*8 arg1d
      REAL*8 result1
      REAL*8 abs0d
      REAL*8 abs0
C
C
C
C
C
C
C
      dx(1) = x2(1) - x1(1)
      dx(2) = x3(1) - x2(1)
      dx(3) = x4(1) - x3(1)
      dx(4) = x1(1) - x4(1)
C
      dy(1) = x2(2) - x1(2)
      dy(2) = x3(2) - x2(2)
      dy(3) = x4(2) - x3(2)
      dy(4) = x1(2) - x4(2)
C
      rid = -(qd(1)/q(1)**2)
      ri = 1.d0/q(1)
      ud = rid*q(2) + ri*qd(2)
      u = ri*q(2)
      vd = rid*q(3) + ri*qd(3)
      v = ri*q(3)
      pd = gm1*(qd(4)-0.5d0*(rid*(q(2)**2+q(3)**2)+ri*(2*q(2)*qd(2)+2*q(
     +  3)*qd(3))))
      p = gm1*(q(4)-0.5d0*ri*(q(2)**2+q(3)**2))
      arg1d = gam*(rid*p+ri*pd)
      arg1 = gam*ri*p
      IF (arg1 .EQ. 0.0) THEN
        cd = 0.D0
      ELSE
        cd = arg1d/(2.0*SQRT(arg1))
      END IF
      c = SQRT(arg1)
C
      adt = 0
      adtd = 0.0
C
      DO i=1,4
        IF (u*dy(i) - v*dx(i) .GE. 0.) THEN
          abs0d = dy(i)*ud - dx(i)*vd
          abs0 = u*dy(i) - v*dx(i)
        ELSE
          abs0d = -(dy(i)*ud-dx(i)*vd)
          abs0 = -(u*dy(i)-v*dx(i))
        END IF
        arg1 = dx(i)**2 + dy(i)**2
        result1 = SQRT(arg1)
        adtd = adtd + abs0d + result1*cd
        adt = adt + abs0 + c*result1
      ENDDO
C
      adtd = adtd/cfl
      adt = adt/cfl
C
      RETURN
      END