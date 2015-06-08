ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c      Nonlinear routines used in airfoil calculations              c
c                                                                   c
c      linear/adjoint versions are created on-the-fly by Tapenade   c
c      complex version is generated using compiler flag -DCOMPLEX   c
c                                                                   c
c      Copyright Devendra Ghate and Mike Giles, 2005                c
c      but can be freely used with due acknowledgement              c
c                                                                   c
c      Tapenade developed by Laurent Hascoet and others at INRIA    c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
#ifdef COMPLEX
       subroutine Ctime_cell(x1,x2,x3,x4,q,adt)
#else
       subroutine time_cell(x1,x2,x3,x4,q,adt)
#endif
c
c      compute area/timestep for an individual cell
c
       implicit none
#include "const.inc"
c
       integer  i
c
#ifdef COMPLEX
       complex*16
#else
       real*8
#endif
     &            x1(2),x2(2),x3(2),x4(2),q(4),adt,
     &            dx(4),dy(4), ri,u,v,p,c
c
#ifdef COMPLEX
       external   min, max, abs
       complex*16 min, max, abs
#endif
c
       dx(1) = x2(1) - x1(1)
       dx(2) = x3(1) - x2(1)
       dx(3) = x4(1) - x3(1)
       dx(4) = x1(1) - x4(1)
c
       dy(1) = x2(2) - x1(2)
       dy(2) = x3(2) - x2(2)
       dy(3) = x4(2) - x3(2)
       dy(4) = x1(2) - x4(2)
c
       ri = 1.d0/q(1)
       u  =   ri*q(2)
       v  =   ri*q(3)
       p  = gm1*(q(4)-0.5d0*ri*(q(2)**2+q(3)**2))
       c  = sqrt(gam*ri*p)
c
       adt = 0
c
       do i = 1,4
         adt = adt + abs(u*dy(i)-v*dx(i)) + c*sqrt(dx(i)**2+dy(i)**2)
       enddo
c
       adt = adt / cfl
c
       return
       end
