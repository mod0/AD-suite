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
#ifdef COMPLEX
       subroutine Cflux_wall(x1,x2,q,res)
#else
       subroutine flux_wall(x1,x2,q,res)
#endif
c
c      compute momentum flux from an individual wall face
c
       implicit none
#include "const.inc"
c
       integer  n
c
#ifdef COMPLEX
       complex*16
#else
       real*8
#endif
     &            x1(2),x2(2),q(4),res(4),
     &            dx,dy, ri,u,v,p
c
       dx = x1(1) - x2(1)
       dy = x1(2) - x2(2)
c
       ri = 1.d0/q(1)
       u  =   ri*q(2)
       v  =   ri*q(3)
       p  = gm1*(q(4)-0.5d0*ri*(q(2)**2+q(3)**2))
c
       res(2) = res(2) - p*dy
       res(3) = res(3) + p*dx
c
       return
       end
