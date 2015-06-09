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
       subroutine Clift_wall(x1,x2,q,lift)
#else
       subroutine lift_wall(x1,x2,q,lift)
#endif
c
c      compute momentum flux from an individual wall face
c
       implicit none
#include "const.inc"
c
c
#ifdef COMPLEX
       complex*16
#else
       real*8
#endif
     &          x1(2),x2(2),q(4),lift,
     &          dx, ri,u,v,p
c
c----------------------------------------------------------------------
c      OpenAD Directives
c----------------------------------------------------------------------
#ifdef OPENAD_D
c$OPENAD INDEPENDENT(q)
c$OPENAD DEPENDENT(lift)
#elif OPENAD_DX
c$OPENAD INDEPENDENT(x1)
c$OPENAD INDEPENDENT(x2)
c$OPENAD INDEPENDENT(q)
c$OPENAD DEPENDENT(lift)
#endif
c
       dx = x1(1) - x2(1)
c
       ri = 1.d0/q(1)
       u  =   ri*q(2)
       v  =   ri*q(3)
       p  = gm1*(q(4)-0.5d0*ri*(q(2)**2+q(3)**2))
c
       lift = lift + p*dx
c
       return
       end
