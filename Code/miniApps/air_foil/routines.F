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


#ifdef COMPLEX
       subroutine Cflux_face(x1,x2,q1,q2,adt1,adt2,res1,res2)
#else
       subroutine flux_face(x1,x2,q1,q2,adt1,adt2,res1,res2)
#endif
c
c      compute flux through an individual face
c
       implicit none
#include "const.inc" 
c
       integer n,r
c
#ifdef COMPLEX
       complex*16
#else
       real*8
#endif
     &            x1(2),x2(2),q1(4),q2(4),adt1,adt2, res1(4),res2(4),
     &            dx,dy,mu, r1i,u1,v1,p1,vol1, r2i,u2,v2,p2,vol2, f(4)
c
       dx = x1(1) - x2(1)
       dy = x1(2) - x2(2)
       mu = 0.5d0*(adt1+adt2)*eps
c
       r1i  = 1.d0/q1(1)
       u1   =  r1i*q1(2)
       v1   =  r1i*q1(3)
       p1   = gm1*(q1(4)-0.5d0*r1i*(q1(2)**2+q1(3)**2))
       vol1 = u1*dy - v1*dx
c
       r2i  = 1.d0/q2(1)
       u2   =  r2i*q2(2)
       v2   =  r2i*q2(3)
       p2   = gm1*(q2(4)-0.5d0*r2i*(q2(2)**2+q2(3)**2))
       vol2 = u2*dy - v2*dx
c
       f(1) = 0.5d0*( vol1* q1(1)         + vol2* q2(1)         )
       f(2) = 0.5d0*( vol1* q1(2) + p1*dy + vol2* q2(2) + p2*dy )
       f(3) = 0.5d0*( vol1* q1(3) - p1*dx + vol2* q2(3) - p2*dx )
       f(4) = 0.5d0*( vol1*(q1(4)+p1)     + vol2*(q2(4)+p2)     )
c
       do n = 1, 4
         f(n)    = f(n) + mu*(q1(n)-q2(n))
         res1(n) = res1(n) + f(n)
         res2(n) = res2(n) - f(n)
       enddo
c
       return
       end


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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c      complex versions of intrinsic functions with analytic extension c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
#ifdef COMPLEX
       function max(arg1,arg2)
c
       implicit   none
       complex*16 arg1, arg2, max
c
       if(dble(arg1).gt.dble(arg2)) then
         max = arg1
       else
         max = arg2
       endif
c
       return
       end

       function min(arg1,arg2)
c
       implicit   none
       complex*16 arg1, arg2, min
c
       if(dble(arg1).le.dble(arg2)) then
         min = arg1
       else
         min = arg2
       endif
c
       return
       end

       function abs(arg1)
c
       implicit   none
       complex*16 arg1, abs
c
       if(dble(arg1).ge.0.0) then
         abs = arg1
       else
         abs = -arg1
       endif
c
       return
       end
#endif
