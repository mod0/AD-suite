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
