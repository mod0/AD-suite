ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c      Nonlinear airfoil lift calculation                           c
c                                                                   c
c      Copyright Devendra Ghate and Mike Giles, 2005                c
c      but can be freely used with due acknowledgement              c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        program airfoil
        implicit none
c
#include "const.inc"
c
        integer    maxnode,     maxcell,       maxedge
        parameter (maxnode=9900,maxcell=9702+1,maxedge=19502)
c
        integer*4 ecell(2,maxedge), boun(maxedge),
     &            edge(2,maxedge), cell(4,maxcell)
        real*8    x(2,maxnode), q(4,maxcell),qold(4,maxcell),
     &            adt(maxcell),res(4,maxcell)
c
        integer*4 in1,in2,in3,in4,ic,ic1,ic2,ie,ipde,
     &            iter,niter,k, nnode,ncell,nedge
        real*8    lift, rms
c
c------ read in grid and flow data
c
        call input(maxnode,maxcell,maxedge, nnode,ncell,nedge,
     &             x, q,cell, edge,ecell,boun)
c
c------ main time-marching loop
c
        niter = 20000
c
        do iter = 1, niter
c
c-------- save old flow solution
c
          do ic = 1, ncell-1
            do ipde = 1,4
              qold(ipde,ic) = q(ipde,ic)
            enddo
          enddo
c 
c-------- predictor/corrector update loop
c
          do k = 1, 2
            do ic = 1, ncell
              do ipde = 1, 4
                res(ipde,ic) = 0.0d0
              enddo
            enddo
c
c---------- calculate area/timstep
c
            do ic = 1, ncell-1
              in1 = cell(1,ic)
              in2 = cell(2,ic)
              in3 = cell(3,ic)
              in4 = cell(4,ic)
              call time_cell(x(1,in1),x(1,in2),x(1,in3),x(1,in4),
     &                       q(1,ic),adt(ic))
            enddo
c
            adt(ncell) = 0.0d0
c 
c---------- flux evaluation loop 
c
            do ie = 1, nedge
              in1 = edge(1,ie)
              in2 = edge(2,ie)
              ic1 = ecell(1,ie)
              ic2 = ecell(2,ie)
c
              if(boun(ie).eq.0) then
                call flux_face(x(1,in1),x(1,in2),q(1,ic1),q(1,ic2), 
     &                         adt(ic1),adt(ic2),res(1,ic1),res(1,ic2))
              elseif (boun(ie).eq.1) then
                call flux_wall(x(1,in1),x(1,in2),q(1,ic2),res(1,ic2))
	      elseif (boun(ie).eq.2) then
                stop 'Problem'
              endif
            enddo
c
c---------- flow field update
c
            rms = 0.0d0
c
            do ic = 1, ncell-1
              do ipde = 1, 4
                q(ipde,ic) = qold(ipde,ic) - res(ipde,ic)/adt(ic)
                rms = rms + (q(ipde,ic)-qold(ipde,ic))**2
              enddo
            enddo
          enddo
c
          rms = dsqrt(rms/ncell)
c
c-------- print iteration history, including lift calculation
c
          if(mod(iter,100).eq.0) then
            lift = 0.0d0
c
            do ie = 1,nedge
              if(boun(ie).eq.1) then
                in1 = edge(1,ie)
                in2 = edge(2,ie)
                ic2 = ecell(2,ie)
                call lift_wall(x(1,in1),x(1,in2),q(1,ic2),lift)
              endif
            enddo
c
            write(*,'(i5,e13.6,2e23.16)') iter,rms,lift
          endif
        enddo
c
        call output(ncell,q)
c 
        stop
        end
