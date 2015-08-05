ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c      Adjoint airfoil lift perturbation calculation                c
c      using adjoint routines created by Tapenade in reverse mode   c
c                                                                   c
c      Copyright Devendra Ghate and Mike Giles, 2005                c
c      but can be freely used with due acknowledgement              c
c                                                                   c
c      Tapenade developed by Laurent Hascoet and others at INRIA    c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        program air_adj
        implicit none
c
#include "const.inc"
c
        integer    maxnode,     maxcell,       maxedge
        parameter (maxnode=721801,maxcell=720000+1,maxedge=1441400)
c
        integer*4 ecell(2,maxedge), boun(maxedge),
     &            edge(2,maxedge), cell(4,maxcell)
        real*8    x(2,maxnode), xb(2,maxnode), xb2(2,maxnode),
     &            q(4,maxcell), qb(4,maxcell), qb2(4,maxcell),
     &            adt(maxcell), adtb(maxcell),  xd(2,maxnode),
     &            resb(4,maxcell), resbold(4,maxcell)
c
        integer*4 in,in1,in2,in3,in4,ic,ic1,ic2,ie,ipde,
     &            iter,niter,k, nnode,ncell,nedge
        real*8    liftd, liftb, rms, dum1(4), dum2(4), dalpha
c
        do ipde = 1, 4
          dum1(ipde) = 0.0d0
          dum2(ipde) = 0.0d0
        enddo
c
c------ read in grid and flow data
c
        call input(maxnode,maxcell,maxedge, nnode,ncell,nedge,
     &             x, q,cell, edge,ecell,boun)
c
c------ define grid perturbation (simple rotation) for validation
c
        dalpha = datan(1.0d0)/45.0d0
c
        do in = 1, nnode
          xd(1,in) =  dalpha*x(2,in)
          xd(2,in) = -dalpha*x(1,in)
        enddo
c
c------ adjoint lift calculation
c
        do in = 1, nnode
          xb2(1,in) = 0.0d0
          xb2(2,in) = 0.0d0
        enddo
c
        do ic = 1, ncell
          do ipde = 1, 4
            qb2(ipde,ic)  = 0.0d0
            resb(ipde,ic) = 0.0d0
          enddo
        enddo
c
        liftb = 1.0d0
c
c        write(*,*) 'entering lift_wall_bx'
c
        do ie = 1, nedge
          if(boun(ie).eq.1) then
            in1 = edge(1,ie)
            in2 = edge(2,ie)
            ic2 = ecell(2,ie)
            call lift_wall_bx(x(1,in1),xb2(1,in1),
     &                        x(1,in2),xb2(1,in2),
     &                        q(1,ic2),qb2(1,ic2), dum1,liftb)
          endif
        enddo
c
c-------- calculate area/timstep
c
c        write(*,*) 'entering time_cell'
c
        do ic = 1, ncell-1
          in1 = cell(1,ic)
          in2 = cell(2,ic)
          in3 = cell(3,ic)
          in4 = cell(4,ic)
          call time_cell(x(1,in1),x(1,in2),x(1,in3),x(1,in4),
     &                   q(1,ic),adt(ic))
        enddo
c
        adt(ncell)  = 0.0d0
c
c------ main time-marching loop
c
        niter = 20000
c
        do iter = 1, niter
c
c-------- save old adjoint solution
c
          do ic = 1, ncell-1
            do ipde = 1,4
              resbold(ipde,ic) = resb(ipde,ic)
            enddo
          enddo
c
c-------- predictor/corrector update loop
c
          do k = 1, 2
            do ic = 1, ncell
              do ipde = 1, 4
                qb(ipde,ic) = 0.0d0
              enddo
              adtb(ic) = 0.0d0
            enddo
c
c---------- flux evaluation loop
c
c            write(*,*) 'entering flux_fac_b/flux_wall_b'
c
            do ie = 1, nedge
              in1 = edge(1,ie)
              in2 = edge(2,ie)
              ic1 = ecell(1,ie)
              ic2 = ecell(2,ie)
c
              if(boun(ie).eq.0) then
                call flux_face_b(x(1,in1),x(1,in2),
     &                           q(1,ic1),qb(1,ic1),
     &                           q(1,ic2),qb(1,ic2),
     &                           adt(ic1),adtb(ic1),
     &                           adt(ic2),adtb(ic2),
     &                           dum1,resb(1,ic1), dum2,resb(1,ic2))
              elseif (boun(ie).eq.1) then
                call flux_wall_b(x(1,in1),x(1,in2),
     &                           q(1,ic2),qb(1,ic2), dum2,resb(1,ic2))
              endif
            enddo
c
c---------- adjoint calculation of area/timstep
c
c            write(*,*) 'entering time_cell_b'
c
            do ic = 1, ncell-1
              in1 = cell(1,ic)
              in2 = cell(2,ic)
              in3 = cell(3,ic)
              in4 = cell(4,ic)
              call time_cell_b(x(1,in1),x(1,in2),x(1,in3),x(1,in4),
     &                         q(1,ic),qb(1,ic),adt(ic),adtb(ic))
            enddo
c
c---------- adjoint solution update
c
            rms = 0.0d0
c
            do ic = 1, ncell-1
              do ipde = 1, 4
                qb(ipde,ic)   = qb(ipde,ic) + qb2(ipde,ic)
                resb(ipde,ic) = resbold(ipde,ic) - qb(ipde,ic)/adt(ic)
                rms = rms + qb(ipde,ic)**2
              enddo
            enddo
          enddo
c
          rms = dsqrt(rms/ncell)
c
c-------- print iteration history, including adjoint lift calculation
c
          if(mod(iter,100).eq.0) then
c
c---------- adjoint calculation of residual perturbation due to perturbed grid
c
            do in = 1, nnode
              xb(1,in) = xb2(1,in)
              xb(2,in) = xb2(2,in)
            enddo
c
c            write(*,*) 'entering flux_face_bx/flux_wall_bx'
c
            do ie = 1, nedge
              in1 = edge(1,ie)
              in2 = edge(2,ie)
              ic1 = ecell(1,ie)
              ic2 = ecell(2,ie)
c
              if(boun(ie).eq.0) then
                call flux_face_bx(x(1,in1),xb(1,in1),
     &                            x(1,in2),xb(1,in2),
     &                            q(1,ic1),qb(1,ic1),
     &                            q(1,ic2),qb(1,ic2),
     &                            adt(ic1),adtb(ic1),
     &                            adt(ic2),adtb(ic2),
     &                            dum1,resb(1,ic1), dum2,resb(1,ic2))
              elseif (boun(ie).eq.1) then
                call flux_wall_bx(x(1,in1),xb(1,in1),
     &                            x(1,in2),xb(1,in2),
     &                            q(1,ic2),qb(1,ic2), dum2,resb(1,ic2))
              endif
            enddo
c
c            write(*,*) 'entering time_cell_bx'
c
            do ic = 1, ncell-1
              in1 = cell(1,ic)
              in2 = cell(2,ic)
              in3 = cell(3,ic)
              in4 = cell(4,ic)
              call time_cell_bx(x(1,in1),xb(1,in1),
     &                          x(1,in2),xb(1,in2),
     &                          x(1,in3),xb(1,in3),
     &                          x(1,in4),xb(1,in4),
     &                          q(1,ic),qb(1,ic), dum2,adtb(ic))
            enddo
c
            liftd = 0.0d0
c
            do in = 1, nnode
              liftd = liftd + xb(1,in)*xd(1,in) + xb(2,in)*xd(2,in)
            enddo
c
            write(*,'(i5,e13.6,e23.16)') iter,rms,liftd
          endif
        enddo
c
        stop
        end
