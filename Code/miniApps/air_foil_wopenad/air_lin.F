ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c      Linear airfoil lift perturbation calculation                 c
c      using linearised routines created by Tapenade                c
c                                                                   c
c      Copyright Devendra Ghate and Mike Giles, 2005                c
c      but can be freely used with due acknowledgement              c
c                                                                   c
c      Tapenade developed by Laurent Hascoet and others at INRIA    c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        program air_lin
        use OAD_active
        use w2f__types
        implicit none

c
#include "const.inc"
c
        integer    maxnode,     maxcell,       maxedge
        parameter (maxnode=9900,maxcell=9702+1,maxedge=19502)
c
        integer*4 ecell(2,maxedge), boun(maxedge),
     &            edge(2,maxedge), cell(4,maxcell)
        real*8    x(2,maxnode), q(4,maxcell), qdold(4,maxcell)
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c     &           xd(2,maxnode), qd(4,maxcell),
c     &            adt(maxcell), adtd(maxcell),
c     &            resd(4,maxcell), resd2(4,maxcell)
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<

        type(active) x_OAD(2,maxnode), q_OAD(4,maxcell),
     &               adt_OAD(maxcell), res1_OAD(4, maxcell),
     &               res2_OAD(4, maxcell), lift_OAD

        integer*4 in,in1,in2,in3,in4,ic,ic1,ic2,ie,ipde,
     &            iter,niter,k, nnode,ncell,nedge
        real*8     rms, dalpha
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c     &    ,liftd, dum1(4), dum2(4)
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
c
        do ipde = 1, 4
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c          dum1(ipde) = 0.0d0
c          dum2(ipde) = 0.0d0
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
           do ic = 1, maxcell
              res1_OAD(ipde, ic)%v = 0.0d0
              res2_OAD(ipde, ic)%v = 0.0d0
           enddo
c--------- TODO: (1) The original code pases dum1, dum2 everywhere res
c--------- is needed. Need to check that if we just pass res itself,
c--------- there are no unintended consequences.
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
        enddo
c
c------ read in grid and flow data
c
        call input(maxnode,maxcell,maxedge, nnode,ncell,nedge,
     &             x, q,cell, edge,ecell,boun)

c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c       Copy the values from x to x_OAD
        do in = 1, nnode
            x_OAD(1, in)%v = x(1, in)
            x_OAD(2, in)%v = x(2, in)
        enddo

c       Copy the values from q to q_OAD
        do ic = 1, maxcell
            q_OAD(1, ic)%v = q(1, ic)
            q_OAD(2, ic)%v = q(2, ic)
            q_OAD(3, ic)%v = q(3, ic)
            q_OAD(4, ic)%v = q(4, ic)
        enddo
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<

c------ define grid perturbation (simple rotation)
c
        dalpha = datan(1.0d0)/45.0d0
c
        do in = 1, nnode
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c          xd(1,in) =  dalpha*x(2,in)
c          xd(2,in) = -dalpha*x(1,in)
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
          x_OAD(1,in)%d =  dalpha*x(2,in)
          x_OAD(2,in)%d = -dalpha*x(1,in)
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
        enddo
c
c------ calculate residual perturbation due to perturbed grid
c
        do ic = 1, ncell
          do ipde = 1, 4
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c            qd(ipde,ic)    = 0.0d0
c            resd2(ipde,ic) = 0.0d0
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
            q_OAD(ipde,ic)%d    = 0.0d0
            res2_OAD(ipde,ic)%d = 0.0d0
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
          enddo
        enddo
c
        do ic = 1, ncell-1
          in1 = cell(1,ic)
          in2 = cell(2,ic)
          in3 = cell(3,ic)
          in4 = cell(4,ic)
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c          call time_cell_dx(x(1,in1),xd(1,in1),
c     &                      x(1,in2),xd(1,in2),
c     &                      x(1,in3),xd(1,in3),
c     &                      x(1,in4),xd(1,in4),
c     &                      q(1,ic),qd(1,ic), adt(ic),adtd(ic))
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
          call time_cell_dx(x_OAD(1,in1),
     &                      x_OAD(1,in2),
     &                      x_OAD(1,in3),
     &                      x_OAD(1,in4),
     &                      q_OAD(1,ic),
     &                      adt_OAD(ic))
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
        enddo
c

c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c        adt(ncell)  = 0.0d0
c        adtd(ncell) = 0.0d0
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
        adt_OAD(ncell)%v = 0.0d0
        adt_OAD(ncell)%d = 0.0d0
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<

c
        do ie = 1, nedge
          in1 = edge(1,ie)
          in2 = edge(2,ie)
          ic1 = ecell(1,ie)
          ic2 = ecell(2,ie)
c
          if(boun(ie).eq.0) then
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c            call flux_face_dx(x(1,in1),xd(1,in1),
c     &                        x(1,in2),xd(1,in2),
c     &                        q(1,ic1),qd(1,ic1),
c     &                        q(1,ic2),qd(1,ic2),
c     &                        adt(ic1),adtd(ic1),
c     &                        adt(ic2),adtd(ic2),
c     &                        dum1,resd2(1,ic1), dum2,resd2(1,ic2))
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
            call flux_face_dx(x_OAD(1,in1),
     &                        x_OAD(1,in2),
     &                        q_OAD(1,ic1),
     &                        q_OAD(1,ic2),
     &                        adt_OAD(ic1),
     &                        adt_OAD(ic2),
     &                        res2_OAD(1,ic1),
     &                        res2_OAD(1,ic2))
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
          elseif (boun(ie).eq.1) then
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c            call flux_wall_dx(x(1,in1),xd(1,in1),
c     &                        x(1,in2),xd(1,in2),
c     &                        q(1,ic2),qd(1,ic2), dum2,resd2(1,ic2))
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
            call flux_wall_dx(x_OAD(1,in1),
     &                        x_OAD(1,in2),
     &                        q_OAD(1,ic2),
     &                        res2_OAD(1,ic2))
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
          endif
        enddo
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
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c              qdold(ipde,ic) = qd(ipde,ic)
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
              qdold(ipde,ic) = q_OAD(ipde,ic)%d
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
            enddo
          enddo
c 
c-------- predictor/corrector update loop
c
          do k = 1, 2
            do ic = 1, ncell
              do ipde = 1, 4
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c                resd(ipde,ic) = 0.0d0
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
                res1_OAD(ipde,ic)%d = 0.0d0
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
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
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c              call time_cell_d(x(1,in1),x(1,in2),x(1,in3),x(1,in4),
c     &                         q(1,ic),qd(1,ic), adt(ic),adtd(ic))
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
              call time_cell_d( x_OAD(1,in1),
     &                          x_OAD(1,in2),
     &                          x_OAD(1,in3),
     &                          x_OAD(1,in4),
     &                          q_OAD(1,ic),
     &                          adt_OAD(ic))
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
            enddo
c
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c        adt(ncell)  = 0.0d0
c        adtd(ncell) = 0.0d0
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
        adt_OAD(ncell)%v  = 0.0d0
        adt_OAD(ncell)%d = 0.0d0
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
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
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c                call flux_face_d(x(1,in1),x(1,in2),
c     &                           q(1,ic1),qd(1,ic1),
c     &                           q(1,ic2),qd(1,ic2),
c     &                           adt(ic1),adtd(ic1),
c     &                           adt(ic2),adtd(ic2),
c     &                           dum1,resd(1,ic1), dum2,resd(1,ic2))
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
                call flux_face_d(x_OAD(1,in1),
     &                           q_OAD(1,ic1),
     &                           q_OAD(1,ic2),
     &                           adt_OAD(ic1),
     &                           adt_OAD(ic2),
     &                           res1_OAD(1,ic1),
     &                           res1_OAD(1,ic2))
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
              elseif (boun(ie).eq.1) then
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c                call flux_wall_d(x(1,in1),x(1,in2),
c     &                           q(1,ic2),qd(1,ic2), dum2,resd(1,ic2))
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
                call flux_wall_d( x_OAD(1,in1),
     &                            x_OAD(1,in2),
     &                            q_OAD(1,ic2),
     &                            res1_OAD(1,ic2))
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
              endif
            enddo
c
c---------- flow field update
c
            rms = 0.0d0
c
            do ic = 1, ncell-1
              do ipde = 1, 4
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c                resd(ipde,ic) = resd(ipde,ic) + resd2(ipde,ic)
c                qd(ipde,ic)   = qdold(ipde,ic) - resd(ipde,ic)/adt(ic)
c                rms = rms + resd(ipde,ic)**2
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
                res1_OAD(ipde,ic)%d = res1_OAD(ipde,ic)%d
     &                                + res1_OAD(ipde,ic)%d
                q_OAD(ipde,ic)%d = qdold(ipde,ic)
     &                           - (res1_OAD(ipde,ic)%d)/(adt_OAD(ic)%v)
                rms = rms + (res1_OAD(ipde,ic)%d)**2
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
              enddo
            enddo
          enddo
c
          rms = dsqrt(rms/ncell)
c
c-------- print iteration history, including lift calculation
c
          if(mod(iter,100).eq.0) then
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c            liftd = 0.0d0
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
            lift_OAD%d = 0.0d0
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
c
            do ie = 1,nedge
              if(boun(ie).eq.1) then
                in1 = edge(1,ie)
                in2 = edge(2,ie)
                ic2 = ecell(2,ie)
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c                call lift_wall_dx(x(1,in1),xd(1,in1),
c     &                            x(1,in2),xd(1,in2),
c     &                            q(1,ic2),qd(1,ic2), dum1,liftd)
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
                call lift_wall_dx(x_OAD(1,in1),
     &                            x_OAD(1,in2),
     &                            q_OAD(1,ic2),
     &                            lift_OAD)
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
              endif
            enddo
c
c--------- BEGIN CHANGES TO ACCOMMODATE OPENAD >>>>>>>>>>
c--------- BEGIN ORIGINAL CODE COMMENT >>>>>>>>>>
c            write(*,'(i5,e13.6,e23.16)') iter,rms,liftd
c--------- END ORIGINAL CODE COMMENT <<<<<<<<<<
            write(*,'(i5,e13.6,e23.16)') iter,rms,lift_OAD%d
c--------- END CHANGES TO ACCOMMODATE OPENAD <<<<<<<<<<
          endif
        enddo
c 
        stop
        end
