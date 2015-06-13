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
        use print_active
        implicit none
c
#include "const.inc"
c
        integer    maxnode,     maxcell,       maxedge
        parameter (maxnode=9900,maxcell=9702+1,maxedge=19502)
c
        integer*4 ecell(2,maxedge), boun(maxedge),
     &            edge(2,maxedge), cell(4,maxcell)
        real*8    x(2,maxnode), xd(2,maxnode), 
     &            q(4,maxcell), qd(4,maxcell), qdold(4,maxcell),
     &            adt(maxcell), adtd(maxcell),
     &            resd(4,maxcell), resd2(4,maxcell)
c
        integer*4 in,in1,in2,in3,in4,ic,ic1,ic2,ie,ipde,
     &            iter,niter,k, nnode,ncell,nedge
        real*8    liftd, rms, dum1(4), dum2(4), dalpha
c
c--------- NO STATEMENTS before this -----------------

#if DEBUG
c---------- DECLARE Variables to hold filenumbers.
        integer*4 :: x_OAD_v_filenum,
     &            x_OAD_d_filenum, q_OAD_v_filenum,
     &            q_OAD_d_filenum, adt_OAD_v_filenum, adt_OAD_d_filenum,
     &            res1_OAD_d_filenum, res1_OAD_v_filenum,
     &            res2_OAD_d_filenum, res2_OAD_v_filenum,
     &            lift_OAD_d_filenum, lift_OAD_v_filenum, tmp_v1, tmp_v2

        character(len=10) :: filename

c---------- SET filenum values
        x_OAD_v_filenum = 11
        x_OAD_d_filenum = 12
        q_OAD_v_filenum = 13
        q_OAD_d_filenum = 14
        adt_OAD_v_filenum = 15
        adt_OAD_d_filenum = 16
        res1_OAD_v_filenum = 17
        res1_OAD_d_filenum = 18
        res2_OAD_v_filenum = 19
        res2_OAD_d_filenum = 20
        lift_OAD_v_filenum = 21
        lift_OAD_d_filenum = 22
c---------- END SET filenum values

c--------- Create filename from filenumber
c---------- OPEN files
        write(filename, '(A6, I2)') "output", x_OAD_v_filenum
        open(unit = x_OAD_v_filenum, file = filename)
        write(filename, '(A6, I2)') "output", x_OAD_d_filenum
        open(unit = x_OAD_d_filenum, file = filename)
        write(filename, '(A6, I2)') "output", q_OAD_v_filenum
        open(unit = q_OAD_v_filenum, file = filename)
        write(filename, '(A6, I2)') "output", q_OAD_d_filenum
        open(unit = q_OAD_d_filenum, file = filename)
        write(filename, '(A6, I2)') "output", adt_OAD_v_filenum
        open(unit = adt_OAD_v_filenum, file = filename)
        write(filename, '(A6, I2)') "output", adt_OAD_d_filenum
        open(unit = adt_OAD_d_filenum, file = filename)
        write(filename, '(A6, I2)') "output", res1_OAD_v_filenum
        open(unit = res1_OAD_v_filenum, file = filename)
        write(filename, '(A6, I2)') "output", res1_OAD_d_filenum
        open(unit = res1_OAD_d_filenum, file = filename)
        write(filename, '(A6, I2)') "output", res2_OAD_v_filenum
        open(unit = res2_OAD_v_filenum, file = filename)
        write(filename, '(A6, I2)') "output", res2_OAD_d_filenum
        open(unit = res2_OAD_d_filenum, file = filename)
        write(filename, '(A6, I2)') "output", lift_OAD_v_filenum
        open(unit = lift_OAD_v_filenum, file = filename)
        write(filename, '(A6, I2)') "output", lift_OAD_d_filenum
        open(unit = lift_OAD_d_filenum, file = filename)
#endif

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
c------ define grid perturbation (simple rotation)
c
        dalpha = datan(1.0d0)/45.0d0
c
        do in = 1, nnode
          xd(1,in) =  dalpha*x(2,in)
          xd(2,in) = -dalpha*x(1,in)
        enddo
c
c------ calculate residual perturbation due to perturbed grid
c
        do ic = 1, ncell
          do ipde = 1, 4
            qd(ipde,ic)    = 0.0d0
            resd2(ipde,ic) = 0.0d0
          enddo
        enddo
c
        do ic = 1, ncell-1
          in1 = cell(1,ic)
          in2 = cell(2,ic)
          in3 = cell(3,ic)
          in4 = cell(4,ic)
          call time_cell_dx(x(1,in1),xd(1,in1),
     &                      x(1,in2),xd(1,in2),
     &                      x(1,in3),xd(1,in3),
     &                      x(1,in4),xd(1,in4),
     &                      q(1,ic),qd(1,ic), adt(ic),adtd(ic))
        enddo

c
        adt(ncell)  = 0.0d0
        adtd(ncell) = 0.0d0
c
        
#if DEBUG
c--------- Write to files
        tmp_v1 = 0
        call write_array(adt, 1, tmp_v1,
     &     adt_OAD_v_filenum)

        tmp_v1 = 0
        call write_array(adtd, 1, tmp_v1,
     &     adt_OAD_d_filenum)

        tmp_v1 = 0
        tmp_v2 = 0
        call write_array(x, 1, tmp_v1, 1, tmp_v2,
     &     x_OAD_v_filenum)

        tmp_v1 = 0
        tmp_v2 = 0
        call write_array(xd, 1, tmp_v1, 1, tmp_v2,
     &     x_OAD_d_filenum)

        tmp_v1 = 0
        tmp_v2 = 0
        call write_array(q, 1, tmp_v1, 1, tmp_v2,
     &     q_OAD_v_filenum)

        tmp_v1 = 0
        tmp_v2 = 0
        call write_array(qd, 1, tmp_v1, 1, tmp_v2,
     &     q_OAD_d_filenum)
#endif

        do ie = 1, nedge
          in1 = edge(1,ie)
          in2 = edge(2,ie)
          ic1 = ecell(1,ie)
          ic2 = ecell(2,ie)
c
          if(boun(ie).eq.0) then
            call flux_face_dx(x(1,in1),xd(1,in1),
     &                        x(1,in2),xd(1,in2),
     &                        q(1,ic1),qd(1,ic1),
     &                        q(1,ic2),qd(1,ic2), 
     &                        adt(ic1),adtd(ic1),
     &                        adt(ic2),adtd(ic2),
     &                        dum1,resd2(1,ic1), dum2,resd2(1,ic2))
          elseif (boun(ie).eq.1) then
            call flux_wall_dx(x(1,in1),xd(1,in1),
     &                        x(1,in2),xd(1,in2),
     &                        q(1,ic2),qd(1,ic2), dum2,resd2(1,ic2))
          endif
        enddo

#if DEBUG
c--------- Write to files
        tmp_v1 = 0
        call write_array(adt, 1, tmp_v1,
     &     adt_OAD_v_filenum)

        tmp_v1 = 0
        call write_array(adtd, 1, tmp_v1,
     &     adt_OAD_d_filenum)

        tmp_v1 = 0
        tmp_v2 = 0
        call write_array(resd2, 1, tmp_v1, 1, tmp_v2,
     &     res2_OAD_d_filenum)

        tmp_v1 = 0
        tmp_v2 = 0
        call write_array(x, 1, tmp_v1, 1, tmp_v2,
     &     x_OAD_v_filenum)

        tmp_v1 = 0
        tmp_v2 = 0
        call write_array(xd, 1, tmp_v1, 1, tmp_v2,
     &     x_OAD_d_filenum)

        tmp_v1 = 0
        tmp_v2 = 0
        call write_array(q, 1, tmp_v1, 1, tmp_v2,
     &     q_OAD_v_filenum)

        tmp_v1 = 0
        tmp_v2 = 0
        call write_array(qd, 1, tmp_v1, 1, tmp_v2,
     &     q_OAD_d_filenum)
#endif

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
              qdold(ipde,ic) = qd(ipde,ic)
            enddo
          enddo
c 
c-------- predictor/corrector update loop
c
          do k = 1, 2
            do ic = 1, ncell
              do ipde = 1, 4
                resd(ipde,ic) = 0.0d0
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
              call time_cell_d(x(1,in1),x(1,in2),x(1,in3),x(1,in4),
     &                         q(1,ic),qd(1,ic), adt(ic),adtd(ic))
            enddo
c
#if DEBUG
c--------- Write to files
            tmp_v1 = 0
            call write_array(adt, 1, tmp_v1,
     &       adt_OAD_v_filenum)

            tmp_v1 = 0
            call write_array(adtd, 1, tmp_v1,
     &      adt_OAD_d_filenum)

            tmp_v1 = 0
            tmp_v2 = 0
            call write_array(x, 1, tmp_v1, 1, tmp_v2,
     &       x_OAD_v_filenum)

            tmp_v1 = 0
            tmp_v2 = 0
            call write_array(xd, 1, tmp_v1, 1, tmp_v2,
     &       x_OAD_d_filenum)

            tmp_v1 = 0
            tmp_v2 = 0
            call write_array(q, 1, tmp_v1, 1, tmp_v2,
     &       q_OAD_v_filenum)

            tmp_v1 = 0
            tmp_v2 = 0
            call write_array(qd, 1, tmp_v1, 1, tmp_v2,
     &       q_OAD_d_filenum)
#endif

            adt(ncell)  = 0.0d0
            adtd(ncell) = 0.0d0
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
                call flux_face_d(x(1,in1),x(1,in2),
     &                           q(1,ic1),qd(1,ic1),
     &                           q(1,ic2),qd(1,ic2), 
     &                           adt(ic1),adtd(ic1),
     &                           adt(ic2),adtd(ic2),
     &                           dum1,resd(1,ic1), dum2,resd(1,ic2))
              elseif (boun(ie).eq.1) then
                call flux_wall_d(x(1,in1),x(1,in2),
     &                           q(1,ic2),qd(1,ic2), dum2,resd(1,ic2))
              endif
            enddo

#if DEBUG
c--------- Write to files
          tmp_v1 = 0
          call write_array(adt, 1, tmp_v1,
     &     adt_OAD_v_filenum)

          tmp_v1 = 0
          call write_array(adtd, 1, tmp_v1,
     &     adt_OAD_d_filenum)

          tmp_v1 = 0
          tmp_v2 = 0
          call write_array(resd, 1, tmp_v1, 1, tmp_v2,
     &     res1_OAD_d_filenum)

          tmp_v1 = 0
          tmp_v2 = 0
          call write_array(x, 1, tmp_v1, 1, tmp_v2,
     &     x_OAD_v_filenum)

          tmp_v1 = 0
          tmp_v2 = 0
          call write_array(xd, 1, tmp_v1, 1, tmp_v2,
     &     x_OAD_d_filenum)

          tmp_v1 = 0
          tmp_v2 = 0
          call write_array(q, 1, tmp_v1, 1, tmp_v2,
     &     q_OAD_v_filenum)

          tmp_v1 = 0
          tmp_v2 = 0
          call write_array(qd, 1, tmp_v1, 1, tmp_v2,
     &     q_OAD_d_filenum)
#endif

c
c---------- flow field update
c
            rms = 0.0d0

c
            do ic = 1, ncell-1
              do ipde = 1, 4
                resd(ipde,ic) = resd(ipde,ic) + resd2(ipde,ic)
                qd(ipde,ic)   = qdold(ipde,ic) - resd(ipde,ic)/adt(ic)
                rms = rms + resd(ipde,ic)**2
              enddo
            enddo

#if DEBUG
c--------- Write to files
          tmp_v1 = 0
          call write_array(adt, 1, tmp_v1,
     &     adt_OAD_v_filenum)

          tmp_v1 = 0
          call write_array(adtd, 1, tmp_v1,
     &     adt_OAD_d_filenum)

          tmp_v1 = 0
          tmp_v2 = 0
          call write_array(resd, 1, tmp_v1, 1, tmp_v2,
     &     res1_OAD_d_filenum)

          tmp_v1 = 0
          tmp_v2 = 0
          call write_array(x, 1, tmp_v1, 1, tmp_v2,
     &     x_OAD_v_filenum)

          tmp_v1 = 0
          tmp_v2 = 0
          call write_array(xd, 1, tmp_v1, 1, tmp_v2,
     &     x_OAD_d_filenum)

          tmp_v1 = 0
          tmp_v2 = 0
          call write_array(q, 1, tmp_v1, 1, tmp_v2,
     &     q_OAD_v_filenum)

          tmp_v1 = 0
          tmp_v2 = 0
          call write_array(qd, 1, tmp_v1, 1, tmp_v2,
     &     q_OAD_d_filenum)
#endif

          enddo
c
          rms = dsqrt(rms/ncell)
c          print *, rms

c
c-------- print iteration history, including lift calculation
c
          if(mod(iter,100).eq.0) then
            liftd = 0.0d0
c
            do ie = 1,nedge
              if(boun(ie).eq.1) then
                in1 = edge(1,ie)
                in2 = edge(2,ie)
                ic2 = ecell(2,ie)
                call lift_wall_dx(x(1,in1),xd(1,in1),
     &                            x(1,in2),xd(1,in2),
     &                            q(1,ic2),qd(1,ic2), dum1,liftd)
              endif
            enddo
c
            write(*,'(i5,e13.6,e23.16)') iter,rms,liftd
          endif
        enddo
c 

#if DEBUG
        close(x_OAD_v_filenum)
        close(x_OAD_d_filenum)
        close(q_OAD_v_filenum)
        close(q_OAD_d_filenum)
        close(adt_OAD_v_filenum)
        close(adt_OAD_d_filenum)
        close(res1_OAD_v_filenum)
        close(res1_OAD_d_filenum)
        close(res2_OAD_v_filenum)
        close(res2_OAD_d_filenum)
        close(lift_OAD_v_filenum)
        close(lift_OAD_d_filenum)
#endif
        stop
        end
