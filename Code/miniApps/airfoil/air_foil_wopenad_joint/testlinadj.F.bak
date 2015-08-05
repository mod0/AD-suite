c**********************************************************************c
c                                                                      c
c   testlinadj -- linear validation versus complex nonlinear routines  c
c                 and adjoint validation versus linear routines        c
c                                                                      c
c**********************************************************************c
c
       program testlinadj
c
       implicit none
#include "const.inc" 
c
       real*8      x1(2), x2(2), x3(2), x4(2), q1(4), q2(4),
     &             adt1, adt2, res1(4), res2(4), lift,
     &            dx1(2),dx2(2),dx3(2),dx4(2),dq1(4),dq2(4),
     &            dadt1,dadt2,dres1(4),dres2(4),dlift,
     &            ax1(2),ax2(2),ax3(2),ax4(2),aq1(4),aq2(4),
     &            aadt1,aadt2,ares1(4),ares2(4),alift,adt
c
       complex*16 Cx1(2),Cx2(2),Cx3(2),Cx4(2),Cq1(4),Cq2(4),
     &            Cadt1,Cadt2,Cres1(4),Cres2(4),Clift
c
       real*8     a1(4,4), a2(4,4), a3(4,4), scale, emax
       integer*4  m,n,i,j, ipass, nmax
c
       character*80 format2, format4
c
c      write(format,'(a,i1,a,i1,a)') '(',4,'(',4,'F12.6,/)/)'
       format2 = '(4(2F12.6/))'
       format4 = '(4(4F12.6/))'
c
c----- set constants
c
       scale = 1.d-20
c
       gam = 1.4d0
       gm1 = gam - 1.d0
c
       eps = 0.25d0
       cfl = 1.d0
c
       alpha = 0.1
c
c----- define the mesh and the base flow
c
       x1(1) = 0.0d0
       x1(2) = 0.0d0
c
       x2(1) = 0.2d0
       x2(2) = 0.3d0
c
       x3(1) = 2.0d0
       x3(2) = 7.0d0
c
       x4(1) = -5.0d0
       x4(2) =  3.5d0
c
c..... flow in cell
       q1(1) = 1.234d0
       q1(2) = 0.345d0
       q1(3) = 0.456d0
       q1(4) = 4.567d0
c
       q2(1) = 1.123d0
       q2(2) = 0.234d0
       q2(3) = 0.345d0
       q2(4) = 4.456d0
c
       adt1 = 1.234d0
       adt2 = 2.345d0
c
c----- initialise max error to zero
c
       emax = 0.d0
c
c*********************************************************************
c
       do ipass = 1, 2
         nmax = 2 + 2*mod(ipass,2)
c
c------- linear timestep calculation
c
         do n = 1, nmax
           do m = 1, 2
             dx1(m) = 0.d0
             dx2(m) = 0.d0
             dx3(m) = 0.d0
             dx4(m) = 0.d0
           enddo
c
           do m = 1, 4
             dq1(m) = 0.d0
           enddo
c
           if(ipass.eq.1) then
             dq1(n) = 1.0d0
             call TIME_CELL_D(x1,x2,x3,x4,q1,dq1,adt1,dadt1)
           else
             dx1(n) = 1.0d0
             call TIME_CELL_DX(x1,dx1,x2,dx2,x3,dx3,x4,dx4,
     &                         q1,dq1, adt1,dadt1)
           endif
c
           a1(1,n) = dadt1
         enddo
c
c------- complex timestep calculation
c
         do n = 1, nmax
           do m = 1, 2
             Cx1(m) = x1(m)
             Cx2(m) = x2(m)
             Cx3(m) = x3(m)
             Cx4(m) = x4(m)
           enddo
c
           do m = 1, 4
             Cq1(m) = q1(m)
           enddo
c
           if(ipass.eq.1) then
             Cq1(n) = Cq1(n) + dcmplx(0.d0,scale)
           else
             Cx1(n) = Cx1(n) + dcmplx(0.d0,scale)
           endif
c
           call Ctime_cell(Cx1,Cx2,Cx3,Cx4,Cq1,Cadt1)
c
           a2(1,n) = aimag(Cadt1)/scale
           emax = max(emax, abs(a2(1,n)-a1(1,n)))
         enddo
c
c------- adjoint timestep calculation
c
         do m = 1, 4
           aq1(m) = 0.d0
         enddo
c
         aadt1 = 1.0d0
c
         if(ipass.eq.1) then
           call TIME_CELL_B(x1,x2,x3,x4,q1,aq1,adt1,aadt1)
         else
           call TIME_CELL_BX(x1,ax1,x2,ax2,x3,ax3,x4,ax4,
     &                       q1,aq1,adt1,aadt1)
         endif
c
         do n = 1, nmax
           if(ipass.eq.1) then
             a3(1,n) = aq1(n)
           else
             a3(1,n) = ax1(n)
           endif
           emax = max(emax, abs(a3(1,n)-a1(1,n)))
         enddo
c
c------- print out comparison
c
         if(ipass.eq.1) then
           write(*,*) 'linear/complex time_cell difference due to q'
           write(*,format4) ((a1(i,j)-a2(i,j),j=1,4),i=1,1)
c
           write(*,*) 'linear/adjoint time_cell difference due to q'
           write(*,format4) ((a1(i,j)-a3(i,j),j=1,4),i=1,1)
         else
           write(*,*) 'linear/complex time_cell difference due to x1'
           write(*,format2) ((a1(i,j)-a2(i,j),j=1,2),i=1,1)
c
           write(*,*) 'linear/adjoint time_cell difference due to x1'
           write(*,format2) ((a1(i,j)-a3(i,j),j=1,2),i=1,1)
         endif
       enddo
c
c*********************************************************************
c
       do ipass = 1, 4
         nmax = 2 + 2*mod(ipass,2)
c
c------- linear flux calculation
c
         do n = 1, nmax
           do m = 1, 2
             dx1(m) = 0.d0
             dx2(m) = 0.d0
           enddo
c
           do m = 1, 4
             dq1(m) = 0.d0
             dq2(m) = 0.d0
             dres1(m) = 0.d0
             dres2(m) = 0.d0
           enddo
c
           dadt1 = 0.d0
           dadt2 = 0.d0
c
           if (ipass.eq.1) then
             dq1(n) = 1.0d0
             call FLUX_FACE_D(x1,x2,q1,dq1,q2,dq2,adt1,dadt1,adt2,dadt2,
     &                        res1,dres1,res2,dres2)
           elseif (ipass.eq.2) then
             dx1(n) = 1.0d0
             call FLUX_FACE_DX(x1,dx1,x2,dx2,q1,dq1,q2,dq2,adt1,dadt1,
     &                         adt2,dadt2,res1,dres1,res2,dres2)
           elseif (ipass.eq.3) then
             dq1(n) = 1.0d0
             call FLUX_WALL_D(x1,x2,q1,dq1,res1,dres1)
           elseif (ipass.eq.4) then
             dx1(n) = 1.0d0
             call FLUX_WALL_DX(x1,dx1,x2,dx2,q1,dq1,res1,dres1)
           endif
c
           do m = 1, 4
             a1(m,n) = dres1(m)
           enddo
         enddo
c
c------- complex flux calculation
c
         do n = 1, nmax
           do m = 1, 2
             Cx1(m) = x1(m)
             Cx2(m) = x2(m)
           enddo
c
           do m = 1, 4
             Cq1(m) = q1(m)
             Cq2(m) = q2(m)
             Cres1(m) = 0.d0
             Cres2(m) = 0.d0
           enddo
c
           Cadt1 = adt1
           Cadt2 = adt2
c
           if (ipass.eq.1) then
             Cq1(n) = Cq1(n) + dcmplx(0.d0,scale)
             call Cflux_face(Cx1,Cx2,Cq1,Cq2,Cadt1,Cadt2,Cres1,Cres2)
           elseif (ipass.eq.2) then
             Cx1(n) = Cx1(n) + dcmplx(0.d0,scale)
             call Cflux_face(Cx1,Cx2,Cq1,Cq2,Cadt1,Cadt2,Cres1,Cres2)
           elseif (ipass.eq.3) then
             Cq1(n) = Cq1(n) + dcmplx(0.d0,scale)
             call Cflux_wall(Cx1,Cx2,Cq1,Cres1)
           elseif (ipass.eq.4) then
             Cx1(n) = Cx1(n) + dcmplx(0.d0,scale)
             call Cflux_wall(Cx1,Cx2,Cq1,Cres1)
           endif
c
           do m = 1, 4
             a2(m,n) = aimag(Cres1(m))/scale
             emax = max(emax, abs(a2(m,n)-a1(m,n)))
           enddo
         enddo
c
c------- adjoint flux calculation
c
         do n = 1, 4
           do m = 1, 2
             ax1(m) = 0.d0
             ax2(m) = 0.d0
           enddo
c
           do m = 1, 4
             aq1(m) = 0.d0
             aq2(m) = 0.d0
             ares1(m) = 0.d0
             ares2(m) = 0.d0
           enddo
c
           aadt1 = 0.d0
           aadt2 = 0.d0
c
           ares1(n) = 1.d0
c
           if (ipass.eq.1) then
             call FLUX_FACE_B(x1,x2,q1,aq1,q2,aq2,adt1,aadt1,adt2,aadt2,
     &                        res1,ares1,res2,ares2)
           elseif (ipass.eq.2) then
             call FLUX_FACE_BX(x1,ax1,x2,ax2,q1,aq1,q2,aq2,adt1,aadt1,
     &                         adt2,aadt2,res1,ares1,res2,ares2)
           elseif (ipass.eq.3) then
             call FLUX_WALL_B(x1,x2,q1,aq1,res1,ares1)
           else
             call FLUX_WALL_BX(x1,ax1,x2,ax2,q1,aq1,res1,ares1)
           endif
c
           do m = 1, nmax
             if (nmax.eq.4) then
               a3(n,m) = aq1(m)
             else
               a3(n,m) = ax1(m)
             endif
             emax = max(emax, abs(a3(n,m)-a1(n,m)))
           enddo
         enddo
c
c------- print out comparison
c
         if (ipass.eq.1) then
           write(*,*) 'linear/complex flux_face difference due to q'
           write(*,format4) ((a1(i,j)-a2(i,j),j=1,4),i=1,4)
           write(*,*) 'linear/adjoint flux_face difference due to q'
           write(*,format4) ((a1(i,j)-a3(i,j),j=1,4),i=1,4)
         else if(ipass.eq.2) then
           write(*,*) 'linear/complex flux_face difference due to x1'
           write(*,format2) ((a1(i,j)-a2(i,j),j=1,2),i=1,4)
           write(*,*) 'linear/adjoint flux_face difference due to x1'
           write(*,format2) ((a1(i,j)-a3(i,j),j=1,2),i=1,4)
         else if(ipass.eq.3) then
           write(*,*) 'linear/complex flux_wall difference due to q'
           write(*,format4) ((a1(i,j)-a2(i,j),j=1,4),i=1,4)
           write(*,*) 'linear/adjoint flux_wall difference due to q'
           write(*,format4) ((a1(i,j)-a3(i,j),j=1,4),i=1,4)
         else if(ipass.eq.4) then
           write(*,*) 'linear/complex flux_wall difference due to x1'
           write(*,format2) ((a1(i,j)-a2(i,j),j=1,2),i=1,4)
           write(*,*) 'linear/adjoint flux_wall difference due to x1'
           write(*,format2) ((a1(i,j)-a3(i,j),j=1,2),i=1,4)
         endif
       enddo
c
c************************************************************
c
       do ipass = 1, 2
         nmax = 2 + 2*mod(ipass,2)
c
c------- linear lift calculation
c
         do n = 1, nmax
           do m = 1, 2
             dx1(m) = 0.d0
             dx2(m) = 0.d0
           enddo
c
           do m = 1, 4
             dq1(m) = 0.d0
           enddo
c
           dlift = 0.d0
c
           if (ipass.eq.1) then
             dq1(n) = 1.0d0
             call LIFT_WALL_D(x1,x2,q1,dq1,lift,dlift)
           elseif (ipass.eq.2) then
             dx1(n) = 1.0d0
             call LIFT_WALL_DX(x1,dx1,x2,dx2,q1,dq1,lift,dlift)
           endif
c
           a1(1,n) = dlift
         enddo
c
c------- complex lift calculation
c
         do n = 1, nmax
           do m = 1, 2
             Cx1(m) = x1(m)
             Cx2(m) = x2(m)
           enddo
c
           do m = 1, 4
             Cq1(m) = q1(m)
           enddo
c
           Clift = 0.d0
c
           if (ipass.eq.1) then
             Cq1(n) = Cq1(n) + dcmplx(0.d0,scale)
           elseif (ipass.eq.2) then
             Cx1(n) = Cx1(n) + dcmplx(0.d0,scale)
           endif
c
           call Clift_wall(Cx1,Cx2,Cq1,Clift)
c
           a2(1,n) = aimag(Clift)/scale
           emax = max(emax, abs(a2(1,n)-a1(1,n)))
         enddo
c
c------- adjoint lift calculation
c
         do m = 1, 2
           ax1(m) = 0.d0
           ax2(m) = 0.d0
         enddo
c
         do m = 1, 4
           aq1(m) = 0.d0
         enddo
c
         alift = 1.d0
c
         if (ipass.eq.1) then
           call LIFT_WALL_B(x1,x2,q1,aq1,lift,alift)
         elseif (ipass.eq.2) then
           call LIFT_WALL_BX(x1,ax1,x2,ax2,q1,aq1,lift,alift)
         endif
c
         do n = 1, nmax
           if (ipass.eq.1) then
             a3(1,n) = aq1(n)
           elseif (ipass.eq.2) then
             a3(1,n) = ax1(n)
           endif
           emax = max(emax, abs(a3(1,n)-a1(1,n)))
         enddo
c
c------- print out comparison
c
         if (ipass.eq.1) then
           write(*,*) 'linear/complex lift_wall difference due to q'
           write(*,format4) ((a1(i,j)-a2(i,j),j=1,4),i=1,1)
           write(*,*) 'linear/adjoint lift_wall difference due to q'
           write(*,format4) ((a1(i,j)-a3(i,j),j=1,4),i=1,1)
         elseif (ipass.eq.2) then
           write(*,*) 'linear/complex lift_wall difference due to x1'
           write(*,format2) ((a1(i,j)-a2(i,j),j=1,2),i=1,1)
           write(*,*) 'linear/adjoint lift_wall difference due to x1'
           write(*,format2) ((a1(i,j)-a3(i,j),j=1,2),i=1,1)
         endif
       enddo
c
c*********************************************************************
c
c----- report max error
c
       write(*,'('' *********************************'',/)')
       if(emax .lt. 1e-8) then
         write(*,'('' max error ='',e9.2,'' --- seems OK'',/)') emax
       else
         write(*,'('' max error ='',e9.2,'' --- problem?'',/)') emax
       endif
c
       stop
       end
