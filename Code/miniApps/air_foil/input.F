        subroutine input(maxnode,maxcell,maxedge, nnode,ncell,nedge,
     &                   x, q,cell, edge,ecell,boun)
        implicit none
c
        integer*4 maxnode,maxcell,maxedge, nnode,ncell,nedge,
     &            cell(4,maxcell), edge(2,maxedge),
     &            ecell(2,maxedge), boun(maxedge)
        real*8    x(2,maxnode), q(4,maxcell)
c
#include "const.inc"
c
        integer in,ie,ic,ipde
        real*8  p,r,u,e, xt,yt
c
c------ set universal constants
c
        gam = 1.4d0
        gm1 = gam - 1.0d0
        cfl = 0.9d0
        eps = 0.05d0
c
c------ read in data from grid file
c
        open(1,file='grid.dat',status='old')
c
        read(1,*) nnode, ncell, nedge
        ncell = ncell+1  ! far-field dummy cell
        if(nnode .gt. maxnode) stop 'maxnode too small'
        if(ncell .gt. maxcell) stop 'maxcell too small'
        if(nedge .gt. maxedge) stop 'maxedge too small'
c
        do in = 1, nnode
          read(1,*) x(1,in),x(2,in)
        enddo
c
        do ic = 1, ncell-1
          read(1,*) cell(1,ic),cell(2,ic),cell(3,ic),cell(4,ic)
        enddo
c
        do ie = 1, nedge
          read(1,*) edge(1,ie),edge(2,ie),
     &              ecell(1,ie),ecell(2,ie),boun(ie)
          if(boun(ie).eq.2) then
           ecell(2,ie) = ncell
            boun(ie) = 0
          endif
        enddo
c
        close(1)
c
c------ read in data from flow file, initialising if necessary
c
        open(1,file='flow.dat',status='old')
        read(1,*) p, r, mach, alpha
        alpha = alpha*datan(1.0d0)/45.0d0
        p = 1.0d0
        r = 1.0d0
        u = dsqrt(gam*p/r)*mach
        e = p/(r*gm1) + 0.5d0*u**2
c
        do ic = 1, ncell
          q(1,ic) = r
          q(2,ic) = r*u
          q(3,ic) = 0.d0
          q(4,ic) = r*e
        enddo
c
        do ic = 1, ncell
          read(1,*,err=999,end=999) (q(ipde,ic),ipde=1,4)
        enddo
 999    close(1)
c
c------ rotate grid to specified angle of attack
c
        do in = 1, nnode
          xt = x(1,in)
          yt = x(2,in)
          x(1,in) = cos(alpha)*xt + sin(alpha)*yt
          x(2,in) = cos(alpha)*yt - sin(alpha)*xt
        enddo
c
        return
        end


        subroutine output(ncell, q)
        implicit none
c
        integer*4 ncell, ic, ipde
        real*8    q(4,ncell), p, r, mach,alpha
c
        open(1,file='flow.dat',status='old')
        read(1,*) p, r, mach, alpha
        close(1)
c
        open(1,file='flow.dat',status='old')
        write(1,*) p, r, mach, alpha
        do ic = 1, ncell-1
          write(1,*) (q(ipde,ic),ipde=1,4)
        enddo
        close(1)
c
        return
        end
