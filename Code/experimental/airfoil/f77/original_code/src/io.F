      subroutine read_grid(maxnode, maxcell, maxedge, nnode,
     $     ncell, nedge, x, cell, edge, ecell, boun, data_directory)
      implicit none
c     
      integer*4 in, ie, ic
      integer*4   maxnode, maxcell, maxedge, nnode, ncell, nedge,
     &     cell(4,maxcell), edge(2,maxedge),
     &     ecell(2,maxedge), boun(maxedge)
      real*8      x(2,maxnode)
      character(len = 100) data_directory
c     ===============================================
c     INCLUDE VARIABLES IN COMMON BLOCK
c     ===============================================
#include "const.inc"
c           


c     ===============================================
c     read in data from grid file
c     ===============================================
      open(1, file = trim(adjustl(data_directory))//'/grid.dat',
     $     status='old')
c     
      read(1,*) nnode, ncell, nedge
      ncell = ncell+1           ! far-field dummy cell
      if(nnode .gt. maxnode) stop 'maxnode too small'
      if(ncell .gt. maxcell) stop 'maxcell too small'
      if(nedge .gt. maxedge) stop 'maxedge too small'
c     
      do in = 1, nnode
         read(1,*,err=999,end=999) x(1,in), x(2,in)
      enddo
c     
      do ic = 1, ncell-1
         read(1,*,err=999,end=999) cell(1,ic), cell(2,ic),
     $        cell(3,ic), cell(4,ic)
      enddo
c     
      do ie = 1, nedge
         read(1,*,err=999,end=999) edge(1,ie), edge(2,ie),
     &        ecell(1,ie), ecell(2,ie), boun(ie)
         if(boun(ie).eq.2) then
            ecell(2,ie) = ncell
            boun(ie) = 0
         endif
      enddo
c     
 999  close(1)
      end subroutine read_grid
      


      subroutine read_flow(maxcell, ncell, p, r, u, e, q,
     $     data_directory)
      implicit none
      integer*4  ic, ipde
      integer*4  maxcell, ncell
      real*8   p, r, u, e
      real*8   q(4, maxcell)
      character(len = 100) data_directory
c     ===============================================
c     INCLUDE VARIABLES IN COMMON BLOCK
c     ===============================================
#include "const.inc"
c           

      if(ncell .gt. maxcell) stop 'maxcell too small'
    
      do ic = 1, ncell
         q(1,ic) = r
         q(2,ic) = r*u
         q(3,ic) = 0.d0
         q(4,ic) = r*e
      enddo          

c     ======================================================
c     read in data from flow file, initialising if necessary
c     ======================================================
      open(1,file=trim(adjustl(data_directory))//'/flow.dat',
     $     status='old')

      do ic = 1, ncell
         read(1,*,err=999,end=999) (q(ipde,ic),ipde=1,4)
      enddo
 999  close(1)

      return
      end subroutine read_flow




      subroutine write_flow(ncell, q, results_directory)
      implicit none
c     
      integer*4 ncell, ic, ipde
      real*8    q(4,ncell)
      character(len = 100) results_directory
c          
      open(1,file=trim(adjustl(results_directory))//'/updated_flow.dat',
     $     status='old')
      do ic = 1, ncell-1
         write(1,*,err=999) (q(ipde,ic),ipde=1,4)
      enddo
 999  close(1)
c     
      return
      end subroutine write_flow
