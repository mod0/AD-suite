      subroutine Movelib_mov(
     .    coor,coorp,node2d3d,nu,logfac,nubo,nsfac,ns,nt,nfac,nseg,nsp)

      implicit none
c
      include 'param.h'
c
c     logfr is local to the movelib subroutines
c      
      integer kspring, ns, nt, nseg, nfac, is,isp,nsp

      integer nu(4,ntmax), logfr(nsmax), nubo(2,nsegmax)
      integer logfac(nfacmax), nsfac(3,nfacmax),node2d3d(nsp)
      real*8 coor(3,nsmax), xw(3,nsmax),coorp(3,nsp),dimod

      integer dfile
C [llh Nan]      real*4 time

      IF (ns.GT.nsmax) THEN
        WRITE(6,*) 'MOVELIB: ns > nsmax : increase nsmax'
        STOP
      ENDIF
      IF (nt.GT.ntmax) THEN
        WRITE(6,*) 'MOVELIB: nt > ntmax : increase ntmax'
        STOP
      ENDIF
      IF (nfac.GT.nfacmax) THEN
        WRITE(6,*) 'MOVELIB: nfac > nfacmax : increase nfacmax'
        STOP
      ENDIF
      IF (nseg.GT.nsegmax) THEN
        WRITE(6,*) 'MOVELIB: nseg > nsegmax : increase nsegmax'
        STOP
      ENDIF

      do is=1,ns
        xw(1,is) = coor(1,is)
        xw(2,is) = coor(2,is)
        xw(3,is) = coor(3,is)
      enddo

      dimod=0.0d0
      do isp=1,nsp
        is= node2d3d(isp)
        xw(1,is) = coorp(1,isp)
        xw(2,is) = coorp(2,isp)
        xw(3,is) = coorp(3,isp) 
        dimod=dimod
     .    + (xw(1,is)-coor(1,is))*(xw(1,is)-coor(1,is))
     .    + (xw(2,is)-coor(2,is))*(xw(2,is)-coor(2,is))
     .    + (xw(3,is)-coor(3,is))*(xw(3,is)-coor(3,is))
      enddo

      dimod=sqrt(dimod)

      if (dimod.gt.1.0e-10) then
        
        call LOGICFR_mov(nsfac,logfac,logfr,ns,nfac)
        
        call VOLTETRA_mov(coor,nu,ns,nt)
        
        write(6,*) 'Begin with MovMsh_PCG'
        
        kspring=2
c...    kspring= 1, lineal spring / 2, torsional spring
        
        
        call MOVMSH_PCG_mov(coor,nu,logfr,nubo,kspring,ns,nt,
     &    nseg,xw)
        
        call VOLTETRA_mov(xw,nu,ns,nt)
        
        do is=1,ns
          coor(1,is)= xw(1,is) 
          coor(2,is)= xw(2,is) 
          coor(3,is)= xw(3,is) 
        enddo

      else

        write (6,*) 'No skin displacement'
        write (6,*) 'Movelib returned!'
        
      end if
      
      end
