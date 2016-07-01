      subroutine suctrl(volcontrol)
c
c  This function performs a surface integral of the gradient which is used to keep
c  the optimized profile volume constant
c
      
      include 'Param3D.h'
      include 'Paramopt3D.h'

      real*8 volcontrol,proynn,vncoq(3,nnsp),ctrl(nnsp)
      integer isp

      volcontrol=0.d0
      if (isuct.eq.0) return

c      call normcoq(ctrl,vncoq)
      
      proynn=1.0d0
      do isp=1,nsp
c        proynn=
c     .    vnon(1,isp)*vncoq(1,isp)+vnon(2,isp)*vncoq(2,isp)
c     .    + vnon(3,isp)*vncoq(3,isp)
        volcontrol=volcontrol+airesp(isp)*grad(isp)*proynn
      end do
      
      end
