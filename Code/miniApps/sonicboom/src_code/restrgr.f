      subroutine restrgr
c
c     This subroutine restricts the gradient according to a given criteria,
c     in order to limit the modifications to a part of the shape.
c
      include 'Param3D.h'
      include 'Paramopt3D.h'

      real*8 swean,pian,macan,resan
      integer isp,iwhico

c...  CHECK THE SWEEP ANGLE (ONLY FOR VARIABLE SWEEPS WINGS)

c      pian  = acos(-1.0d0)
c      macan = 180.0d0*asin(1.0d0/xmach)/pian
c      resan = macan - 5.0d0                                 !5 deg smaller than the mach angle
c      
c      do isp=1,nsp
c
c        swean= 180.0d0*acos(vnon(1,isp))/pian
c        if (swean.lt.resan) grad(isp)=0.0d0
c
c        write (6,*) swean,macan,resan
c        
c      end do

c...  CHECK DISTANCE TO THE SYMETRY AXE

      iwhico= int(faczz(8))
      resan = 0.4d0
      macan = 0.5d0
      
      do isp=1,nsp

        swean= coorp(iwhico,isp)

        if (swean.lt.macan) then
          if (swean.gt.resan) then
            grad(isp)=grad(isp)*(swean-resan)/(macan-resan)
          else
            grad(isp)=0.0d0
          end if
        end if

      end do
      

      end
