      subroutine voltetra_mov(coor,nu,ns,nt)

      implicit none

      include 'param.h'

      integer ns, nt
      integer nu(4, ntmax)
      integer k, nut(4), jt, is
      real*8  coor(3,nsmax), x(4), y(4), z(4)
      real*8 z12, z13, z14, z23, z24, z34, b1, b2, b3, b4, vol6
      real*8 volt, voltmin, voltmax


      voltmin = 1.e+10
      voltmax = -1.e+10

      do jt=1,nt
        do k=1,4
            nut(k)                = nu(k,jt)

            is                     = nut(k)
c
            x(k)                  = coor(1,is)
            y(k)                  = coor(2,is)
            z(k)                  = coor(3,is)

        enddo
c
         z12                       = z(2) - z(1)
         z13                       = z(3) - z(1)
         z14                       = z(4) - z(1)
         z23                       = z(3) - z(2)
         z24                       = z(4) - z(2)
         z34                       = z(4) - z(3)
         b1                        = y(2)*z34 - y(3)*z24 + y(4)*z23
         b2                        =-y(1)*z34 + y(3)*z14 - y(4)*z13
         b3                        = y(1)*z24 - y(2)*z14 + y(4)*z12
         b4                        =-y(1)*z23 + y(2)*z13 - y(3)*z12
c
         vol6                      = x(1)*b1 + x(2)*b2 +
     &                               x(3)*b3 + x(4)*b4

         volt = vol6/6.

         voltmin = min(voltmin,volt)
         voltmax = max(voltmax,volt)

       enddo


      write(6,*) 'voltmin: ',voltmin
      write(6,*) 'voltmax: ',voltmax


      return
      end
