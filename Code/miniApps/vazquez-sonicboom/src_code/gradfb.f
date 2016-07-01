      SUBROUTINE GRADFB(x, y, z, b, c, d, vol6)
c---------------------------------------------------------------------   
      IMPLICIT NONE
c---------------------------------------------------------------------
c     Local variables definition
      REAL*8 b(4), c(4), d(4), x(4), y(4), z(4), vol6
      REAL*8 x1  , y1  , z1  , x2  , y2  , z2  
      REAL*8 x3  , y3  , z3  , x4  , y4  , z4
      REAL*8 z12 , z13 , z14 , z23 , z24 , z34
      REAL*8 y12 , y13 , y14 , y23 , y24 , y34
c
      x1                           = x(1)
      y1                           = y(1)
      z1                           = z(1)
      x2                           = x(2)
      y2                           = y(2)
      z2                           = z(2)
      x3                           = x(3)
      y3                           = y(3)
      z3                           = z(3)
      x4                           = x(4)
      y4                           = y(4)
      z4                           = z(4)
c
      z12                          = z2 - z1
      z13                          = z3 - z1
      z14                          = z4 - z1
      z23                          = z3 - z2
      z24                          = z4 - z2
      z34                          = z4 - z3
c
      b(1)                         = y2*z34 - y3*z24 + y4*z23
      b(2)                         =-y1*z34 + y3*z14 - y4*z13
      b(3)                         = y1*z24 - y2*z14 + y4*z12
      b(4)                         =-y1*z23 + y2*z13 - y3*z12
c
      c(1)                         =-x2*z34 + x3*z24 - x4*z23
      c(2)                         = x1*z34 - x3*z14 + x4*z13
      c(3)                         =-x1*z24 + x2*z14 - x4*z12
      c(4)                         = x1*z23 - x2*z13 + x3*z12
c
      y12                          = y2 - y1
      y13                          = y3 - y1
      y14                          = y4 - y1
      y23                          = y3 - y2
      y24                          = y4 - y2
      y34                          = y4 - y3
c
      d(1)                         = x2*y34 - x3*y24 + x4*y23
      d(2)                         =-x1*y34 + x3*y14 - x4*y13
      d(3)                         = x1*y24 - x2*y14 + x4*y12
      d(4)                         =-x1*y23 + x2*y13 - x3*y12
c
c     Six times the volume
c
      vol6                         = x1*b(1) + x2*b(2) + 
     &                               x3*b(3) + x4*b(4)
c
      RETURN
      END
