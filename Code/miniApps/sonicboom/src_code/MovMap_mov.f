c
      subroutine MovMap_mov(nufaces,nucycles,lowlambda,highlambda,
     &                  map_pointer)
c
      implicit none
c
      real*8 lowlambda,highlambda
      integer nufaces,nucycles,map_pointer(4,3,4)
c
c     contains mapping for the tetra computation
c     in improved MovMsh
c
      nufaces = 4

      nucycles = 3
c
      lowlambda  = 0.0

      highlambda = 1.0
c
c
c     ------> mapping needed for local tetra computations
c
c
c   ---> face 1, cycle 1
c
      map_pointer(1,1,1) = 1
      map_pointer(1,1,2) = 2
      map_pointer(1,1,3) = 3
      map_pointer(1,1,4) = 4
c
c   ---> face 1, cycle 2
c
      map_pointer(1,2,1) = 3
      map_pointer(1,2,2) = 1
      map_pointer(1,2,3) = 2
      map_pointer(1,2,4) = 4
c
c   ---> face 1, cycle 3
c
      map_pointer(1,3,1) = 2
      map_pointer(1,3,2) = 3
      map_pointer(1,3,3) = 1
      map_pointer(1,3,4) = 4
c
c   ---> face 2, cycle 1
c
      map_pointer(2,1,1) = 2
      map_pointer(2,1,2) = 1
      map_pointer(2,1,3) = 4
      map_pointer(2,1,4) = 3
c
c   ---> face 2, cycle 2
c
      map_pointer(2,2,1) = 4
      map_pointer(2,2,2) = 2
      map_pointer(2,2,3) = 1
      map_pointer(2,2,4) = 3
c
c   ---> face 2, cycle 3
c
      map_pointer(2,3,1) = 1
      map_pointer(2,3,2) = 4
      map_pointer(2,3,3) = 2
      map_pointer(2,3,4) = 3
c
c   ---> face 3, cycle 1
c
      map_pointer(3,1,1) = 4
      map_pointer(3,1,2) = 1
      map_pointer(3,1,3) = 3
      map_pointer(3,1,4) = 2
c
c   ---> face 3, cycle 2
c
      map_pointer(3,2,1) = 3
      map_pointer(3,2,2) = 4
      map_pointer(3,2,3) = 1
      map_pointer(3,2,4) = 2
c
c   ---> face 3, cycle 3
c
      map_pointer(3,3,1) = 1
      map_pointer(3,3,2) = 3
      map_pointer(3,3,3) = 4
      map_pointer(3,3,4) = 2
c
c   ---> face 4, cycle 1
c
      map_pointer(4,1,1) = 2
      map_pointer(4,1,2) = 4
      map_pointer(4,1,3) = 3
      map_pointer(4,1,4) = 1
c
c   ---> face 4, cycle 2
c
      map_pointer(4,2,1) = 4
      map_pointer(4,2,2) = 3
      map_pointer(4,2,3) = 2
      map_pointer(4,2,4) = 1
c
c   ---> face 4, cycle 3
c
      map_pointer(4,3,1) = 3
      map_pointer(4,3,2) = 2
      map_pointer(4,3,3) = 4
      map_pointer(4,3,4) = 1
c
c
c
      return
      end
