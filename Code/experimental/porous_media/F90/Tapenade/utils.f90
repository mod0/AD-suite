module utils
  public :: sort
  
  interface sort
    module procedure sort1
    module procedure sort2
  end interface sort
contains
! from: http://www.cs.mtu.edu/~shene/courses/cs201/notes/chap08/sorting.f90
! --------------------------------------------------------------------
! integer function  findminimum():
!    this function returns the location of the minimum in the section
! between start and end.
! --------------------------------------------------------------------

   integer function  findminimum(x, start, end)
      implicit  none
      integer, dimension(1:), intent(in) :: x
      integer, intent(in)                :: start, end
      integer                            :: minimum
      integer                            :: location
      integer                            :: i

      minimum  = x(start)               ! assume the first is the min
      location = start                  ! record its position
      do i = start+1, end               ! start with next elements
         if (x(i) < minimum) then       !   if x(i) less than the min?
            minimum  = x(i)             !      yes, a new minimum found
            location = i                !      record its position
         end if
      end do
      findminimum = location            ! return the position
   end function  findminimum

! --------------------------------------------------------------------
! subroutine  swap():
!    this subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   subroutine  swap1(a, b)
      implicit  none
      integer, intent(inout) :: a, b
      integer                :: temp

      temp = a
      a    = b
      b    = temp
   end subroutine  swap1

! --------------------------------------------------------------------
! subroutine  sort():
!    this subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   subroutine  sort1(x, size)
      implicit  none
      integer, dimension(1:), intent(inout) :: x
      integer, intent(in)                   :: size
      integer                               :: i
      integer                               :: location

      do i = 1, size-1                  ! except for the last
         location = findminimum(x, i, size)     ! find min from this to last
         call  swap1(x(i), x(location))  ! swap this and the minimum
      end do
   end subroutine  sort1
   
! --------------------------------------------------------------------
! subroutine  swap(): [Key, Value]'
!    this subroutine swaps the [key, value]' of its two formal arguments.
! --------------------------------------------------------------------

   subroutine  swap2(a, b)
      implicit  none
      integer, dimension(2), intent(inout) :: a, b
      integer, dimension(2)                :: temp

      temp = a
      a    = b
      b    = temp
   end subroutine  swap2

! --------------------------------------------------------------------
! subroutine  sort(): [Key, Value]'
!    this subroutine receives a 2xN array x() and sorts it into ascending
! order based on the key stored in the first row
! This is currently doing insertion sort, need to do quick/merge sort
! --------------------------------------------------------------------

   subroutine  sort2(x, size)
      implicit  none
      integer, dimension(2,size), intent(inout) :: x
      integer, intent(in)                       :: size
      integer                                   :: i
      integer                                   :: location

      do i = 1, size-1                           ! except for the last
         location = findminimum(x(1,:), i, size) ! find min key from this to last
         call  swap2(x(:,i), x(:,location))      ! swap this and the minimum
      end do
   end subroutine  sort2
end module utils