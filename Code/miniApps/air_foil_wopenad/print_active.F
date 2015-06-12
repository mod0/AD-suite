      module print_active
        use OAD_active;
        public :: print_array, print_active_array, write_array
     &            write_active_array

        interface print_array
          module procedure print_array1
          module procedure print_array2
        end interface

        interface print_active_array
           module procedure print_active_array1
           module procedure print_active_array2
        end interface

        interface write_array
           module procedure write_array1
           module procedure write_array2
        end interface

        interface write_active_array
           module procedure write_active_array1
           module procedure write_active_array2
        end interface
      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Print real 1d array to the screen
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine print_array1(a, minrow, maxrow)
          real*8 :: a(:)
          integer*4 :: i, l, minrow, maxrow
          l = size(a, 1)

          if (maxrow.eq.0) then
            maxrow = l
          endif

          do i = 1,maxrow
c-------- print the values
              print *, a(i)
          enddo
        end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Print real 2d array to the screen
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine print_array2(a, minrow, maxrow, mincol,
     &   maxcol)
          real*8 :: a(:,:)
          integer*4 :: i, j, l1, l2, minrow, maxrow, mincol, maxcol
          l1 = size(a, 1)
          l2 = size(a, 2)

          if (maxrow.eq.0) then
            maxrow = l1
          endif

          if (maxcol.eq.0) then
            maxcol = l2
          endif

          do i = 1,maxrow
            do j = 1, maxcol
c-------- print the values
              write (*,'(e25.16, a,$)') a(i,j), ' '
            enddo
            print *,
          enddo
        end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Print active 1d array to the screen
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine print_active_array1(a, c, minrow, maxrow)
          type(active) :: a(:)
          integer*4 :: c, i, l, minrow, maxrow
          l = size(a, 1)

          if (maxrow.eq.0) then
            maxrow = l
          endif

          do i = 1,maxrow
            if(c .eq. 0) then
c-------- print the values
              write (*,'(e25.16)'), a(i)%v
            elseif (c .eq. 1) then
c-------- print the derivatives
              write (*,'(e25.16)'), a(i)%d
            endif
          enddo
        end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Print active 2d array to the screen
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine print_active_array2(a, c, minrow, maxrow, mincol,
     &   maxcol)
          type(active) :: a(:,:)
          integer*4 :: c, i, j, l1, l2, minrow, maxrow, mincol, maxcol
          l1 = size(a, 1)
          l2 = size(a, 2)

          if (maxrow.eq.0) then
            maxrow = l1
          endif

          if (maxcol.eq.0) then
            maxcol = l2
          endif

          do i = 1,maxrow
            do j = 1, maxcol
              if(c .eq. 0) then
c-------- print the values
                write (*,'(e25.16, a,$)') a(i,j)%v, ' '
              elseif (c .eq. 1) then
c-------- print the derivatives
                write (*,'(e25.16, a,$)') a(i,j)%d, ' '
              endif
            enddo
            print *,
          enddo
        end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Write real 1d array to the specified filenumber
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine write_array1(a, minrow, maxrow, filenum)
          real*8  :: a(:)
          integer*4 :: i, l, minrow, maxrow, filenum
!          logical :: exists
!          character(len=10) :: filename
!
!c---------Create filename from filenum
!          write(filename, '(A6, I2)') "output", filenum
!
!c---------Check if file exists
!          inquire(file=filename, exist=exists)
!
!          if (exists) then
!            open(filenum, file=filename, status="old", position="append"
!     &                                , action="write")
!          else
!            open(filenum, file=filename, status="new", action="write")
!          endif

          l = size(a, 1)

          if (maxrow.eq.0) then
            maxrow = l
          endif

          print *, minrow, maxrow, filenum, a(1), a(maxrow), a(9703)

          do i = minrow, maxrow
c-------- print the values
            write (filenum,'(e25.16)'), a(i)
          enddo

          close(filenum)
        end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Write real 2d array to the specified filenumber
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine write_array2(a, minrow, maxrow, mincol,
     &   maxcol, filenum)
          real*8  :: a(:,:)
          integer*4 :: i, j, l1, l2, minrow, maxrow, mincol, maxcol,
     &                 filenum
!          logical :: exists
!          character(len=10) :: filename
!
!c---------Create filename from filenum
!          write(filename, '(A6, I2)') "output", filenum
!
!c---------Check if file exists
!          inquire(file=filename, exist=exists)
!
!          if (exists) then
!            open(filenum, file=filename, status="old", position="append"
!     &                                , action="write")
!          else
!            open(filenum, file=filename, status="new", action="write")
!          endif

          l1 = size(a, 1)
          l2 = size(a, 2)

          if (maxrow.eq.0) then
            maxrow = l1
          endif

          if (maxcol.eq.0) then
            maxcol = l2
          endif

          do i = minrow, maxrow
            do j = mincol, maxcol
c-------- print the values
              write (filenum,'(e25.16, a,$)') a(i,j), ' '
            enddo
            write (filenum,*)
          enddo

          close(filenum)
        end subroutine


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Write active 1d array to the specified filenumber
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine write_active_array1(a, c, minrow, maxrow, filenum)
          use OAD_active;
          type(active), dimension(:), intent(in) :: a
          integer :: c, i, l, minrow, maxrow, filenum
!          logical :: exists
!          character(len=10) :: filename
!
!c---------Create filename from filenum
!          write(filename, '(A6, I2)') "output", filenum
!
!c---------Check if file exists
!          inquire(file=filename, exist=exists)
!
!          if (exists) then
!            open(filenum, file=filename, status="old", position="append"
!     &                                , action="write")
!          else
!            open(filenum, file=filename, status="new", action="write")
!          endif

          l = size(a, 1)

          if (maxrow .eq. 0) then
            maxrow = l
          endif

          print *, minrow, maxrow, filenum, a(1), a(maxrow), a(9703)

          do i = minrow, maxrow
            if(c .eq. 0) then
c-------- print the values
              write (filenum,'(e25.16)'), a(i)%v
            elseif (c .eq. 1) then
c-------- print the derivatives
              write (filenum,'(e25.16)'), a(i)%d
            else
              write (filenum,*) , "This is nonsense"
            endif
          enddo
!
!          close(filenum)
        end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Write active 2d array to the specified filenumber
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine write_active_array2(a, c, minrow, maxrow, mincol,
     &   maxcol, filenum)
          type(active) :: a(:,:)
          integer*4 :: c, i, j, l1, l2, minrow, maxrow, mincol, maxcol,
     &                 filenum
          l1 = size(a, 1)
          l2 = size(a, 2)

          if (maxrow.eq.0) then
            maxrow = l1
          endif

          if (maxcol.eq.0) then
            maxcol = l2
          endif

          do i = minrow, maxrow
            do j = mincol, maxcol
              if(c .eq. 0) then
c-------- print the values
                write (filenum,'(e25.16, a,$)') a(i,j)%v, ' '
              elseif (c .eq. 1) then
c-------- print the derivatives
                write (filenum,'(e25.16, a,$)') a(i,j)%d, ' '
              endif
            enddo
            write (filenum,*)
          enddo
        end subroutine

      end module print_active
