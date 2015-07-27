      module print_active
*        use OAD_active;
        public :: print_array, print_active_array, write_array,
     &            write_active_array

        interface print_array
          module procedure print_array1
          module procedure print_array2
          module procedure print_array3
          module procedure print_array4
        end interface

*        interface print_active_array
*           module procedure print_active_array1
*           module procedure print_active_array2
*        end interface

        interface write_array
           module procedure write_array1
           module procedure write_array2
        end interface

*        interface write_active_array
*           module procedure write_active_array1
*           module procedure write_active_array2
*        end interface
      contains

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       Print real 1d array to the screen
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine print_array1(a, minrow, maxrow)
          real*8 :: a(:)
          integer*4 :: i, l, minrow, maxrow, maxrow_temp

          l = size(a, 1)

          if (maxrow.eq.0) then
            maxrow_temp = l
          else
            maxrow_temp = maxrow
          endif

          do i = minrow, maxrow_temp
!-------- print the values
              print *, a(i)
          enddo
        end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       Print real 2d array to the screen
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine print_array2(a, minrow, maxrow, mincol,
     &   maxcol)
          real*8 :: a(:,:)
          integer*4 :: i, j, l1, l2, minrow, maxrow, mincol, maxcol,
     &                 maxrow_temp, maxcol_temp

          l1 = size(a, 1)
          l2 = size(a, 2)

          if (maxrow.eq.0) then
            maxrow_temp = l1
          else
            maxrow_temp = maxrow
          endif

          if (maxcol.eq.0) then
            maxcol_temp = l2
          else
            maxcol_temp = maxcol
          endif

          do i = minrow, maxrow_temp
            do j = mincol, maxcol_temp
!-------- print the values
              write (*,'(e23.16, a,$)') a(i,j), ' '
            enddo
            print *, ""
          enddo
        end subroutine


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       Print real 3d array to the screen
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine print_array3(a, minrow, maxrow, mincol,
     &   maxcol, minstack, maxstack)
          real*8 :: a(:,:,:)
          integer*4 :: i, j, l1, l2, minrow, maxrow, mincol, maxcol,
     &                 k, l3, minstack, maxstack,
     &                 maxrow_temp, maxcol_temp, maxstack_temp

          l1 = size(a, 1)
          l2 = size(a, 2)
          l3 = size(a, 3)

          if (maxrow.eq.0) then
            maxrow_temp = l1
          else
            maxrow_temp = maxrow
          endif

          if (maxcol.eq.0) then
            maxcol_temp = l2
          else
            maxcol_temp = maxcol
          endif

          if (maxstack.eq.0) then
            maxstack_temp = l3
          else
            maxstack_temp = maxstack
          endif

          do k = minstack, maxstack_temp
              print *, "a(:,:," , k , "):"
              do i = minrow, maxrow_temp
                do j = mincol, maxcol_temp
    !-------- print the values
                  write (*,'(e23.16, a,$)') a(i,j,k), ' '
                enddo
                print *, ""
              enddo
          end do
        end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       Print real 4d array to the screen
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine print_array4(a, minrow, maxrow, mincol,
     &   maxcol, minstack, maxstack, minbox, maxbox)
          real*8 :: a(:,:,:,:)
          integer*4 :: i, j, l1, l2, minrow, maxrow, mincol, maxcol,
     &                 k, l3, minstack, maxstack, maxstack_temp,
     &                 l, l4, minbox, maxbox, maxbox_temp,
     &                 maxrow_temp, maxcol_temp

          l1 = size(a, 1)
          l2 = size(a, 2)
          l3 = size(a, 3)
          l4 = size(a, 4)

          if (maxrow.eq.0) then
            maxrow_temp = l1
          else
            maxrow_temp = maxrow
          endif

          if (maxcol.eq.0) then
            maxcol_temp = l2
          else
            maxcol_temp = maxcol
          endif

          if (maxstack.eq.0) then
            maxstack_temp = l3
          else
            maxstack_temp = maxstack
          endif

          if (maxbox.eq.0) then
            maxbox_temp = l4
          else
            maxbox_temp = maxbox
          endif

          do l = minbox, maxbox_temp
              do k = minstack, maxstack_temp
                  print *, "a(:,:," , k , ",", l, "):"
                  do i = minrow, maxrow_temp
                    do j = mincol, maxcol_temp
!-------- print the values
                      write (*,'(e23.16, a,$)') a(i,j,k,l), ' '
                    enddo
                    print *, ""
                  enddo
              end do
          end do
        end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       Print active 1d array to the screen
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*        subroutine print_active_array1(a, c, minrow, maxrow)
*          use OAD_active
*          type(active), dimension(:), intent(in) :: a
*          integer*4 :: c, i, l, minrow, maxrow
*          l = size(a, 1)
*
*          if (maxrow.eq.0) then
*            maxrow = l
*          endif
*
*          do i = 1,maxrow
*            if(c .eq. 0) then
*!-------- print the values
*              write (*,'(e23.16)'), a(i)%v
*            elseif (c .eq. 1) then
*!-------- print the derivatives
*              write (*,'(e23.16)'), a(i)%d
*            endif
*          enddo
*        end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       Print active 2d array to the screen
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*        subroutine print_active_array2(a, c, minrow, maxrow, mincol,
*     &   maxcol)
*          use OAD_active
*          type(active), dimension(:,:), intent(in) :: a
*          integer*4 :: c, i, j, l1, l2, minrow, maxrow, mincol, maxcol
*          l1 = size(a, 1)
*          l2 = size(a, 2)
*
*          if (maxrow.eq.0) then
*            maxrow = l1
*          endif
*
*          if (maxcol.eq.0) then
*            maxcol = l2
*          endif
*
*          do i = 1,maxrow
*            do j = 1, maxcol
*              if(c .eq. 0) then
*!-------- print the values
*                write (*,'(e23.16, a,$)') a(i,j)%v, ' '
*              elseif (c .eq. 1) then
*!-------- print the derivatives
*                write (*,'(e23.16, a,$)') a(i,j)%d, ' '
*              endif
*            enddo
*            print *,
*          enddo
*        end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       Write real 1d array to the specified filenumber
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine write_array1(a, minrow, maxrow, filenum)
          real*8  :: a(:)
          integer*4 :: i, l, minrow, maxrow, filenum, maxrow_temp

          l = size(a, 1)

          if (maxrow.eq.0) then
            maxrow_temp = l
          else
            maxrow_temp = maxrow
          endif

          do i = minrow, maxrow_temp
!-------- print the values
            write (filenum,'(e23.16)'), a(i)
          enddo
        end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       Write real 2d array to the specified filenumber
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine write_array2(a, minrow, maxrow, mincol,
     &   maxcol, filenum)
          real*8  :: a(:,:)
          integer*4 :: i, j, l1, l2, minrow, maxrow, mincol, maxcol,
     &                 filenum, maxrow_temp, maxcol_temp

          l1 = size(a, 1)
          l2 = size(a, 2)

          if (maxrow.eq.0) then
            maxrow_temp = l1
          else
            maxrow_temp = maxrow
          endif

          if (maxcol.eq.0) then
            maxcol_temp = l2
          else
            maxcol_temp = maxcol
          endif

          do i = minrow, maxrow_temp
            do j = mincol, maxcol_temp
!-------- print the values
              write (filenum,'(e23.16,a,$)') a(i,j), ","
            enddo
            write (filenum, *) ""
          enddo
        end subroutine


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       Write active 1d array to the specified filenumber
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*        subroutine write_active_array1(a, c, minrow, maxrow, filenum,
*     &   dim1)
*          use OAD_active
*          integer :: c, i, l, minrow, maxrow, filenum, dim1
*          type(active), dimension(dim1), intent(in) :: a
*
*
*          l = size(a, 1)
*
*          if (maxrow .eq. 0) then
*            maxrow = l
*          endif
*
*          do i = minrow, maxrow
*            if(c .eq. 0) then
*!-------- print the values
*              write (filenum,'(e23.16)'), a(i)%v
*            elseif (c .eq. 1) then
*!-------- print the derivatives
*              write (filenum,'(e23.16)'), a(i)%d
*            endif
*          enddo
*        end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       Write active 2d array to the specified filenumber
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*        subroutine write_active_array2(a, c, minrow, maxrow, mincol,
*     &   maxcol, filenum, dim1, dim2)
*          use OAD_active
*          integer*4 :: c, i, j, l1, l2, minrow, maxrow, mincol, maxcol,
*     &                 filenum, dim1, dim2
*          type(active), dimension(dim1,dim2), intent(in) :: a
*
*          l1 = size(a, 1)
*          l2 = size(a, 2)
*
*          if (maxrow.eq.0) then
*            maxrow = l1
*          endif
*
*          if (maxcol.eq.0) then
*            maxcol = l2
*          endif
*
*          do i = minrow, maxrow
*            do j = mincol, maxcol
*              if(c .eq. 0) then
*!-------- print the values
*                write (filenum,'(e23.16)') a(i,j)%v
*              elseif (c .eq. 1) then
*!-------- print the derivatives
*                write (filenum,'(e23.16)') a(i,j)%d
*              endif
*            enddo
*          enddo
*        end subroutine

      end module print_active
