      subroutine template()
      use OAD_tape
      use OAD_rev
      use OAD_cp
      use w2f__types
      !use OAD_cp
      !use OAD_tape
      !use OAD_rev
      use conj_grad_mod
      use conj_grad_ad_mod
      use stream_vel_variables
!$TEMPLATE_PRAGMA_DECLARATIONS

      type(modeType) :: our_orig_mode

      real(8), dimension(n) :: x_p
      real(8), dimension(n) :: b_p
      real(8), dimension(n,3) :: A_p
      real(8), dimension(n) :: x_d
      real(8), dimension(n) :: b_d
      real(8), dimension(n,3) :: A_d
      integer :: j

      if (our_rev_mode%plain) then
! set up for plain execution
         our_orig_mode=our_rev_mode
         our_rev_mode%arg_store=.FALSE.
         our_rev_mode%arg_restore=.FALSE.
         our_rev_mode%plain=.TRUE.
         our_rev_mode%tape=.FALSE.
         our_rev_mode%adjoint=.FALSE.
         b_p  = b%v
         x_p  = x%v
         A_p  = A%v
         call solve(x_p,b_p,A_p)
         x%v = x_p 
! reset the mode
         our_rev_mode=our_orig_mode
      end if
      if (our_rev_mode%tape) then
! set up for plain execution
         our_orig_mode=our_rev_mode
         our_rev_mode%arg_store=.FALSE.
         our_rev_mode%arg_restore=.FALSE.
         our_rev_mode%plain=.TRUE.
         our_rev_mode%tape=.FALSE.
         our_rev_mode%adjoint=.FALSE.
         b_p  = b%v
         x_p  = x%v
         A_p  = A%v
         !call oad_tape_push(b_p)      
         !call oad_tape_push(A_p)
         !call push_s1(b_p)
         !call push_s2(A_p) 
         !WRITE(*,*) 'Size(b,1) = ', size(b,1)
         !do i=1,n
         !   WRITE(*,*) 'x(', i, ') = ', b_p(i)
         !end do 
         do i=1,n
            call push_s0(b_p(i))
         end do 
         do i=1,n
           do j=1,3
            call push_s0(A_p(i,j))
           end do
         end do     
         call solve(x_p,b_p,A_p)
         !call oad_tape_push(x_p)     
         !call push_s1(x_p)
         !WRITE(*,*) 'Size(x,1) = ', size(x,1)  
         do i=1,n
            call push_s0(x_p(i))
         end do   
         x%v = x_p
! reset the mode
         our_rev_mode=our_orig_mode
      end if
      if (our_rev_mode%adjoint) then
         !call oad_tape_pop(x_p)     
         !call oad_tape_pop(A_p)     
         !call oad_tape_pop(b_p)
         !call pop_s1(x_p)     
         !call pop_s2(A_p)     
         !call pop_s1(b_p)
         do i=n, 1, -1
            call pop_s0(x_p(i))
         end do
         do i=n, 1, -1
           do j=3,1, -1
            call pop_s0(A_p(i,j))
           end do
         end do 
         do i=n,1,-1
            call pop_s0(b_p(i))
         end do   
! set up for plain execution
         our_orig_mode=our_rev_mode
         our_rev_mode%arg_store=.FALSE.
         our_rev_mode%arg_restore=.FALSE.
         our_rev_mode%plain=.TRUE.
         our_rev_mode%tape=.FALSE.
         our_rev_mode%adjoint=.FALSE.
         b_d  = b%d
         x_d  = x%d
         A_d  = A%d
      call adsolve( x_p, x_d, b_p, b_d, a_p, a_d )
! reset the mode
         b%d  = b_d
         x%d  = x_d
         A%d  = A_d
         our_rev_mode=our_orig_mode
      end if 
      end subroutine template
