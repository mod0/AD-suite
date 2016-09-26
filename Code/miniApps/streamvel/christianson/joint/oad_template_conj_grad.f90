      subroutine template()
      use OAD_tape
      use OAD_rev
      use OAD_cp
      use w2f__types
      use conj_grad_mod
      use conj_grad_ad_mod
      use stream_vel_variables
!$TEMPLATE_PRAGMA_DECLARATIONS
      type(modeType) :: our_orig_mode
! checkpointing stacks and offsets
      integer :: cp_loop_variable_1
! floats 'F'
      double precision, dimension(:), allocatable, save :: theArgFStack
      integer, save :: theArgFStackoffset=0, theArgFStackSize=0

      real(8), dimension(n) :: x_p
      real(8), dimension(n) :: b_p
      real(8), dimension(n,3) :: A_p
      real(8), dimension(n) :: x_d
      real(8), dimension(n) :: b_d
      real(8), dimension(n,3) :: A_d
      integer :: j

      if (our_rev_mode%arg_store) then
! store arguments
        do i=1,n
          call cp_store_real_scalar(x(i)%v,theArgFStack,theArgFStackoffset,theArgFStackSize)
        end do
        do i=1,n
          call cp_store_real_scalar(b(i)%v,theArgFStack,theArgFStackoffset,theArgFStackSize)
        end do
        do i=1,n
          do j=1,3
            call cp_store_real_scalar(A(i,j)%v,theArgFStack,theArgFStackoffset,theArgFStackSize)
          end do
        end do
      end if
      if (our_rev_mode%arg_restore) then
! restore arguments
        do i=n, 1, -1
          do j=3, 1, -1
            A(i,j)%v = theArgFStack(theArgFStackoffset)
            theArgFStackoffset = theArgFStackoffset-1
          end do
        end do
        do cp_loop_variable_1 = n,1,-1
          b(cp_loop_variable_1)%v = theArgFStack(theArgFStackoffset)
          theArgFStackoffset = theArgFStackoffset-1
        end do
        do cp_loop_variable_1 = n,1,-1
          x(cp_loop_variable_1)%v = theArgFStack(theArgFStackoffset)
          theArgFStackoffset = theArgFStackoffset-1
        end do
      end if
      if (our_rev_mode%res_store) then
        WRITE(*,*) 'UNHANDLED template: our_rev_mode%res_store called'
      end if
      if (our_rev_mode%res_restore) then
        WRITE(*,*) 'UNHANDLED template: our_rev_mode%res_restore called'
      end if
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
        our_orig_mode = our_rev_mode
        our_rev_mode%arg_store=.TRUE.
        our_rev_mode%arg_restore=.FALSE.
        our_rev_mode%plain=.TRUE.
        our_rev_mode%tape=.FALSE.
        our_rev_mode%adjoint=.FALSE.
        b_p  = b%v
        x_p  = x%v
        A_p  = A%v
        b_d  = b%d
        x_d  = x%d
        A_d  = A%d
        call adsolve( x_p, x_d, b_p, b_d, a_p, a_d )
        b%d  = b_d
        x%d  = x_d
        A%d  = A_d
! reset the mode
        our_rev_mode=our_orig_mode
!adjoint end
      end if 
      end subroutine template
