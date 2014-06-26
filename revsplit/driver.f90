program driver
  use OAD_active
  use OAD_rev
  use stream_vel_variables
  implicit none 
  external head

  TYPE (active), dimension(79) :: bb 
  TYPE (active), dimension(80) :: u
  TYPE (active) :: fc, fc2
  integer :: ii, jj
  real(8) :: fdfc, fc_0, accuracyAD
  
  external stream_vel
  accuracyAD=0.0
!  n = 79
!  initialize
  do ii = 1, n
    bb(ii)%v    = 0.
    bb(ii)%d = 0.0D0
    u(ii)%d =  0.0D0
  end do
  u(n+1)%d =  0.0D0
  fc%v = 0.0
  fc%d = 1.0D0
  fc2%v = 0.0
! call adjoint model
our_rev_mode%arg_store=.FALSE.
our_rev_mode%arg_restore=.FALSE.
our_rev_mode%plain=.FALSE.
our_rev_mode%adjoint=.FALSE.
  our_rev_mode%tape=.TRUE.
         call stream_vel( u, bb, fc )
our_rev_mode%arg_store=.FALSE.
our_rev_mode%arg_restore=.FALSE.
our_rev_mode%plain=.FALSE.
  our_rev_mode%tape=.FALSE.
  our_rev_mode%adjoint=.TRUE.
!  do ii = 1, n
!    bb(ii)%v    = 0.
!    bb(ii)%d = 0.0D0
!    u(ii)%d =  0.0D0
!  end do
!  u(n+1)%d =  0.0D0
!  fc%v = 0.	
!  fc%d = 1.0D0
        call stream_vel( u, bb, fc )


       our_rev_mode%arg_store=.FALSE.
       our_rev_mode%arg_restore=.FALSE.
       our_rev_mode%plain=.TRUE.
       our_rev_mode%tape=.FALSE.
       our_rev_mode%adjoint=.FALSE.
       

 print *, '    position       velocity'
 print *, '----------------------------------------------------------------'
        do ii=1,n+1
         print *,ii,u(ii)%v
        enddo
        print *, '                                    position       ', &
         '            dfc/dM [ADM]           dfc/dM [f.d.]            relative accuracy'
        print *, '----------------------------------------------------------------', &
         '-----------------------------------------------------------------------------'
        do ii=1,n
         our_rev_mode%arg_store=.FALSE.
         our_rev_mode%arg_restore=.FALSE.
         our_rev_mode%plain=.TRUE.
         our_rev_mode%tape=.FALSE.
         our_rev_mode%adjoint=.FALSE.

         do jj=1,n
          bb(jj)%v = 0.0
         enddo
         bb(ii)%v = eps
         call stream_vel (u, bb, fc2)
         fdfc = (fc2%v - fc%v)/eps
         if ( fdfc .NE. 0. ) then
              accuracyAD = 1.d0 - bb(ii)%d/fdfc
         else
              accuracyAD = 0.
         end if

         
         
!        enddo



    
         



! loop over each directional derivative
!        do ii = 1, n
! write output
           print *, 'ph: ii, adfc, fdfc, acc = ', &
                     ii, bb(ii)%d, fdfc, accuracyAD
!                    ii, adbb(ii), fdfc, accuracyAD
        end do

	end program

