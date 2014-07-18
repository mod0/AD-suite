program driver
  use OAD_active
  use OAD_rev
  implicit none 
  external head

  TYPE (active), dimension(79) :: bb 
  TYPE (active), dimension(79) :: h
  TYPE (active), dimension(80) :: u
  TYPE (active) :: fc
  integer :: ii, jj, n
  real(8) :: fc_0, accuracyAD, fdfc
  real(8), parameter :: ep = 1.d-7

  external stream_vel
  n = 79
!  initialize
  do ii = 1, n
    bb(ii)%v    = 0.
    bb(ii)%d = 0.0D0
    u(ii)%d =  0.0D0
  end do
  u(n+1)%d =  0.0D0
  fc%v = 0.	
  fc%d = 1.0D0

           our_rev_mode%arg_store=.FALSE.
           our_rev_mode%arg_restore=.FALSE.
           our_rev_mode%arg_look=.FALSE.
           our_rev_mode%res_store=.FALSE.
           our_rev_mode%res_restore=.FALSE.
           our_rev_mode%plain=.TRUE.
           our_rev_mode%tape=.FALSE.
           our_rev_mode%adjoint=.FALSE.
      call stream_vel_timedep( h, u, bb, fc )
 print *, '    position       velocity'
 print *, '----------------------------------------------------------------'
        do ii=1,n
         print *,ii,u(ii+1)%v, h(ii)%v
        enddo

! call adjoint model
     our_rev_mode%arg_store=.FALSE.
     our_rev_mode%arg_restore=.FALSE.
     !our_rev_mode%arg_look=.FALSE.
     our_rev_mode%res_store=.FALSE.
     our_rev_mode%res_restore=.FALSE.
     our_rev_mode%plain=.FALSE.
     our_rev_mode%tape=.TRUE.
     our_rev_mode%adjoint=.TRUE.
        call stream_vel_timedep( h, u, bb, fc )



         

print *, '                                    position       ', &
         '            dfc/dM [ADM]           dfc/dM [f.d.]            relative accuracy'
print *, '----------------------------------------------------------------', &
         '-----------------------------------------------------------------------------'


! loop over each directional derivative
        do ii = 1, n

           do jj=1,n
            bb(jj)%v=0.0
           enddo
           bb(ii)%v=ep
           our_rev_mode%arg_store=.FALSE.
           our_rev_mode%arg_restore=.FALSE.
           our_rev_mode%arg_look=.FALSE.
           our_rev_mode%res_store=.FALSE.
           our_rev_mode%res_restore=.FALSE.
           our_rev_mode%plain=.TRUE.
           our_rev_mode%tape=.FALSE.
           our_rev_mode%adjoint=.FALSE.

           call stream_vel_timedep (h,u,bb,fc)
           fdfc = (fc%v-fc_0)/ep      
           accuracyAD = 1-bb(ii)%d/fdfc

 
! write output
           print *, 'ph: ii, adfc, fdfc, acc = ', &
                     ii, bb(ii)%d, fdfc, accuracyAD
        end do

	end program

