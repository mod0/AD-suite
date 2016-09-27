program driver
  implicit none 
  external head

  real(8), dimension(79) :: bb 
  real(8), dimension(79) :: h
  real(8), dimension(80) :: u
  real(8) :: fc
  integer :: ii, jj, n
  real(8) :: fc_0, accuracyAD, fdfc
  real(8), parameter :: ep = 1.d-4

  external stream_vel
  n = 79
!  initialize
  do ii = 1, n
    bb(ii)    = 0.
  end do
  u(n+1) =  0.0D0
  fc = 0.	
      call stream_vel_timedep( h, u, bb, fc )
      fc_0 = fc
 print *, '    position       velocity     thickness'
 print *, '----------------------------------------------------------------'
        do ii=1,n
         print *,ii,u(ii+1), h(ii)
        enddo



         

print *, '                                    position       ', &
         '            dfc/dM [ADM]           dfc/dM [f.d.]            relative accuracy'
print *, '----------------------------------------------------------------', &
         '-----------------------------------------------------------------------------'


! loop over each directional derivative
        do ii = 1, n
!        do ii = 1, 0

           do jj=1,n
            bb(jj)=0.0
           enddo
           bb(ii)=ep

           call stream_vel_timedep (h,u,bb,fc)
           fdfc = (fc-fc_0)/ep      

 
! write output
          print *, '                          ',ii, fdfc
        end do

	end program

