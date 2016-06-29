module mathutil
implicit none
contains

subroutine scalar_max(scalarin1, scalarin2, scalarout)
  double precision :: scalarin1, scalarin2, scalarout

  if (scalarin2 <= scalarin1) then
    scalarout = scalarin1
  else
    scalarout = scalarin2
  end if
end subroutine scalar_max

!
! Subroutine computes the 2 norm of the vector
! adapted from the original function version
! of the corresponding blas routine from
! NETLIB
!
subroutine dnrm2(v, len_v, n)
    implicit none
    integer :: i, len_v
    double precision :: n
    double precision, dimension(len_v) :: v
    double precision :: scalein, scaleout

    n = 0.0d0
    scalein = 0.0d0
    scaleout = 0.0d0

    do  i = 1, len_v
        call scalar_max(scalein, abs(v(i)), scaleout)
        scalein = scaleout
    enddo

    if (scaleout .eq. 0.0d0) then
        n = 0.0d0
    else
        do i = 1, len_v
            n = n + (v(i)/scaleout)**2
        enddo

        n = scaleout*sqrt(n)
    endif
end subroutine dnrm2
end module mathutil