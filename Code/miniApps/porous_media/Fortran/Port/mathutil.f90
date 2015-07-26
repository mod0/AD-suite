module mathutil
implicit none
contains
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
    double precision :: scale

    n = 0.0d0
    scale = 0.0d0

    do  i = 1, len_v
        scale = max(scale, abs(v(i)))
    enddo

    if (scale .eq. 0.0d0) then
        n = 0.0d0
    else
        do i = 1, len_v
            n = n + (v(i)/scale)**2
        enddo

        n = scale*sqrt(n)
    endif
end subroutine dnrm2
end module mathutil
