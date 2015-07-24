module linsolve
use matrix

implicit none

! export solve
public :: solve

! interface to linear solvers
interface solve
    module procedure sparse_solve
end interface solve

contains

!
! Calls a specific solver - here the jacobi method
!
subroutine sparse_solve(A, b, x)
    implicit none
    type(spmat) :: A
    double precision, dimension(:) :: b
    double precision, dimension(:) :: x

    call sparse_jacobi_method(A, b, x)
end subroutine sparse_solve

!
! The subroutine assumes the diagonals are non-zero
!
subroutine sparse_jacobi_method(A, b, x)
    implicit none
    integer :: i, j
    type(spmat) :: A
    double precision :: nrm
    double precision, dimension(A%rows) :: b, x, x_old, main_diag

    do i = 1, 10000       ! Max iterations
        ! start with old x
        x_old = x

        ! start x with right hand side
        x = b

        do j = 1, A%nnz
            ! Not the main diagonal
            if (A%row_index(j) /= A%col_index(j)) then
                ! update the new iterate at the corresponding row
                x(A%row_index(j)) = x(A%row_index(j)) &
                                    - A%values(j) * x_old(A%col_index(j))
            else
                ! place the entry in the corresponding diagonal.
                main_diag(A%row_index(j)) = A%values(j)
            end if
        end do

        ! divide by the diagonal entries
        x = x/main_diag

        ! get the norm
        call dnrm2((x - x_old), A%rows, nrm)

        ! Exit when converged.
        if (nrm < 1.0d-8) then
            print *, "Converged in norm"
            exit
        end if
    end do

    print *, "Iteration count: ", i
end subroutine sparse_jacobi_method

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
end module linsolve
