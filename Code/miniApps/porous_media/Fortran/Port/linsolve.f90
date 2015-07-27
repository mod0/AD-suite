module linsolve
use matrix
use mathutil

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

    x = b

    do i = 1, 5000      ! Max iterations
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
            !print *, "Converged in norm"
            exit
        end if
    end do

    if (i > 5000) then
        print *, "Norm: ", nrm
    end if
end subroutine sparse_jacobi_method


!
! The subroutine assumes the diagonals are non-zero
!
subroutine sparse_gauss_seidel_method(A, b, x)
    implicit none
    integer :: i, j, k
    type(spmat) :: A
    double precision :: nrm
    double precision, dimension(A%rows) :: b, x, x_old, main_diag

    x = b

    do k = 1, 5000      ! Max iterations
        ! start with old x
        x_old = x

        ! start x with right hand side
        x = b

        do j = 1, A%nnz
            ! Not the main diagonal
            ! This is same as Jacobi (upper triangle part of A)
            if (A%row_index(j) < A%col_index(j)) then
                ! update the new iterate at the corresponding row
                x(A%row_index(j)) = x(A%row_index(j)) &
                                    - A%values(j) * x_old(A%col_index(j))
            else
                ! place the entry in the corresponding diagonal.
                main_diag(A%row_index(j)) = A%values(j)
            end if
        end do

        do i = 1,A%rows                 ! An artificial bound for the row
            do j = 1,A%nnz
                if(A%row_index(j) == i .and. A%col_index(j) < i) then ! 1st part not required.
                    ! Use the new solution upto row < i
                    x(A%row_index(j)) = x(A%row_index(j)) &
                                    - A%values(j) * x(A%col_index(j))
                else if(A%col_index(j) >= i) then
                    ! exit when the column is more than or equal i
                    exit
                end if
            end do

            ! Divide the ith component of soln by the ith component of maindiag
            x(i) = x(i)/main_diag(i)
        end do

        ! get the norm
        call dnrm2((x - x_old), A%rows, nrm)

        ! Exit when converged.
        if (nrm < 1.0d-8) then
            !print *, "Converged in norm"
            exit
        end if
    end do

    if (k > 5000) then
        print *, "Norm: ", nrm
    end if
end subroutine sparse_gauss_seidel_method

end module linsolve
