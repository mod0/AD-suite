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
    if(1 .ne. 0) then
        x = 0.0d0
    end if
    call sparse_mgmres_method(A, b, x)
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


    do i = 1, 100000      ! Max iterations

        x_old = x
        ! start x with right hand side
        x = b

!$omp parallel
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
!$omp end parallel

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

    do k = 1, 100000      ! Max iterations
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
!            print *, "Iterations: ", k
!            print *, "Converged in norm"
            exit
        end if
    end do

    if (k > 5000) then
        print *, "Norm: ", nrm
    end if
end subroutine sparse_gauss_seidel_method

!
! A wrapper for mgmres
!
subroutine sparse_mgmres_method(A, b, x)
    use mgmres
    implicit none
    type(spmat) :: A
    integer :: itr_max, mr
    double precision :: tol_abs, tol_rel
    double precision, dimension(:) :: b, x


    tol_abs = 1.0d-6
    tol_rel = 1.0d-6
    itr_max = 5
    mr = (A%rows)

    call mgmres_st ( A%rows, A%nnz, A%row_index, A%col_index, A%values, x, b, &
                    itr_max, mr, tol_abs, tol_rel )

end subroutine sparse_mgmres_method


end module linsolve
