module linsolve
use matrix
use mathutil

implicit none

! export solve
public :: solve

! interface to linear solvers
interface solve
    module procedure sparse_solve2
end interface solve

contains

!
! Calls a specific solver - here the jacobi method
!
subroutine sparse_solve2(arows, acols, annz, alen, arow_index, acol_index, avalues, b, x)
    implicit none
    integer :: arows, acols, annz, alen
    integer, dimension(alen) :: arow_index
    integer, dimension(alen) :: acol_index
    double precision, dimension(alen) :: avalues
    double precision, dimension(arows) :: b
    double precision, dimension(arows) :: x
    if(1 .ne. 0) then
        x = 0.0d0
    end if
    call sparse_gauss_seidel_method2(arows, acols, annz, alen, arow_index, acol_index, avalues, b, x)
end subroutine sparse_solve2

!
! The subroutine assumes the diagonals are non-zero
!
subroutine sparse_jacobi_method2(arows, acols, annz, alen, arow_index, acol_index, avalues, b, x)
    implicit none
    integer :: i, j
    double precision :: nrm
    integer :: arows, acols, annz, alen
    integer, dimension(alen) :: arow_index
    integer, dimension(alen) :: acol_index
    double precision, dimension(alen) :: avalues
    double precision, dimension(arows) :: b
    double precision, dimension(arows) :: x
    double precision, dimension(arows) :: x_old
    double precision, dimension(arows) :: residual
    double precision, dimension(arows) :: main_diag


    do i = 1, 100000      ! Max iterations

        x_old = x
        ! start x with right hand side
        x = b

!$omp parallel
        do j = 1, annz
            ! Not the main diagonal
            if (arow_index(j) /= acol_index(j)) then
                ! update the new iterate at the corresponding row
                x(arow_index(j)) = x(arow_index(j)) &
                                    - avalues(j) * x_old(acol_index(j))
            else
                ! place the entry in the corresponding diagonal.
                main_diag(arow_index(j)) = avalues(j)
            end if
        end do
!$omp end parallel

        ! divide by the diagonal entries
        x = x/main_diag

        ! Get the residual between the old and new values
        residual = (x - x_old)

        ! get the norm
        call dnrm2(residual, arows, nrm)

        ! Exit when converged.
        if (nrm < 1.0d-8) then
            !print *, "Converged in norm"
            exit
        end if
    end do

    if (i > 5000) then
        !print *, "Norm: ", nrm
    end if
end subroutine sparse_jacobi_method2

!
! The subroutine assumes the diagonals are non-zero
!
subroutine sparse_gauss_seidel_method2(arows, acols, annz, alen, arow_index, &
                                       acol_index, avalues, b, x)
    implicit none
    integer :: i, j, k
    double precision :: nrm
    integer :: arows, acols, annz, alen
    integer, dimension(alen) :: arow_index
    integer, dimension(alen) :: acol_index
    double precision, dimension(alen) :: avalues
    double precision, dimension(arows) :: b
    double precision, dimension(arows) :: x
    double precision, dimension(arows) :: x_old
    double precision, dimension(arows) :: residual
    double precision, dimension(arows) :: main_diag

    do k = 1, 100000      ! Max iterations
        ! start with old x
        x_old = x

        ! start x with right hand side
        x = b

        do j = 1, annz
            ! Not the main diagonal
            ! This is same as Jacobi (upper triangle part of A)
            if (arow_index(j) < acol_index(j)) then
                ! update the new iterate at the corresponding row
                x(arow_index(j)) = x(arow_index(j)) &
                                    - avalues(j) * x_old(acol_index(j))
            else
                ! place the entry in the corresponding diagonal.
                main_diag(arow_index(j)) = avalues(j)
            end if
        end do

        do i = 1,arows                 ! An artificial bound for the row
            do j = 1,annz
                if(arow_index(j) == i .and. acol_index(j) < i) then ! 1st part not required.
                    ! Use the new solution upto row < i
                    x(arow_index(j)) = x(arow_index(j)) &
                                    - avalues(j) * x(acol_index(j))
                else if(acol_index(j) >= i) then
                    ! exit when the column is more than or equal i
                    exit
                end if
            end do

            ! Divide the ith component of soln by the ith component of maindiag
            x(i) = x(i)/main_diag(i)
        end do


        ! Get the residual between the old and new values
        residual = (x - x_old)

        ! get the norm
        call dnrm2(residual, arows, nrm)

        ! Exit when converged.
        if (nrm < 1.0d-8) then
!            print *, "Iterations: ", k
!            print *, "Converged in norm"
            exit
        end if
    end do

    if (k > 5000) then
        !print *, "Norm: ", nrm
    end if
end subroutine sparse_gauss_seidel_method2

! !
! ! A wrapper for mgmres
! !
! subroutine sparse_mgmres_method2(arows, acols, annz, alen, arow_index, &
!                                  acol_index, avalues, b, x)
!     use mgmres
!     implicit none
!     integer :: itr_max, mr
!     double precision :: tol_abs, tol_rel
!     integer :: arows, acols, annz, alen
!     integer, dimension(alen) :: arow_index
!     integer, dimension(alen) :: acol_index
!     double precision, dimension(alen) :: avalues
!     double precision, dimension(arows) :: b
!     double precision, dimension(arows) :: x
!
!     tol_abs = 1.0d-6
!     tol_rel = 1.0d-6
!     itr_max = 5
!     mr = (arows)
!
!     call mgmres_st ( arows, annz, arow_index, acol_index, avalues, x, b, &
!                     itr_max, mr, tol_abs, tol_rel )
! end subroutine sparse_mgmres_method2

!
! A method to test OpenAD without solver
!
subroutine sparse_dummy_method2(arows, acols, annz, alen, arow_index, &
                                acol_index, avalues, b, x)
        integer :: i
        double precision :: sum
        integer :: arows, acols, annz, alen
        integer, dimension(alen) :: arow_index
        integer, dimension(alen) :: acol_index
        double precision, dimension(alen) :: avalues
        double precision, dimension(arows) :: b
        double precision, dimension(arows) :: x

        sum = 0.0d0
        do i = 1, annz
            sum = sum + avalues(i)
        end do

        x = b/sum
end subroutine sparse_dummy_method2

end module linsolve
