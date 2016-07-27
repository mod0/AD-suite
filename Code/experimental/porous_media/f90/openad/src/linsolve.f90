module linsolve
use parameters
use mathutil
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
subroutine sparse_solve(n, annz, arow_index, arow_compressed, &
                        acol_index, avalues, b, x)
  implicit none
 
  integer :: n, annz
  integer, dimension(7 * n) :: arow_index
  integer, dimension(7 * n) :: acol_index
  double precision, dimension(7 * n) :: avalues
  integer, dimension(n + 1) :: arow_compressed

  double precision, dimension(n) :: b
  double precision, dimension(n) :: x

  call sparse_dummy_method(n, annz, 7*n, arow_index, arow_compressed,&
                           acol_index, avalues, b, x, solver_inner, solver_outer, verbose)
end subroutine sparse_solve

subroutine sparse_dummy_method(n, annz, alen, arow_index, arow_compressed, &
                               acol_index, avalues, b, x, solver_inner, solver_outer, verbose)
  integer :: i
  logical :: verbose
  integer :: solver_inner, solver_outer
  integer :: n, annz, alen
  integer, dimension(alen) :: arow_index
  integer, dimension(alen) :: acol_index
  double precision, dimension(alen) :: avalues
  integer, dimension(n + 1) :: arow_compressed
  
  double precision, dimension(n) :: b
  double precision, dimension(n) :: x

  call spmat_multiply_vector(n, annz, arow_index, arow_compressed, &
                              acol_index, avalues, b, x, "PRE")
end subroutine sparse_dummy_method
end module linsolve
