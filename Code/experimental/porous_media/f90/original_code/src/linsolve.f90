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
subroutine sparse_solve(annz, arow_index, arow_compressed, &
                        acol_index, avalues, b, x)
  implicit none
 
  integer :: annz
  integer, parameter :: matdim = N_
  integer, parameter :: maxlen = 7 * matdim
  integer, dimension(7 * N_) :: arow_index
  integer, dimension(7 * N_) :: acol_index
  double precision, dimension(7 * N_) :: avalues
  integer, dimension(N_ + 1) :: arow_compressed

  double precision, dimension(N_) :: b
  double precision, dimension(N_) :: x

  call sparse_dummy_method(matdim, annz, maxlen, arow_index, arow_compressed,&
                           acol_index, avalues, b, x, solver_inner, solver_outer, verbose)
end subroutine sparse_solve

subroutine sparse_dummy_method(n, annz, alen, arow_index, arow_compressed, &
                               acol_index, avalues, b, x, solver_inner, solver_outer, verbose)
  integer :: i
  logical :: verbose
  integer :: solver_inner, solver_outer
  double precision :: sum
  integer :: n, annz, alen
  integer, dimension(alen) :: arow_index
  integer, dimension(alen) :: acol_index
  double precision, dimension(alen) :: avalues
  integer, dimension(n + 1) :: arow_compressed
  
  double precision, dimension(n) :: b
  double precision, dimension(n) :: x

  x = 0.0d0
  sum = 0.0d0
  do i = 1, annz
      sum = sum + avalues(i)
  end do

  x = b/sum
end subroutine sparse_dummy_method
end module linsolve
