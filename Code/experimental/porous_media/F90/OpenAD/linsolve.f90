module linsolve
use grid
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
subroutine sparse_solve(annz, arow_index, arow_compressed, &
                        acol_index, avalues, b, x, solver_inner, solver_outer, verbose)
  implicit none
  logical :: verbose
  integer :: solver_inner, solver_outer
  
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

!
! A method to test OpenAD without solver
!
subroutine sparse_dummy_method(n, annz, alen, arow_index, arow_compressed, &
                               acol_index, avalues, b, x, solver_inner, solver_outer, verbose)
  integer :: i
  logical :: verbose
  integer :: solver_inner, solver_outer
  double precision :: summe
  integer :: n, annz, alen
  integer, dimension(alen) :: arow_index
  integer, dimension(alen) :: acol_index
  double precision, dimension(alen) :: avalues
  integer, dimension(n + 1) :: arow_compressed
  
  double precision, dimension(n) :: b
  double precision, dimension(n) :: x
! 
!   x(:) = 0.0d0
!   summe = 0.0d0
!   do i = 1, annz
!       summe = summe + avalues(i)
!   end do
! 
!   x(:) = b(:)*summe
  call spmat_multiply_vector2(annz, arow_index, arow_compressed, &
                              acol_index, avalues, b, x, "PRE")
end subroutine sparse_dummy_method


! a wrapper for mgmres
!
! subroutine sparse_mgmres_method(annz, arow_index, arow_compressed, &
!                                  acol_index, avalues, b, x)
!   use mgmres
!   implicit none
!   integer :: itr_max, mr
!   double precision :: tol_abs, tol_rel
!   integer :: annz
!   integer, dimension(7 * N_) :: arow_index
!   integer, dimension(7 * N_) :: acol_index
!   double precision, dimension(7 * N_) :: avalues
!   integer, dimension(N_ + 1) :: arow_compressed
!   double precision, dimension(N_) :: b
!   double precision, dimension(N_) :: x
! 
!   tol_abs = 1.0d-8
!   tol_rel = 1.0d-8
!   itr_max = int(sqrt(N_ * 1.0)) + 1
!   mr = min(N_, 100)
! 
!   x = 0.0d0
! 
!   call mgmres_st (N_, annz, arow_index, acol_index, avalues, x, b, &
!                   itr_max, mr, tol_abs, tol_rel )
! end subroutine sparse_mgmres_method
end module linsolve