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

    integer :: annz
    integer :: n
    integer, dimension(7 * n) :: arow_index
    integer, dimension(7 * n) :: acol_index
    double precision, dimension(7 * n) :: avalues
    integer, dimension(n + 1) :: arow_compressed

    double precision, dimension(n) :: b
    double precision, dimension(n) :: x

    call sparse_dummy_method(n, annz, 7 * n, arow_index, arow_compressed,&
         acol_index, avalues, b, x, solver_inner, solver_outer, verbose)
  end subroutine sparse_solve
end module linsolve
