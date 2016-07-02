!
! A wrapper for mgmres
!
SUBROUTINE SPARSE_MGMRES_METHOD(n, annz, alen, arow_index, &
                                acol_index, avalues, b, x, &
                                solver_inner, solver_outer, verbose)
    USE MGMRES
    IMPLICIT NONE
    INTEGER :: itr_max, mr
    LOGICAL :: verbose
    INTEGER :: solver_inner, solver_outer
    DOUBLE PRECISION :: tol_abs, tol_rel
    INTEGER :: n, annz, alen
    INTEGER, DIMENSION(alen) :: arow_index
    INTEGER, DIMENSION(alen) :: acol_index
    DOUBLE PRECISION, DIMENSION(alen) :: avalues
    DOUBLE PRECISION, DIMENSION(n) :: b
    DOUBLE PRECISION, DIMENSION(n) :: x

    tol_abs = 1.0d-8
    tol_rel = 1.0d-8
    itr_max = solver_outer
    mr = solver_inner

    x = 0.0d0

    CALL MGMRES_ST ( n, annz, arow_index, acol_index, avalues, x, b, &
                    itr_max, mr, tol_abs, tol_rel, verbose )
END SUBROUTINE SPARSE_MGMRES_METHOD