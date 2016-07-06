!
! A wrapper for mgmres
!
SUBROUTINE SPARSE_MGMRES_METHOD_D(n annz, alen, arow_index, &
                                acol_index, avalues, avaluesd, b, bd, x, xd, &
                                solver_inner, solver_outer, verbose)
    USE MGMRES
    IMPLICIT NONE
    INTEGER :: i, itr_max, mr
    INTEGER :: solver_inner, solver_outer
    LOGICAL :: verbose
    DOUBLE PRECISION :: tol_abs, tol_rel
    INTEGER :: n, annz, alen
    INTEGER, DIMENSION(alen) :: arow_index
    INTEGER, DIMENSION(alen) :: acol_index
    DOUBLE PRECISION, DIMENSION(alen) :: avalues
    DOUBLE PRECISION, DIMENSION(alen) :: avaluesd
    DOUBLE PRECISION, DIMENSION(n) :: b
    DOUBLE PRECISION, DIMENSION(n) :: bd
    DOUBLE PRECISION, DIMENSION(n) :: x
    DOUBLE PRECISION, DIMENSION(n) :: xd
    DOUBLE PRECISION, DIMENSION(n) :: RHSd

    tol_abs = 1.0d-8
    tol_rel = 1.0d-8
    itr_max = solver_outer
    mr = solver_inner

    x = 0.0d0
    xd = 0.0d0
    
    CALL MGMRES_ST ( n, annz, arow_index, acol_index, avalues, x, b, &
                    itr_max, mr, tol_abs, tol_rel, verbose)

    RHSd = bd
    
    DO i = 1, annz
      RHSd(arow_index(i)) = RHSd(arow_index(i)) - avaluesd(i) * x(acol_index(i))
    END DO

    CALL MGMRES_ST ( n, annz, arow_index, acol_index, avalues, xd, RHSd, &
                    itr_max, mr, tol_abs, tol_rel, verbose)
END SUBROUTINE SPARSE_MGMRES_METHOD_D