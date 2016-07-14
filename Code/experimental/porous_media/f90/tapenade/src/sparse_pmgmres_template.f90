!
! A wrapper for pmgmres_ilu_cr
!
SUBROUTINE SPARSE_PMGMRES_METHOD(n, annz, alen, arow_index, arow_compressed, &
                                  acol_index, avalues, b, x, solver_inner, &
                                  solver_outer, verbose)
    USE MGMRES
    USE MATHUTIL
    IMPLICIT NONE
    INTEGER :: itr_max, mr
    LOGICAL :: verbose
    INTEGER :: solver_inner, solver_outer
    DOUBLE PRECISION :: tol_abs, tol_rel, nrm
    INTEGER :: n,  annz, alen
    INTEGER, DIMENSION(alen) :: arow_index
    INTEGER, DIMENSION(alen) :: acol_index
    DOUBLE PRECISION, DIMENSION(alen) :: avalues
    INTEGER, DIMENSION(n + 1) :: arow_compressed
    DOUBLE PRECISION, DIMENSION(n) :: b
    DOUBLE PRECISION, DIMENSION(n) :: x

    tol_abs = 1.0d-8
    tol_rel = 1.0d-8
    itr_max = solver_outer
    mr = solver_inner
    
    CALL DNRM2(b, n, nrm)
    x = 0.0d0   
    IF(nrm /= 0.0d0) THEN
       CALL PMGMRES_ILU_CR (n, annz, arow_compressed, acol_index, avalues, &
                  x, b, itr_max, mr, tol_abs, tol_rel, verbose)
    ENDIF
    RETURN
END SUBROUTINE SPARSE_PMGMRES_METHOD

