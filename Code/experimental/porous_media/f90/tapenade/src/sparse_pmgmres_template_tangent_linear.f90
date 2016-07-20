!
! A wrapper for pmgmres_ilu_cr
!
SUBROUTINE SPARSE_PMGMRES_METHOD(n, annz, alen, arow_index, arow_compressed, &
     acol_index, avalues, b, x, solver_inner, &
     solver_outer, verbose)
  USE MGMRES
  USE MATHUTIL_D
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


!
! A wrapper for pmgmres_ilu_cr - tangent linear version
!
SUBROUTINE SPARSE_PMGMRES_METHOD_D(n, annz, alen, arow_index, arow_compressed, &
                                acol_index, avalues, avaluesd, b, bd, x, xd, solver_inner, &
                                solver_outer, verbose)
    USE MGMRES
    USE MATHUTIL_D
    IMPLICIT NONE
    INTEGER :: i, itr_max, mr
    INTEGER :: solver_inner, solver_outer
    LOGICAL :: verbose
    DOUBLE PRECISION :: tol_abs, tol_rel, nrm
    INTEGER :: n, annz, alen
    INTEGER, DIMENSION(alen) :: arow_index
    INTEGER, DIMENSION(alen) :: acol_index
    DOUBLE PRECISION, DIMENSION(alen) :: avalues
    DOUBLE PRECISION, DIMENSION(alen) :: avaluesd
    INTEGER, DIMENSION(n + 1) :: arow_compressed
    DOUBLE PRECISION, DIMENSION(n) :: b
    DOUBLE PRECISION, DIMENSION(n) :: bd
    DOUBLE PRECISION, DIMENSION(n) :: x
    DOUBLE PRECISION, DIMENSION(n) :: xd
    DOUBLE PRECISION, DIMENSION(n) :: RHSd

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
    
    RHSd = bd
    
    DO i = 1, annz
      RHSd(arow_index(i)) = RHSd(arow_index(i)) - avaluesd(i) * x(acol_index(i))
    END DO


    CALL DNRM2(RHSd, n, nrm)
    xd = 0.0d0  
    IF(nrm /= 0.0d0) THEN  
      CALL PMGMRES_ILU_CR (n, annz, arow_compressed, acol_index, avalues, &
                          xd, RHSd, itr_max, mr, tol_abs, tol_rel, verbose)
    ENDIF
    RETURN
END SUBROUTINE SPARSE_PMGMRES_METHOD_D

