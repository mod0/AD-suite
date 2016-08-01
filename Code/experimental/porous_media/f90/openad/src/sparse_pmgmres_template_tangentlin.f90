SUBROUTINE SPARSE_DUMMY_METHOD(N, ANNZ, ALEN, AROW_INDEX, AROW_COMPRESSED, ACO&
     &L_INDEX, AVALUES, B, X, SOLVER_INNER, SOLVER_OUTER, VERBOSE)

  use w2f__types
  use mgmres
  use mathutil
  IMPLICIT NONE
  !
  !       **** Parameters and Result ****
  !
  INTEGER(w2f__i4) N
  INTEGER(w2f__i4) ANNZ
  INTEGER(w2f__i4) ALEN
  INTEGER(w2f__i4) AROW_INDEX(1 : ALEN)
  INTEGER(w2f__i4) AROW_COMPRESSED(1 : INT((N + 1)))
  INTEGER(w2f__i4) ACOL_INDEX(1 : ALEN)
  type(active) :: AVALUES(1:ALEN)
  type(active) :: B(1:N)
  type(active) :: X(1:N)
  DOUBLE PRECISION :: AV(1:ALEN)
  DOUBLE PRECISION :: AD(1:ALEN)
  DOUBLE PRECISION :: BV(1:N)
  DOUBLE PRECISION :: BD(1:N)
  DOUBLE PRECISION :: XV(1:N)
  DOUBLE PRECISION :: XD(1:N)
  INTEGER(w2f__i4) SOLVER_INNER
  INTEGER(w2f__i4) SOLVER_OUTER
  LOGICAL(w2f__i4) VERBOSE

  INTEGER :: I, ITR_MAX, MR
  DOUBLE PRECISION :: TOL_ABS, TOL_REL, NRM
  DOUBLE PRECISION, DIMENSION(N) :: RHSD

  tol_abs = 1.0d-8
  tol_rel = 1.0d-8
  itr_max = solver_outer
  mr = solver_inner

  ! Copy active into passive
  do i = 1, N
     BV(i) = B(i)%v
     BD(i) = B(i)%d
  end do

  do i = 1, ANNZ
     AV(i) = AVALUES(i)%v
     AD(i) = AVALUES(i)%d
  end do

  CALL DNRM2(BV, n, nrm)

  XV = 0.0d0
  IF(nrm /= 0.0d0) THEN  
     CALL PMGMRES_ILU_CR (n, annz, arow_compressed, acol_index, AV, &
          XV, BV, itr_max, mr, tol_abs, tol_rel, verbose)
  ENDIF

  RHSd = BD

  DO i = 1, annz
     RHSd(arow_index(i)) = RHSd(arow_index(i)) - AD(i) * XV(acol_index(i))
  END DO


  CALL DNRM2(RHSd, n, nrm)
  XD = 0.0d0  
  IF(nrm /= 0.0d0) THEN  
     CALL PMGMRES_ILU_CR (n, annz, arow_compressed, acol_index, AV, &
          XD, RHSd, itr_max, mr, tol_abs, tol_rel, verbose)
  ENDIF

  ! Copy passive into active
  do i = 1, N
     X(i)%v = XV(i)
     X(i)%d = XD(i)
  end do

  RETURN
END SUBROUTINE SPARSE_DUMMY_METHOD

