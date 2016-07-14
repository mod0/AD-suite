!
! A wrapper for pmgmres_ilu_cr
!
SUBROUTINE SPARSE_PMGMRES_METHOD_B(n, annz, alen, arow_index,   &
                                  &   arow_compressed, acol_index, avalues,&
                                  &   avaluesb, b, bb, x, xb,              &
                                  &   solver_inner, solver_outer, verbose)
                                  
  USE MGMRES
  USE MATHUTIL_B
  USE UTILS
  IMPLICIT NONE
  INTEGER :: i, curr_nnz, curr_row
  INTEGER :: itr_max, mr
  LOGICAL :: verbose
  INTEGER :: solver_inner, solver_outer
  DOUBLE PRECISION :: tol_abs, tol_rel, nrm
  INTEGER :: n, annz, alen
  INTEGER, DIMENSION(alen) :: arow_index
  INTEGER, DIMENSION(alen) :: acol_index
  DOUBLE PRECISION, DIMENSION(alen) :: avalues
  DOUBLE PRECISION, DIMENSION(alen) :: avaluesb
  INTEGER, DIMENSION(n + 1) :: arow_compressed
  DOUBLE PRECISION, DIMENSION(n) :: b
  DOUBLE PRECISION, DIMENSION(n) :: bb
  DOUBLE PRECISION, DIMENSION(n) :: x
  DOUBLE PRECISION, DIMENSION(n) :: xb

  DOUBLE PRECISION, DIMENSION(n) :: incrbb
  
  INTEGER, DIMENSION(alen) :: arow_index_transposed
  INTEGER, DIMENSION(alen) :: acol_index_transposed
  DOUBLE PRECISION, DIMENSION(alen) :: avalues_transposed
  INTEGER, DIMENSION(n + 1) :: arow_compressed_transposed
  
  INTEGER, DIMENSION(2,annz)  :: column_index_key_value 
  
  tol_abs = 1.0d-8
  tol_rel = 1.0d-8
  itr_max = solver_outer
  mr = solver_inner
 
  ! Transpose the matrix A
  ! First get the columns and its indices
  DO i = 1,annz
    column_index_key_value(1, i) = acol_index(i)        ! key
    column_index_key_value(2, i) = i                    ! value
  ENDDO
  
  ! Sort the indices based on the column indices
  CALL SORT(column_index_key_value, annz)
  
  ! Now use the transposed entries to create A^T compressed row
  DO i = 1,annz
    arow_index_transposed(i) = acol_index(column_index_key_value(2, i))
    acol_index_transposed(i) = arow_index(column_index_key_value(2, i))
    avalues_transposed(i) = avalues(column_index_key_value(2, i))
  ENDDO
       
  ! Now create row compressed for A^T
  i = 1
  curr_nnz = 0
  curr_row = 1
  arow_compressed_transposed(:) = 0
  arow_compressed_transposed(1) = 1

  DO WHILE (i <= annz .and. curr_row <= n) 
    IF(arow_index_transposed(i) == curr_row) THEN
      curr_nnz = curr_nnz + 1
      i = i + 1
    ELSE
      curr_row = curr_row + 1
      arow_compressed_transposed(curr_row) = &
                    arow_compressed_transposed(curr_row - 1) + curr_nnz
      curr_nnz = 0  
    ENDIF
  ENDDO
  
  IF(curr_row == n) THEN 
    curr_row = curr_row + 1
    arow_compressed_transposed(curr_row) = &
          arow_compressed_transposed(curr_row - 1) + curr_nnz
  ELSEIF(curr_row < n .or. i <= annz) THEN
    stop "The matrix is singular"
  ENDIF
  
  call dnrm2(xb,n, nrm)
  incrbb = 0.0d0
  if(nrm /= 0.0d0) then
    CALL PMGMRES_ILU_CR (n, annz, arow_compressed_transposed, &
                        acol_index_transposed, avalues_transposed, &
                        incrbb, xb, itr_max, mr, tol_abs, tol_rel, verbose)
    bb = bb + incrbb
  endif
    
  CALL DNRM2(b, n, nrm)
  x = 0.0d0
  IF(nrm /= 0.0d0) THEN  
    CALL PMGMRES_ILU_CR (n, annz, arow_compressed, acol_index, avalues, &
                x, b, itr_max, mr, tol_abs, tol_rel, verbose)
  ENDIF
  
  DO i = 1,annz
    avaluesb(i) = avaluesb(i) - x(acol_index(i)) * incrbb(arow_index(i))
  ENDDO
  
  xb = 0.0d0

  RETURN
END SUBROUTINE SPARSE_PMGMRES_METHOD_B

