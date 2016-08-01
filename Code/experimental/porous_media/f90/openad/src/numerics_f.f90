
MODULE oad_intrinsics
use OAD_active
use w2f__types
IMPLICIT NONE
SAVE
!
!     **** Statements ****
!
END MODULE

MODULE parameters
use OAD_active
use w2f__types
IMPLICIT NONE
SAVE
!
!     **** Global Variables & Derived Type Definitions ****
!
REAL(w2f__8) HX
REAL(w2f__8) HY
REAL(w2f__8) HZ
REAL(w2f__8) IR
REAL(w2f__8) PERM(:, :, :, :)
ALLOCATABLE PERM
REAL(w2f__8) POR(:)
ALLOCATABLE POR
INTEGER(w2f__i4) SCENARIO_ID
INTEGER(w2f__i4) SOLVER_INNER
INTEGER(w2f__i4) SOLVER_OUTER
REAL(w2f__8) SOR
REAL(w2f__8) SWC
LOGICAL(w2f__i4) VERBOSE
REAL(w2f__8) VO
REAL(w2f__8) VOL
REAL(w2f__8) VW
!
!     **** Statements ****
!
END MODULE

MODULE mathutil
use OAD_active
use w2f__types
use oad_intrinsics
IMPLICIT NONE
SAVE
!
!     **** Statements ****
!
CONTAINS

  SUBROUTINE SCALAR_MAX(SCALARIN1, SCALARIN2, SCALAROUT)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  REAL(w2f__8) SCALARIN1
  REAL(w2f__8) SCALARIN2
  REAL(w2f__8) SCALAROUT
!
!       **** Statements ****
!
  IF(SCALARIN1 .GE. SCALARIN2) THEN
    SCALAROUT = SCALARIN1
  ELSE
    SCALAROUT = SCALARIN2
  ENDIF
  END SUBROUTINE

  SUBROUTINE DNRM2(V, LEN_V, N)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) LEN_V
  REAL(w2f__8) N
  REAL(w2f__8) V(1 : LEN_V)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  REAL(w2f__8) OAD_CTMP0
  REAL(w2f__8) SCALEIN
  REAL(w2f__8) SCALEOUT
!
!       **** Statements ****
!
  N = 0.0D00
  SCALEIN = 0.0D00
  SCALEOUT = 0.0D00
  DO I = 1, LEN_V, 1
    OAD_CTMP0 = ABS(V(I))
    CALL SCALAR_MAX(SCALEIN,OAD_CTMP0,SCALEOUT)
    SCALEIN = SCALEOUT
  END DO
  IF(SCALEOUT .eq. 0.0D00) THEN
    N = 0.0D00
  ELSE
    DO I = 1, LEN_V, 1
      N = (N +((V(I) / SCALEOUT) ** 2))
    END DO
    N = (SCALEOUT * SQRT(N))
  ENDIF
  END SUBROUTINE
END

MODULE matrix
use OAD_active
use w2f__types
use oad_intrinsics
use parameters
IMPLICIT NONE
SAVE
!
!     **** Top Level Pragmas ****
!
interface ADD_X
  module procedure ADDX_DIAGONAL
  module procedure ADDX_ELEM

end interface
      
interface MYMAX
  module procedure MYMAX_1_1_DOUBLE
  module procedure MYMAX_1_0_DOUBLE
  module procedure MYMAX_0_0_DOUBLE

end interface
      
interface MYMIN
  module procedure MYMIN_1_1_DOUBLE
  module procedure MYMIN_1_0_DOUBLE
  module procedure MYMIN_0_0_DOUBLE

end interface
      
interface MYRESHAPE
  module procedure MYRESHAPE_4_1
  module procedure MYRESHAPE_1_4
  module procedure MYRESHAPE_3_1
  module procedure MYRESHAPE_1_3
  module procedure MYRESHAPE_2_1
  module procedure MYRESHAPE_1_2

end interface
      
interface SPMAT_MULTIPLY
  module procedure SCALAR_MULTIPLY_SPMAT
  module procedure SPMAT_MULTIPLY_VECTOR
  module procedure SPMAT_MULTIPLY_DIAGONAL

end interface

!
!     **** Statements ****
!
CONTAINS

  SUBROUTINE ADDX_ELEM(N, INNZ, IROW_INDEX, IROW_COMPRESSED, ICOL_INDEX, IVALUES&
     &, X, ROW, COL)

  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) N
  INTEGER(w2f__i4) INNZ
  INTEGER(w2f__i4) IROW_INDEX(1 : INT((N * 7)))
  INTEGER(w2f__i4) IROW_COMPRESSED(1 : INT((N + 1)))
  INTEGER(w2f__i4) ICOL_INDEX(1 : INT((N * 7)))
  REAL(w2f__8) IVALUES(1 : INT((N * 7)))
  REAL(w2f__8) X
  INTEGER(w2f__i4) ROW
  INTEGER(w2f__i4) COL
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
!
!       **** Statements ****
!
 2      CONTINUE
  GO TO 3
 3      CONTINUE
  GO TO 4
 4      CONTINUE
  I = 1
  GO TO 13
 5      CONTINUE
  I = I + 1
 13     CONTINUE
  IF(I .LE. INNZ) THEN
    GO TO 6
  ELSE
    GO TO 11
  ENDIF
 6      CONTINUE
  IF((IROW_INDEX(I) .eq. ROW) .AND.(ICOL_INDEX(I) .eq. COL)) THEN
    GO TO 10
  ELSE
    GO TO 7
  ENDIF
 7      CONTINUE
  GO TO 8
 8      CONTINUE
  GO TO 9
 9      CONTINUE
  GO TO 5
 10     CONTINUE
  IVALUES(INT(I)) = (IVALUES(I) + X)
  GO TO 11
 11     CONTINUE
  GO TO 12
 12     CONTINUE
  GO TO 1
 1      CONTINUE
  END SUBROUTINE

  SUBROUTINE ADDX_DIAGONAL(N, INNZ, IROW_INDEX, IROW_COMPRESSED, ICOL_INDEX, IVA&
     &LUES, X, DIAG)

  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) N
  INTEGER(w2f__i4) INNZ
  INTEGER(w2f__i4) IROW_INDEX(1 : INT((N * 7)))
  INTEGER(w2f__i4) IROW_COMPRESSED(1 : INT((N + 1)))
  INTEGER(w2f__i4) ICOL_INDEX(1 : INT((N * 7)))
  type(active) :: IVALUES(1:INT((N*7)))
  REAL(w2f__8) X
  INTEGER(w2f__i4) DIAG
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  type(active) :: OpenAD_prp_0
!
!       **** Statements ****
!
  DO I = 1, INNZ, 1
    IF(DIAG .eq.(ICOL_INDEX(I) - IROW_INDEX(I))) THEN
      IVALUES(INT(I))%v = (IVALUES(I)%v+X)
      CALL setderiv(OpenAD_prp_0,IVALUES(I))
      CALL setderiv(IVALUES(I),OpenAD_prp_0)
    ENDIF
  END DO
  END SUBROUTINE

  SUBROUTINE GETDIAG(ROW, COL, DIAG)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) ROW
  INTEGER(w2f__i4) COL
  INTEGER(w2f__i4) DIAG
!
!       **** Statements ****
!
  DIAG = (COL - ROW)
  END SUBROUTINE

  SUBROUTINE FIRSTELM(DIAG, ROW, COL)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) DIAG
  INTEGER(w2f__i4) ROW
  INTEGER(w2f__i4) COL
!
!       **** Statements ****
!
  IF(DIAG .LT. 0) THEN
    ROW = (1 - DIAG)
    COL = 1
  ELSE
    IF(DIAG .eq. 0) THEN
      ROW = 1
      COL = 1
    ELSE
      ROW = 1
      COL = (DIAG + 1)
    ENDIF
  ENDIF
  END SUBROUTINE

  SUBROUTINE NOELEMS(DIAG, OROWS, OCOLS, N)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) DIAG
  INTEGER(w2f__i4) OROWS
  INTEGER(w2f__i4) OCOLS
  INTEGER(w2f__i4) N
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) COL
  INTEGER(w2f__i4) ROW
!
!       **** Statements ****
!
  CALL FIRSTELM(DIAG,ROW,COL)
  IF(DIAG .LT. 0) THEN
    N = (OROWS - ROW + 1)
  ELSE
    IF(DIAG .eq. 0) THEN
      N = MIN(OROWS, OCOLS)
    ELSE
      N = (OCOLS - COL + 1)
    ENDIF
  ENDIF
  END SUBROUTINE

  SUBROUTINE MYRESHAPE_2_1(AMATRIX, BMATRIX)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  REAL(w2f__8) AMATRIX(1 :, 1 :)
  REAL(w2f__8) BMATRIX(1 :)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) J
  INTEGER(w2f__i4) K
!
!       **** Statements ****
!
  K = 0
  DO J = 1, SIZE(AMATRIX, 2), 1
    DO I = 1, SIZE(AMATRIX, 1), 1
      K = (K + 1)
      BMATRIX(INT(K)) = AMATRIX(I, J)
    END DO
  END DO
  END SUBROUTINE

  SUBROUTINE MYRESHAPE_1_2(AMATRIX, BMATRIX)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  REAL(w2f__8) AMATRIX(1 :)
  REAL(w2f__8) BMATRIX(1 :, 1 :)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) J
  INTEGER(w2f__i4) K
!
!       **** Statements ****
!
  K = 0
  DO J = 1, SIZE(BMATRIX, 2), 1
    DO I = 1, SIZE(BMATRIX, 1), 1
      K = (K + 1)
      BMATRIX(INT(I), INT(J)) = AMATRIX(K)
    END DO
  END DO
  END SUBROUTINE

  SUBROUTINE MYRESHAPE_3_1(AMATRIX, BMATRIX)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  type(active) :: AMATRIX(1:,1:,1:)
  type(active) :: BMATRIX(1:)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) J
  INTEGER(w2f__i4) K
  INTEGER(w2f__i4) L
!
!       **** Statements ****
!
  L = 0
  DO K = 1, SIZE(AMATRIX, 3), 1
    DO J = 1, SIZE(AMATRIX, 2), 1
      DO I = 1, SIZE(AMATRIX, 1), 1
        L = (L + 1)
        BMATRIX(INT(L))%v = AMATRIX(I,J,K)%v
        CALL setderiv(BMATRIX(L),AMATRIX(I,J,K))
      END DO
    END DO
  END DO
  END SUBROUTINE

  SUBROUTINE MYRESHAPE_1_3(AMATRIX, BMATRIX)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  type(active) :: AMATRIX(1:)
  type(active) :: BMATRIX(1:,1:,1:)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) J
  INTEGER(w2f__i4) K
  INTEGER(w2f__i4) L
!
!       **** Statements ****
!
  L = 0
  DO K = 1, SIZE(BMATRIX, 3), 1
    DO J = 1, SIZE(BMATRIX, 2), 1
      DO I = 1, SIZE(BMATRIX, 1), 1
        L = (L + 1)
        BMATRIX(INT(I),INT(J),INT(K))%v = AMATRIX(L)%v
        CALL setderiv(BMATRIX(I,J,K),AMATRIX(L))
      END DO
    END DO
  END DO
  END SUBROUTINE

  SUBROUTINE MYRESHAPE_4_1(AMATRIX, BMATRIX)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  REAL(w2f__8) AMATRIX(1 :, 1 :, 1 :, 1 :)
  REAL(w2f__8) BMATRIX(1 :)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) J
  INTEGER(w2f__i4) K
  INTEGER(w2f__i4) L
  INTEGER(w2f__i4) M
!
!       **** Statements ****
!
  M = 0
  DO L = 1, SIZE(AMATRIX, 4), 1
    DO K = 1, SIZE(AMATRIX, 3), 1
      DO J = 1, SIZE(AMATRIX, 2), 1
        DO I = 1, SIZE(AMATRIX, 1), 1
          M = (M + 1)
          BMATRIX(INT(M)) = AMATRIX(I, J, K, L)
        END DO
      END DO
    END DO
  END DO
  END SUBROUTINE

  SUBROUTINE MYRESHAPE_1_4(AMATRIX, BMATRIX)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  type(active) :: AMATRIX(1:)
  type(active) :: BMATRIX(1:,1:,1:,1:)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) J
  INTEGER(w2f__i4) K
  INTEGER(w2f__i4) L
  INTEGER(w2f__i4) M
!
!       **** Statements ****
!
  M = 0
  DO L = 1, SIZE(BMATRIX, 4), 1
    DO K = 1, SIZE(BMATRIX, 3), 1
      DO J = 1, SIZE(BMATRIX, 2), 1
        DO I = 1, SIZE(BMATRIX, 1), 1
          M = (M + 1)
          BMATRIX(INT(I),INT(J),INT(K),INT(L))%v = AMATRIX(M)%v
          CALL setderiv(BMATRIX(I,J,K,L),AMATRIX(M))
        END DO
      END DO
    END DO
  END DO
  END SUBROUTINE

  SUBROUTINE SPMAT_MULTIPLY_DIAGONAL(N, ANNZ, AROW_INDEX, AROW_COMPRESSED, ACOL_&
     &INDEX, AVALUES, DMATRIX, RNNZ, RROW_INDEX , RROW_COMPRESSED, RCOL_INDEX, R&
     &VALUES, ORDER)

  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) N
  INTEGER(w2f__i4) ANNZ
  INTEGER(w2f__i4) AROW_INDEX(1 : INT((N * 7)))
  INTEGER(w2f__i4) AROW_COMPRESSED(1 : INT((N + 1)))
  INTEGER(w2f__i4) ACOL_INDEX(1 : INT((N * 7)))
  type(active) :: AVALUES(1:INT((N*7)))
  type(active) :: DMATRIX(1:N)
  INTEGER(w2f__i4) RNNZ
  INTEGER(w2f__i4) RROW_INDEX(1 : INT((N * 7)))
  INTEGER(w2f__i4) RROW_COMPRESSED(1 : INT((N + 1)))
  INTEGER(w2f__i4) RCOL_INDEX(1 : INT((N * 7)))
  type(active) :: RVALUES(1:INT((N*7)))
  CHARACTER(3) ORDER
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  REAL(w2f__8) OpenAD_lin_0
  REAL(w2f__8) OpenAD_lin_1
  REAL(w2f__8) OpenAD_lin_2
  REAL(w2f__8) OpenAD_lin_3
  type(active) :: OpenAD_prp_1
  type(active) :: OpenAD_prp_2
!
!       **** Statements ****
!
  RNNZ = ANNZ
  RROW_COMPRESSED(1 : INT(N + 1)) = AROW_COMPRESSED(1 : INT(N + 1 ))
  IF(ORDER .EQ. 'PRE') THEN
    DO I = 1, ANNZ, 1
      RROW_INDEX(I) = AROW_INDEX(I)
      RCOL_INDEX(I) = ACOL_INDEX(I)
      OpenAD_lin_0 = DMATRIX(ACOL_INDEX(I))%v
      OpenAD_lin_1 = AVALUES(I)%v
      RVALUES(INT(I))%v = (AVALUES(I)%v*DMATRIX(ACOL_INDEX(I))%v)
      CALL setderiv(OpenAD_prp_1,AVALUES(I))
      CALL sax(OpenAD_lin_0,OpenAD_prp_1,RVALUES(I))
      CALL saxpy(OpenAD_lin_1,DMATRIX(ACOL_INDEX(I)),RVALUES(I))
    END DO
  ELSE
    IF (ORDER.EQ.'POS') THEN
      DO I = 1,ANNZ,1
        RROW_INDEX(I) = AROW_INDEX(I)
        RCOL_INDEX(I) = ACOL_INDEX(I)
        OpenAD_lin_2 = DMATRIX(AROW_INDEX(I))%v
        OpenAD_lin_3 = AVALUES(I)%v
        RVALUES(INT(I))%v = (AVALUES(I)%v*DMATRIX(AROW_INDEX(I))%v)
        CALL setderiv(OpenAD_prp_2,AVALUES(I))
        CALL sax(OpenAD_lin_2,OpenAD_prp_2,RVALUES(I))
        CALL saxpy(OpenAD_lin_3,DMATRIX(AROW_INDEX(I)),RVALUES(I))
      END DO
    ENDIF
  ENDIF
  END SUBROUTINE

  SUBROUTINE SPMAT_MULTIPLY_VECTOR(N, ANNZ, AROW_INDEX, AROW_COMPRESSED, ACOL_IN&
     &DEX, AVALUES, BVECTOR, CVECTOR, ORDER)

  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) N
  INTEGER(w2f__i4) ANNZ
  INTEGER(w2f__i4) AROW_INDEX(1 : INT((N * 7)))
  INTEGER(w2f__i4) AROW_COMPRESSED(1 : INT((N + 1)))
  INTEGER(w2f__i4) ACOL_INDEX(1 : INT((N * 7)))
  type(active) :: AVALUES(1:INT((N*7)))
  type(active) :: BVECTOR(1:N)
  type(active) :: CVECTOR(1:N)
  CHARACTER(3) ORDER
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  REAL(w2f__8) OpenAD_lin_4
  REAL(w2f__8) OpenAD_lin_5
  REAL(w2f__8) OpenAD_lin_6
  REAL(w2f__8) OpenAD_lin_7
  type(active) :: OpenAD_prp_3
  type(active) :: OpenAD_prp_4
!
!       **** Statements ****
!
  CVECTOR(1:INT(N))%v = 0.0D00
  CALL zero_deriv(CVECTOR(1:INT(N)))
  IF (ORDER.EQ.'PRE') THEN
    DO I = 1,ANNZ,1
      OpenAD_lin_4 = BVECTOR(ACOL_INDEX(I))%v
      OpenAD_lin_5 = AVALUES(I)%v
      CVECTOR(AROW_INDEX(INT(I)))%v = (CVECTOR(AROW_INDEX(I))%v+AVALUES(I)%v*BVE&
     &CTOR(ACOL_INDEX(I))%v)

      CALL setderiv(OpenAD_prp_3,CVECTOR(AROW_INDEX(I)))
      CALL setderiv(CVECTOR(AROW_INDEX(I)),OpenAD_prp_3)
      CALL saxpy(OpenAD_lin_4,AVALUES(I),CVECTOR(AROW_INDEX(I)))
      CALL saxpy(OpenAD_lin_5,BVECTOR(ACOL_INDEX(I)),CVECTOR(AROW_INDEX(I)))
    END DO
  ELSE
    IF (ORDER.EQ.'POS') THEN
      DO I = 1,ANNZ,1
        OpenAD_lin_6 = BVECTOR(AROW_INDEX(I))%v
        OpenAD_lin_7 = AVALUES(I)%v
        CVECTOR(ACOL_INDEX(INT(I)))%v = (CVECTOR(ACOL_INDEX(I))%v+AVALUES(I)%v*B&
     &VECTOR(AROW_INDEX(I))%v)

        CALL setderiv(OpenAD_prp_4,CVECTOR(ACOL_INDEX(I)))
        CALL setderiv(CVECTOR(ACOL_INDEX(I)),OpenAD_prp_4)
        CALL saxpy(OpenAD_lin_6,AVALUES(I),CVECTOR(ACOL_INDEX(I)))
        CALL saxpy(OpenAD_lin_7,BVECTOR(AROW_INDEX(I)),CVECTOR(ACOL_INDEX(I)))
      END DO
    ENDIF
  ENDIF
  END SUBROUTINE

  SUBROUTINE SCALAR_MULTIPLY_SPMAT(N, ANNZ, AROW_INDEX, AROW_COMPRESSED, ACOL_IN&
     &DEX, AVALUES, SCALAR, RNNZ, RROW_INDEX, RROW_COMPRESSED, RCOL_INDEX, RVALU&
     &ES)

  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) N
  INTEGER(w2f__i4) ANNZ
  INTEGER(w2f__i4) AROW_INDEX(1 : INT((N * 7)))
  INTEGER(w2f__i4) AROW_COMPRESSED(1 : INT((N + 1)))
  INTEGER(w2f__i4) ACOL_INDEX(1 : INT((N * 7)))
  REAL(w2f__8) AVALUES(1 : INT((N * 7)))
  REAL(w2f__8) SCALAR
  INTEGER(w2f__i4) RNNZ
  INTEGER(w2f__i4) RROW_INDEX(1 : INT((N * 7)))
  INTEGER(w2f__i4) RROW_COMPRESSED(1 : INT((N + 1)))
  INTEGER(w2f__i4) RCOL_INDEX(1 : INT((N * 7)))
  REAL(w2f__8) RVALUES(1 : INT((N * 7)))
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
!
!       **** Statements ****
!
  RNNZ = ANNZ
  RROW_COMPRESSED(1 : INT(N + 1)) = AROW_COMPRESSED(1 : INT(N + 1 ))
  DO I = 1, ANNZ, 1
    RROW_INDEX(I) = AROW_INDEX(I)
    RCOL_INDEX(I) = ACOL_INDEX(I)
    RVALUES(INT(I)) = (AVALUES(I) * SCALAR)
  END DO
  END SUBROUTINE

  SUBROUTINE MYMIN_0_0_DOUBLE(SCALARIN1, SCALARIN2, SCALAROUT)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  REAL(w2f__8) SCALARIN1
  REAL(w2f__8) SCALARIN2
  REAL(w2f__8) SCALAROUT
!
!       **** Statements ****
!
  IF(SCALARIN1 .LE. SCALARIN2) THEN
    SCALAROUT = SCALARIN1
  ELSE
    SCALAROUT = SCALARIN2
  ENDIF
  END SUBROUTINE

  SUBROUTINE MYMIN_1_0_DOUBLE(VECTORIN, SCALARIN, VECTOROUT)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  type(active) :: VECTORIN(1:)
  REAL(w2f__8) SCALARIN
  type(active) :: VECTOROUT(1:)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
!
!       **** Statements ****
!
  DO I = 1, SIZE(VECTORIN, 1), 1
    IF (VECTORIN(I)%v.GE.SCALARIN) THEN
      VECTOROUT(INT(I))%v = SCALARIN
      CALL zero_deriv(VECTOROUT(INT(I)))
    ELSE
      VECTOROUT(INT(I))%v = VECTORIN(I)%v
      CALL setderiv(VECTOROUT(I),VECTORIN(I))
    ENDIF
  END DO
  END SUBROUTINE

  SUBROUTINE MYMIN_1_1_DOUBLE(VECTORIN1, VECTORIN2, VECTOROUT)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  REAL(w2f__8) VECTORIN1(1 :)
  REAL(w2f__8) VECTORIN2(1 :)
  REAL(w2f__8) VECTOROUT(1 :)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
!
!       **** Statements ****
!
  DO I = 1, SIZE(VECTORIN1, 1), 1
    IF(VECTORIN1(I) .GE. VECTORIN2(I)) THEN
      VECTOROUT(INT(I)) = VECTORIN2(I)
    ELSE
      VECTOROUT(INT(I)) = VECTORIN1(I)
    ENDIF
  END DO
  END SUBROUTINE

  SUBROUTINE MYMAX_0_0_DOUBLE(SCALARIN1, SCALARIN2, SCALAROUT)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  REAL(w2f__8) SCALARIN1
  REAL(w2f__8) SCALARIN2
  REAL(w2f__8) SCALAROUT
!
!       **** Statements ****
!
  IF(SCALARIN1 .GE. SCALARIN2) THEN
    SCALAROUT = SCALARIN1
  ELSE
    SCALAROUT = SCALARIN2
  ENDIF
  END SUBROUTINE

  SUBROUTINE MYMAX_1_0_DOUBLE(VECTORIN, SCALARIN, VECTOROUT)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  type(active) :: VECTORIN(1:)
  REAL(w2f__8) SCALARIN
  type(active) :: VECTOROUT(1:)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
!
!       **** Statements ****
!
  DO I = 1, SIZE(VECTORIN, 1), 1
    IF (VECTORIN(I)%v.LE.SCALARIN) THEN
      VECTOROUT(INT(I))%v = SCALARIN
      CALL zero_deriv(VECTOROUT(INT(I)))
    ELSE
      VECTOROUT(INT(I))%v = VECTORIN(I)%v
      CALL setderiv(VECTOROUT(I),VECTORIN(I))
    ENDIF
  END DO
  END SUBROUTINE

  SUBROUTINE MYMAX_1_1_DOUBLE(VECTORIN1, VECTORIN2, VECTOROUT)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  REAL(w2f__8) VECTORIN1(1 :)
  REAL(w2f__8) VECTORIN2(1 :)
  REAL(w2f__8) VECTOROUT(1 :)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
!
!       **** Statements ****
!
  DO I = 1, SIZE(VECTORIN1, 1), 1
    IF(VECTORIN1(I) .LE. VECTORIN2(I)) THEN
      VECTOROUT(INT(I)) = VECTORIN2(I)
    ELSE
      VECTOROUT(INT(I)) = VECTORIN1(I)
    ENDIF
  END DO
  END SUBROUTINE
END

MODULE linsolve
use OAD_active
use w2f__types
use oad_intrinsics
use parameters
use mathutil
use matrix
IMPLICIT NONE
SAVE
!
!     **** Top Level Pragmas ****
!
interface SOLVE
  module procedure SPARSE_SOLVE

end interface

!
!     **** Statements ****
!
CONTAINS

  SUBROUTINE SPARSE_SOLVE(N, ANNZ, AROW_INDEX, AROW_COMPRESSED, ACOL_INDEX, AVAL&
     &UES, B, X)

  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) N
  INTEGER(w2f__i4) ANNZ
  INTEGER(w2f__i4) AROW_INDEX(1 : INT((N * 7)))
  INTEGER(w2f__i4) AROW_COMPRESSED(1 : INT((N + 1)))
  INTEGER(w2f__i4) ACOL_INDEX(1 : INT((N * 7)))
  type(active) :: AVALUES(1:INT((N*7)))
  type(active) :: B(1:N)
  type(active) :: X(1:N)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) OAD_CTMP0
!
!       **** Statements ****
!
  OAD_CTMP0 = (N * 7)
  call sparse_pmgmres_method(N,ANNZ,OAD_CTMP0,AROW_INDEX,AROW_COMPRESSED,ACOL_INDE&
     &X,AVALUES,B,X,SOLVER_INNER,SOLVER_OUTER,VERBOSE)

  END SUBROUTINE

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

END

MODULE finitevolume
use OAD_active
use w2f__types
use oad_intrinsics
use parameters
use mathutil
use matrix
use linsolve
IMPLICIT NONE
SAVE
!
!     **** Top Level Pragmas ****
!
interface RELPERM
  module procedure RELPERM_VECTOR
  module procedure RELPERM_SCALAR

end interface

!
!     **** Statements ****
!
CONTAINS

  SUBROUTINE NEWTRAPH(NX, NY, NZ, ND, PT, ST, Q, V, S)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  INTEGER(w2f__i4) ND
  INTEGER(w2f__i4) PT
  INTEGER(w2f__i4) ST
  type(active) :: Q(1:INT((NZ*NX*NY)))
  type(active) :: V(1:3,1:INT((NX+1)),1:INT((NY+1)),1:INT((NZ+1)))
  type(active) :: S(1:INT((NZ*NX*NY)))
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) ACOL_INDEX(1 : INT((NZ * NX * NY * 7)))
  INTEGER(w2f__i4) ANNZ
  INTEGER(w2f__i4) AROW_COMPRESSED(1 : INT((NZ * NX * NY + 1)))
  INTEGER(w2f__i4) AROW_INDEX(1 : INT((NZ * NX * NY * 7)))
  type(active) :: AVALUES(1:INT((NZ*NX*NY*7)))
  INTEGER(w2f__i4) BCOL_INDEX(1 : INT((NZ * NX * NY * 7)))
  type(active) :: BFW(1:INT((NZ*NX*NY)))
  INTEGER(w2f__i4) BNNZ
  INTEGER(w2f__i4) BROW_COMPRESSED(1 : INT((NZ * NX * NY + 1)))
  INTEGER(w2f__i4) BROW_INDEX(1 : INT((NZ * NX * NY * 7)))
  type(active) :: BVALUES(1:INT((NZ*NX*NY*7)))
  LOGICAL(w2f__i4) CONVERGED
  type(active) :: DF(1:INT((NZ*NX*NY)))
  INTEGER(w2f__i4) DGCOL_INDEX(1 : INT((NZ * NX * NY * 7)))
  INTEGER(w2f__i4) DGNNZ
  INTEGER(w2f__i4) DGROW_COMPRESSED(1 : INT((NZ * NX * NY + 1)))
  INTEGER(w2f__i4) DGROW_INDEX(1 : INT((NZ * NX * NY * 7)))
  type(active) :: DGVALUES(1:INT((NZ*NX*NY*7)))
  type(active) :: DMO(1:INT((NZ*NX*NY)))
  type(active) :: DMW(1:INT((NZ*NX*NY)))
  type(active) :: DS(1:INT((NZ*NX*NY)))
  REAL(w2f__8) DSN
  REAL(w2f__8) DT
  REAL(w2f__8) DTX(1 : INT((NZ * NX * NY)))
  type(active) :: FI(1:INT((NZ*NX*NY)))
  type(active) :: FW(1:INT((NZ*NX*NY)))
  type(active) :: G(1:INT((NZ*NX*NY)))
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) IT
  INTEGER(w2f__i4) J
  type(active) :: MO(1:INT((NZ*NX*NY)))
  type(active) :: MW(1:INT((NZ*NX*NY)))
  INTEGER(w2f__i4) N
  type(active) :: S_COPY(1:INT((NZ*NX*NY)))
  type(active) :: S_ITER_COPY(1:INT((NZ*NX*NY)))
  REAL(w2f__8) OpenAD_acc_0(:)
  ALLOCATABLE OpenAD_acc_0
  REAL(w2f__8) OpenAD_acc_1(:)
  ALLOCATABLE OpenAD_acc_1
  REAL(w2f__8) OpenAD_acc_2(:)
  ALLOCATABLE OpenAD_acc_2
  REAL(w2f__8) OpenAD_aux_0(:)
  ALLOCATABLE OpenAD_aux_0
  REAL(w2f__8) OpenAD_aux_1(:)
  ALLOCATABLE OpenAD_aux_1
  REAL(w2f__8) OpenAD_aux_2(:)
  ALLOCATABLE OpenAD_aux_2
  REAL(w2f__8) OpenAD_aux_3(:)
  ALLOCATABLE OpenAD_aux_3
  REAL(w2f__8) OpenAD_aux_4(:)
  ALLOCATABLE OpenAD_aux_4
  REAL(w2f__8) OpenAD_aux_5(:)
  ALLOCATABLE OpenAD_aux_5
  REAL(w2f__8) OpenAD_lin_10(:)
  ALLOCATABLE OpenAD_lin_10
  REAL(w2f__8) OpenAD_lin_11(:)
  ALLOCATABLE OpenAD_lin_11
  REAL(w2f__8) OpenAD_lin_12(:)
  ALLOCATABLE OpenAD_lin_12
  REAL(w2f__8) OpenAD_lin_13(:)
  ALLOCATABLE OpenAD_lin_13
  REAL(w2f__8) OpenAD_lin_14(:)
  ALLOCATABLE OpenAD_lin_14
  REAL(w2f__8) OpenAD_lin_15(:)
  ALLOCATABLE OpenAD_lin_15
  REAL(w2f__8) OpenAD_lin_16(:)
  ALLOCATABLE OpenAD_lin_16
  REAL(w2f__8) OpenAD_lin_17(:)
  ALLOCATABLE OpenAD_lin_17
  REAL(w2f__8) OpenAD_lin_8(:)
  ALLOCATABLE OpenAD_lin_8
  REAL(w2f__8) OpenAD_lin_9(:)
  ALLOCATABLE OpenAD_lin_9
  type(active) :: OpenAD_prp_10(:)
  ALLOCATABLE OpenAD_prp_10
  type(active) :: OpenAD_prp_5(:)
  ALLOCATABLE OpenAD_prp_5
  type(active) :: OpenAD_prp_6(:)
  ALLOCATABLE OpenAD_prp_6
  type(active) :: OpenAD_prp_7(:)
  ALLOCATABLE OpenAD_prp_7
  type(active) :: OpenAD_prp_8(:)
  ALLOCATABLE OpenAD_prp_8
  type(active) :: OpenAD_prp_9(:)
  ALLOCATABLE OpenAD_prp_9
  type(active) :: OpenAD_tyc_0(:)
  ALLOCATABLE OpenAD_tyc_0
  REAL(w2f__8) OpenAD_tyc_1(:)
  ALLOCATABLE OpenAD_tyc_1
!
!       **** Statements ****
!
  CONVERGED = .FALSE.
  N = (NZ * NX * NY)
  CALL GENA(NX,NY,NZ,V,Q,ANNZ,AROW_INDEX,AROW_COMPRESSED,ACOL_INDEX,AVALUES)
  S_COPY(1:NZ*NX*NY)%v = S(1:NZ*NX*NY)%v
  IT = 0
  CALL setderiv(S_COPY(1:NZ*NX*NY),S(1:NZ*NX*NY))
  do while (.not. CONVERGED)
    DT = (ST/(2**IT))
    DTX(1:NZ*NX*NY) = (DT/(POR(:)*VOL))
    CALL MYMAX_1_0_DOUBLE(Q,0.0D00,FI)
!         $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
    CALL oad_AllocateMatching(OpenAD_lin_8,DTX(1:NZ*NX*NY))
    OpenAD_lin_8 = DTX(1:NZ*NX*NY)
    FI(1:NZ*NX*NY)%v = (DTX(1:NZ*NX*NY)*FI(1:NZ*NX*NY)%v)
!         $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
    CALL oad_AllocateMatching(OpenAD_prp_5,FI(1:NZ*NX*NY))
    CALL setderiv(OpenAD_prp_5,FI(1:NZ*NX*NY))
    CALL sax(OpenAD_lin_8,OpenAD_prp_5,FI(1:NZ*NX*NY))
!         $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
    CALL oad_AllocateMatching(OpenAD_tyc_0,DTX)
!         $OpenAD$ INLINE oad_convert(subst,subst)
    CALL oad_convert(OpenAD_tyc_0,DTX)
    CALL SPMAT_MULTIPLY_DIAGONAL(N,ANNZ,AROW_INDEX,AROW_COMPRESSED,ACOL_INDEX,AV&
     &ALUES,OpenAD_tyc_0,BNNZ,BROW_INDEX,BROW_COMPRESSED,BCOL_INDEX,BVALUES,'POS&
     &')

!         $OpenAD$ INLINE oad_ShapeTest(subst,subst)
    CALL oad_ShapeTest(OpenAD_tyc_0,DTX)
!         $OpenAD$ INLINE oad_convert(subst,subst)
    CALL oad_convert(DTX,OpenAD_tyc_0)
    I = 0
    do while (I.LT.(2**IT))
      J = 0
      I = (I+1)
      DSN = 1.0D00
      S_ITER_COPY(1:NZ*NX*NY)%v = S(1:NZ*NX*NY)%v
      CALL setderiv(S_ITER_COPY(1:NZ*NX*NY),S(1:NZ*NX*NY))
      do while ((J.LT.10).AND.(DSN.GT.1.00000000000000002082D-03))
        CALL RELPERM_VECTOR(NX,NY,NZ,S,MW,MO,DMW,DMO)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_aux_0,MW(1:NZ*NX*NY))
        OpenAD_aux_0 = (MO(1:NZ*NX*NY)%v+MW(1:NZ*NX*NY)%v)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_aux_4,MW(1:NZ*NX*NY))
        OpenAD_aux_4 = (MO(1:NZ*NX*NY)%v+MW(1:NZ*NX*NY)%v)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_aux_3,MW(1:NZ*NX*NY))
        OpenAD_aux_3 = (OpenAD_aux_4**2)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_aux_1,MW(1:NZ*NX*NY))
        OpenAD_aux_1 = (MW(1:NZ*NX*NY)%v/OpenAD_aux_3)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_aux_2,DMW(1:NZ*NX*NY))
        OpenAD_aux_2 = (DMO(1:NZ*NX*NY)%v+DMW(1:NZ*NX*NY)%v)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_lin_9,OpenAD_aux_0)
        OpenAD_lin_9 = (INT(1_w2f__i8)/OpenAD_aux_0)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_lin_10,OpenAD_aux_0)
        OpenAD_lin_10 = (-(DMW(1:NZ*NX*NY)%v/(OpenAD_aux_0*OpenAD_aux_0)))
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_lin_13,OpenAD_aux_3)
        OpenAD_lin_13 = (INT(1_w2f__i8)/OpenAD_aux_3)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_lin_15,OpenAD_aux_4)
        OpenAD_lin_15 = (2*(OpenAD_aux_4**(2-INT(1_w2f__i8))))
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_lin_14,OpenAD_aux_3)
        OpenAD_lin_14 = (-(MW(1:NZ*NX*NY)%v/(OpenAD_aux_3*OpenAD_aux_3)))
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_lin_11,OpenAD_aux_2)
        OpenAD_lin_11 = OpenAD_aux_2
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_lin_12,OpenAD_aux_1)
        OpenAD_lin_12 = OpenAD_aux_1
        DF(1:NZ*NX*NY)%v = ((DMW(1:NZ*NX*NY)%v/OpenAD_aux_0)-OpenAD_aux_1*OpenAD&
     &_aux_2)

!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_acc_0,OpenAD_lin_12)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_acc_1,OpenAD_lin_11)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_acc_2,OpenAD_lin_11)
        OpenAD_acc_0 = (OpenAD_lin_12*INT((-1_w2f__i8)))
        OpenAD_acc_1 = (OpenAD_lin_13*OpenAD_lin_11*INT((-1_w2f__i8)))
        OpenAD_acc_2 = (OpenAD_lin_15*OpenAD_lin_14*OpenAD_lin_11*INT((-1_w2f__i&
     &8)))

!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_prp_6,MO(1:NZ*NX*NY))
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_prp_7,MO(1:NZ*NX*NY))
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_prp_8,DMO(1:NZ*NX*NY))
        CALL setderiv(OpenAD_prp_6,MO(1:NZ*NX*NY))
        CALL inc_deriv(OpenAD_prp_6,MW(1:NZ*NX*NY))
        CALL setderiv(OpenAD_prp_7,MO(1:NZ*NX*NY))
        CALL inc_deriv(OpenAD_prp_7,MW(1:NZ*NX*NY))
        CALL setderiv(OpenAD_prp_8,DMO(1:NZ*NX*NY))
        CALL inc_deriv(OpenAD_prp_8,DMW(1:NZ*NX*NY))
        CALL sax(OpenAD_lin_9,DMW(1:NZ*NX*NY),DF(1:NZ*NX*NY))
        CALL saxpy(OpenAD_lin_10,OpenAD_prp_6,DF(1:NZ*NX*NY))
        CALL saxpy(OpenAD_acc_0,OpenAD_prp_8,DF(1:NZ*NX*NY))
        CALL saxpy(OpenAD_acc_1,MW(1:NZ*NX*NY),DF(1:NZ*NX*NY))
        CALL saxpy(OpenAD_acc_2,OpenAD_prp_7,DF(1:NZ*NX*NY))
        CALL SPMAT_MULTIPLY_DIAGONAL(N,BNNZ,BROW_INDEX,BROW_COMPRESSED,BCOL_INDE&
     &X,BVALUES,DF,DGNNZ,DGROW_INDEX,DGROW_COMPRESSED,DGCOL_INDEX,DGVALUES,'PRE'&
     &)

        CALL ADDX_DIAGONAL(N,DGNNZ,DGROW_INDEX,DGROW_COMPRESSED,DGCOL_INDEX,DGVA&
     &LUES,-1.0D00,0)

!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_aux_5,MW(1:NZ*NX*NY))
        OpenAD_aux_5 = (MO(1:NZ*NX*NY)%v+MW(1:NZ*NX*NY)%v)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_lin_16,OpenAD_aux_5)
        OpenAD_lin_16 = (INT(1_w2f__i8)/OpenAD_aux_5)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_lin_17,OpenAD_aux_5)
        OpenAD_lin_17 = (-(MW(1:NZ*NX*NY)%v/(OpenAD_aux_5*OpenAD_aux_5)))
        FW(1:NZ*NX*NY)%v = (MW(1:NZ*NX*NY)%v/OpenAD_aux_5)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_prp_9,MO(1:NZ*NX*NY))
        CALL setderiv(OpenAD_prp_9,MO(1:NZ*NX*NY))
        CALL inc_deriv(OpenAD_prp_9,MW(1:NZ*NX*NY))
        CALL sax(OpenAD_lin_16,MW(1:NZ*NX*NY),FW(1:NZ*NX*NY))
        CALL saxpy(OpenAD_lin_17,OpenAD_prp_9,FW(1:NZ*NX*NY))
        CALL SPMAT_MULTIPLY_VECTOR(N,BNNZ,BROW_INDEX,BROW_COMPRESSED,BCOL_INDEX,&
     &BVALUES,FW,BFW,'PRE')

        G(1:NZ*NX*NY)%v = (S(1:NZ*NX*NY)%v-S_ITER_COPY(1:NZ*NX*NY)%v-BFW(1:NZ*NX&
     &*NY)%v-FI(1:NZ*NX*NY)%v)

        CALL setderiv(G(1:NZ*NX*NY),S(1:NZ*NX*NY))
        CALL dec_deriv(G(1:NZ*NX*NY),FI(1:NZ*NX*NY))
        CALL dec_deriv(G(1:NZ*NX*NY),BFW(1:NZ*NX*NY))
        CALL dec_deriv(G(1:NZ*NX*NY),S_ITER_COPY(1:NZ*NX*NY))
        CALL SPARSE_SOLVE(N,DGNNZ,DGROW_INDEX,DGROW_COMPRESSED,DGCOL_INDEX,DGVAL&
     &UES,G,DS)

        S(1:NZ*NX*NY)%v = (S(1:NZ*NX*NY)%v+DS(1:NZ*NX*NY)%v)
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_prp_10,S(1:NZ*NX*NY))
        CALL setderiv(OpenAD_prp_10,S(1:NZ*NX*NY))
        CALL setderiv(S(1:NZ*NX*NY),OpenAD_prp_10)
        CALL inc_deriv(S(1:NZ*NX*NY),DS(1:NZ*NX*NY))
!             $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
        CALL oad_AllocateMatching(OpenAD_tyc_1,DS)
!             $OpenAD$ INLINE oad_convert(subst,subst)
        CALL oad_convert(OpenAD_tyc_1,DS)
        CALL DNRM2(OpenAD_tyc_1,N,DSN)
!             $OpenAD$ INLINE oad_ShapeTest(subst,subst)
        CALL oad_ShapeTest(OpenAD_tyc_1,DS)
!             $OpenAD$ INLINE oad_convert(subst,subst)
        CALL oad_convert(DS,OpenAD_tyc_1)
        J = (J+1)
      END DO
      IF (DSN.GT.1.00000000000000002082D-03) THEN
        I = (2**IT)
        S(1:NZ*NX*NY)%v = S_COPY(1:NZ*NX*NY)%v
        CALL setderiv(S(1:NZ*NX*NY),S_COPY(1:NZ*NX*NY))
      ENDIF
    END DO
    IF (DSN.LT.1.00000000000000002082D-03) THEN
      CONVERGED = .TRUE.
    ELSE
      IT = (IT+1)
    ENDIF
  END DO
  END SUBROUTINE

  SUBROUTINE PRES(NX, NY, NZ, Q, S, P, V)
  use w2f__types
  use parameters
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  type(active) :: Q(1:INT((NZ*NX*NY)))
  type(active) :: S(1:INT((NZ*NX*NY)))
  type(active) :: P(1:NX,1:NY,1:NZ)
  type(active) :: V(1:3,1:INT((NX+1)),1:INT((NY+1)),1:INT((NZ+1)))
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) H
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) J
  type(active) :: KM(1:3,1:NX,1:NY,1:NZ)
  type(active) :: M(1:INT((NZ*NX*NY*3)))
  type(active) :: MO(1:INT((NZ*NX*NY)))
  type(active) :: MW(1:INT((NZ*NX*NY)))
  INTEGER(w2f__i4) N
  REAL(w2f__8) OpenAD_lin_27
  type(active) :: OpenAD_prp_17
  type(active) :: OpenAD_prp_18
  type(active) :: OpenAD_prp_19
!
!       **** Statements ****
!
  CALL RELPERM_VECTOR(NX,NY,NZ,S,MW,MO)
  DO I = 1,(NZ*NX*NY),1
    M(I*3+(-2))%v = (MO(I)%v+MW(I)%v)
    CALL setderiv(M(I*3+(-2)),MO(I))
    CALL inc_deriv(M(I*3+(-2)),MW(I))
    M(I*3+(-1))%v = M(I*3+(-2))%v
    CALL setderiv(OpenAD_prp_17,M(I*3+(-2)))
    CALL setderiv(M(I*3+(-1)),OpenAD_prp_17)
    M(I*3)%v = M(I*3+(-2))%v
    CALL setderiv(OpenAD_prp_18,M(I*3+(-2)))
    CALL setderiv(M(I*3),OpenAD_prp_18)
  END DO
  CALL MYRESHAPE_1_4(M,KM)
  DO H = 1,3,1
    DO I = 1,NX,1
      DO J = 1,NY,1
        DO N = 1,NZ,1
          OpenAD_lin_27 = PERM(H,I,J,N)
          KM(INT(H),INT(I),INT(J),INT(N))%v = (KM(H,I,J,N)%v*PERM(H,I,J,N))
          CALL setderiv(OpenAD_prp_19,KM(H,I,J,N))
          CALL sax(OpenAD_lin_27,OpenAD_prp_19,KM(H,I,J,N))
        END DO
      END DO
    END DO
  END DO
  CALL TPFA(NX,NY,NZ,KM,Q,P,V)
  END SUBROUTINE

  SUBROUTINE RELPERM_VECTOR(NX, NY, NZ, S, MW, MO, DMW, DMO)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  type(active) :: S(1:INT((NZ*NX*NY)))
  type(active) :: MW(1:INT((NZ*NX*NY)))
  type(active) :: MO(1:INT((NZ*NX*NY)))
  type(active) :: DMW(1:INT((NZ*NX*NY)))
  OPTIONAL DMW
  type(active) :: DMO(1:INT((NZ*NX*NY)))
  OPTIONAL DMO
!
!       **** Local Variables and Functions ****
!
  type(active) :: S_TEMP(1:INT((NZ*NX*NY)))
  LOGICAL(w2f__i4) t__196
  REAL(w2f__8) OpenAD_acc_3(:)
  ALLOCATABLE OpenAD_acc_3
  REAL(w2f__8) OpenAD_acc_4(:)
  ALLOCATABLE OpenAD_acc_4
  REAL(w2f__8) OpenAD_acc_5
  REAL(w2f__8) OpenAD_acc_6
  REAL(w2f__8) OpenAD_aux_10(:)
  ALLOCATABLE OpenAD_aux_10
  REAL(w2f__8) OpenAD_aux_11(:)
  ALLOCATABLE OpenAD_aux_11
  REAL(w2f__8) OpenAD_aux_12
  REAL(w2f__8) OpenAD_aux_13(:)
  ALLOCATABLE OpenAD_aux_13
  REAL(w2f__8) OpenAD_aux_14(:)
  ALLOCATABLE OpenAD_aux_14
  REAL(w2f__8) OpenAD_aux_15
  REAL(w2f__8) OpenAD_aux_16(:)
  ALLOCATABLE OpenAD_aux_16
  REAL(w2f__8) OpenAD_aux_17(:)
  ALLOCATABLE OpenAD_aux_17
  REAL(w2f__8) OpenAD_aux_6(:)
  ALLOCATABLE OpenAD_aux_6
  REAL(w2f__8) OpenAD_aux_7
  REAL(w2f__8) OpenAD_aux_8(:)
  ALLOCATABLE OpenAD_aux_8
  REAL(w2f__8) OpenAD_aux_9(:)
  ALLOCATABLE OpenAD_aux_9
  REAL(w2f__8) OpenAD_lin_18
  REAL(w2f__8) OpenAD_lin_19
  REAL(w2f__8) OpenAD_lin_20(:)
  ALLOCATABLE OpenAD_lin_20
  REAL(w2f__8) OpenAD_lin_21
  REAL(w2f__8) OpenAD_lin_22(:)
  ALLOCATABLE OpenAD_lin_22
  REAL(w2f__8) OpenAD_lin_23
  REAL(w2f__8) OpenAD_lin_24
  REAL(w2f__8) OpenAD_lin_25
  REAL(w2f__8) OpenAD_lin_26
!
!       **** Statements ****
!
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_aux_6,S(1:NZ*NX*NY))
  OpenAD_aux_6 = (S(1:NZ*NX*NY)%v-SWC)
  OpenAD_aux_7 = (1.0D00-SWC-SOR)
  OpenAD_lin_18 = (INT(1_w2f__i8)/OpenAD_aux_7)
  S_TEMP(1:NZ*NX*NY)%v = (OpenAD_aux_6/OpenAD_aux_7)
  CALL sax(OpenAD_lin_18,S(1:NZ*NX*NY),S_TEMP(1:NZ*NX*NY))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_aux_8,S_TEMP(1:NZ*NX*NY))
  OpenAD_aux_8 = (S_TEMP(1:NZ*NX*NY)%v**2)
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_lin_20,S_TEMP(1:NZ*NX*NY))
  OpenAD_lin_20 = (2*(S_TEMP(1:NZ*NX*NY)%v**(2-INT(1_w2f__i8))))
  OpenAD_lin_19 = (INT(1_w2f__i8)/VW)
  MW(1:NZ*NX*NY)%v = (OpenAD_aux_8/VW)
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_aux_10,S_TEMP(1:NZ*NX*NY))
  OpenAD_aux_10 = (1.0D00-S_TEMP(1:NZ*NX*NY)%v)
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_aux_9,S_TEMP(1:NZ*NX*NY))
  OpenAD_aux_9 = (OpenAD_aux_10**2)
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_lin_22,OpenAD_aux_10)
  OpenAD_lin_22 = (2*(OpenAD_aux_10**(2-INT(1_w2f__i8))))
  OpenAD_lin_21 = (INT(1_w2f__i8)/VO)
  MO(1:NZ*NX*NY)%v = (OpenAD_aux_9/VO)
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_acc_3,OpenAD_lin_20)
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_acc_4,OpenAD_lin_22)
  OpenAD_acc_3 = (OpenAD_lin_20*OpenAD_lin_19)
  OpenAD_acc_4 = (INT((-1_w2f__i8))*OpenAD_lin_22*OpenAD_lin_21)
  CALL sax(OpenAD_acc_3,S_TEMP(1:NZ*NX*NY),MW(1:NZ*NX*NY))
  CALL sax(OpenAD_acc_4,S_TEMP(1:NZ*NX*NY),MO(1:NZ*NX*NY))
  t__196 = .TRUE.
  IF (.not. PRESENT(DMO)) THEN
    t__196 = .FALSE.
  ELSE
    IF (.not. PRESENT(DMW)) THEN
      t__196 = .FALSE.
    ENDIF
  ENDIF
  IF(t__196) THEN
!         $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
    CALL oad_AllocateMatching(OpenAD_aux_13,S_TEMP(1:NZ*NX*NY))
    OpenAD_aux_13 = (S_TEMP(1:NZ*NX*NY)%v*2.0D00)
!         $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
    CALL oad_AllocateMatching(OpenAD_aux_11,S_TEMP(1:NZ*NX*NY))
    OpenAD_aux_11 = (OpenAD_aux_13/VW)
    OpenAD_aux_12 = (1.0D00-SWC-SOR)
    OpenAD_lin_24 = (INT(1_w2f__i8)/VW)
    OpenAD_lin_23 = (INT(1_w2f__i8)/OpenAD_aux_12)
    DMW(1:NZ*NX*NY)%v = (OpenAD_aux_11/OpenAD_aux_12)
!         $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
    CALL oad_AllocateMatching(OpenAD_aux_17,S_TEMP(1:NZ*NX*NY))
    OpenAD_aux_17 = (1.0D00-S_TEMP(1:NZ*NX*NY)%v)
!         $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
    CALL oad_AllocateMatching(OpenAD_aux_16,S_TEMP(1:NZ*NX*NY))
    OpenAD_aux_16 = (OpenAD_aux_17*2.0D00)
!         $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
    CALL oad_AllocateMatching(OpenAD_aux_14,S_TEMP(1:NZ*NX*NY))
    OpenAD_aux_14 = (OpenAD_aux_16/VO)
    OpenAD_aux_15 = (1.0D00-SWC-SOR)
    OpenAD_lin_26 = (INT(1_w2f__i8)/VO)
    OpenAD_lin_25 = (INT(1_w2f__i8)/OpenAD_aux_15)
    DMO(1:NZ*NX*NY)%v = (-(OpenAD_aux_14/OpenAD_aux_15))
    OpenAD_acc_5 = (2.0D00*OpenAD_lin_24*OpenAD_lin_23)
    OpenAD_acc_6 = (INT((-1_w2f__i8))*2.0D00*OpenAD_lin_26*OpenAD_lin_25*INT((-1&
     &_w2f__i8)))

    CALL sax(OpenAD_acc_5,S_TEMP(1:NZ*NX*NY),DMW(1:NZ*NX*NY))
    CALL sax(OpenAD_acc_6,S_TEMP(1:NZ*NX*NY),DMO(1:NZ*NX*NY))
  ENDIF
  END SUBROUTINE

  SUBROUTINE RELPERM_SCALAR(S, MW, MO, DMW, DMO)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  type(active) :: S
  type(active) :: MW
  type(active) :: MO
  REAL(w2f__8) DMW
  OPTIONAL DMW
  REAL(w2f__8) DMO
  OPTIONAL DMO
!
!       **** Local Variables and Functions ****
!
  type(active) :: S_TEMP
  LOGICAL(w2f__i4) t__197
  REAL(w2f__8) OpenAD_acc_7
  REAL(w2f__8) OpenAD_acc_8
  REAL(w2f__8) OpenAD_aux_24
  REAL(w2f__8) OpenAD_aux_25
  REAL(w2f__8) OpenAD_aux_26
  REAL(w2f__8) OpenAD_aux_27
  REAL(w2f__8) OpenAD_aux_28
  REAL(w2f__8) OpenAD_lin_38
  REAL(w2f__8) OpenAD_lin_39
  REAL(w2f__8) OpenAD_lin_40
  REAL(w2f__8) OpenAD_lin_41
  REAL(w2f__8) OpenAD_lin_42
!
!       **** Statements ****
!
  OpenAD_aux_24 = (S%v-SWC)
  OpenAD_aux_25 = (1.0D00-SWC-SOR)
  OpenAD_lin_38 = (INT(1_w2f__i8)/OpenAD_aux_25)
  S_TEMP%v = (OpenAD_aux_24/OpenAD_aux_25)
  OpenAD_aux_26 = (S_TEMP%v**2)
  OpenAD_lin_40 = (2*(S_TEMP%v**(2-INT(1_w2f__i8))))
  OpenAD_lin_39 = (INT(1_w2f__i8)/VW)
  MW%v = (OpenAD_aux_26/VW)
  OpenAD_aux_28 = (1.0D00-S_TEMP%v)
  OpenAD_aux_27 = (OpenAD_aux_28**2)
  OpenAD_lin_42 = (2*(OpenAD_aux_28**(2-INT(1_w2f__i8))))
  OpenAD_lin_41 = (INT(1_w2f__i8)/VO)
  MO%v = (OpenAD_aux_27/VO)
  OpenAD_acc_7 = (OpenAD_lin_40*OpenAD_lin_39)
  OpenAD_acc_8 = (INT((-1_w2f__i8))*OpenAD_lin_42*OpenAD_lin_41)
  CALL sax(OpenAD_lin_38,S,S_TEMP)
  CALL sax(OpenAD_acc_7,S_TEMP,MW)
  CALL sax(OpenAD_acc_8,S_TEMP,MO)
  t__197 = .TRUE.
  IF (.not. PRESENT(DMO)) THEN
    t__197 = .FALSE.
  ELSE
    IF (.not. PRESENT(DMW)) THEN
      t__197 = .FALSE.
    ENDIF
  ENDIF
  IF(t__197) THEN
    DMW = (((S_TEMP%v*2.0D00)/VW)/(1.0D00-SWC-SOR))
    DMO = (-((((1.0D00-S_TEMP%v)*2.0D00)/VO)/(1.0D00-SWC-SOR)))
  ENDIF
  END SUBROUTINE

  SUBROUTINE GENA(NX, NY, NZ, V, Q, ANNZ, AROW_INDEX, AROW_COMPRESSED, ACOL_INDE&
     &X, AVALUES)

  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  type(active) :: V(1:3,1:INT((NX+1)),1:INT((NY+1)),1:INT((NZ+1)))
  type(active) :: Q(1:INT((NZ*NX*NY)))
  INTEGER(w2f__i4) ANNZ
  INTEGER(w2f__i4) AROW_INDEX(1 : INT((NZ * NX * NY * 7)))
  INTEGER(w2f__i4) AROW_COMPRESSED(1 : INT((NZ * NX * NY + 1)))
  INTEGER(w2f__i4) ACOL_INDEX(1 : INT((NZ * NX * NY * 7)))
  type(active) :: AVALUES(1:INT((NZ*NX*NY*7)))
!
!       **** Local Variables and Functions ****
!
  type(active) :: DIAGS(1:INT((NZ*NX*NY)),1:7)
  type(active) :: DIAG_TMP(1:INT((NZ*NX*NY)))
  type(active) :: VXYZ(1:NX,1:NY,1:NZ)
  type(active) :: OpenAD_prp_11(:)
  ALLOCATABLE OpenAD_prp_11
  type(active) :: OpenAD_prp_12(:)
  ALLOCATABLE OpenAD_prp_12
  type(active) :: OpenAD_prp_13(:)
  ALLOCATABLE OpenAD_prp_13
  type(active) :: OpenAD_prp_14(:)
  ALLOCATABLE OpenAD_prp_14
  type(active) :: OpenAD_prp_15(:)
  ALLOCATABLE OpenAD_prp_15
  type(active) :: OpenAD_prp_16(:)
  ALLOCATABLE OpenAD_prp_16
!
!       **** Statements ****
!
  DIAGS(1:NZ*NX*NY,1:7)%v = 0.0D00
  VXYZ(1:INT(NX),1:INT(NY),1:INT(NZ))%v = V(3,1:NX,1:NY,2:NZ+1)%v
  CALL zero_deriv(DIAGS(1:NZ*NX*NY,1:7))
  CALL setderiv(VXYZ(1:NX,1:NY,1:NZ),V(3,1:NX,1:NY,2:NZ+1))
  CALL MYRESHAPE_3_1(VXYZ,DIAGS(1:NZ*NX*NY,1))
  VXYZ(1:INT(NX),1:INT(NY),1:INT(NZ))%v = V(2,1:NX,2:NY+1,1:NZ)%v
  CALL setderiv(VXYZ(1:NX,1:NY,1:NZ),V(2,1:NX,2:NY+1,1:NZ))
  CALL MYRESHAPE_3_1(VXYZ,DIAGS(1:NZ*NX*NY,2))
  VXYZ(1:INT(NX),1:INT(NY),1:INT(NZ))%v = V(1,2:NX+1,1:NY,1:NZ)%v
  CALL setderiv(VXYZ(1:NX,1:NY,1:NZ),V(1,2:NX+1,1:NY,1:NZ))
  CALL MYRESHAPE_3_1(VXYZ,DIAGS(1:NZ*NX*NY,3))
  VXYZ(1:INT(NX),1:INT(NY),1:INT(NZ))%v = V(1,1:NX,1:NY,1:NZ)%v
  CALL setderiv(VXYZ(1:NX,1:NY,1:NZ),V(1,1:NX,1:NY,1:NZ))
  CALL MYRESHAPE_3_1(VXYZ,DIAGS(1:NZ*NX*NY,5))
  VXYZ(1:INT(NX),1:INT(NY),1:INT(NZ))%v = V(2,1:NX,1:NY,1:NZ)%v
  CALL setderiv(VXYZ(1:NX,1:NY,1:NZ),V(2,1:NX,1:NY,1:NZ))
  CALL MYRESHAPE_3_1(VXYZ,DIAGS(1:NZ*NX*NY,6))
  VXYZ(1:INT(NX),1:INT(NY),1:INT(NZ))%v = V(3,1:NX,1:NY,1:NZ)%v
  CALL setderiv(VXYZ(1:NX,1:NY,1:NZ),V(3,1:NX,1:NY,1:NZ))
  CALL MYRESHAPE_3_1(VXYZ,DIAGS(1:NZ*NX*NY,7))
  DIAG_TMP(1:NZ*NX*NY)%v = 0.0D00
  CALL zero_deriv(DIAG_TMP(1:NZ*NX*NY))
  CALL MYMAX_1_0_DOUBLE(DIAGS(1:NZ*NX*NY,1),0.0D00,DIAG_TMP)
  DIAGS(1:NZ*NX*NY,1)%v = DIAG_TMP(1:NZ*NX*NY)%v
  DIAG_TMP(1:NZ*NX*NY)%v = 0.0D00
  CALL zero_deriv(DIAG_TMP(1:NZ*NX*NY))
  CALL setderiv(DIAGS(1:NZ*NX*NY,1),DIAG_TMP(1:NZ*NX*NY))
  CALL MYMAX_1_0_DOUBLE(DIAGS(1:NZ*NX*NY,2),0.0D00,DIAG_TMP)
  DIAGS(1:NZ*NX*NY,2)%v = DIAG_TMP(1:NZ*NX*NY)%v
  DIAG_TMP(1:NZ*NX*NY)%v = 0.0D00
  CALL zero_deriv(DIAG_TMP(1:NZ*NX*NY))
  CALL setderiv(DIAGS(1:NZ*NX*NY,2),DIAG_TMP(1:NZ*NX*NY))
  CALL MYMAX_1_0_DOUBLE(DIAGS(1:NZ*NX*NY,3),0.0D00,DIAG_TMP)
  DIAGS(1:NZ*NX*NY,3)%v = DIAG_TMP(1:NZ*NX*NY)%v
  DIAG_TMP(1:NZ*NX*NY)%v = 0.0D00
  CALL zero_deriv(DIAG_TMP(1:NZ*NX*NY))
  CALL setderiv(DIAGS(1:NZ*NX*NY,3),DIAG_TMP(1:NZ*NX*NY))
  CALL MYMIN_1_0_DOUBLE(DIAGS(1:NZ*NX*NY,5),0.0D00,DIAG_TMP)
  DIAGS(1:NZ*NX*NY,5)%v = (-DIAG_TMP(1:NZ*NX*NY)%v)
  DIAG_TMP(1:NZ*NX*NY)%v = 0.0D00
  CALL zero_deriv(DIAG_TMP(1:NZ*NX*NY))
  CALL set_neg_deriv(DIAGS(1:NZ*NX*NY,5),DIAG_TMP(1:NZ*NX*NY))
  CALL MYMIN_1_0_DOUBLE(DIAGS(1:NZ*NX*NY,6),0.0D00,DIAG_TMP)
  DIAGS(1:NZ*NX*NY,6)%v = (-DIAG_TMP(1:NZ*NX*NY)%v)
  DIAG_TMP(1:NZ*NX*NY)%v = 0.0D00
  CALL zero_deriv(DIAG_TMP(1:NZ*NX*NY))
  CALL set_neg_deriv(DIAGS(1:NZ*NX*NY,6),DIAG_TMP(1:NZ*NX*NY))
  CALL MYMIN_1_0_DOUBLE(DIAGS(1:NZ*NX*NY,7),0.0D00,DIAG_TMP)
  DIAGS(1:NZ*NX*NY,7)%v = (-DIAG_TMP(1:NZ*NX*NY)%v)
  DIAG_TMP(1:NZ*NX*NY)%v = 0.0D00
  CALL zero_deriv(DIAG_TMP(1:NZ*NX*NY))
  CALL set_neg_deriv(DIAGS(1:NZ*NX*NY,7),DIAG_TMP(1:NZ*NX*NY))
  CALL MYMIN_1_0_DOUBLE(Q,0.0D00,DIAG_TMP)
  DIAGS(1:NZ*NX*NY,4)%v = (DIAG_TMP(1:NZ*NX*NY)%v-DIAGS(1:NZ*NX*NY,5)%v-DIAGS(1:&
     &NZ*NX*NY,3)%v-DIAGS(1:NZ*NX*NY,6)%v-DIAGS(1:NZ*NX*NY,2)%v-DIAGS(1:NZ*NX*NY&
     &,7)%v-DIAGS(1:NZ*NX*NY,1)%v)

!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_11,DIAGS(1:NZ*NX*NY,5))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_12,DIAGS(1:NZ*NX*NY,3))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_13,DIAGS(1:NZ*NX*NY,6))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_14,DIAGS(1:NZ*NX*NY,2))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_15,DIAGS(1:NZ*NX*NY,7))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_16,DIAGS(1:NZ*NX*NY,1))
  CALL setderiv(OpenAD_prp_11,DIAGS(1:NZ*NX*NY,5))
  CALL setderiv(OpenAD_prp_12,DIAGS(1:NZ*NX*NY,3))
  CALL setderiv(OpenAD_prp_13,DIAGS(1:NZ*NX*NY,6))
  CALL setderiv(OpenAD_prp_14,DIAGS(1:NZ*NX*NY,2))
  CALL setderiv(OpenAD_prp_15,DIAGS(1:NZ*NX*NY,7))
  CALL setderiv(OpenAD_prp_16,DIAGS(1:NZ*NX*NY,1))
  CALL setderiv(DIAGS(1:NZ*NX*NY,4),DIAG_TMP(1:NZ*NX*NY))
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_16)
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_15)
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_14)
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_13)
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_12)
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_11)
  CALL SPDIAGS_FVM_CSR(NX,NY,NZ,DIAGS,ANNZ,AROW_INDEX,AROW_COMPRESSED,ACOL_INDEX&
     &,AVALUES)

  END SUBROUTINE

  SUBROUTINE TPFA(NX, NY, NZ, K, Q, P, V)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  type(active) :: K(1:3,1:NX,1:NY,1:NZ)
  type(active) :: Q(1:INT((NZ*NX*NY)))
  type(active) :: P(1:NX,1:NY,1:NZ)
  type(active) :: V(1:3,1:INT((NX+1)),1:INT((NY+1)),1:INT((NZ+1)))
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) ACOL_INDEX(1 : INT((NZ * NX * NY * 7)))
  INTEGER(w2f__i4) ANNZ
  INTEGER(w2f__i4) AROW_COMPRESSED(1 : INT((NZ * NX * NY + 1)))
  INTEGER(w2f__i4) AROW_INDEX(1 : INT((NZ * NX * NY * 7)))
  type(active) :: AVALUES(1:INT((NZ*NX*NY*7)))
  type(active) :: DIAGS(1:INT((NZ*NX*NY)),1:7)
  INTEGER(w2f__i4) H
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) J
  type(active) :: L(1:3,1:NX,1:NY,1:NZ)
  INTEGER(w2f__i4) M
  INTEGER(w2f__i4) N
  type(active) :: TX(1:INT((NX+1)),1:NY,1:NZ)
  type(active) :: TXYZ(1:NX,1:NY,1:NZ)
  REAL(w2f__8) TX_
  type(active) :: TY(1:NX,1:INT((NY+1)),1:NZ)
  REAL(w2f__8) TY_
  type(active) :: TZ(1:NX,1:NY,1:INT((NZ+1)))
  REAL(w2f__8) TZ_
  type(active) :: U(1:INT((NZ*NX*NY)))
  REAL(w2f__8) OpenAD_aux_18
  REAL(w2f__8) OpenAD_aux_19
  REAL(w2f__8) OpenAD_aux_20
  REAL(w2f__8) OpenAD_aux_21
  REAL(w2f__8) OpenAD_aux_22
  REAL(w2f__8) OpenAD_aux_23
  REAL(w2f__8) OpenAD_lin_28
  REAL(w2f__8) OpenAD_lin_29
  REAL(w2f__8) OpenAD_lin_30
  REAL(w2f__8) OpenAD_lin_31
  REAL(w2f__8) OpenAD_lin_32
  REAL(w2f__8) OpenAD_lin_33
  REAL(w2f__8) OpenAD_lin_34
  REAL(w2f__8) OpenAD_lin_35
  REAL(w2f__8) OpenAD_lin_36
  REAL(w2f__8) OpenAD_lin_37
  type(active) :: OpenAD_prp_20
  type(active) :: OpenAD_prp_21
  type(active) :: OpenAD_prp_22
  type(active) :: OpenAD_prp_23(:)
  ALLOCATABLE OpenAD_prp_23
  type(active) :: OpenAD_prp_24(:)
  ALLOCATABLE OpenAD_prp_24
  type(active) :: OpenAD_prp_25(:)
  ALLOCATABLE OpenAD_prp_25
  type(active) :: OpenAD_prp_26(:)
  ALLOCATABLE OpenAD_prp_26
  type(active) :: OpenAD_prp_27(:)
  ALLOCATABLE OpenAD_prp_27
  type(active) :: OpenAD_prp_28(:)
  ALLOCATABLE OpenAD_prp_28
  type(active) :: OpenAD_prp_29
  type(active) :: OpenAD_prp_30
  type(active) :: OpenAD_prp_31
!
!       **** Statements ****
!
  N = (NZ * NX * NY)
  DO H = 1, 3, 1
    DO I = 1, NX, 1
      DO J = 1, NY, 1
        DO M = 1, NZ, 1
          OpenAD_lin_28 = (-(1.0D00/(K(H,I,J,M)%v*K(H,I,J,M)%v)))
          L(INT(H),INT(I),INT(J),INT(M))%v = (1.0D00/K(H,I,J,M)%v)
          CALL sax(OpenAD_lin_28,K(H,I,J,M),L(H,I,J,M))
        END DO
      END DO
    END DO
  END DO
  TX_ = ((HZ*HY*2.0D00)/HX)
  TY_ = ((HZ*HX*2.0D00)/HY)
  TZ_ = ((HX*HY*2.0D00)/HZ)
  TX(1:NX+1,1:INT(NY),1:INT(NZ))%v = 0.0D00
  TY(1:INT(NX),1:NY+1,1:INT(NZ))%v = 0.0D00
  TZ(1:INT(NX),1:INT(NY),1:NZ+1)%v = 0.0D00
  CALL zero_deriv(TX(1:NX+1,1:INT(NY),1:INT(NZ)))
  CALL zero_deriv(TY(1:INT(NX),1:NY+1,1:INT(NZ)))
  CALL zero_deriv(TZ(1:INT(NX),1:INT(NY),1:NZ+1))
  DO I = 2,NX,1
    DO J = 1,NY,1
      DO M = 1,NZ,1
        OpenAD_aux_18 = (L(1,I,J,M)%v+L(1,I+(-1),J,M)%v)
        OpenAD_lin_29 = (-(TX_/(OpenAD_aux_18*OpenAD_aux_18)))
        TX(INT(I),INT(J),INT(M))%v = (TX_/OpenAD_aux_18)
        CALL setderiv(OpenAD_prp_20,L(1,I,J,M))
        CALL inc_deriv(OpenAD_prp_20,L(1,I+(-1),J,M))
        CALL sax(OpenAD_lin_29,OpenAD_prp_20,TX(I,J,M))
      END DO
    END DO
  END DO
  DO I = 1,NX,1
    DO J = 2,NY,1
      DO M = 1,NZ,1
        OpenAD_aux_19 = (L(2,I,J,M)%v+L(2,I,J+(-1),M)%v)
        OpenAD_lin_30 = (-(TY_/(OpenAD_aux_19*OpenAD_aux_19)))
        TY(INT(I),INT(J),INT(M))%v = (TY_/OpenAD_aux_19)
        CALL setderiv(OpenAD_prp_21,L(2,I,J,M))
        CALL inc_deriv(OpenAD_prp_21,L(2,I,J+(-1),M))
        CALL sax(OpenAD_lin_30,OpenAD_prp_21,TY(I,J,M))
      END DO
    END DO
  END DO
  DO I = 1,NX,1
    DO J = 1,NY,1
      DO M = 2,NZ,1
        OpenAD_aux_20 = (L(3,I,J,M)%v+L(3,I,J,M+(-1))%v)
        OpenAD_lin_31 = (-(TZ_/(OpenAD_aux_20*OpenAD_aux_20)))
        TZ(INT(I),INT(J),INT(M))%v = (TZ_/OpenAD_aux_20)
        CALL setderiv(OpenAD_prp_22,L(3,I,J,M))
        CALL inc_deriv(OpenAD_prp_22,L(3,I,J,M+(-1)))
        CALL sax(OpenAD_lin_31,OpenAD_prp_22,TZ(I,J,M))
      END DO
    END DO
  END DO
  DIAGS(1:NZ*NX*NY,1:7)%v = 0.0D00
  CALL zero_deriv(DIAGS(1:NZ*NX*NY,1:7))
  DO I = 1,NX,1
    DO J = 1,NY,1
      DO M = 1,NZ,1
        TXYZ(INT(I),INT(J),INT(M))%v = (-TX(I,J,M)%v)
        CALL set_neg_deriv(TXYZ(I,J,M),TX(I,J,M))
      END DO
    END DO
  END DO
  CALL MYRESHAPE_3_1(TXYZ,DIAGS(1:NZ*NX*NY,5))
  DO I = 1,NX,1
    DO J = 1,NY,1
      DO M = 1,NZ,1
        TXYZ(INT(I),INT(J),INT(M))%v = (-TY(I,J,M)%v)
        CALL set_neg_deriv(TXYZ(I,J,M),TY(I,J,M))
      END DO
    END DO
  END DO
  CALL MYRESHAPE_3_1(TXYZ,DIAGS(1:NZ*NX*NY,6))
  DO I = 1,NX,1
    DO J = 1,NY,1
      DO M = 1,NZ,1
        TXYZ(INT(I),INT(J),INT(M))%v = (-TZ(I,J,M)%v)
        CALL set_neg_deriv(TXYZ(I,J,M),TZ(I,J,M))
      END DO
    END DO
  END DO
  CALL MYRESHAPE_3_1(TXYZ,DIAGS(1:NZ*NX*NY,7))
  DO I = 2,(NX+1),1
    DO J = 1,NY,1
      DO M = 1,NZ,1
        TXYZ(INT(I),INT(J),INT(M))%v = (-TX(I,J,M)%v)
        CALL set_neg_deriv(TXYZ(I,J,M),TX(I,J,M))
      END DO
    END DO
  END DO
  CALL MYRESHAPE_3_1(TXYZ,DIAGS(1:NZ*NX*NY,3))
  DO I = 1,NX,1
    DO J = 2,(NY+1),1
      DO M = 1,NZ,1
        TXYZ(INT(I),INT(J),INT(M))%v = (-TY(I,J,M)%v)
        CALL set_neg_deriv(TXYZ(I,J,M),TY(I,J,M))
      END DO
    END DO
  END DO
  CALL MYRESHAPE_3_1(TXYZ,DIAGS(1:NZ*NX*NY,2))
  DO I = 1,NX,1
    DO J = 1,NY,1
      DO M = 2,(NZ+1),1
        TXYZ(INT(I),INT(J),INT(M))%v = (-TZ(I,J,M)%v)
        CALL set_neg_deriv(TXYZ(I,J,M),TZ(I,J,M))
      END DO
    END DO
  END DO
  CALL MYRESHAPE_3_1(TXYZ,DIAGS(1:NZ*NX*NY,1))
  DIAGS(1:NZ*NX*NY,4)%v = (-(DIAGS(1:NZ*NX*NY,7)%v+DIAGS(1:NZ*NX*NY,6)%v+DIAGS(1&
     &:NZ*NX*NY,5)%v+DIAGS(1:NZ*NX*NY,3)%v+DIAGS(1:NZ*NX*NY,1)%v+DIAGS(1:NZ*NX*N&
     &Y,2)%v))

!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_23,DIAGS(1:NZ*NX*NY,7))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_24,DIAGS(1:NZ*NX*NY,6))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_25,DIAGS(1:NZ*NX*NY,5))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_26,DIAGS(1:NZ*NX*NY,3))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_27,DIAGS(1:NZ*NX*NY,1))
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_prp_28,DIAGS(1:NZ*NX*NY,2))
  CALL setderiv(OpenAD_prp_23,DIAGS(1:NZ*NX*NY,7))
  CALL setderiv(OpenAD_prp_24,DIAGS(1:NZ*NX*NY,6))
  CALL setderiv(OpenAD_prp_25,DIAGS(1:NZ*NX*NY,5))
  CALL setderiv(OpenAD_prp_26,DIAGS(1:NZ*NX*NY,3))
  CALL setderiv(OpenAD_prp_27,DIAGS(1:NZ*NX*NY,1))
  CALL setderiv(OpenAD_prp_28,DIAGS(1:NZ*NX*NY,2))
  CALL set_neg_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_23)
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_24)
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_25)
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_26)
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_27)
  CALL dec_deriv(DIAGS(1:NZ*NX*NY,4),OpenAD_prp_28)
  CALL SPDIAGS_FVM_CSR(NX,NY,NZ,DIAGS,ANNZ,AROW_INDEX,AROW_COMPRESSED,ACOL_INDEX&
     &,AVALUES)

  DO I = 1,ANNZ,1
    IF ((MOD(AROW_INDEX(I),NY).eq.1).AND.(AROW_INDEX(I).LT.(NX*NY))) THEN
      IF (ACOL_INDEX(I).eq.AROW_INDEX(I)) THEN
        AVALUES(INT(I))%v = 1
        CALL zero_deriv(AVALUES(INT(I)))
      ELSE
        AVALUES(INT(I))%v = 0
        CALL zero_deriv(AVALUES(INT(I)))
      ENDIF
    ENDIF
  END DO
  CALL SPARSE_SOLVE(N,ANNZ,AROW_INDEX,AROW_COMPRESSED,ACOL_INDEX,AVALUES,Q,U)
  CALL MYRESHAPE_1_3(U,P)
  DO I = 2,NX,1
    DO J = 1,NY,1
      DO M = 1,NZ,1
        OpenAD_aux_21 = (P(I+(-1),J,M)%v-P(I,J,M)%v)
        OpenAD_lin_32 = OpenAD_aux_21
        OpenAD_lin_33 = TX(I,J,M)%v
        V(1,INT(I),INT(J),INT(M))%v = (TX(I,J,M)%v*OpenAD_aux_21)
        CALL setderiv(OpenAD_prp_29,P(I+(-1),J,M))
        CALL dec_deriv(OpenAD_prp_29,P(I,J,M))
        CALL sax(OpenAD_lin_32,TX(I,J,M),V(1,I,J,M))
        CALL saxpy(OpenAD_lin_33,OpenAD_prp_29,V(1,I,J,M))
      END DO
    END DO
  END DO
  DO I = 1,NX,1
    DO J = 2,NY,1
      DO M = 1,NZ,1
        OpenAD_aux_22 = (P(I,J+(-1),M)%v-P(I,J,M)%v)
        OpenAD_lin_34 = OpenAD_aux_22
        OpenAD_lin_35 = TY(I,J,M)%v
        V(2,INT(I),INT(J),INT(M))%v = (TY(I,J,M)%v*OpenAD_aux_22)
        CALL setderiv(OpenAD_prp_30,P(I,J+(-1),M))
        CALL dec_deriv(OpenAD_prp_30,P(I,J,M))
        CALL sax(OpenAD_lin_34,TY(I,J,M),V(2,I,J,M))
        CALL saxpy(OpenAD_lin_35,OpenAD_prp_30,V(2,I,J,M))
      END DO
    END DO
  END DO
  DO I = 1,NX,1
    DO J = 1,NY,1
      DO M = 2,NZ,1
        OpenAD_aux_23 = (P(I,J,M+(-1))%v-P(I,J,M)%v)
        OpenAD_lin_36 = OpenAD_aux_23
        OpenAD_lin_37 = TZ(I,J,M)%v
        V(3,INT(I),INT(J),INT(M))%v = (TZ(I,J,M)%v*OpenAD_aux_23)
        CALL setderiv(OpenAD_prp_31,P(I,J,M+(-1)))
        CALL dec_deriv(OpenAD_prp_31,P(I,J,M))
        CALL sax(OpenAD_lin_36,TZ(I,J,M),V(3,I,J,M))
        CALL saxpy(OpenAD_lin_37,OpenAD_prp_31,V(3,I,J,M))
      END DO
    END DO
  END DO
  END SUBROUTINE

  SUBROUTINE SPDIAGS_FVM(NX, NY, NZ, IMATRIX, ONNZ, OROW_INDEX, OROW_COMPRESSED,&
     & OCOL_INDEX, OVALUES)

  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  REAL(w2f__8) IMATRIX(1 : INT((NZ * NX * NY)), 1 : 7)
  INTEGER(w2f__i4) ONNZ
  INTEGER(w2f__i4) OROW_INDEX(1 : INT((NZ * NX * NY * 7)))
  INTEGER(w2f__i4) OROW_COMPRESSED(1 : INT((NZ * NX * NY + 1)))
  INTEGER(w2f__i4) OCOL_INDEX(1 : INT((NZ * NX * NY * 7)))
  REAL(w2f__8) OVALUES(1 : INT((NZ * NX * NY * 7)))
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) COL
  INTEGER(w2f__i4) END_ROW_IMATRIX
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) IDIAGS(1 : 7)
  INTEGER(w2f__i4) J
  INTEGER(w2f__i4) ROW
  INTEGER(w2f__i4) START_ROW_IMATRIX
!
!       **** Statements ****
!
  IDIAGS(1) = (-(NX * NY))
  IDIAGS(2) = (- NX)
  IDIAGS(3) = (-1)
  IDIAGS(4) = 0
  IDIAGS(5) = 1
  IDIAGS(6) = NX
  IDIAGS(7) = (NX * NY)
  ONNZ = 0
  OROW_COMPRESSED(1 : INT(NZ * NX * NY + 1)) = 0
  DO I = 1, 7, 1
    IF(IDIAGS(I) .GT. 0) THEN
      START_ROW_IMATRIX = (IDIAGS(I) + 1)
      END_ROW_IMATRIX = (NZ * NX * NY)
    ELSE
      IF(IDIAGS(I) .LE. 0) THEN
        START_ROW_IMATRIX = 1
        END_ROW_IMATRIX = (IDIAGS(I) + NZ * NX * NY)
      ENDIF
    ENDIF
    CALL FIRSTELM(IDIAGS(I),ROW,COL)
    DO J = START_ROW_IMATRIX, END_ROW_IMATRIX, 1
      IF((COL .eq. ROW) .OR.(IMATRIX(J, I) .ne. 0.0D00)) THEN
        ONNZ = (ONNZ + 1)
        OVALUES(INT(ONNZ)) = IMATRIX(J, I)
        OROW_INDEX(ONNZ) = ROW
        OCOL_INDEX(ONNZ) = COL
      ENDIF
      ROW = (ROW + 1)
      COL = (COL + 1)
    END DO
  END DO
  END SUBROUTINE

  SUBROUTINE SPDIAGS_FVM_CSR(NX, NY, NZ, IMATRIX, ONNZ, OROW_INDEX, OROW_COMPRES&
     &SED, OCOL_INDEX, OVALUES)

  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  type(active) :: IMATRIX(1:INT((NZ*NX*NY)),1:7)
  INTEGER(w2f__i4) ONNZ
  INTEGER(w2f__i4) OROW_INDEX(1 : INT((NZ * NX * NY * 7)))
  INTEGER(w2f__i4) OROW_COMPRESSED(1 : INT((NZ * NX * NY + 1)))
  INTEGER(w2f__i4) OCOL_INDEX(1 : INT((NZ * NX * NY * 7)))
  type(active) :: OVALUES(1:INT((NZ*NX*NY*7)))
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) COL_DIAG(1 : 7)
  INTEGER(w2f__i4) END_ROW_IMATRIX(1 : 7)
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) IDIAGS(1 : 7)
  INTEGER(w2f__i4) J
  INTEGER(w2f__i4) ROWNNZ
  INTEGER(w2f__i4) ROW_DIAG(1 : 7)
  INTEGER(w2f__i4) START_ROW_IMATRIX(1 : 7)
!
!       **** Statements ****
!
  IDIAGS(1) = (-(NX * NY))
  IDIAGS(2) = (- NX)
  IDIAGS(3) = (-1)
  IDIAGS(4) = 0
  IDIAGS(5) = 1
  IDIAGS(6) = NX
  IDIAGS(7) = (NX * NY)
  ONNZ = 0
  OROW_COMPRESSED(1) = 1
  DO I = 1, 7, 1
    IF(IDIAGS(I) .GT. 0) THEN
      START_ROW_IMATRIX(I) = (IDIAGS(I) + 1)
      END_ROW_IMATRIX(I) = (NZ * NX * NY)
    ELSE
      IF(IDIAGS(I) .LE. 0) THEN
        START_ROW_IMATRIX(I) = 1
        END_ROW_IMATRIX(I) = (IDIAGS(I) + NZ * NX * NY)
      ENDIF
    ENDIF
    CALL FIRSTELM(IDIAGS(I),ROW_DIAG(I),COL_DIAG(I))
  END DO
  DO I = 1, (NZ * NX * NY), 1
    ROWNNZ = 0
    DO J = 1, 7, 1
      IF((END_ROW_IMATRIX(J) .GE. START_ROW_IMATRIX(J)) .AND.( ROW_DIAG(J) .LE. &
     &I)) THEN

        IF ((IDIAGS(J).eq.0).OR.(IMATRIX(START_ROW_IMATRIX(J),J)%v.ne.0.0D00)) T&
     &HEN

          ONNZ = (ONNZ+1)
          ROWNNZ = (ROWNNZ+1)
          OVALUES(INT(ONNZ))%v = IMATRIX(START_ROW_IMATRIX(J),J)%v
          CALL setderiv(OVALUES(ONNZ),IMATRIX(START_ROW_IMATRIX(J),J))
          OROW_INDEX(ONNZ) = I
          OCOL_INDEX(ONNZ) = COL_DIAG(J)
        ENDIF
        START_ROW_IMATRIX(J) = (START_ROW_IMATRIX(J)+1)
        COL_DIAG(J) = (COL_DIAG(J)+1)
      ENDIF
    END DO
    OROW_COMPRESSED(I+1) = (OROW_COMPRESSED(I)+ROWNNZ)
  END DO
  END SUBROUTINE
END

MODULE simulation
use OAD_active
use w2f__types
use oad_intrinsics
use parameters
use matrix
use finitevolume
IMPLICIT NONE
SAVE
!
!     **** Statements ****
!
CONTAINS

  SUBROUTINE INIT_FLW_TRNC_NORM_XIN_PT_OUT(NX, NY, NZ, MU, SIGMA, Q)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  type(active) :: MU
  type(active) :: SIGMA
  type(active) :: Q(1:INT((NZ*NX*NY)))
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  REAL(w2f__8) IDX(1 : NX)
  INTEGER(w2f__i4) J
  type(active) :: MASS
  type(active) :: PDF
  REAL(w2f__8) PI
  type(active) :: Q_X(1:NX)
  REAL(w2f__8) X
  REAL(w2f__8) OpenAD_acc_10
  REAL(w2f__8) OpenAD_acc_11
  REAL(w2f__8) OpenAD_acc_12
  REAL(w2f__8) OpenAD_acc_13
  REAL(w2f__8) OpenAD_acc_14
  REAL(w2f__8) OpenAD_acc_9
  REAL(w2f__8) OpenAD_aux_29
  REAL(w2f__8) OpenAD_aux_30
  REAL(w2f__8) OpenAD_aux_31
  REAL(w2f__8) OpenAD_aux_32
  REAL(w2f__8) OpenAD_aux_33
  REAL(w2f__8) OpenAD_aux_34
  REAL(w2f__8) OpenAD_aux_35
  REAL(w2f__8) OpenAD_aux_36
  REAL(w2f__8) OpenAD_lin_43
  REAL(w2f__8) OpenAD_lin_44
  REAL(w2f__8) OpenAD_lin_45
  REAL(w2f__8) OpenAD_lin_46
  REAL(w2f__8) OpenAD_lin_47
  REAL(w2f__8) OpenAD_lin_48
  REAL(w2f__8) OpenAD_lin_49
  REAL(w2f__8) OpenAD_lin_50
  REAL(w2f__8) OpenAD_lin_51
  REAL(w2f__8) OpenAD_lin_52
  REAL(w2f__8) OpenAD_lin_53
  type(active) :: OpenAD_prp_32
  type(active) :: OpenAD_prp_33
!
!       **** Statements ****
!
  PI = 3.141592653589793116D00
  MASS%v = 0.0D00
  Q_X(1:INT(NX))%v = 0.0D00
  CALL zero_deriv(MASS)
  CALL zero_deriv(Q_X(1:INT(NX)))
  DO I = 1,NX,1
    X = ((((I+(-1))*3.0D00)/(NX+(-1)))+(-1.5D00))
    OpenAD_aux_33 = (X-MU%v)
    OpenAD_aux_32 = (OpenAD_aux_33/SIGMA%v)
    OpenAD_aux_31 = (OpenAD_aux_32**2.0D00)
    OpenAD_aux_29 = EXP(OpenAD_aux_31*(-5.0D-01))
    OpenAD_aux_35 = SQRT(PI*2.0D00)
    OpenAD_aux_34 = (SIGMA%v*OpenAD_aux_35)
    OpenAD_aux_30 = (1.0D00/OpenAD_aux_34)
    OpenAD_lin_47 = (INT(1_w2f__i8)/SIGMA%v)
    OpenAD_lin_48 = (-(OpenAD_aux_33/(SIGMA%v*SIGMA%v)))
    OpenAD_lin_46 = (2.0D00*(OpenAD_aux_32**(2.0D00-INT(1_w2f__i8))))
    OpenAD_lin_45 = OpenAD_aux_29
    OpenAD_lin_43 = OpenAD_aux_30
    OpenAD_lin_50 = OpenAD_aux_35
    OpenAD_lin_49 = (-(1.0D00/(OpenAD_aux_34*OpenAD_aux_34)))
    OpenAD_lin_44 = OpenAD_aux_29
    PDF%v = (OpenAD_aux_29*OpenAD_aux_30)
    Q_X(INT(I))%v = PDF%v
    MASS%v = (MASS%v+PDF%v)
    IDX(INT(I)) = I
    OpenAD_acc_9 = (OpenAD_lin_46*(-5.0D-01)*OpenAD_lin_45*OpenAD_lin_43)
    OpenAD_acc_10 = (OpenAD_lin_50*OpenAD_lin_49*OpenAD_lin_44)
    OpenAD_acc_11 = (OpenAD_lin_48*OpenAD_acc_9)
    OpenAD_acc_12 = (INT((-1_w2f__i8))*OpenAD_lin_47*OpenAD_acc_9)
    CALL setderiv(OpenAD_prp_32,MASS)
    CALL sax(OpenAD_acc_10,SIGMA,PDF)
    CALL saxpy(OpenAD_acc_11,SIGMA,PDF)
    CALL saxpy(OpenAD_acc_12,MU,PDF)
    CALL setderiv(Q_X(I),PDF)
    CALL setderiv(MASS,OpenAD_prp_32)
    CALL inc_deriv(MASS,PDF)
  END DO
  DO I = 1,NX,1
    OpenAD_aux_36 = (Q_X(I)%v/MASS%v)
    OpenAD_lin_52 = (INT(1_w2f__i8)/MASS%v)
    OpenAD_lin_53 = (-(Q_X(I)%v/(MASS%v*MASS%v)))
    OpenAD_lin_51 = IR
    Q_X(INT(I))%v = (IR*OpenAD_aux_36)
    OpenAD_acc_13 = (OpenAD_lin_52*OpenAD_lin_51)
    OpenAD_acc_14 = (OpenAD_lin_53*OpenAD_lin_51)
    CALL setderiv(OpenAD_prp_33,Q_X(I))
    CALL sax(OpenAD_acc_13,OpenAD_prp_33,Q_X(I))
    CALL saxpy(OpenAD_acc_14,MASS,Q_X(I))
  END DO
  J = 1
  DO I = 1,(NX*NY),NY
    Q(INT(I))%v = Q_X(J)%v
    CALL setderiv(Q(I),Q_X(J))
    J = (J+1)
  END DO
  Q(NZ*NX*NY)%v = (-IR)
  CALL zero_deriv(Q(NZ*NX*NY))
  END SUBROUTINE

  SUBROUTINE SIMULATE_RESERVOIR(NX, NY, NZ, ND, PT, ST, Q, S, P, V, TT, PC, OIL)
  use w2f__types
  use parameters
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  INTEGER(w2f__i4) ND
  INTEGER(w2f__i4) PT
  INTEGER(w2f__i4) ST
  type(active) :: Q(1:INT((NZ*NX*NY)))
  type(active) :: S(1:INT((NZ*NX*NY)))
  type(active) :: P(1:NX,1:NY,1:NZ)
  type(active) :: V(1:3,1:INT((NX+1)),1:INT((NY+1)),1:INT((NZ+1)))
  REAL(w2f__8) TT(1 : INT(((ND / ST) + 1)))
  type(active) :: PC(1:2,1:INT(((ND/ST)+1)))
  type(active) :: OIL
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) J
  INTEGER(w2f__i4) K
  type(active) :: MO
  type(active) :: MT
  type(active) :: MW
  type(active) :: TEMPOIL1
  type(active) :: TEMPOIL2
  REAL(w2f__8) OpenAD_lin_54
  REAL(w2f__8) OpenAD_lin_55
  REAL(w2f__8) OpenAD_lin_56
  REAL(w2f__8) OpenAD_lin_57
!
!       **** Statements ****
!
  S(1:NZ*NX*NY)%v = SWC
  PC(1,1)%v = 0.0D00
  PC(2,1)%v = 1.0D00
  TT(1) = 0.0D00
  TEMPOIL1%v = 0.0D00
  TEMPOIL2%v = 0.0D00
  K = 1
  CALL zero_deriv(S(1:NZ*NX*NY))
  CALL zero_deriv(PC(1,1))
  CALL zero_deriv(PC(2,1))
  CALL zero_deriv(TEMPOIL1)
  CALL zero_deriv(TEMPOIL2)
  DO I = 1,(ND/PT),1
    DO J = 1,(PT/ST),1
      K = (K+1)
      IF (J.eq.1) THEN
        CALL STEPFORWARD(NX,NY,NZ,ND,PT,ST,1,Q,S,P,V,MW,MO)
      ELSE
        CALL STEPFORWARD(NX,NY,NZ,ND,PT,ST,0,Q,S,P,V,MW,MO)
      ENDIF
      MT%v = (MO%v+MW%v)
      TT(INT(K)) = (ST*K)
      OpenAD_lin_54 = (INT(1_w2f__i8)/MT%v)
      OpenAD_lin_55 = (-(MW%v/(MT%v*MT%v)))
      PC(1,INT(K))%v = (MW%v/MT%v)
      CALL setderiv(MT,MO)
      CALL inc_deriv(MT,MW)
      CALL sax(OpenAD_lin_54,MW,PC(1,K))
      CALL saxpy(OpenAD_lin_55,MT,PC(1,K))
      OpenAD_lin_56 = (INT(1_w2f__i8)/MT%v)
      OpenAD_lin_57 = (-(MO%v/(MT%v*MT%v)))
      PC(2,INT(K))%v = (MO%v/MT%v)
      CALL sax(OpenAD_lin_56,MO,PC(2,K))
      CALL saxpy(OpenAD_lin_57,MT,PC(2,K))
      CALL UPDATE_OIL(ND,PT,ST,PC,K,TEMPOIL1,TEMPOIL2)
      TEMPOIL1%v = TEMPOIL2%v
      CALL setderiv(TEMPOIL1,TEMPOIL2)
    END DO
  END DO
  OIL%v = TEMPOIL2%v
  CALL setderiv(OIL,TEMPOIL2)
  END SUBROUTINE

  SUBROUTINE STEPFORWARD(NX, NY, NZ, ND, PT, ST, PRESSURE_STEP, Q , S, P, V, MW,&
     & MO)

  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  INTEGER(w2f__i4) ND
  INTEGER(w2f__i4) PT
  INTEGER(w2f__i4) ST
  INTEGER(w2f__i4) PRESSURE_STEP
  type(active) :: Q(1:INT((NZ*NX*NY)))
  type(active) :: S(1:INT((NZ*NX*NY)))
  type(active) :: P(1:NX,1:NY,1:NZ)
  type(active) :: V(1:3,1:INT((NX+1)),1:INT((NY+1)),1:INT((NZ+1)))
  type(active) :: MW
  type(active) :: MO
!
!       **** Statements ****
!
  IF(PRESSURE_STEP .eq. 1) THEN
    CALL PRES(NX,NY,NZ,Q,S,P,V)
  ENDIF
  CALL NEWTRAPH(NX,NY,NZ,ND,PT,ST,Q,V,S)
  CALL RELPERM_SCALAR(S(NZ*NX*NY),MW,MO)
  END SUBROUTINE

  SUBROUTINE UPDATE_OIL(ND, PT, ST, PC, K, OILIN, OILOUT)
  use w2f__types
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) ND
  INTEGER(w2f__i4) PT
  INTEGER(w2f__i4) ST
  type(active) :: PC(1:2,1:INT(((ND/ST)+1)))
  INTEGER(w2f__i4) K
  type(active) :: OILIN
  type(active) :: OILOUT
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) OpenAD_lin_58
!
!       **** Statements ****
!
  OpenAD_lin_58 = ST
  OILOUT%v = (OILIN%v+ST*PC(2,K)%v)
  CALL setderiv(OILOUT,OILIN)
  CALL saxpy(OpenAD_lin_58,PC(2,K),OILOUT)
  END SUBROUTINE

  SUBROUTINE WRAPPER(NX, NY, NZ, ND, PT, ST, MU, SIGMA, Q, S, P, V, TT, PC, OIL)
  use w2f__types
  use parameters
  IMPLICIT NONE
!
!       **** Parameters and Result ****
!
  INTEGER(w2f__i4) NX
  INTEGER(w2f__i4) NY
  INTEGER(w2f__i4) NZ
  INTEGER(w2f__i4) ND
  INTEGER(w2f__i4) PT
  INTEGER(w2f__i4) ST
  type(active) :: MU
  type(active) :: SIGMA
  type(active) :: Q(1:INT((NZ*NX*NY)))
  type(active) :: S(1:INT((NZ*NX*NY)))
  REAL(w2f__8) P(1 : NX, 1 : NY, 1 : NZ)
  REAL(w2f__8) V(1 : 3, 1 : INT((NX + 1)), 1 : INT((NY + 1)), 1 : INT((NZ + 1)))
  REAL(w2f__8) TT(1 : INT(((ND / ST) + 1)))
  REAL(w2f__8) PC(1 : 2, 1 : INT(((ND / ST) + 1)))
  type(active) :: OIL
!
!       **** Local Variables and Functions ****
!
  type(active) :: OpenAD_tyc_2(:,:,:)
  ALLOCATABLE OpenAD_tyc_2
  type(active) :: OpenAD_tyc_3(:,:,:,:)
  ALLOCATABLE OpenAD_tyc_3
  type(active) :: OpenAD_tyc_4(:,:)
  ALLOCATABLE OpenAD_tyc_4
!
!       **** Top Level Pragmas ****
!
!$OPENAD INDEPENDENT(MU)
!$OPENAD INDEPENDENT(SIGMA)
!$OPENAD DEPENDENT(OIL)
!
!       **** Statements ****
!
  CALL INIT_FLW_TRNC_NORM_XIN_PT_OUT(NX,NY,NZ,MU,SIGMA,Q)
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_tyc_2,P)
!       $OpenAD$ INLINE oad_convert(subst,subst)
  CALL oad_convert(OpenAD_tyc_2,P)
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_tyc_3,V)
!       $OpenAD$ INLINE oad_convert(subst,subst)
  CALL oad_convert(OpenAD_tyc_3,V)
!       $OpenAD$ INLINE oad_AllocateMatching(subst,subst)
  CALL oad_AllocateMatching(OpenAD_tyc_4,PC)
!       $OpenAD$ INLINE oad_convert(subst,subst)
  CALL oad_convert(OpenAD_tyc_4,PC)
  CALL SIMULATE_RESERVOIR(NX,NY,NZ,ND,PT,ST,Q,S,OpenAD_tyc_2,OpenAD_tyc_3,TT,Ope&
     &nAD_tyc_4,OIL)

!       $OpenAD$ INLINE oad_ShapeTest(subst,subst)
  CALL oad_ShapeTest(OpenAD_tyc_2,P)
!       $OpenAD$ INLINE oad_convert(subst,subst)
  CALL oad_convert(P,OpenAD_tyc_2)
!       $OpenAD$ INLINE oad_ShapeTest(subst,subst)
  CALL oad_ShapeTest(OpenAD_tyc_3,V)
!       $OpenAD$ INLINE oad_convert(subst,subst)
  CALL oad_convert(V,OpenAD_tyc_3)
!       $OpenAD$ INLINE oad_ShapeTest(subst,subst)
  CALL oad_ShapeTest(OpenAD_tyc_4,PC)
!       $OpenAD$ INLINE oad_convert(subst,subst)
  CALL oad_convert(PC,OpenAD_tyc_4)
  END SUBROUTINE
END
