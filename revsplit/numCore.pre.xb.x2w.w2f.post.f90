
MODULE oad_intrinsics
use OAD_active
use w2f__types
IMPLICIT NONE
SAVE
!
!     **** Statements ****
!
END MODULE

MODULE stream_vel_variables
use OAD_active
use w2f__types
IMPLICIT NONE
SAVE
!
!     **** Global Variables & Derived Type Definitions ****
!
REAL(w2f__8) AX(1 : 79)
REAL(w2f__8) DX
REAL(w2f__8) P(1 : 79)
REAL(w2f__8) P_OLD(1 : 79)
REAL(w2f__8) R(1 : 79)
REAL(w2f__8) R_OLD(1 : 79)
REAL(w2f__8) TRIDIAG_0(1 : 79, 1 : 3)
REAL(w2f__8) X_OLD(1 : 79)
!
!     **** Local Variables and Functions ****
!
REAL(w2f__8) AGLEN
PARAMETER ( AGLEN = 5.00019993114791044342D-17)
REAL(w2f__8) BETA_CONST
PARAMETER ( BETA_CONST = 5.0D00)
REAL(w2f__8) EPS
PARAMETER ( EPS = 9.99999974737875163555D-06)
REAL(w2f__8) EP_GLEN
PARAMETER ( EP_GLEN = 1.00000001168609742308D-07)
REAL(w2f__8) G
PARAMETER ( G = 9.81000041961669921875D00)
REAL(w2f__8) H_LEFT
PARAMETER ( H_LEFT = 1.05D+03)
REAL(w2f__8) H_RIGHT
PARAMETER ( H_RIGHT = 1.05D+03)
REAL(w2f__8) LX
PARAMETER ( LX = 7.9D+04)
INTEGER(w2f__i4) N
PARAMETER ( N = 79)
REAL(w2f__8) NGLEN
PARAMETER ( NGLEN = 3.0D00)
INTEGER(w2f__i4) N_NL
PARAMETER ( N_NL = 20)
REAL(w2f__8) PI
PARAMETER ( PI = 3.141592653589793116D00)
REAL(w2f__8) RHOI
PARAMETER ( RHOI = 9.1D+02)
REAL(w2f__8) RHOW
PARAMETER ( RHOW = 1.035D+03)
REAL(w2f__8) R_BED
PARAMETER ( R_BED = -9.0D+02)
!
!     **** Statements ****
!
END MODULE

MODULE conj_gradstub_mod
use OAD_active
use w2f__types
use oad_intrinsics
IMPLICIT NONE
SAVE
!
!     **** Top Level Pragmas ****
!
interface SOLVE
  module procedure CONJ_GRADSTUB
end interface

!
!     **** Statements ****
!
CONTAINS

  SUBROUTINE CONJ_GRADSTUB(X, B, A)
      use OAD_tape
      use OAD_rev
      use OAD_cp
      use w2f__types
!use OAD_cp
!     !use OAD_tape
!     !use OAD_rev
      use conj_grad_mod
      use conj_grad_ad_mod
      use stream_vel_variables
   use w2f__types
  use stream_vel_variables
  use stream_vel_variables
  use stream_vel_variables
  IMPLICIT NONE
!
!       **** Global Variables & Derived Type Definitions ****
!
  INTEGER(w2f__i8) OpenAD_Symbol_0
  INTEGER(w2f__i8) OpenAD_Symbol_1
  INTEGER(w2f__i8) OpenAD_Symbol_2
  INTEGER(w2f__i8) OpenAD_Symbol_3
  INTEGER(w2f__i8) OpenAD_Symbol_4
  INTEGER(w2f__i8) OpenAD_Symbol_5
!
!       **** Parameters and Result ****
!
  type(active) :: X(1:79)
  type(active) :: B(1:79)
  type(active) :: A(1:79,1:3)
!
!       **** Local Variables and Functions ****
!
  INTEGER(w2f__i4) I
  INTEGER(w2f__i4) OpenAD_Symbol_99


      type(modeType) :: our_orig_mode

      real(8), dimension(n) :: x_p
      real(8), dimension(n) :: b_p
      real(8), dimension(n,3) :: A_p
      real(8), dimension(n) :: x_d
      real(8), dimension(n) :: b_d
      real(8), dimension(n,3) :: A_d
      integer :: j
!
!       **** Statements ****
!

      if (our_rev_mode%plain) then
! set up for plain execution
         our_orig_mode=our_rev_mode
         our_rev_mode%arg_store=.FALSE.
         our_rev_mode%arg_restore=.FALSE.
         our_rev_mode%plain=.TRUE.
         our_rev_mode%tape=.FALSE.
         our_rev_mode%adjoint=.FALSE.
         b_p = b%v
         x_p = x%v
         A_p = A%v
         call solve(x_p,b_p,A_p)
         x%v = x_p
! reset the mode
         our_rev_mode=our_orig_mode
      end if
      if (our_rev_mode%tape) then
! set up for plain execution
         our_orig_mode=our_rev_mode
         our_rev_mode%arg_store=.FALSE.
         our_rev_mode%arg_restore=.FALSE.
         our_rev_mode%plain=.TRUE.
         our_rev_mode%tape=.FALSE.
         our_rev_mode%adjoint=.FALSE.
         b_p = b%v
         x_p = x%v
         A_p = A%v
!call oad_tape_push(b_p)      
!        !call oad_tape_push(A_p)
!        !call push_s1(b_p)
!        !call push_s2(A_p) 
!        !WRITE(*,*) 'Size(b,1) = ', size(b,1)
!        !do i=1,n
!        !   WRITE(*,*) 'x(', i, ') = ', b_p(i)
!        !end do

         do i=1,n
            call push_s0(b_p(i))
         end do
         do i=1,n
           do j=1,3
            call push_s0(A_p(i,j))
           end do
         end do
         call solve(x_p,b_p,A_p)
!call oad_tape_push(x_p)     
!        !call push_s1(x_p)
!        !WRITE(*,*) 'Size(x,1) = ', size(x,1)

         do i=1,n
            call push_s0(x_p(i))
         end do
         x%v = x_p
! reset the mode
         our_rev_mode=our_orig_mode
      end if
      if (our_rev_mode%adjoint) then
!call oad_tape_pop(x_p)     
!        !call oad_tape_pop(A_p)     
!        !call oad_tape_pop(b_p)
!        !call pop_s1(x_p)     
!        !call pop_s2(A_p)     
!        !call pop_s1(b_p)
         do i=n, 1, -1
            call pop_s0(x_p(i))
         end do
         do i=n, 1, -1
           do j=3,1, -1
            call pop_s0(A_p(i,j))
           end do
         end do
         do i=n,1,-1
            call pop_s0(b_p(i))
         end do
! set up for plain execution
         our_orig_mode=our_rev_mode
         our_rev_mode%arg_store=.FALSE.
         our_rev_mode%arg_restore=.FALSE.
         our_rev_mode%plain=.TRUE.
         our_rev_mode%tape=.FALSE.
         our_rev_mode%adjoint=.FALSE.
         b_d = b%d
         x_d = x%d
         A_d = A%d
      call adsolve( x_p, x_d, b_p, b_d, a_p, a_d )
! reset the mode
         b%d = b_d
         x%d = x_d
         A%d = A_d
         our_rev_mode=our_orig_mode
      end if
      end subroutine CONJ_GRADSTUB
END
!#########################################################
! This file is part of OpenAD released under the LGPL.   #
! The full COPYRIGHT notice can be found in the top      #
! level directory of the OpenAD distribution             #
!#########################################################

      SUBROUTINE stream_vel_init(H, BETA)
    use OAD_tape
    use OAD_rev

! original arguments get inserted before version
!         ! and declared here together with all local variables
!         ! generated by xaifBooster

      use OAD_active
      use w2f__types
      use oad_intrinsics
      use stream_vel_variables
      use oad_intrinsics
      use stream_vel_variables
      use oad_intrinsics
      use stream_vel_variables
      IMPLICIT NONE
!
!     **** Global Variables & Derived Type Definitions ****
!
      INTEGER(w2f__i8) OpenAD_Symbol_10
      INTEGER(w2f__i8) OpenAD_Symbol_11
      INTEGER(w2f__i8) OpenAD_Symbol_6
      INTEGER(w2f__i8) OpenAD_Symbol_7
      INTEGER(w2f__i8) OpenAD_Symbol_8
      INTEGER(w2f__i8) OpenAD_Symbol_9
!
!     **** Parameters and Result ****
!
      REAL(w2f__8) H(1 : 79)
      REAL(w2f__8) BETA(1 : 79)
!
!     **** Local Variables and Functions ****
!
      INTEGER(w2f__i4) I


    integer iaddr
    external iaddr
!
!     **** Statements ****
!

   if (our_rev_mode%plain) then
! original function
      DX = (7.9D+04 / REAL(79))
      DO I = 1, 79, 1
        BETA(INT(I)) = 5.0D00
        H(INT(I)) = (DX *(I + DBLE((-5.0E-01))) * 0.0D00 + 1.05D+03)
      END DO
    end if
    if (our_rev_mode%tape) then
! taping
      DX = (7.9D+04 / REAL(79))
      OpenAD_Symbol_8 = 0_w2f__i8
      DO I = 1, 79, 1
        BETA(INT(I)) = 5.0D00
        H(INT(I)) = (DX *(I + DBLE((-5.0E-01))) * 0.0D00 + 1.05D+03)
        OpenAD_Symbol_8 = (INT(OpenAD_Symbol_8) + INT(1_w2f__i8))
      END DO
      integer_tape(integer_tape_pointer) = OpenAD_Symbol_8
      integer_tape_pointer = integer_tape_pointer+1
    end if
    if (our_rev_mode%adjoint) then
! adjoint
      integer_tape_pointer = integer_tape_pointer-1
      OpenAD_Symbol_6 = integer_tape(integer_tape_pointer)
      OpenAD_Symbol_7 = 1
      DO WHILE(INT(OpenAD_Symbol_7) .LE. INT(OpenAD_Symbol_6))
        OpenAD_Symbol_7 = INT(OpenAD_Symbol_7) + 1
      END DO
    end if
  end subroutine stream_vel_init
!#########################################################
! This file is part of OpenAD released under the LGPL.   #
! The full COPYRIGHT notice can be found in the top      #
! level directory of the OpenAD distribution             #
!#########################################################

SUBROUTINE stream_vel(U, BB, FC)
    use OAD_tape
    use OAD_rev

! original arguments get inserted before version
!         ! and declared here together with all local variables
!         ! generated by xaifBooster

use OAD_active
use w2f__types
use oad_intrinsics
use stream_vel_variables
use conj_gradstub_mod
use oad_intrinsics
use stream_vel_variables
use conj_gradstub_mod
use oad_intrinsics
use stream_vel_variables
use conj_gradstub_mod
IMPLICIT NONE
!
!     **** Global Variables & Derived Type Definitions ****
!
INTEGER(w2f__i8) OpenAD_Symbol_12
INTEGER(w2f__i8) OpenAD_Symbol_13
INTEGER(w2f__i8) OpenAD_Symbol_14
INTEGER(w2f__i8) OpenAD_Symbol_15
INTEGER(w2f__i8) OpenAD_Symbol_16
INTEGER(w2f__i8) OpenAD_Symbol_17
INTEGER(w2f__i8) OpenAD_Symbol_18
INTEGER(w2f__i8) OpenAD_Symbol_19
INTEGER(w2f__i8) OpenAD_Symbol_20
INTEGER(w2f__i8) OpenAD_Symbol_21
INTEGER(w2f__i8) OpenAD_Symbol_22
INTEGER(w2f__i8) OpenAD_Symbol_23
INTEGER(w2f__i8) OpenAD_Symbol_24
INTEGER(w2f__i8) OpenAD_Symbol_25
INTEGER(w2f__i8) OpenAD_Symbol_26
INTEGER(w2f__i8) OpenAD_Symbol_27
INTEGER(w2f__i8) OpenAD_Symbol_28
INTEGER(w2f__i8) OpenAD_Symbol_29
INTEGER(w2f__i8) OpenAD_Symbol_30
INTEGER(w2f__i8) OpenAD_Symbol_31
INTEGER(w2f__i8) OpenAD_Symbol_32
INTEGER(w2f__i8) OpenAD_Symbol_33
INTEGER(w2f__i8) OpenAD_Symbol_34
INTEGER(w2f__i8) OpenAD_Symbol_35
INTEGER(w2f__i8) OpenAD_Symbol_36
INTEGER(w2f__i8) OpenAD_Symbol_37
INTEGER(w2f__i8) OpenAD_Symbol_38
INTEGER(w2f__i8) OpenAD_Symbol_39
INTEGER(w2f__i8) OpenAD_Symbol_40
INTEGER(w2f__i8) OpenAD_Symbol_41
!
!     **** Parameters and Result ****
!
type(active) :: U(1:80)
type(active) :: BB(1:79)
type(active) :: FC
!
!     **** Local Variables and Functions ****
!
type(active) :: A(1:79,1:3)
type(active) :: B(1:79)
REAL(w2f__8) BETA_0(1 : 79)
REAL(w2f__8) BETA_FRIC(1 : 79)
type(active) :: F(1:79)
type(active) :: FEND
type(active) :: H(1:79)
REAL(w2f__8) H0(1 : 79)
INTEGER(w2f__i4) I
INTEGER(w2f__i4) J
type(active) :: NU(1:79)
EXTERNAL stream_assemble
EXTERNAL stream_vel_init
EXTERNAL stream_vel_taud
EXTERNAL stream_vel_visc
type(active) :: UTMP(1:79)
INTEGER(w2f__i4) OpenAD_Symbol_100
REAL(w2f__8) OpenAD_Symbol_101
INTEGER(w2f__i4) OpenAD_Symbol_102
REAL(w2f__8) OpenAD_Symbol_103
INTEGER(w2f__i4) OpenAD_Symbol_104
INTEGER(w2f__i4) OpenAD_Symbol_105
INTEGER(w2f__i4) OpenAD_Symbol_106
INTEGER(w2f__i4) OpenAD_Symbol_107
REAL(w2f__8) OpenAD_Symbol_108
REAL(w2f__8) OpenAD_Symbol_109
INTEGER(w2f__i4) OpenAD_Symbol_42
INTEGER(w2f__i4) OpenAD_Symbol_43
REAL(w2f__8) OpenAD_acc_0
REAL(w2f__8) OpenAD_acc_1
REAL(w2f__8) OpenAD_lin_0
REAL(w2f__8) OpenAD_lin_1
REAL(w2f__8) OpenAD_lin_2
REAL(w2f__8) OpenAD_lin_3
type(active) :: OpenAD_prp_0
type(active) :: OpenAD_prp_1
type(active) :: OpenAD_prp_2


    integer iaddr
    external iaddr
!
!     **** Top Level Pragmas ****
!
!$OPENAD INDEPENDENT(U)
!$OPENAD INDEPENDENT(BB)
!$OPENAD DEPENDENT(FC)
!
!     **** Statements ****
!

   if (our_rev_mode%plain) then
! original function
CALL stream_vel_init(H0,BETA_0)
BETA_FRIC(1 : 79) = BETA_0(1 : 79)
H(1:79)%v = (BB(1:79)%v+H0(1:79))
CALL stream_vel_taud(H,F,FEND)
U(1:80)%v = 0.0
DO I = 1,79,1
  B(INT(I))%v = (-(F(I)%v*DX))
  IF (I.LT.79) THEN
    B(INT(I))%v = (B(I)%v-F(I+1)%v*DX)
  ENDIF
END DO
B(79)%v = (B(79)%v+FEND%v)
DO I = 1,20,1
  CALL stream_vel_visc(H,U,NU)
  CALL stream_assemble(NU,BETA_FRIC,A)
  UTMP(1:79)%v = 0.0
  CALL CONJ_GRADSTUB(UTMP,B,A)
  DO J = 1,79,1
    U(J+1)%v = UTMP(J)%v
  END DO
END DO
FC%v = 0.0D00
DO I = 2,80,1
  FC%v = (FC%v+U(I)%v*U(I)%v)
END DO
    end if
    if (our_rev_mode%tape) then
! taping
CALL stream_vel_init(H0,BETA_0)
BETA_FRIC(1:79) = BETA_0(1:79)
H(1:79)%v = (BB(1:79)%v+H0(1:79))
CALL stream_vel_taud(H,F,FEND)
U(1:80)%v = 0.0
OpenAD_Symbol_21 = 0_w2f__i8
DO I = 1,79,1
  OpenAD_lin_0 = DX
  B(INT(I))%v = (-(F(I)%v*DX))
  OpenAD_acc_0 = (OpenAD_lin_0*INT((-1_w2f__i8)))
  double_tape(double_tape_pointer) = OpenAD_acc_0
  double_tape_pointer = double_tape_pointer+1
  integer_tape(integer_tape_pointer) = I
  integer_tape_pointer = integer_tape_pointer+1
  IF (I.LT.79) THEN
    OpenAD_lin_1 = DX
    B(INT(I))%v = (B(I)%v-F(I+1)%v*DX)
    OpenAD_acc_1 = (OpenAD_lin_1*INT((-1_w2f__i8)))
    double_tape(double_tape_pointer) = OpenAD_acc_1
    double_tape_pointer = double_tape_pointer+1
    OpenAD_Symbol_42 = (I+1)
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_42
    integer_tape_pointer = integer_tape_pointer+1
    integer_tape(integer_tape_pointer) = I
    integer_tape_pointer = integer_tape_pointer+1
    OpenAD_Symbol_22 = 1_w2f__i8
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_22
    integer_tape_pointer = integer_tape_pointer+1
  ELSE
    OpenAD_Symbol_23 = 0_w2f__i8
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_23
    integer_tape_pointer = integer_tape_pointer+1
  ENDIF
  OpenAD_Symbol_21 = (INT(OpenAD_Symbol_21)+INT(1_w2f__i8))
END DO
integer_tape(integer_tape_pointer) = OpenAD_Symbol_21
integer_tape_pointer = integer_tape_pointer+1
B(79)%v = (B(79)%v+FEND%v)
OpenAD_Symbol_24 = 0_w2f__i8
DO I = 1,20,1
  CALL stream_vel_visc(H,U,NU)
  CALL stream_assemble(NU,BETA_FRIC,A)
  UTMP(1:79)%v = 0.0
  CALL CONJ_GRADSTUB(UTMP,B,A)
  OpenAD_Symbol_25 = 0_w2f__i8
  DO J = 1,79,1
    U(J+1)%v = UTMP(J)%v
    OpenAD_Symbol_43 = (J+1)
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_43
    integer_tape_pointer = integer_tape_pointer+1
    integer_tape(integer_tape_pointer) = J
    integer_tape_pointer = integer_tape_pointer+1
    OpenAD_Symbol_25 = (INT(OpenAD_Symbol_25)+INT(1_w2f__i8))
  END DO
  integer_tape(integer_tape_pointer) = OpenAD_Symbol_25
  integer_tape_pointer = integer_tape_pointer+1
  OpenAD_Symbol_24 = (INT(OpenAD_Symbol_24)+INT(1_w2f__i8))
END DO
integer_tape(integer_tape_pointer) = OpenAD_Symbol_24
integer_tape_pointer = integer_tape_pointer+1
FC%v = 0.0D00
OpenAD_Symbol_26 = 0_w2f__i8
DO I = 2,80,1
  OpenAD_lin_2 = U(I)%v
  OpenAD_lin_3 = U(I)%v
  FC%v = (FC%v+U(I)%v*U(I)%v)
  double_tape(double_tape_pointer) = OpenAD_lin_2
  double_tape_pointer = double_tape_pointer+1
  double_tape(double_tape_pointer) = OpenAD_lin_3
  double_tape_pointer = double_tape_pointer+1
  integer_tape(integer_tape_pointer) = I
  integer_tape_pointer = integer_tape_pointer+1
  OpenAD_Symbol_26 = (INT(OpenAD_Symbol_26)+INT(1_w2f__i8))
END DO
integer_tape(integer_tape_pointer) = OpenAD_Symbol_26
integer_tape_pointer = integer_tape_pointer+1
    end if
    if (our_rev_mode%adjoint) then
! adjoint
integer_tape_pointer = integer_tape_pointer-1
OpenAD_Symbol_12 = integer_tape(integer_tape_pointer)
OpenAD_Symbol_13 = 1
do while (INT(OpenAD_Symbol_13).LE.INT(OpenAD_Symbol_12))
  integer_tape_pointer = integer_tape_pointer-1
  OpenAD_Symbol_107 = integer_tape(integer_tape_pointer)
  double_tape_pointer = double_tape_pointer-1
  OpenAD_Symbol_108 = double_tape(double_tape_pointer)
  double_tape_pointer = double_tape_pointer-1
  OpenAD_Symbol_109 = double_tape(double_tape_pointer)
  U(OpenAD_Symbol_107)%d = U(OpenAD_Symbol_107)%d+FC%d*(OpenAD_Symbol_108)
  U(OpenAD_Symbol_107)%d = U(OpenAD_Symbol_107)%d+FC%d*(OpenAD_Symbol_109)
  OpenAD_prp_2%d = OpenAD_prp_2%d+FC%d
  FC%d = 0.0d0
  FC%d = FC%d+OpenAD_prp_2%d
  OpenAD_prp_2%d = 0.0d0
  OpenAD_Symbol_13 = INT(OpenAD_Symbol_13)+1
END DO
FC%d = 0.0d0
integer_tape_pointer = integer_tape_pointer-1
OpenAD_Symbol_14 = integer_tape(integer_tape_pointer)
OpenAD_Symbol_15 = 1
do while (INT(OpenAD_Symbol_15).LE.INT(OpenAD_Symbol_14))
  integer_tape_pointer = integer_tape_pointer-1
  OpenAD_Symbol_16 = integer_tape(integer_tape_pointer)
  OpenAD_Symbol_17 = 1
  do while (INT(OpenAD_Symbol_17).LE.INT(OpenAD_Symbol_16))
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_105 = integer_tape(integer_tape_pointer)
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_106 = integer_tape(integer_tape_pointer)
    UTMP(OpenAD_Symbol_105)%d = UTMP(OpenAD_Symbol_105)%d+U(OpenAD_Symbol_106)%d
    U(OpenAD_Symbol_106)%d = 0.0d0
    OpenAD_Symbol_17 = INT(OpenAD_Symbol_17)+1
  END DO
  CALL CONJ_GRADSTUB(UTMP,B,A)
  UTMP(1:79)%d = 0.0d0
  CALL stream_assemble(NU,BETA_FRIC,A)
  CALL stream_vel_visc(H,U,NU)
  OpenAD_Symbol_15 = INT(OpenAD_Symbol_15)+1
END DO
FEND%d = FEND%d+B(79)%d
OpenAD_prp_1%d = OpenAD_prp_1%d+B(79)%d
B(79)%d = 0.0d0
B(79)%d = B(79)%d+OpenAD_prp_1%d
OpenAD_prp_1%d = 0.0d0
integer_tape_pointer = integer_tape_pointer-1
OpenAD_Symbol_18 = integer_tape(integer_tape_pointer)
OpenAD_Symbol_19 = 1
do while (INT(OpenAD_Symbol_19).LE.INT(OpenAD_Symbol_18))
  integer_tape_pointer = integer_tape_pointer-1
  OpenAD_Symbol_20 = integer_tape(integer_tape_pointer)
  IF (OpenAD_Symbol_20.ne.0) THEN
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_102 = integer_tape(integer_tape_pointer)
    double_tape_pointer = double_tape_pointer-1
    OpenAD_Symbol_103 = double_tape(double_tape_pointer)
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_104 = integer_tape(integer_tape_pointer)
    F(OpenAD_Symbol_104)%d = F(OpenAD_Symbol_104)%d+B(OpenAD_Symbol_102)%d*(Open&
     &AD_Symbol_103)

    OpenAD_prp_0%d = OpenAD_prp_0%d+B(OpenAD_Symbol_102)%d
    B(OpenAD_Symbol_102)%d = 0.0d0
    B(OpenAD_Symbol_102)%d = B(OpenAD_Symbol_102)%d+OpenAD_prp_0%d
    OpenAD_prp_0%d = 0.0d0
  ENDIF
  integer_tape_pointer = integer_tape_pointer-1
  OpenAD_Symbol_100 = integer_tape(integer_tape_pointer)
  double_tape_pointer = double_tape_pointer-1
  OpenAD_Symbol_101 = double_tape(double_tape_pointer)
  F(OpenAD_Symbol_100)%d = F(OpenAD_Symbol_100)%d+B(OpenAD_Symbol_100)%d*(OpenAD&
     &_Symbol_101)

  B(OpenAD_Symbol_100)%d = 0.0d0
  OpenAD_Symbol_19 = INT(OpenAD_Symbol_19)+1
END DO
U(1:80)%d = 0.0d0
CALL stream_vel_taud(H,F,FEND)
BB(1:79)%d = BB(1:79)%d+H(1:79)%d
H(1:79)%d = 0.0d0
CALL stream_vel_init(H0,BETA_0)
    end if
  end subroutine stream_vel
!#########################################################
! This file is part of OpenAD released under the LGPL.   #
! The full COPYRIGHT notice can be found in the top      #
! level directory of the OpenAD distribution             #
!#########################################################

SUBROUTINE stream_vel_taud(H, F, FEND)
    use OAD_tape
    use OAD_rev

! original arguments get inserted before version
!         ! and declared here together with all local variables
!         ! generated by xaifBooster

use OAD_active
use w2f__types
use oad_intrinsics
use stream_vel_variables
use oad_intrinsics
use stream_vel_variables
use oad_intrinsics
use stream_vel_variables
IMPLICIT NONE
!
!     **** Global Variables & Derived Type Definitions ****
!
INTEGER(w2f__i8) OpenAD_Symbol_44
INTEGER(w2f__i8) OpenAD_Symbol_45
INTEGER(w2f__i8) OpenAD_Symbol_46
INTEGER(w2f__i8) OpenAD_Symbol_47
INTEGER(w2f__i8) OpenAD_Symbol_48
INTEGER(w2f__i8) OpenAD_Symbol_49
INTEGER(w2f__i8) OpenAD_Symbol_50
INTEGER(w2f__i8) OpenAD_Symbol_51
INTEGER(w2f__i8) OpenAD_Symbol_52
INTEGER(w2f__i8) OpenAD_Symbol_53
INTEGER(w2f__i8) OpenAD_Symbol_54
INTEGER(w2f__i8) OpenAD_Symbol_55
INTEGER(w2f__i8) OpenAD_Symbol_56
INTEGER(w2f__i8) OpenAD_Symbol_57
INTEGER(w2f__i8) OpenAD_Symbol_58
INTEGER(w2f__i8) OpenAD_Symbol_59
INTEGER(w2f__i8) OpenAD_Symbol_60
INTEGER(w2f__i8) OpenAD_Symbol_61
INTEGER(w2f__i8) OpenAD_Symbol_62
INTEGER(w2f__i8) OpenAD_Symbol_63
INTEGER(w2f__i8) OpenAD_Symbol_64
INTEGER(w2f__i8) OpenAD_Symbol_65
INTEGER(w2f__i8) OpenAD_Symbol_66
INTEGER(w2f__i8) OpenAD_Symbol_67
!
!     **** Parameters and Result ****
!
type(active) :: H(1:79)
type(active) :: F(1:79)
type(active) :: FEND
!
!     **** Local Variables and Functions ****
!
INTEGER(w2f__i4) I
INTEGER(w2f__i4) OpenAD_Symbol_110
REAL(w2f__8) OpenAD_Symbol_111
REAL(w2f__8) OpenAD_Symbol_112
INTEGER(w2f__i4) OpenAD_Symbol_113
INTEGER(w2f__i4) OpenAD_Symbol_114
INTEGER(w2f__i4) OpenAD_Symbol_115
REAL(w2f__8) OpenAD_Symbol_116
REAL(w2f__8) OpenAD_Symbol_117
INTEGER(w2f__i4) OpenAD_Symbol_118
INTEGER(w2f__i4) OpenAD_Symbol_119
REAL(w2f__8) OpenAD_Symbol_120
REAL(w2f__8) OpenAD_Symbol_121
INTEGER(w2f__i4) OpenAD_Symbol_122
REAL(w2f__8) OpenAD_Symbol_123
INTEGER(w2f__i4) OpenAD_Symbol_68
INTEGER(w2f__i4) OpenAD_Symbol_69
INTEGER(w2f__i4) OpenAD_Symbol_70
INTEGER(w2f__i4) OpenAD_Symbol_71
REAL(w2f__8) OpenAD_acc_2
REAL(w2f__8) OpenAD_acc_3
REAL(w2f__8) OpenAD_acc_4
REAL(w2f__8) OpenAD_acc_5
REAL(w2f__8) OpenAD_acc_6
REAL(w2f__8) OpenAD_acc_7
REAL(w2f__8) OpenAD_acc_8
REAL(w2f__8) OpenAD_acc_9
REAL(w2f__8) OpenAD_aux_0
REAL(w2f__8) OpenAD_aux_1
REAL(w2f__8) OpenAD_aux_10
REAL(w2f__8) OpenAD_aux_11
REAL(w2f__8) OpenAD_aux_2
REAL(w2f__8) OpenAD_aux_3
REAL(w2f__8) OpenAD_aux_4
REAL(w2f__8) OpenAD_aux_5
REAL(w2f__8) OpenAD_aux_6
REAL(w2f__8) OpenAD_aux_7
REAL(w2f__8) OpenAD_aux_8
REAL(w2f__8) OpenAD_aux_9
REAL(w2f__8) OpenAD_lin_10
REAL(w2f__8) OpenAD_lin_11
REAL(w2f__8) OpenAD_lin_12
REAL(w2f__8) OpenAD_lin_13
REAL(w2f__8) OpenAD_lin_4
REAL(w2f__8) OpenAD_lin_5
REAL(w2f__8) OpenAD_lin_6
REAL(w2f__8) OpenAD_lin_7
REAL(w2f__8) OpenAD_lin_8
REAL(w2f__8) OpenAD_lin_9
type(active) :: OpenAD_prp_3
type(active) :: OpenAD_prp_4
type(active) :: OpenAD_prp_5


    integer iaddr
    external iaddr
!
!     **** Statements ****
!

   if (our_rev_mode%plain) then
! original function
DO I = 1, 79, 1
  IF((I .GT. 1) .AND.(I .LT. 79)) THEN
    F(INT(I))%v = ((H(I)%v*8.92710038185119628906D+03*(H(I+1)%v-H(I+(-1))%v)*5.0&
     &D-01)/DX)

  ELSE
    IF (I.eq.1) THEN
      F(INT(I))%v = ((H(I)%v*8.92710038185119628906D+03*(H(I+1)%v-H(I)%v))/DX)
    ELSE
      IF (I.eq.79) THEN
        F(INT(I))%v = ((H(I)%v*8.92710038185119628906D+03*(H(I)%v-H(I+(-1))%v))/&
     &DX)

      ENDIF
    ENDIF
  ENDIF
END DO
FEND%v = (((H(79)%v**2)*8.92710038185119628906D+03+(-8.22421385178565979004D+09)&
     &)*5.0D-01)

    end if
    if (our_rev_mode%tape) then
! taping
OpenAD_Symbol_49 = 0_w2f__i8
DO I = 1,79,1
  IF ((I.GT.1).AND.(I.LT.79)) THEN
    OpenAD_aux_2 = (H(I)%v*8.92710038185119628906D+03)
    OpenAD_aux_3 = (H(I+1)%v-H(I+(-1))%v)
    OpenAD_aux_1 = (OpenAD_aux_2*OpenAD_aux_3)
    OpenAD_aux_0 = (OpenAD_aux_1*5.0D-01)
    OpenAD_lin_5 = OpenAD_aux_3
    OpenAD_lin_6 = OpenAD_aux_2
    OpenAD_lin_4 = (INT(1_w2f__i8)/DX)
    F(INT(I))%v = (OpenAD_aux_0/DX)
    OpenAD_acc_2 = (5.0D-01*OpenAD_lin_4)
    OpenAD_acc_3 = (OpenAD_lin_6*OpenAD_acc_2)
    OpenAD_acc_4 = (8.92710038185119628906D+03*OpenAD_lin_5*OpenAD_acc_2)
    OpenAD_Symbol_68 = (I+1)
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_68
    integer_tape_pointer = integer_tape_pointer+1
    OpenAD_Symbol_69 = (I+(-1))
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_69
    integer_tape_pointer = integer_tape_pointer+1
    double_tape(double_tape_pointer) = OpenAD_acc_3
    double_tape_pointer = double_tape_pointer+1
    double_tape(double_tape_pointer) = OpenAD_acc_4
    double_tape_pointer = double_tape_pointer+1
    integer_tape(integer_tape_pointer) = I
    integer_tape_pointer = integer_tape_pointer+1
    OpenAD_Symbol_54 = 1_w2f__i8
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_54
    integer_tape_pointer = integer_tape_pointer+1
  ELSE
    IF (I.eq.1) THEN
      OpenAD_aux_5 = (H(I)%v*8.92710038185119628906D+03)
      OpenAD_aux_6 = (H(I+1)%v-H(I)%v)
      OpenAD_aux_4 = (OpenAD_aux_5*OpenAD_aux_6)
      OpenAD_lin_8 = OpenAD_aux_6
      OpenAD_lin_9 = OpenAD_aux_5
      OpenAD_lin_7 = (INT(1_w2f__i8)/DX)
      F(INT(I))%v = (OpenAD_aux_4/DX)
      OpenAD_acc_5 = (OpenAD_lin_9*OpenAD_lin_7)
      OpenAD_acc_6 = (8.92710038185119628906D+03*OpenAD_lin_8*OpenAD_lin_7)
      OpenAD_Symbol_70 = (I+1)
      integer_tape(integer_tape_pointer) = OpenAD_Symbol_70
      integer_tape_pointer = integer_tape_pointer+1
      double_tape(double_tape_pointer) = OpenAD_acc_5
      double_tape_pointer = double_tape_pointer+1
      double_tape(double_tape_pointer) = OpenAD_acc_6
      double_tape_pointer = double_tape_pointer+1
      integer_tape(integer_tape_pointer) = I
      integer_tape_pointer = integer_tape_pointer+1
      OpenAD_Symbol_52 = 1_w2f__i8
      integer_tape(integer_tape_pointer) = OpenAD_Symbol_52
      integer_tape_pointer = integer_tape_pointer+1
    ELSE
      IF (I.eq.79) THEN
        OpenAD_aux_8 = (H(I)%v*8.92710038185119628906D+03)
        OpenAD_aux_9 = (H(I)%v-H(I+(-1))%v)
        OpenAD_aux_7 = (OpenAD_aux_8*OpenAD_aux_9)
        OpenAD_lin_11 = OpenAD_aux_9
        OpenAD_lin_12 = OpenAD_aux_8
        OpenAD_lin_10 = (INT(1_w2f__i8)/DX)
        F(INT(I))%v = (OpenAD_aux_7/DX)
        OpenAD_acc_7 = (OpenAD_lin_12*OpenAD_lin_10)
        OpenAD_acc_8 = (8.92710038185119628906D+03*OpenAD_lin_11*OpenAD_lin_10)
        OpenAD_Symbol_71 = (I+(-1))
        integer_tape(integer_tape_pointer) = OpenAD_Symbol_71
        integer_tape_pointer = integer_tape_pointer+1
        double_tape(double_tape_pointer) = OpenAD_acc_7
        double_tape_pointer = double_tape_pointer+1
        double_tape(double_tape_pointer) = OpenAD_acc_8
        double_tape_pointer = double_tape_pointer+1
        integer_tape(integer_tape_pointer) = I
        integer_tape_pointer = integer_tape_pointer+1
        OpenAD_Symbol_50 = 1_w2f__i8
        integer_tape(integer_tape_pointer) = OpenAD_Symbol_50
        integer_tape_pointer = integer_tape_pointer+1
      ELSE
        OpenAD_Symbol_51 = 0_w2f__i8
        integer_tape(integer_tape_pointer) = OpenAD_Symbol_51
        integer_tape_pointer = integer_tape_pointer+1
      ENDIF
      OpenAD_Symbol_53 = 0_w2f__i8
      integer_tape(integer_tape_pointer) = OpenAD_Symbol_53
      integer_tape_pointer = integer_tape_pointer+1
    ENDIF
    OpenAD_Symbol_55 = 0_w2f__i8
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_55
    integer_tape_pointer = integer_tape_pointer+1
  ENDIF
  OpenAD_Symbol_49 = (INT(OpenAD_Symbol_49)+INT(1_w2f__i8))
END DO
integer_tape(integer_tape_pointer) = OpenAD_Symbol_49
integer_tape_pointer = integer_tape_pointer+1
OpenAD_aux_11 = (H(79)%v**2)
OpenAD_aux_10 = (OpenAD_aux_11*8.92710038185119628906D+03+(-8.224213851785659790&
     &04D+09))

OpenAD_lin_13 = (2*(H(79)%v**(2-INT(1_w2f__i8))))
FEND%v = (OpenAD_aux_10*5.0D-01)
OpenAD_acc_9 = (OpenAD_lin_13*8.92710038185119628906D+03*5.0D-01)
double_tape(double_tape_pointer) = OpenAD_acc_9
double_tape_pointer = double_tape_pointer+1
    end if
    if (our_rev_mode%adjoint) then
! adjoint
double_tape_pointer = double_tape_pointer-1
OpenAD_Symbol_123 = double_tape(double_tape_pointer)
H(79)%d = H(79)%d+FEND%d*(OpenAD_Symbol_123)
FEND%d = 0.0d0
integer_tape_pointer = integer_tape_pointer-1
OpenAD_Symbol_44 = integer_tape(integer_tape_pointer)
OpenAD_Symbol_45 = 1
do while (INT(OpenAD_Symbol_45).LE.INT(OpenAD_Symbol_44))
  integer_tape_pointer = integer_tape_pointer-1
  OpenAD_Symbol_46 = integer_tape(integer_tape_pointer)
  IF (OpenAD_Symbol_46.ne.0) THEN
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_110 = integer_tape(integer_tape_pointer)
    double_tape_pointer = double_tape_pointer-1
    OpenAD_Symbol_111 = double_tape(double_tape_pointer)
    double_tape_pointer = double_tape_pointer-1
    OpenAD_Symbol_112 = double_tape(double_tape_pointer)
    H(OpenAD_Symbol_110)%d = H(OpenAD_Symbol_110)%d+F(OpenAD_Symbol_110)%d*(Open&
     &AD_Symbol_111)

    OpenAD_prp_3%d = OpenAD_prp_3%d+F(OpenAD_Symbol_110)%d*(OpenAD_Symbol_112)
    F(OpenAD_Symbol_110)%d = 0.0d0
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_113 = integer_tape(integer_tape_pointer)
    H(OpenAD_Symbol_113)%d = H(OpenAD_Symbol_113)%d-OpenAD_prp_3%d
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_114 = integer_tape(integer_tape_pointer)
    H(OpenAD_Symbol_114)%d = H(OpenAD_Symbol_114)%d+OpenAD_prp_3%d
    OpenAD_prp_3%d = 0.0d0
  ELSE
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_47 = integer_tape(integer_tape_pointer)
    IF (OpenAD_Symbol_47.ne.0) THEN
      integer_tape_pointer = integer_tape_pointer-1
      OpenAD_Symbol_115 = integer_tape(integer_tape_pointer)
      double_tape_pointer = double_tape_pointer-1
      OpenAD_Symbol_116 = double_tape(double_tape_pointer)
      double_tape_pointer = double_tape_pointer-1
      OpenAD_Symbol_117 = double_tape(double_tape_pointer)
      H(OpenAD_Symbol_115)%d = H(OpenAD_Symbol_115)%d+F(OpenAD_Symbol_115)%d*(Op&
     &enAD_Symbol_116)

      OpenAD_prp_4%d = OpenAD_prp_4%d+F(OpenAD_Symbol_115)%d*(OpenAD_Symbol_117)
      F(OpenAD_Symbol_115)%d = 0.0d0
      H(OpenAD_Symbol_115)%d = H(OpenAD_Symbol_115)%d-OpenAD_prp_4%d
      integer_tape_pointer = integer_tape_pointer-1
      OpenAD_Symbol_118 = integer_tape(integer_tape_pointer)
      H(OpenAD_Symbol_118)%d = H(OpenAD_Symbol_118)%d+OpenAD_prp_4%d
      OpenAD_prp_4%d = 0.0d0
    ELSE
      integer_tape_pointer = integer_tape_pointer-1
      OpenAD_Symbol_48 = integer_tape(integer_tape_pointer)
      IF (OpenAD_Symbol_48.ne.0) THEN
        integer_tape_pointer = integer_tape_pointer-1
        OpenAD_Symbol_119 = integer_tape(integer_tape_pointer)
        double_tape_pointer = double_tape_pointer-1
        OpenAD_Symbol_120 = double_tape(double_tape_pointer)
        double_tape_pointer = double_tape_pointer-1
        OpenAD_Symbol_121 = double_tape(double_tape_pointer)
        H(OpenAD_Symbol_119)%d = H(OpenAD_Symbol_119)%d+F(OpenAD_Symbol_119)%d*(&
     &OpenAD_Symbol_120)

        OpenAD_prp_5%d = OpenAD_prp_5%d+F(OpenAD_Symbol_119)%d*(OpenAD_Symbol_12&
     &1)

        F(OpenAD_Symbol_119)%d = 0.0d0
        integer_tape_pointer = integer_tape_pointer-1
        OpenAD_Symbol_122 = integer_tape(integer_tape_pointer)
        H(OpenAD_Symbol_122)%d = H(OpenAD_Symbol_122)%d-OpenAD_prp_5%d
        H(OpenAD_Symbol_119)%d = H(OpenAD_Symbol_119)%d+OpenAD_prp_5%d
        OpenAD_prp_5%d = 0.0d0
      ENDIF
    ENDIF
  ENDIF
  OpenAD_Symbol_45 = INT(OpenAD_Symbol_45)+1
END DO
    end if
  end subroutine stream_vel_taud
!#########################################################
! This file is part of OpenAD released under the LGPL.   #
! The full COPYRIGHT notice can be found in the top      #
! level directory of the OpenAD distribution             #
!#########################################################

SUBROUTINE stream_vel_visc(H, U, NU)
    use OAD_tape
    use OAD_rev

! original arguments get inserted before version
!         ! and declared here together with all local variables
!         ! generated by xaifBooster

use OAD_active
use w2f__types
use oad_intrinsics
use stream_vel_variables
use oad_intrinsics
use stream_vel_variables
use oad_intrinsics
use stream_vel_variables
IMPLICIT NONE
!
!     **** Global Variables & Derived Type Definitions ****
!
INTEGER(w2f__i8) OpenAD_Symbol_72
INTEGER(w2f__i8) OpenAD_Symbol_73
INTEGER(w2f__i8) OpenAD_Symbol_74
INTEGER(w2f__i8) OpenAD_Symbol_75
INTEGER(w2f__i8) OpenAD_Symbol_76
INTEGER(w2f__i8) OpenAD_Symbol_77
!
!     **** Parameters and Result ****
!
type(active) :: H(1:79)
type(active) :: U(1:80)
type(active) :: NU(1:79)
!
!     **** Local Variables and Functions ****
!
INTEGER(w2f__i4) I
type(active) :: TMP
type(active) :: UX
INTEGER(w2f__i4) OpenAD_Symbol_124
REAL(w2f__8) OpenAD_Symbol_125
REAL(w2f__8) OpenAD_Symbol_126
INTEGER(w2f__i4) OpenAD_Symbol_127
INTEGER(w2f__i4) OpenAD_Symbol_78
REAL(w2f__8) OpenAD_acc_10
REAL(w2f__8) OpenAD_acc_11
REAL(w2f__8) OpenAD_aux_12
REAL(w2f__8) OpenAD_aux_13
REAL(w2f__8) OpenAD_aux_14
REAL(w2f__8) OpenAD_aux_15
REAL(w2f__8) OpenAD_lin_14
REAL(w2f__8) OpenAD_lin_15
REAL(w2f__8) OpenAD_lin_16
REAL(w2f__8) OpenAD_lin_17
REAL(w2f__8) OpenAD_lin_18
type(active) :: OpenAD_prp_6


    integer iaddr
    external iaddr
!
!     **** Statements ****
!

   if (our_rev_mode%plain) then
! original function
DO I = 1, 79, 1
  UX%v = ((U(I+1)%v-U(I)%v)/DX)
  TMP%v = ((UX%v**2)+1.00000002337219498273D-14)
  NU(INT(I))%v = ((TMP%v**(-3.3333333333333331483D-01))*H(I)%v*5.0D-01*2.7143814&
     &3778400379233D+05)

END DO
    end if
    if (our_rev_mode%tape) then
! taping
OpenAD_Symbol_74 = 0_w2f__i8
DO I = 1,79,1
  OpenAD_aux_12 = (U(I+1)%v-U(I)%v)
  OpenAD_lin_14 = (INT(1_w2f__i8)/DX)
  UX%v = (OpenAD_aux_12/DX)
  OpenAD_lin_15 = (2*(UX%v**(2-INT(1_w2f__i8))))
  TMP%v = ((UX%v**2)+1.00000002337219498273D-14)
  OpenAD_aux_13 = (TMP%v**(-3.3333333333333331483D-01))
  OpenAD_aux_15 = (H(I)%v*5.0D-01)
  OpenAD_aux_14 = (OpenAD_aux_15*2.71438143778400379233D+05)
  OpenAD_lin_18 = ((-3.3333333333333331483D-01)*(TMP%v**((-3.3333333333333331483&
     &D-01)-INT(1_w2f__i8))))

  OpenAD_lin_16 = OpenAD_aux_14
  OpenAD_lin_17 = OpenAD_aux_13
  NU(INT(I))%v = (OpenAD_aux_13*OpenAD_aux_14)
  OpenAD_acc_10 = (5.0D-01*2.71438143778400379233D+05*OpenAD_lin_17)
  OpenAD_acc_11 = (OpenAD_lin_14*OpenAD_lin_15*OpenAD_lin_18*OpenAD_lin_16)
  OpenAD_Symbol_78 = (I+1)
  integer_tape(integer_tape_pointer) = OpenAD_Symbol_78
  integer_tape_pointer = integer_tape_pointer+1
  double_tape(double_tape_pointer) = OpenAD_acc_10
  double_tape_pointer = double_tape_pointer+1
  double_tape(double_tape_pointer) = OpenAD_acc_11
  double_tape_pointer = double_tape_pointer+1
  integer_tape(integer_tape_pointer) = I
  integer_tape_pointer = integer_tape_pointer+1
  OpenAD_Symbol_74 = (INT(OpenAD_Symbol_74)+INT(1_w2f__i8))
END DO
integer_tape(integer_tape_pointer) = OpenAD_Symbol_74
integer_tape_pointer = integer_tape_pointer+1
    end if
    if (our_rev_mode%adjoint) then
! adjoint
integer_tape_pointer = integer_tape_pointer-1
OpenAD_Symbol_72 = integer_tape(integer_tape_pointer)
OpenAD_Symbol_73 = 1
do while (INT(OpenAD_Symbol_73).LE.INT(OpenAD_Symbol_72))
  integer_tape_pointer = integer_tape_pointer-1
  OpenAD_Symbol_124 = integer_tape(integer_tape_pointer)
  double_tape_pointer = double_tape_pointer-1
  OpenAD_Symbol_125 = double_tape(double_tape_pointer)
  double_tape_pointer = double_tape_pointer-1
  OpenAD_Symbol_126 = double_tape(double_tape_pointer)
  OpenAD_prp_6%d = OpenAD_prp_6%d+NU(OpenAD_Symbol_124)%d*(OpenAD_Symbol_125)
  H(OpenAD_Symbol_124)%d = H(OpenAD_Symbol_124)%d+NU(OpenAD_Symbol_124)%d*(OpenA&
     &D_Symbol_126)

  NU(OpenAD_Symbol_124)%d = 0.0d0
  U(OpenAD_Symbol_124)%d = U(OpenAD_Symbol_124)%d-OpenAD_prp_6%d
  integer_tape_pointer = integer_tape_pointer-1
  OpenAD_Symbol_127 = integer_tape(integer_tape_pointer)
  U(OpenAD_Symbol_127)%d = U(OpenAD_Symbol_127)%d+OpenAD_prp_6%d
  OpenAD_prp_6%d = 0.0d0
  OpenAD_Symbol_73 = INT(OpenAD_Symbol_73)+1
END DO
    end if
  end subroutine stream_vel_visc
!#########################################################
! This file is part of OpenAD released under the LGPL.   #
! The full COPYRIGHT notice can be found in the top      #
! level directory of the OpenAD distribution             #
!#########################################################

SUBROUTINE stream_assemble(NU, BETA_FRIC, A)
    use OAD_tape
    use OAD_rev

! original arguments get inserted before version
!         ! and declared here together with all local variables
!         ! generated by xaifBooster

use OAD_active
use w2f__types
use oad_intrinsics
use stream_vel_variables
use oad_intrinsics
use stream_vel_variables
use oad_intrinsics
use stream_vel_variables
IMPLICIT NONE
!
!     **** Global Variables & Derived Type Definitions ****
!
INTEGER(w2f__i8) OpenAD_Symbol_79
INTEGER(w2f__i8) OpenAD_Symbol_80
INTEGER(w2f__i8) OpenAD_Symbol_81
INTEGER(w2f__i8) OpenAD_Symbol_82
INTEGER(w2f__i8) OpenAD_Symbol_83
INTEGER(w2f__i8) OpenAD_Symbol_84
INTEGER(w2f__i8) OpenAD_Symbol_85
INTEGER(w2f__i8) OpenAD_Symbol_86
INTEGER(w2f__i8) OpenAD_Symbol_87
INTEGER(w2f__i8) OpenAD_Symbol_88
INTEGER(w2f__i8) OpenAD_Symbol_89
INTEGER(w2f__i8) OpenAD_Symbol_90
INTEGER(w2f__i8) OpenAD_Symbol_91
INTEGER(w2f__i8) OpenAD_Symbol_92
INTEGER(w2f__i8) OpenAD_Symbol_93
INTEGER(w2f__i8) OpenAD_Symbol_94
INTEGER(w2f__i8) OpenAD_Symbol_95
INTEGER(w2f__i8) OpenAD_Symbol_96
!
!     **** Parameters and Result ****
!
type(active) :: NU(1:79)
REAL(w2f__8) BETA_FRIC(1 : 79)
INTENT(IN) BETA_FRIC
type(active) :: A(1:79,1:3)
!
!     **** Local Variables and Functions ****
!
INTEGER(w2f__i4) I
INTEGER(w2f__i4) OpenAD_Symbol_128
REAL(w2f__8) OpenAD_Symbol_129
INTEGER(w2f__i4) OpenAD_Symbol_130
REAL(w2f__8) OpenAD_Symbol_131
INTEGER(w2f__i4) OpenAD_Symbol_132
REAL(w2f__8) OpenAD_Symbol_133
INTEGER(w2f__i4) OpenAD_Symbol_134
INTEGER(w2f__i4) OpenAD_Symbol_135
REAL(w2f__8) OpenAD_Symbol_136
INTEGER(w2f__i4) OpenAD_Symbol_137
INTEGER(w2f__i4) OpenAD_Symbol_97
INTEGER(w2f__i4) OpenAD_Symbol_98
REAL(w2f__8) OpenAD_acc_12
REAL(w2f__8) OpenAD_acc_13
REAL(w2f__8) OpenAD_acc_14
REAL(w2f__8) OpenAD_acc_15
REAL(w2f__8) OpenAD_aux_16
REAL(w2f__8) OpenAD_aux_17
REAL(w2f__8) OpenAD_aux_18
REAL(w2f__8) OpenAD_aux_19
REAL(w2f__8) OpenAD_lin_19
REAL(w2f__8) OpenAD_lin_20
REAL(w2f__8) OpenAD_lin_21
REAL(w2f__8) OpenAD_lin_22
type(active) :: OpenAD_prp_7


    integer iaddr
    external iaddr
!
!     **** Statements ****
!

   if (our_rev_mode%plain) then
! original function
DO I = 1, 79, 1
  A(INT(I),2)%v = (((NU(I)%v*4.0D00)/DX)+(DX/3.0D00)*(BETA_FRIC(I)**2))
  IF (I.GT.1) THEN
    A(INT(I),1)%v = ((DX/6.0D00)*(BETA_FRIC(I)**2)-((NU(I)%v*4.0D00)/DX))
  ENDIF
  IF (I.LT.79) THEN
    A(INT(I),2)%v = (A(I,2)%v+((NU(I+1)%v*4.0D00)/DX)+(DX/3.0D00)*(BETA_FRIC(I+1&
     &)**2))

    A(INT(I),3)%v = ((DX/6.0D00)*(BETA_FRIC(I+1)**2)-((NU(I+1)%v*4.0D00)/DX))
  ENDIF
END DO
    end if
    if (our_rev_mode%tape) then
! taping
OpenAD_Symbol_83 = 0_w2f__i8
DO I = 1,79,1
  OpenAD_aux_16 = (NU(I)%v*4.0D00)
  OpenAD_lin_19 = (INT(1_w2f__i8)/DX)
  A(INT(I),2)%v = ((OpenAD_aux_16/DX)+(DX/3.0D00)*(BETA_FRIC(I)**2))
  OpenAD_acc_12 = (4.0D00*OpenAD_lin_19)
  double_tape(double_tape_pointer) = OpenAD_acc_12
  double_tape_pointer = double_tape_pointer+1
  integer_tape(integer_tape_pointer) = I
  integer_tape_pointer = integer_tape_pointer+1
  IF (I.GT.1) THEN
    OpenAD_aux_17 = (NU(I)%v*4.0D00)
    OpenAD_lin_20 = (INT(1_w2f__i8)/DX)
    A(INT(I),1)%v = ((DX/6.0D00)*(BETA_FRIC(I)**2)-(OpenAD_aux_17/DX))
    OpenAD_acc_13 = (4.0D00*OpenAD_lin_20*INT((-1_w2f__i8)))
    double_tape(double_tape_pointer) = OpenAD_acc_13
    double_tape_pointer = double_tape_pointer+1
    integer_tape(integer_tape_pointer) = I
    integer_tape_pointer = integer_tape_pointer+1
    OpenAD_Symbol_84 = 1_w2f__i8
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_84
    integer_tape_pointer = integer_tape_pointer+1
  ELSE
    OpenAD_Symbol_85 = 0_w2f__i8
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_85
    integer_tape_pointer = integer_tape_pointer+1
  ENDIF
  IF (I.LT.79) THEN
    OpenAD_aux_18 = (NU(I+1)%v*4.0D00)
    OpenAD_lin_21 = (INT(1_w2f__i8)/DX)
    A(INT(I),2)%v = (A(I,2)%v+(OpenAD_aux_18/DX)+(DX/3.0D00)*(BETA_FRIC(I+1)**2)&
     &)

    OpenAD_acc_14 = (4.0D00*OpenAD_lin_21)
    double_tape(double_tape_pointer) = OpenAD_acc_14
    double_tape_pointer = double_tape_pointer+1
    OpenAD_Symbol_97 = (I+1)
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_97
    integer_tape_pointer = integer_tape_pointer+1
    integer_tape(integer_tape_pointer) = I
    integer_tape_pointer = integer_tape_pointer+1
    OpenAD_aux_19 = (NU(I+1)%v*4.0D00)
    OpenAD_lin_22 = (INT(1_w2f__i8)/DX)
    A(INT(I),3)%v = ((DX/6.0D00)*(BETA_FRIC(I+1)**2)-(OpenAD_aux_19/DX))
    OpenAD_acc_15 = (4.0D00*OpenAD_lin_22*INT((-1_w2f__i8)))
    double_tape(double_tape_pointer) = OpenAD_acc_15
    double_tape_pointer = double_tape_pointer+1
    OpenAD_Symbol_98 = (I+1)
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_98
    integer_tape_pointer = integer_tape_pointer+1
    integer_tape(integer_tape_pointer) = I
    integer_tape_pointer = integer_tape_pointer+1
    OpenAD_Symbol_86 = 1_w2f__i8
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_86
    integer_tape_pointer = integer_tape_pointer+1
  ELSE
    OpenAD_Symbol_87 = 0_w2f__i8
    integer_tape(integer_tape_pointer) = OpenAD_Symbol_87
    integer_tape_pointer = integer_tape_pointer+1
  ENDIF
  OpenAD_Symbol_83 = (INT(OpenAD_Symbol_83)+INT(1_w2f__i8))
END DO
integer_tape(integer_tape_pointer) = OpenAD_Symbol_83
integer_tape_pointer = integer_tape_pointer+1
    end if
    if (our_rev_mode%adjoint) then
! adjoint
integer_tape_pointer = integer_tape_pointer-1
OpenAD_Symbol_79 = integer_tape(integer_tape_pointer)
OpenAD_Symbol_80 = 1
do while (INT(OpenAD_Symbol_80).LE.INT(OpenAD_Symbol_79))
  integer_tape_pointer = integer_tape_pointer-1
  OpenAD_Symbol_81 = integer_tape(integer_tape_pointer)
  IF (OpenAD_Symbol_81.ne.0) THEN
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_132 = integer_tape(integer_tape_pointer)
    double_tape_pointer = double_tape_pointer-1
    OpenAD_Symbol_133 = double_tape(double_tape_pointer)
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_134 = integer_tape(integer_tape_pointer)
    NU(OpenAD_Symbol_134)%d = NU(OpenAD_Symbol_134)%d+A(OpenAD_Symbol_132,3)%d*(&
     &OpenAD_Symbol_133)

    A(OpenAD_Symbol_132,3)%d = 0.0d0
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_135 = integer_tape(integer_tape_pointer)
    double_tape_pointer = double_tape_pointer-1
    OpenAD_Symbol_136 = double_tape(double_tape_pointer)
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_137 = integer_tape(integer_tape_pointer)
    NU(OpenAD_Symbol_137)%d = NU(OpenAD_Symbol_137)%d+A(OpenAD_Symbol_135,2)%d*(&
     &OpenAD_Symbol_136)

    OpenAD_prp_7%d = OpenAD_prp_7%d+A(OpenAD_Symbol_135,2)%d
    A(OpenAD_Symbol_135,2)%d = 0.0d0
    A(OpenAD_Symbol_135,2)%d = A(OpenAD_Symbol_135,2)%d+OpenAD_prp_7%d
    OpenAD_prp_7%d = 0.0d0
  ENDIF
  integer_tape_pointer = integer_tape_pointer-1
  OpenAD_Symbol_82 = integer_tape(integer_tape_pointer)
  IF (OpenAD_Symbol_82.ne.0) THEN
    integer_tape_pointer = integer_tape_pointer-1
    OpenAD_Symbol_130 = integer_tape(integer_tape_pointer)
    double_tape_pointer = double_tape_pointer-1
    OpenAD_Symbol_131 = double_tape(double_tape_pointer)
    NU(OpenAD_Symbol_130)%d = NU(OpenAD_Symbol_130)%d+A(OpenAD_Symbol_130,1)%d*(&
     &OpenAD_Symbol_131)

    A(OpenAD_Symbol_130,1)%d = 0.0d0
  ENDIF
  integer_tape_pointer = integer_tape_pointer-1
  OpenAD_Symbol_128 = integer_tape(integer_tape_pointer)
  double_tape_pointer = double_tape_pointer-1
  OpenAD_Symbol_129 = double_tape(double_tape_pointer)
  NU(OpenAD_Symbol_128)%d = NU(OpenAD_Symbol_128)%d+A(OpenAD_Symbol_128,2)%d*(Op&
     &enAD_Symbol_129)

  A(OpenAD_Symbol_128,2)%d = 0.0d0
  OpenAD_Symbol_80 = INT(OpenAD_Symbol_80)+1
END DO
    end if
  end subroutine stream_assemble
