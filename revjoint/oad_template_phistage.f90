      subroutine template()
      use OAD_tape
      use OAD_rev
      use OAD_cp
!$TEMPLATE_PRAGMA_DECLARATIONS
    type(modeType) :: our_orig_mode

! dummy variables for calling adjoint computation without side effects, storing
! adjoint variable iterates
    type(active) :: B_DUMMY(1:79)
    type(active) :: H_DUMMY(1:79)
    type(active) :: U_DUMMY(1:80)

! checkpointing stacks and offsets
    integer :: cp_loop_variable_1,cp_loop_variable_2,cp_loop_variable_3,cp_loop_&
     &variable_4,cp_loop_variable_5,cp_loop_variable_6


! floats 'F'
    double precision, dimension(:), allocatable, save :: theArgFStack
    integer, save :: theArgFStackoffset=0, theArgFStackSize=0

    ! added to roll back file offset after "look" into stack
    integer, save :: theArgFStackoffsetTemp=0

! integers 'I'
    integer, dimension(:), allocatable, save :: theArgIStack
    integer, save :: theArgIStackoffset=0, theArgIStackSize=0
! booleans 'B'
    logical, dimension(:), allocatable, save :: theArgBStack
    integer, save :: theArgBStackoffset=0, theArgBStackSize=0
! strings 'S'
    character*(80), dimension(:), allocatable, save :: theArgSStack
    integer, save :: theArgSStackoffset=0, theArgSStackSize=0


! external C function used in inlined code
    integer iaddr
    external iaddr
!
!     **** Statements ****
!

    if (our_rev_mode%arg_store) then
! store arguments
call cp_store_real_scalar(DX,theArgFStack,theArgFStackoffset,theArgFStackSize)
call cp_store_real_vector(U,size(U),theArgFStack,theArgFStackoffset,theArgFStack&
     &Size)

call cp_store_real_vector(U_IP1,size(U_IP1),theArgFStack,theArgFStackoffset,theA&
     &rgFStackSize)

call cp_store_real_vector(B,size(B),theArgFStack,theArgFStackoffset,theArgFStack&
     &Size)

call cp_store_real_vector(H,size(H),theArgFStack,theArgFStackoffset,theArgFStack&
     &Size)

call cp_store_p_real_vector(BETA_FRIC,size(BETA_FRIC),theArgFStack,theArgFStacko&
     &ffset,theArgFStackSize)

    end if
    if (our_rev_mode%arg_restore) then
! restore arguments

! current top of stack -- for look functionality
theArgFStackoffsetTemp=theArgFStackoffset

do cp_loop_variable_1 = ubound(BETA_FRIC,1),lbound(BETA_FRIC,1),-1
BETA_FRIC(cp_loop_variable_1) = theArgFStack(theArgFStackoffset)
theArgFStackoffset = theArgFStackoffset-1
end do
do cp_loop_variable_1 = ubound(H,1),lbound(H,1),-1
H(cp_loop_variable_1)%v = theArgFStack(theArgFStackoffset)
theArgFStackoffset = theArgFStackoffset-1
end do
do cp_loop_variable_1 = ubound(B,1),lbound(B,1),-1
B(cp_loop_variable_1)%v = theArgFStack(theArgFStackoffset)
theArgFStackoffset = theArgFStackoffset-1
end do
do cp_loop_variable_1 = ubound(U_IP1,1),lbound(U_IP1,1),-1
U_IP1(cp_loop_variable_1)%v = theArgFStack(theArgFStackoffset)
theArgFStackoffset = theArgFStackoffset-1
end do
do cp_loop_variable_1 = ubound(U,1),lbound(U,1),-1
U(cp_loop_variable_1)%v = theArgFStack(theArgFStackoffset)
theArgFStackoffset = theArgFStackoffset-1
end do
DX = theArgFStack(theArgFStackoffset)
theArgFStackoffset = theArgFStackoffset-1
    end if
    if (our_rev_mode%plain) then
      our_orig_mode=our_rev_mode
      our_rev_mode%arg_store=.FALSE.
! original function
CALL phi(U,U_IP1,B,H,BETA_FRIC)
U(1:80)%v = U_IP1(1:80)%v

! original function end
      our_rev_mode=our_orig_mode
    end if
    if (our_rev_mode%tape) then
!            print*, " tape       ", our_rev_mode
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%plain=.TRUE.
      our_rev_mode%tape=.FALSE.
      our_rev_mode%adjoint=.FALSE.
! taping

 our_rev_mode%arg_store=.true.
 if (isinloop.eq.1) then
 our_rev_mode%arg_store=.false.
 endif

 CALL phi(U,U_IP1,B,H,BETA_FRIC)
 U(1:80)%v = U_IP1(1:80)%v

 if (isinloop.eq.1) then
 our_rev_mode%arg_store=.true.
 endif

! taping end
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%plain=.FALSE.
      our_rev_mode%tape=.FALSE.
      our_rev_mode%adjoint=.TRUE.
    end if
    if (our_rev_mode%adjoint) then
!            print*, " adjoint    ", our_rev_mode
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.TRUE.
      our_rev_mode%plain=.FALSE.
      our_rev_mode%tape=.TRUE.
      our_rev_mode%adjoint=.FALSE.
! adjoint

 select case (isinloop) 

  case (1)

  B_DUMMY(1:79)%d=B(1:79)%d
  H_DUMMY(1:79)%d=H(1:79)%d
  U_DUMMY(1:80)%d = U(1:80)%d

!  U_IP1(1:80)%d = U_IP1(1:80)%d+U(1:80)%d
!  U(1:80)%d = 0.0d0

  our_rev_mode%arg_look = .true.
  CALL phi(U_DUMMY,U_IP1,B_DUMMY,H_DUMMY,BETA_FRIC)
  our_rev_mode%arg_look = .false.

  U_IP1(1:80)%d = U_DUMMY(1:80)%d

 case (2)

  U_IP1(1:80)%d = U_IP1(1:80)%d+U(1:80)%d
  U(1:80)%d = 0.0d0
  CALL phi(U,U_IP1,B,H,BETA_FRIC)

 case (0)

  CALL phi(U,U_IP1,B,H,BETA_FRIC)

end select

! adjoint end
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.TRUE.
      our_rev_mode%plain=.FALSE.
      our_rev_mode%tape=.TRUE.
      our_rev_mode%adjoint=.FALSE.
    end if

    ! added -- if called in "look" mode this 
    !          will prevent stack offset changes
    if (our_rev_mode%arg_look) then
     theArgFStackoffset=theArgFStackoffsetTemp
    end if


      end subroutine template
