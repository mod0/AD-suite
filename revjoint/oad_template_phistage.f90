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
    integer, save :: theArgFStackoffset=0, theArgFStackSize=0, theArgFStackoffsetTemp=0
! integers 'I'
    integer, dimension(:), allocatable, save :: theArgIStack
    integer, save :: theArgIStackoffset=0, theArgIStackSize=0
! booleans 'B'
    logical, dimension(:), allocatable, save :: theArgBStack
    integer, save :: theArgBStackoffset=0, theArgBStackSize=0
! strings 'S'
    character*(80), dimension(:), allocatable, save :: theArgSStack
    integer, save :: theArgSStackoffset=0, theArgSStackSize=0
    integer :: i

! external C function used in inlined code
    integer iaddr
    external iaddr
!
!     **** Statements ****
!

    if (our_rev_mode%arg_store) then
if(isinloop.eq.2) then
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
    end if
    if (our_rev_mode%arg_restore) then
!Set the Temporary offset to the actual offset
theArgFStackoffsetTemp = theArgFStackoffset

! restore arguments
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
if(isinloop.ne.0) then
! The call in phistage_0 has no effect on the computation
  CALL phi(U,U_IP1,B,H,BETA_FRIC)
end if
! original function end
      our_rev_mode=our_orig_mode
    end if
    if (our_rev_mode%tape) then
!            print*, " tape       ", our_rev_mode
      our_rev_mode%arg_store=.TRUE.
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%plain=.TRUE.
      our_rev_mode%tape=.FALSE.
      our_rev_mode%adjoint=.FALSE.
! taping
CALL phi(U,U_IP1,B,H,BETA_FRIC)

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
    if(isinloop.eq.0) then
      do i=n+1,1,-1
        call pop_s0(U_DUMMY(i)%d)
      end do
! adjoint
      CALL phi(U,U_IP1,B,H,BETA_FRIC)
    end if 
    if(isinloop.eq.1) then
      do i=n+1,1,-1
        call pop_s0(U_DUMMY(i)%d)
      end do
      do i=1,n+1
        call push_s0(U_DUMMY(i)%d)
      end do
      B_DUMMY = B
      H_DUMMY = H
      U_DUMMY%v = U%v
! adjoint
     CALL phi(U_DUMMY,U_IP1,B_DUMMY,H_DUMMY,BETA_FRIC)
     U_IP1%d = U_DUMMY%d
    end if
!We have to push the value of U%d as it is overwritten in the caller (stream_vel in this case)
    if(isinloop.eq.2) then
! adjoint
      CALL phi(U,U_IP1,B,H,BETA_FRIC)
      do i=1,n+1
        call push_s0(U(i)%d)
      end do
    end if
! adjoint end
    if(isinloop.ne.0) then
      !Set the Temporary offset to the actual offset
      theArgFStackoffset = theArgFStackoffsetTemp
    end if
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.TRUE.
      our_rev_mode%plain=.FALSE.
      our_rev_mode%tape=.TRUE.
      our_rev_mode%adjoint=.FALSE.
    end if
      end subroutine template
