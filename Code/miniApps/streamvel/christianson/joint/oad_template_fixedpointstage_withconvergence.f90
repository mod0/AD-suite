! Given an fixed point iteration loop, this template is used to orchestrate 
! the appropriate use of the differentiated loop body. It is assumed that the 
! user is already able to pass the existing code through OpenAD and obtain
! correct output. The intent of these notes is not to explain *why* the 
! template works. Rather, the notes will help a user to modify the existing
! application and edit the template appropritely. 
! The template is inspired by the paper:
!
!@TECHREPORT{Christianson1992RAa,
!       AUTHOR = "Christianson, D. Bruce",
!       TITLE = "Reverse Accumulation and Attractive Fixed Points",
!       TYPE = "Technical Report",
!       NUMBER = "NOC TR 258",
!       INSTITUTION = "The Numerical Optimisation Centre, University of Hertfordshire",
!       ADDRESS = "Hatfield, UK",
!       MONTH = "March",
!       YEAR = "1992"
!} 
!
! The following steps must be performed:
! 1. Encapsulate the contents of the loop into a function. In the code below,
!    this function is assumed to be phi(). It is vital that phi() is written
!    such that U_i is purely an input variable and U_ip1 is the output variable
!    of the fixed point iteratation.
!
!subroutine phi (u_i, u_ip1, b, h, beta_fric)
!  use stream_vel_variables
!  use conj_grad_mod
!  real(8), dimension(n) :: b, h
!  real(8), intent(inout), dimension(n) :: beta_fric
!  real(8), intent(in), dimension(n+1) :: u_i
!  real(8), intent(out), dimension(n+1) :: u_ip1
!  real(8), dimension(n) :: nu,utmp
!  real(8), dimension(n,3) :: A
!  integer :: j
!  call stream_vel_visc (h, u_i, nu)              ! update viscosities              
!  call stream_assemble (nu, beta_fric, A)        ! assemble tridiag matrix
!                                                 ! this represents discretization of
!                                                 !  (nu^(i-1) u^(i)_x)_x - \beta^2 u^(i) = f
!  utmp = 0.
!  call solve (utmp, b, A)                        ! solve linear system for new u
!  do j=1,n
!   u_ip1(j+1) = utmp(j)                          ! effectively apply boundary condition u(1)==0
!  enddo
!end subroutine phi 
!
!2. Instead of calling phi() from within the fixedpoint loop, we instead call phistage()
!   which is a function whose body will be replaced by the contents of this template
!   in the differentiated code. phistage() will simply call phi(). phistage()
!   also contains a directive informing OpenAD that the differentiated contents of phistage()
!   should be replaced by the contents of oad_template_fixedpointstage.f90 (this file)
!
!subroutine phistage (u, u_ip1, b, h, beta_fric, isinloop)
!$openad xxx template oad_template_fixedpointstage.f90
!  use stream_vel_variables
!  real(8), dimension(n) :: b, h
!  real(8), intent(inout), dimension(n) :: beta_fric
!  real(8), intent(in), dimension(n+1) :: u
!  real(8), intent(out), dimension(n+1) :: u_ip1
!  real(8), dimension(n+1) :: u_ip1
!  integer, intent(in) :: isinloop 
!  call phi (u, u_ip1, b, h, beta_fric)
!end subroutine phistage 
!
!3. The template must be changed as follows:
!   a)As a user of the template, one must add code to store the arguments and global variables
!     that should normally be checkpointed for the reverse mode to function correctly. If one is 
!     unsure of what these variables are, then perform the following steps
!     i) From withing phistage() remove the line: $openad xxx template oad_template_fixedpointstage.f90
!     ii) Invoke OpenAD as usual to create differentiated code for phistage()
!     iii) copy the contents of the store and restore cases within the differentiated code into the
!        store and restore portions respectively of this template. They are itdentified by the comments
!        Begin user checkpoint store and Begin user checkpoint restore.
!   b) All portions itedentified by the comment Begin user function call should contain the call to phi(), 
!      just as it was present in the original function phistage()
!   c) Declarations:
!     i) A dummy variables is required for the output of the fixed point iteration which in this case is U. 
!     ii) Dummy variables are required for the active parameters of the fixed point iteration (B and H here). 
!     iii) A variable needed to hold the adjoint of the output of the fixed point iteration which in this case is U.
!          This variable has been called U_d_0 and required the save attibute.
!   d) In the section identified by the comment Begin user initialize dummy args, the following changes should be made:
!     i) The dummy parameters should be intialized by their non-dummy counterparts. 
!     ii)The dummy fixed point loop output's adjoint should be initialized by the adjoint 
!         of the output of the fixed point iteration. In this case, U_d_0.
!   e) The section identified by the comment Begin Begin user function call WITH DUMMY ARGS, 
!      should contain the call to phi(), with dummy arguments replacing the original arguments.
!   f) The section identified by the comment Begin user copy dummy value should copy the adjoint of the,
!      input to the call to phi() (U_DUMMY%d) to the adjoint of the output of the call to phi() (U_IP1%d). 
!   g) The section identified by the comment Begin user save value should contain an assignment of the adjoint 
!         of the output of the fixed point iteration, i.e,. U%d to the static variable U_0_d
!   h) No other changes are necessary at this moment
!
!4. Optional changes for convergence criteria of the original loop and adjoint loop
!   a) Modify phistage to contain declarations and 'stub-statements' used in the actual convergence criteria given in the template
!      The modified subroutine is given below
!
!subroutine phistage (u, u_ip1, b, h, beta_fric, isinloop)
!$openad xxx template oad_template_fixedpointstage.f90
!  use stream_vel_variables
!  real(8), dimension(n) :: b, h
!  real(8), intent(inout), dimension(n) :: beta_fric
!  real(8), intent(in), dimension(n+1) :: u
!  real(8), intent(out), dimension(n+1) :: u_ip1
!  real(8), dimension(n+1) :: u_ip1
!  integer, intent(in) :: isinloop 
!  integer, save :: iter=0, adj_iter=0
!  logical, save :: conv_flag=.FALSE., adj_conv_flag=.FALSE.
!  integer :: k
!  real(8) :: normdiff, normZ, diff
!  if (conv_flag .eqv. .false.) then
!    iter = iter + 1
!     call phi (u, u_ip1, b, h, beta_fric)
!  endif
!  normdiff=0.
!  normZ=0.
!  diff = 0.
!  k=0
!end subroutine phistage 
!   b) Encapsulate the 'plain' mode call to phi() with an if statement that checks if conv_flag remains FALSE. If yes, then
!      call phi() as usual and pefrom convergene criteria computation that will set conv_flag to TRUE upon convergence.
!   c) Encapsulate the 'adjoint' mode call to phi() when isinloop is 1 with 
!      an if statement that checks if conv_flag remains FALSE. If yes, then
!      call phi() as usual for this caseand pefrom convergene criteria computation
!      that will set adj_conv_flag to TRUE upon convergence.

subroutine template()
  use OAD_tape
  use OAD_rev
  use OAD_cp
!$TEMPLATE_PRAGMA_DECLARATIONS
  type(modeType) :: our_orig_mode

! checkpointing stacks and offsets
  integer :: cp_loop_variable_1,cp_loop_variable_2,cp_loop_variable_3,cp_loop_&
   &variable_4,cp_loop_variable_5,cp_loop_variable_6

! floats 'F'
  double precision, dimension(:), allocatable, save :: theArgFStack
  integer, save :: theArgFStackoffset=0, theArgFStackSize=0, theArgFStackoffsetTemp=0
! integers 'I'
  integer, dimension(:), allocatable, save :: theArgIStack
  integer, save :: theArgIStackoffset=0, theArgIStackSize=0, theArgIStackoffsetTemp=0
! booleans 'B'
  logical, dimension(:), allocatable, save :: theArgBStack
  integer, save :: theArgBStackoffset=0, theArgBStackSize=0, theArgBStackoffsetTemp=0
! strings 'S'
  character*(80), dimension(:), allocatable, save :: theArgSStack
  integer, save :: theArgSStackoffset=0, theArgSStackSize=0, theArgSStackoffsetTemp=0

! external C function used in inlined code
  integer iaddr
  external iaddr
!<------------------Begin user declarations ---------------------->!
! Insert declarations of dummy variables for calling adjoint computation 
! without side effects, and storing adjoint variable iterates
  type(active) :: B_DUMMY(1:79)
  type(active) :: H_DUMMY(1:79)
  type(active) :: U_DUMMY(1:80)
! Insert a declaration to hold the adjoint of U, between successive calls
  double precision, save :: U_0_d(1:80)
! Optionally insert a declaration to hold the value of U_IP1 before for convergence testing
  type(active) :: U_IP1_PRE(1:80)
!<------------------End user declarations ------------------------>!
  if (our_rev_mode%arg_store) then
    if(isinloop.eq.2) then
!<------------------Begin user checkpoint store ------------------>!
! Use the checkpointing methods given in OAD_cp to checkpoint the arguments
! to the staging function
      call cp_store_real_scalar(DX,theArgFStack,theArgFStackoffset,theArgFStackSize)
      call cp_store_real_vector(U,size(U),theArgFStack,theArgFStackoffset,theArgFStackSize)
      call cp_store_real_vector(U_IP1,size(U_IP1),theArgFStack,theArgFStackoffset,theArgFStackSize)
      call cp_store_real_vector(B,size(B),theArgFStack,theArgFStackoffset,theArgFStackSize)
      call cp_store_real_vector(H,size(H),theArgFStack,theArgFStackoffset,theArgFStackSize)
      call cp_store_p_real_vector(BETA_FRIC,size(BETA_FRIC),theArgFStack,theArgFStackoffset,theArgFStackSize)
    end if
!<------------------End user checkpoint store -------------------->!
  end if
  if (our_rev_mode%arg_restore) then
!<------------------Begin DO NOT CHANGE THIS CODE----------------->!
!Set the Temporary offset to the actual offset
    theArgFStackoffsetTemp = theArgFStackoffset
    theArgIStackoffsetTemp = theArgIStackoffset
    theArgBStackoffsetTemp = theArgBStackoffset
    theArgSStackoffsetTemp = theArgSStackoffset
!<------------------End DO NOT CHANGE THIS CODE------------------->!
!<------------------Begin user checkpoint restore ---------------->!
! Use the checkpointing methods given in OAD_cp to restore the arguments
! of the staging function
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
!<------------------End user checkpoint restore ------------------>!
    end if
    if (our_rev_mode%plain) then
      if(isinloop.eq.0) then
        CONV_FLAG = .false.
        iter = 0
      endif
      our_orig_mode=our_rev_mode
      our_rev_mode%arg_store=.FALSE.
      if(isinloop.eq.1) then
!!<------------------Begin user function call --------------------->!
!! Comment this call out if convergence criteria is used
!!          CALL phi(U,U_IP1,B,H,BETA_FRIC)
!!<------------------End user function call ----------------------->!
!<------------------Begin USER FORWARD CONVERGENCE CRITERIA------->!
! Comment out the conditional statement and its body if convergence
! criteria is not used
        if (CONV_FLAG.eqv..FALSE.) THEN
!<------------------Begin user function call --------------------->!
          CALL phi(U,U_IP1,B,H,BETA_FRIC)
!<------------------End user function call ----------------------->!
          iter = iter + 1
          normDiff = 0.0
          normZ = 0.0
          do k=1,n+1
            diff = abs(u_ip1(k)%v-u(k)%v)
            if (diff.gt.normDiff) then
              normdiff=diff
            endif
            if (abs(u(k)%v).gt.normZ) then
              normZ=abs(u(k)%v)
            endif
          enddo
          if (normDIff/(normZ + 1.0).le.tol) then
            conv_flag = .TRUE.
            !print *, "forward converged i=", iter
          endif     
        end if
!<------------------End USER FORWARD CONVERGENCE CRITERIA--------->!
      end if
      our_rev_mode=our_orig_mode
    end if
    if (our_rev_mode%tape) then
      our_rev_mode%arg_store=.TRUE.
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%plain=.TRUE.
      our_rev_mode%tape=.FALSE.
      our_rev_mode%adjoint=.FALSE.
!<------------------Begin user function call --------------------->!
      CALL phi(U,U_IP1,B,H,BETA_FRIC)
!<------------------End user function call ----------------------->!
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%plain=.FALSE.
      our_rev_mode%tape=.FALSE.
      our_rev_mode%adjoint=.TRUE.
    end if
    if (our_rev_mode%adjoint) then
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.TRUE.
      our_rev_mode%plain=.FALSE.
      our_rev_mode%tape=.TRUE.
      our_rev_mode%adjoint=.FALSE.
      if(isinloop.eq.0) then
        ADJ_CONV_FLAG = .false.
        adj_iter = 0
!<------------------Begin user function call --------------------->!
        CALL phi(U,U_IP1,B,H,BETA_FRIC)
!<------------------End user function call ----------------------->!
      end if 
      if(isinloop.eq.1) then
!!<------------------Begin USER ADJOINT W/O CONVERGENCE CRITERIA--->!
!! Comment the following code out if convergence criteria is used
!!<------------------Begin user initialize dummy args ------------->!
!        B_DUMMY = B
!        H_DUMMY = H
!        U_DUMMY%d = U_0_d
!!<------------------End user initialize dummy args --------------->!
!        U_IP1_PRE = U_IP1
!!<------------------Begin user function call WITH DUMMY ARGS ----->!
!        CALL phi(U_DUMMY,U_IP1,B_DUMMY,H_DUMMY,BETA_FRIC)
!!<------------------End user function call WITH DUMMY ARGS ------->!
!!<------------------Begin user copy dummy value ------------------>!
!        U_IP1%d = U_DUMMY%d
!!<------------------End user copy dummy value -------------------->!         
!!<------------------End USER ADJOINT W/O CONVERGENCE CRITERIA----->!

!<------------------Begin USER ADJOINT CONVERGENCE CRITERIA------->!
! Comment out the conditional statement and its body if convergence
! criteria is not used
        if (ADJ_CONV_FLAG.eqv..FALSE.) THEN
!<------------------Begin user initialize dummy args ------------->!
          B_DUMMY = B
          H_DUMMY = H
          U_DUMMY%d = U_0_d
!<------------------End user initialize dummy args --------------->!
          U_IP1_PRE = U_IP1
!<------------------Begin user function call WITH DUMMY ARGS ----->!
          CALL phi(U_DUMMY,U_IP1,B_DUMMY,H_DUMMY,BETA_FRIC)
!<------------------End user function call WITH DUMMY ARGS ------->!
!<------------------Begin user copy dummy value ------------------>!
          U_IP1%d = U_DUMMY%d
!<------------------End user copy dummy value -------------------->!
          adj_iter = adj_iter + 1
          normDiff = 0.0
          normZ = 0.0
          do k=1,n+1
            diff = abs(u_ip1(k)%d-u_ip1_pre(k)%d)
            if (diff.gt.normDiff) then
              normdiff=diff
            end if
            if (abs(u_ip1_pre(k)%d).gt.normZ) then
              normZ=abs(u_ip1_pre(k)%d)
            end if
          enddo
          if (normDIff/(normZ + 1.0).le.adjtol) then
            adj_conv_flag = .TRUE.
            !print *, "adjoint converged i=", adj_iter
          end if          
        end if
!<------------------End USER ADJOINT CONVERGENCE CRITERIA--------->!
      end if
!We have to push the value of U%d as it is overwritten in the caller (stream_vel in this case)
      if(isinloop.eq.2) then
!<------------------Begin user function call --------------------->!
        CALL phi(U,U_IP1,B,H,BETA_FRIC)
!<------------------End user function call ----------------------->!
!<------------------Begin user save value ------------------------>!
        U_0_d = U%d
!<------------------End user save value -------------------------->!
      end if
      if(isinloop.ne.0) then
        !Set the Temporary offset to the actual offset
        theArgFStackoffset = theArgFStackoffsetTemp
        theArgIStackoffset = theArgIStackoffsetTemp
        theArgBStackoffset = theArgBStackoffsetTemp
        theArgSStackoffset = theArgSStackoffsetTemp
      end if
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.TRUE.
      our_rev_mode%plain=.FALSE.
      our_rev_mode%tape=.TRUE.
      our_rev_mode%adjoint=.FALSE.
    end if
end subroutine template
