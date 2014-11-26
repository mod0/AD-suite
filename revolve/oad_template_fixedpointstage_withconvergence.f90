subroutine template()
  use OAD_tape
  use OAD_rev
  use OAD_cp
!$TEMPLATE_PRAGMA_DECLARATIONS
  type(modeType) :: our_orig_mode

!
!     **** Global Variables & Derived Type Definitions ****
!
integer myi
! Temporaries to hold the stack pointers
  integer temp_double_tape_pointer,     &
&         temp_integer_tape_pointer,    &
&         temp_logical_tape_pointer,    &
&         temp_character_tape_pointer,  & 
&         temp_stringlength_tape_pointer 
! external C function used in inlined code
  integer iaddr
  external iaddr
!<------------------Begin user declarations ---------------------->!
! Insert declarations of dummy variables for calling adjoint computation 
! without side effects, and storing adjoint variable iterates
  type(active) :: B_DUMMY(1:79)
  type(active) :: H_DUMMY(1:79)
  type(active) :: U_DUMMY(1:80)
  type(active) :: U_TEMP(1:80)
  double precision, save :: U_0_d(1:80)
  type(active) :: U_IP1_PRE(1:80)
!<------------------End user declarations ------------------------>!
    if (our_rev_mode%plain) then
      our_orig_mode=our_rev_mode
! original function
      if(isinloop.eq.0) then
        CONV_FLAG = .false.
        ADJ_CONV_FLAG = .false.
        iter = 0
        adj_iter = 0
      endif
      if(isinloop.ne.0) then
        if (CONV_FLAG.eqv..FALSE.) THEN
          CALL phi(U,U_IP1,B,H,BETA_FRIC)
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
          endif
        end if
      end if
      our_rev_mode=our_orig_mode
    end if
    if (our_rev_mode%tape) then
      !If isinloop==0, then call phi(..) in plain mode once
      !If isinloop==1, then call phi(..) in plain mode till convergence, then call phi(..) once in tape mode to inject partials onto the stack correctly
      !If isinloop==2, then call phi(..) in tape mode once 
      our_orig_mode=our_rev_mode
      if(isinloop.eq.0) then
        CONV_FLAG = .false.
        ADJ_CONV_FLAG = .false.
        iter = 0
        adj_iter = 0      
        !print *, "FW:0  ", double_tape_pointer, integer_tape_pointer, our_rev_mode%plain, our_rev_mode%tape, our_rev_mode%adjoint
      end if
      if(isinloop.eq.1) then
        IF (ITER .EQ. 0) then
            !print *, "FW:1-1", double_tape_pointer, integer_tape_pointer, our_rev_mode%plain, our_rev_mode%tape, our_rev_mode%adjoint
        ENDIF
        IF (.not. (CONV_FLAG)) THEN
          ITER = (ITER+1)
          our_rev_mode%plain=.TRUE.
          our_rev_mode%tape=.FALSE.
          our_rev_mode%adjoint=.FALSE.
          CALL phi(U,U_IP1,B,H,BETA_FRIC)
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
          endif 
          if (conv_flag.eqv..true. .OR. iter.eq.n_nl) then
            our_rev_mode%plain=.FALSE.
            our_rev_mode%tape=.TRUE.
            our_rev_mode%adjoint=.FALSE.
            CALL phi(U,U_IP1,B,H,BETA_FRIC)
            !print *, "FW:1-2", double_tape_pointer, integer_tape_pointer, our_rev_mode%plain, our_rev_mode%tape, our_rev_mode%adjoint
          endif 
        ENDIF
      endif
      if(isinloop.eq.0) then
        our_rev_mode%plain=.TRUE.
        our_rev_mode%tape=.FALSE.
        our_rev_mode%adjoint=.FALSE.
        CALL phi(U,U_IP1,B,H,BETA_FRIC)
      endif 
      if(isinloop.eq.2 ) then
        !print *, "FW:2-1", double_tape_pointer, integer_tape_pointer, our_rev_mode%plain, our_rev_mode%tape, our_rev_mode%adjoint
        CALL phi(U,U_IP1,B,H,BETA_FRIC)
        !print *, "FW:2-2", double_tape_pointer, integer_tape_pointer, our_rev_mode%plain, our_rev_mode%tape, our_rev_mode%adjoint
      endif 
      our_rev_mode=our_orig_mode
    end if
    if (our_rev_mode%adjoint) then
      !If isinloop==0, then call phi(..) in once
      !If isinloop==1, then call phi(..) in till convergence, using the same partials on the stack each time
      !If isinloop==2, then call phi(..) in once 
      our_orig_mode=our_rev_mode
      if(isinloop.eq.2) then
        CONV_FLAG = .false.
        ADJ_CONV_FLAG = .false.
        iter = 0
        adj_iter = 0      
      end if
      if(isinloop.eq.0) then
        do myi=n+1,1,-1
          call pop_s0(U_DUMMY(myi)%d)
        end do
! adjoint
        !print *, "BW:1  ", double_tape_pointer, integer_tape_pointer, our_rev_mode%plain, our_rev_mode%tape, our_rev_mode%adjoint
        CALL phi(U,U_IP1,B,H,BETA_FRIC)
        !print *, "BW:0  ", double_tape_pointer, integer_tape_pointer, our_rev_mode%plain, our_rev_mode%tape, our_rev_mode%adjoint
      end if
      if(isinloop.eq.1) then
        if(ADJ_CONV_FLAG.eqv..false.) then
          adj_iter = adj_iter + 1
          do myi=n+1,1,-1
            call pop_s0(U_dummy(myi)%d)
            U_temp(myi)%d = U_dummy(myi)%d
          end do
          B_DUMMY = B
          H_DUMMY = H
          U_IP1_PRE = U_IP1
          !Store the stack pointers
          temp_double_tape_pointer = double_tape_pointer   
          temp_integer_tape_pointer = integer_tape_pointer     
          temp_logical_tape_pointer = logical_tape_pointer    
          temp_character_tape_pointer = character_tape_pointer 
          temp_stringlength_tape_pointer = stringlength_tape_pointer
          CALL phi(U_dummy,U_IP1,B_DUMMY,H_DUMMY,BETA_FRIC)
          U_IP1%d = U_dummy%d
          normDiff = 0.0
          normZ = 0.0
          do k=1,n+1
            diff = abs(u_ip1(k)%d-u_ip1_pre(k)%d)
            if (diff.gt.normDiff) then
              normdiff=diff
            endif
            if (abs(u_ip1_pre(k)%d).gt.normZ) then
              normZ=abs(u_ip1_pre(k)%d)
            endif
          enddo
          if (normDIff/(normZ + 1.0).le.adjtol) then
            adj_conv_flag = .true.
          endif 
          !Retore the stack pointers
          double_tape_pointer = temp_double_tape_pointer   
          integer_tape_pointer = temp_integer_tape_pointer     
          logical_tape_pointer = temp_logical_tape_pointer    
          character_tape_pointer = temp_character_tape_pointer 
          stringlength_tape_pointer = temp_stringlength_tape_pointer
          do myi=1,n+1
            call push_s0(U_temp(myi)%d)
          end do
        end if
      end if
      if(isinloop.eq.2) then
        !print *, "BW:2-2", double_tape_pointer, integer_tape_pointer, our_rev_mode%plain, our_rev_mode%tape, our_rev_mode%adjoint
! adjoint
        CALL phi(U,U_IP1,B,H,BETA_FRIC)
        !print *, "BW:2-1", double_tape_pointer, integer_tape_pointer, our_rev_mode%plain, our_rev_mode%tape, our_rev_mode%adjoint
        do myi=1,n+1
          call push_s0(U(myi)%d)
        end do
      end if
      our_rev_mode=our_orig_mode
    end if
end subroutine template

