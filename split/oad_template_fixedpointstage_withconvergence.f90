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
      our_orig_mode=our_rev_mode
      if(isinloop.eq.0) then
        CONV_FLAG = .false.
        ADJ_CONV_FLAG = .false.
        iter = 0
        adj_iter = 0       
      end if
      if(isinloop.eq.1) then
        IF (.not. (CONV_FLAG)) THEN
          ITER = (ITER+1)
          !Store the stack pointers
          temp_double_tape_pointer = double_tape_pointer   
          temp_integer_tape_pointer = integer_tape_pointer     
          temp_logical_tape_pointer = logical_tape_pointer    
          temp_character_tape_pointer = character_tape_pointer 
          temp_stringlength_tape_pointer = stringlength_tape_pointer
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
          if (conv_flag.neqv..true. .AND. iter.ne.n_nl) then
            !Restore the stack pointers
            double_tape_pointer = temp_double_tape_pointer   
            integer_tape_pointer = temp_integer_tape_pointer     
            logical_tape_pointer = temp_logical_tape_pointer    
            character_tape_pointer = temp_character_tape_pointer 
            stringlength_tape_pointer = temp_stringlength_tape_pointer
          endif 
        ENDIF
      endif
      if(isinloop.eq.0) then
        !Store the stack pointers
        temp_double_tape_pointer = double_tape_pointer   
        temp_integer_tape_pointer = integer_tape_pointer     
        temp_logical_tape_pointer = logical_tape_pointer    
        temp_character_tape_pointer = character_tape_pointer 
        temp_stringlength_tape_pointer = stringlength_tape_pointer
        CALL phi(U,U_IP1,B,H,BETA_FRIC)
        !Restore the stack pointers
        double_tape_pointer = temp_double_tape_pointer   
        integer_tape_pointer = temp_integer_tape_pointer     
        logical_tape_pointer = temp_logical_tape_pointer    
        character_tape_pointer = temp_character_tape_pointer 
        stringlength_tape_pointer = temp_stringlength_tape_pointer
      endif 
      if(isinloop.eq.2 ) then
        CALL phi(U,U_IP1,B,H,BETA_FRIC)
      endif 
      our_rev_mode=our_orig_mode
    end if
    if (our_rev_mode%adjoint) then
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
        CALL phi(U,U_IP1,B,H,BETA_FRIC)
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
! adjoint
        CALL phi(U,U_IP1,B,H,BETA_FRIC)
        do myi=1,n+1
          call push_s0(U(myi)%d)
        end do
      end if
      our_rev_mode=our_orig_mode
    end if
end subroutine template

