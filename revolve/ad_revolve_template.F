#include "PACKAGES_CONFIG.h"

      subroutine template()
      use OAD_revolve_cp
      use OAD_cp
      use OAD_tape
      use OAD_rev
      use revolve

      
!$TEMPLATE_PRAGMA_DECLARATIONS
      LOGICAL :: initialized=.FALSE.
      TYPE(rvAction),save :: theAction
      CHARACTER(80) :: errorMsg
      integer, save :: jointCPCount
      integer, save :: currIter
      integer :: cp_loop_variable_1,cp_loop_variable_2,
     +     cp_loop_variable_3,cp_loop_variable_4,cp_loop_variable_5

      double precision, dimension(:), allocatable, save :: theArgFStack
      integer, save :: theArgFStackoffset=0, theArgFStackSize=0
      integer, dimension(:), allocatable, save :: theArgIStack
      integer, save :: theArgIStackoffset=0, theArgIStackSize=0
      type(modeType) :: our_orig_mode

      integer iaddr
      external iaddr
!      nIter0 = NINT( (startTime-baseTime)/deltaTClock )
      if (our_rev_mode%arg_store) then 
	call cp_write_open()
!$PLACEHOLDER_PRAGMA$ id=8
       call cp_close()
      end if 
      if (our_rev_mode%arg_restore) then
	call cp_read_open()
!$PLACEHOLDER_PRAGMA$ id=9
	call cp_close()
      end if
      if (our_rev_mode%plain) then
        do nt=1,n_timesteps
        CALL OpenAD_forward_step (h, u, beta_fric)
      enddo
      end if
      if (our_rev_mode%tape) then
         currIter=0
         jointCPcount=cp_fNumber()
         initialized=rvInit(n_timesteps,120,
     +                      errorMsg,theAction)
         IF (.NOT.initialized) WRITE(*,'(A,A)') 'Error: ', errorMsg
         do while (theAction%actionFlag/=rvDone)
           theAction=rvNextAction()
           select case (theAction%actionFlag) 
             case (rvStore)
                call cp_write_open(theAction%cpNum+jointCPCount)
!$PLACEHOLDER_PRAGMA$ id=8
                call cp_close
             case (rvForward)
                call OAD_revPlain
                do currIter=currIter,theAction%iteration-1
                   CALL OpenAD_forward_step (h, u, beta_fric)
                end do
                call OAD_revTape
             case (rvFirstUTurn)
                CALL OpenAD_forward_step (h, u, beta_fric)
! get out now ...
                exit                
             end select   
          end do
      end if 
      if (our_rev_mode%adjoint) then
        IF (.NOT.initialized) WRITE(*,'(A)') 'Error: not initialized'
        do while (theAction%actionFlag/=rvDone)
           select case (theAction%actionFlag) 
             case (rvFirstUTurn)
!we taped already ... see above
                CALL OpenAD_forward_step (h, u, beta_fric)
             case (rvStore)
                call cp_write_open(theAction%cpNum+jointCPCount)
!$PLACEHOLDER_PRAGMA$ id=8
                call cp_close
             case (rvRestore)
                call cp_read_open(theAction%cpNum+jointCPCount)
!$PLACEHOLDER_PRAGMA$ id=9
                currIter=theAction%iteration
                call cp_close
             case (rvForward)
                call OAD_revPlain
                do currIter=currIter,theAction%iteration-1
                   CALL OpenAD_forward_step (h, u, beta_fric)
                end do
                call OAD_revAdjoint
             case (rvUTurn)
                call OAD_revTape
                CALL OpenAD_forward_step (h, u, beta_fric)
                call OAD_revAdjoint
                CALL OpenAD_forward_step (h, u, beta_fric)
           end select   
           theAction=rvNextAction()
        end do
      end if 
      end subroutine template
