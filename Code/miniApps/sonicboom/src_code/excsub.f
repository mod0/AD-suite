      subroutine cputim(rtime)
c***********************************************************************
c
c**** This routine finds out the CPU time in seconds
c
c***********************************************************************      

      implicit none
      real*8 rtime
      real*4 ETIME,stime(2),result
      
      result = ETIME(stime)
      rtime  = stime(1)

      end

cccccccccccccccccccccc
cccccccccccccccccccccc
cccccccccccccccccccccc

      
      subroutine flunow(lunit)

      integer*4 lunit

      call FLUSH(lunit)
      

      end

cccccccccccccccccccccc
cccccccccccccccccccccc
cccccccccccccccccccccc


      
c            subroutine admemo(ioptn,iaddr,lbyts)
c***********************************************************************
c
c**** This routine allocates and releases virtual memory
c
c     If it is compiled with -64 option it is mandatory to use also 
c     -i8 as compilation option. If it is compiled without -64 option
c     any size of integer is valid, but then iaddr must be stored in
c     an integer*4, since `free' needs a 32-bits integer.       
c
c***********************************************************************
c      implicit none
c      integer    iaddr                                      ! pointer
c
c     Modified for LINUX arch
c
c      integer*4  iaddr
c      integer*4  memall
c      integer*4  ioptn,lbyts
c
c      if (ioptn.eq.0) then                                  ! Memory allocation
c        iaddr=memall(lbyts)                                   
c        if (iaddr.eq.0. or .lbyts.lt.0)                      
c     .    call runend('SADMEM: ERROR WHEN CALLING MALLOC  ')
c      else if (ioptn.eq.1) then                             ! To be used for 
c        call runend('SADMEM: TASK NOT READY ')              ! re-allocation
c      else if (ioptn.eq.2) then                             ! Release of memory
c        call memdea(iaddr)                                 
c      end if                                                
c   
c      end


      

