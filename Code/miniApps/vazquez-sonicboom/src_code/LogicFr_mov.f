      subroutine LOGICFR_mov(nsfac,logfac,logfr,ns,nfac)

      implicit none

      include 'param.h'

      integer is, ns, nfac, is1, is2, is3, ifac
      integer nsfac(3,nfacmax), logfac(nfacmax), logfr(nsmax)

      
      do is=1,ns
        logfr(is)=0
      enddo


      do ifac=1,nfac
c
         is1                      = nsfac(1,ifac)
         is2                      = nsfac(2,ifac)
         is3                      = nsfac(3,ifac)
c
	 IF (logfr(is1) .GE. 0) 
     &   logfr(is1)            = logfac(ifac)
         IF (logfr(is2) .GE. 0)
     &      logfr(is2)            = logfac(ifac)
         IF (logfr(is3) .GE. 0) 
     &      logfr(is3)            = logfac(ifac)
	 IF(ABS(logfac(ifac)).EQ.3
     .     .or. ABS(logfac(ifac)).EQ.2) THEN                !ALWAYS SLIPPING NODES
           IF (logfr(is1) .LT. 0) 
     &        logfr(is1)          = -3
           IF (logfr(is2) .LT. 0) 
     &        logfr(is2)          = -3
           IF (logfr(is3) .LT. 0) 
     &        logfr(is3)          = -3
         ENDIF

      enddo

      return
      end
