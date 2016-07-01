      SUBROUTINE SEG3D_mov(nu,ns,nt,nubo,nseg)
c---------------------------------------------------------------------   
c Constructs the mesh edges
c---------------------------------------------------------------------   
      implicit none

      include 'param.h'
c---------------------------------------------------------------------
c     Global variables definition   
      INTEGER itrs, itrt, ifixs, ifixt
      COMMON/comtri/itrs, itrt, ifixs, ifixt
c
c     Loval variables definition
      INTEGER is, jt, kv, k, is1, is2, kv1
      INTEGER i1, iseg0, nub1  , nub2  , kv2, i2
      INTEGER nor(6)   , nex(6)
      INTEGER ierr
      INTEGER nu(4,ntmax), nt, ns
      INTEGER nubo(2,nsegmax), nseg, ndeg(nsmax)
      INTEGER jaret(nsmax,ndegmax)
c
c     Initializations
c
      ierr                         = 0
c
      nex(1)                       = 1
      nex(2)                       = 1
      nex(3)                       = 1
      nex(4)                       = 2
      nex(5)                       = 2
      nex(6)                       = 3
      nor(1)                       = 2
      nor(2)                       = 3
      nor(3)                       = 4
      nor(4)                       = 3
      nor(5)                       = 4
      nor(6)                       = 4
c
      DO 10 is=1,ns
c
	 ndeg(is)                  = 0
c
         DO 5 kv=1,ndegmax
            jaret(is,kv)           = 0
5        CONTINUE
c
10    CONTINUE
c
      nseg                         = 0
c
      DO 1000 jt=1,nt
c
         DO 500 k=1,6
c
            is1                    = nu(nor(k),jt)
            is2                    = nu(nex(k),jt)
c
            DO 100 kv1=1,ndegmax
c
	       IF (kv1 .EQ. ndegmax) THEN 
                  ierr             = 1
                  GOTO 2000
               ENDIF
c
               IF (jaret(is1,kv1) .EQ. 0) THEN
		  i1               = 0
		  GOTO 110
	       ENDIF
c
               iseg0               = jaret(is1,kv1)
	       nub1                = nubo(1,iseg0)
	       nub2                = nubo(2,iseg0)
c
               IF ((is1 .EQ. nub1) .AND. (is2 .EQ. nub2)) THEN
		  i1               = 1
		  GOTO 200
	       ENDIF
c
	       IF ((is2 .EQ. nub1) .AND. (is1 .EQ. nub2)) THEN
		  i1               = 1
		  GOTO 200
	       ENDIF
c
100         CONTINUE
c
110         CONTINUE
c
            DO 120 kv2=1,ndegmax
c
               IF (kv2 .EQ. ndegmax) THEN
                  ierr             = 1
                  GOTO 2000
               ENDIF
c
               IF (jaret(is2,kv2) .EQ.0) THEN
		  i2               = 0
		  GOTO 130
	       ENDIF
c
	       iseg0               = jaret(is2,kv2)
	       nub1                = nubo(1,iseg0)
 	       nub2                = nubo(2,iseg0)
c
               IF ((is1 .EQ. nub1) .AND. (is2 .EQ. nub2)) THEN
		  i2               = 1
		  GOTO 200
	       ENDIF
c
	       IF ((is2 .EQ. nub1) .AND. (is1 .EQ. nub2)) THEN
		  i2               = 1
		  GOTO 200
	       ENDIF
c
120         CONTINUE
c
130         CONTINUE
c
            IF ((i1 .EQ. 0) .AND. (i2 .EQ. 0)) THEN
c
               nseg                = nseg + 1
c
               IF (nseg .GT. nsegmax) THEN
                  ierr             = 1
                  GOTO 2000
               ENDIF
c
	       jaret(is1,kv1)      = nseg
               jaret(is2,kv2)      = nseg
               nubo(1,nseg)        = is1
               nubo(2,nseg)        = is2
c
	    ENDIF
c
200         CONTINUE
c
500      CONTINUE
c
1000  CONTINUE
c
2000  CONTINUE
c  
      IF (ierr .NE. 0) THEN 
c
            WRITE(6, *) 'Error in SEG3D'
            WRITE(6, *) ' Ndegmax  = ', ndegmax, 
     &                  ' Kv1    = ', kv1,
     &                  ' Kv2    = ', kv2
            WRITE(6, *) ' Nsegmax = ', nsegmax, 
     &                  ' Nseg   = ', nseg
c
      ENDIF
c 
      RETURN
      END





