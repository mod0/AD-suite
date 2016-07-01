
      SUBROUTINE SEG3D
c---------------------------------------------------------------------   
c Constructs the mesh edges
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------
c     Global variables definition   
      INTEGER itrs, itrt
      COMMON/comtri/itrs, itrt
c
c     Loval variables definition
      INTEGER is, iseg , jt, kv, k, is1, is2, kv1
      INTEGER i1, iseg0, nub1  , nub2  , kv2, i2, isg, isegis
      INTEGER nor(6)   , nex(6)
      INTEGER nobu(2,nsgmax)
c
c     Initializations
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
      DO 2 iseg=1,nsgmax
	 mark(iseg)                = 0
2     CONTINUE
c
      DO 10 is=1,ns
c
	 ndeg(is)                  = 0
c
         DO 5 kv=1,ndmax
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
            DO 100 kv1=1,ndmax
c
	       IF (kv1 .EQ. ndmax) THEN
		  PRINT*,'Increase NDMAX'
		  CALL TILT
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
            DO 120 kv2=1,ndmax
c
               IF (kv2 .EQ. ndmax) THEN
		  PRINT*,'Increase NDMAX'
		  CALL TILT
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
               IF (nseg .GT. nsgmax) THEN
                  PRINT*,'Increase NSGMAX'
                  CALL TILT
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
C [llh]      WRITE(6, *) 'Nombre total de segments: ', nseg
C [llh]      WRITE(6, *) ' '
C [llh]      WRITE(6, *) ' '
c
      IF (itrs .EQ. 1) THEN
c
         DO 3000 iseg=1,nseg
c
            is1                    = nubo(1,iseg)
            is2                    = nubo(2,iseg)
            ndeg(is1)              = ndeg(is1) + 1
            ndeg(is2)              = ndeg(is2) + 1
c
3000     CONTINUE
c
         CALL TRINS(ns, inew, ndeg, nsmax)
c
         isg                       = 0
c
         DO 5000 is=1,ns
c
            DO 5001 kv=1,ndeg(inew(is))
c
               isegis                   = jaret(inew(is),kv)
c
               IF (mark(isegis) .EQ. 0) THEN
c
                  isg                   = isg + 1
c
                  nobu(1,isg)           = inew(is)
c
                  IF (nubo(1,isegis) .EQ. inew(is)) THEN
	             nobu(2,isg)        = nubo(2,isegis)
		  ELSE
                     IF (nubo(2,isegis) .EQ. inew(is)) THEN
     		        nobu(2,isg)     = nubo(1,isegis)
                     ELSE
                        WRITE(6, *) 'Error in SEG3D'
                        CALL TILT
                     ENDIF
	          ENDIF
c
                  mark(isegis)          = 1
c
               ENDIF
c
5001        CONTINUE
c
5000     CONTINUE
c
         DO 6000 iseg=1,nseg
            nubo(1,iseg)           = nobu(1,iseg)
	    nubo(2,iseg)           = nobu(2,iseg)
6000     CONTINUE
c
         IF (isg .NE. nseg) THEN
            PRINT*,'Error in SEG3D'
            CALL TILT
         ENDIF
c
      ENDIF
c
      RETURN
      END
