
      SUBROUTINE PARTIT(ires, il, ir)
c---------------------------------------------------------------------   
c Reorders inew(il:ir) such that 
c          inew(il:ires-1) : ndeg(inew(i)) > ndeg(inew(ires))
c          inew(ires+1:ir) : ndeg(inew(i)) < ndeg(inew(ires))
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
      INTEGER i, il  , ir, iaux, j, ires, iaux1
c
      iaux                         = ndeg(inew(ir))
      i                            = il - 1
      j                            = ir
c
10    CONTINUE
c
      i                            = i + 1
c
      IF (ndeg(inew(i)) .GT. iaux) GOTO 10
c
20    CONTINUE
c
      IF (j .EQ. il) GOTO 30
c
      j                            = j - 1
c
      IF (ndeg(inew(j)) .LT. iaux) GOTO 20 
c
30    CONTINUE
c
      iaux1                        = inew(i)
      inew(i)                      = inew(j)
      inew(j)                      = iaux1
c
      IF (j .GT. i) GOTO 10
c
      inew(j)                      = inew(i)
      inew(i)                      = inew(ir)
      inew(ir)                     = iaux1
      ires                         = i       
c
      RETURN
      END
