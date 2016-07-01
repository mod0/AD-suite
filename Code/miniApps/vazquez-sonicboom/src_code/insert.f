      SUBROUTINE INSERT(il, ir)
c---------------------------------------------------------------------   
c Sorts array ndeg in the decreasing order of the degrees of 
c the mesh vertices i.e. ndeg(inew(il)) > ... > ndeg(inew(ir)
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
      INTEGER iaux , iaux1
      INTEGER i, il, ir, j
C
      DO 1 i=il+1,ir
c
         iaux                      = ndeg(inew(i))
         iaux1                     = inew(i)
         j                         = i
c
50       IF (ndeg(inew(j-1)) .GE. iaux) GOTO 100
c
         inew(j)                   = inew(j-1)
         j                         = j - 1
c
         IF (j .EQ. il) GOTO 100
c
         GOTO 50
c
100      CONTINUE
c
         inew(j)                   = iaux1
c
1     CONTINUE
c
      RETURN
      END
