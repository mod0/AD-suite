      SUBROUTINE TRINS
c---------------------------------------------------------------------   
c Reorders the mesh vertices
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
      INTEGER i
c
      DO 1 i=1,ns
         inew(i)                   = i
1     CONTINUE
c
      CALL TRIRAP
c
      RETURN
      END
