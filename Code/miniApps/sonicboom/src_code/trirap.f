

      SUBROUTINE TRIRAP
c---------------------------------------------------------------------   
c Reorders the mesh vertices
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
      INTEGER p, il, ir, m, i
      INTEGER stack(0:50)
c
      il                           = 1
      ir                           = ns
      p                            = 2
      m                            = 20
c
100   CONTINUE
c         
      IF ((ir - il) .GT. m) THEN
c
         CALL PARTIT(i, il, ir)
c
         IF ((i - il) .GT. (ir - i)) THEN
            stack(p)               = il
            stack(p+1)             = i - 1 
            il                     = i + 1
         ELSE
            stack(p)               = i + 1
            stack(p+1)             = ir
            ir                     = i - 1
         ENDIF
c
         p                         = p + 2
c
      ELSE
c
         CALL INSERT(il, ir)
c
         p                         = p - 2
         il                        = stack(p)
         ir                        = stack(p+1)
c
      ENDIF
c
      IF (p .NE. 0) GOTO 100
c
      RETURN
      END
