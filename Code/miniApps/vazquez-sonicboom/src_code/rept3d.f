
      SUBROUTINE REPT3D
c---------------------------------------------------------------------   
c 
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
c     Local variables definition
      INTEGER jt, is, in, mi
c
      DO is=1,ns
         nbvoi(is)                 = 0
      ENDDO
c
      DO jt=1,nt
c
         DO in=1,4
c
            mi                     = nu(in,jt)
c
            IF (mi .EQ. 0) THEN
               WRITE(6, *) 'Error in REPT3D'
               CALL TILT
            ENDIF
c
            nbvoi(mi)              = nbvoi(mi) + 1
c
            ivoi(mi,nbvoi(mi))     = jt
c
            IF (nbvoi(mi) .GT. nvmax) THEN 
               WRITE(6, *) 'Error in REPT3D'
               WRITE(6, *) 'Increase NVMAX'
               CALL TILT
            ENDIF
c
         ENDDO
c
      ENDDO
c
      RETURN
      END
