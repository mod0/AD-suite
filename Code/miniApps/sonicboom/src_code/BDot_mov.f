      FUNCTION BDOT_mov(v1x,v1y,v1z,v2x,v2y,v2z,logfr,mvlgfr,ns)
c
      implicit none
c
      include 'param.h'
c
      REAL*8 v1x(*), v1y(*), v1z(*), v2x(*), v2y(*),v2z(*)
      INTEGER is,ns
      INTEGER logfr(nsmax), mvlgfr(nsmax)
      REAL*8 BDOT_mov
C
      bdot_mov = 0.
c
      DO is=1,ns
c
         IF ((logfr(is) .EQ. 0) .OR. (mvlgfr(is) .GT. 0)) THEN
c
            bdot_mov = bdot_mov + v1x(is)*v2x(is)
            bdot_mov = bdot_mov + v1y(is)*v2y(is)
            bdot_mov = bdot_mov + v1z(is)*v2z(is)
c
         ENDIF
c
      ENDDO
c
      RETURN
      END


