      SUBROUTINE CALDEG
c---------------------------------------------------------------------   
c Sorts the mesh vertices in order to optimize the vector 
c performances
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c--------------------------------------------------------------------- 
c     Local variables definition
      INTEGER jt, is, nuk, k, j, ktv
      INTEGER iel(4)
      INTEGER nut(4,ntmax)
c
c     Initializations
c
      DO 2 is=1,ns
c
         ndeg(is)                  = 0
c
         DO 1 k=1,ndmax
            jaret(is,k)            = 0
1        CONTINUE
c
2     CONTINUE
c
      DO 5 jt=1,nt
         mark(jt)                  = 0
         DO 4 k=1,4
            nut(k,jt)              = 0
            nuk                    = nu(k,jt)
            ndeg(nuk)              = ndeg(nuk) + 1
            jaret(nuk,ndeg(nuk))   = jt
4        CONTINUE
5     CONTINUE
c
      CALL TRINS(ns, inew, ndeg, nsmax)
c
c     Verification
c
      IF (ndeg(inew(1)) .GT. ndmax) THEN
         WRITE(6, *) 'Increase NDMAX : ', ndmax
         CALL TILT
      ENDIF
c
      jt                           = 0
c
      DO 1000 is=1,ns
         DO 500 j=1,ndeg(is)
            ktv                    = jaret(is,j)
            IF (ktv .NE. 0) THEN
               IF (mark(ktv) .EQ. 0) THEN
                  jt               = jt + 1
                  DO 100 k=1,4
                     nut(k,jt)     = inew(nu(k,ktv))
100               CONTINUE       
                  mark(ktv)        = 1
               ENDIF
            ENDIF
500      CONTINUE
1000  CONTINUE
c
      IF (jt .NE. nt) THEN
         WRITE(6, *) 'Error while renumeroting the tetraedras'
         CALL TILT
      ENDIF
c
      DO 2000 jt=1,nt
         DO 1500 k=1,4
            iel(k)                 = nut(k,jt)
            nut(k,jt)              = nu(k,jt)
            nu(k,jt)               = nut(k,jt)
1500     CONTINUE
2000  CONTINUE
c
      RETURN
      END
