      SUBROUTINE GRADNOD
c---------------------------------------------------------------------   
c Computes the hermitian nodal gradients
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------
c     Local variables definition
      INTEGER ivar  , k     , jt   , ia  , ib
      INTEGER is    , is1   , is2  , is3 , is4  , id1 , id2
      INTEGER iseg , nub1, nub2 , js1 , js2
      REAL*8    usro  
      REAL*8    vol6  , ait   , ais  
      REAL*8    ds3   , us4  , us6
      REAL*8    x(4)  , y(4)  , z(4) , b(4), c(4) , d(4)
      REAL*8    dbx(5), dby(5), dbz(5)
      REAL*8    dxt(5), dyt(5), dzt(5)
      REAL*8    uph(5,4) 
      REAL*8    pentel(nsmax)
c
c     Initialisations
c
      ds3                          = 2.0d0/3.0d0
      us4                          = 1.0d0/4.0d0
      us6                          = 1.0d0/6.0d0
c
c
      DO 10 is=1,ns
         pentel(is)             = 0.0d0
10    CONTINUE
c
c
c
C$AD II-LOOP
      DO 1000 jt=1,nt
c     ===============
c
c
c
         DO 30 k=1,4
c        -----------
            is                     = nu(k,jt)
c
            uph(1,k)               = ua(1,is)
            usro                   = 1.0d0/uph(1,k)
c
            uph(2,k)               = ua(2,is)*usro
            uph(3,k)               = ua(3,is)*usro
            uph(4,k)               = ua(4,is)*usro
c
            uph(5,k)               = gam1*(ua(5,is) -
     &           0.5d0*uph(1,k)*(
     &           uph(2,k)*uph(2,k) +
     &           uph(3,k)*uph(3,k) +
     &           uph(4,k)*uph(4,k)))
c
30       CONTINUE
c        --------
c
c
         is1                       = nu(1,jt)
         is2                       = nu(2,jt)
         is3                       = nu(3,jt)
         is4                       = nu(4,jt)
c
         x(1)                      = coor(1,is1)
         y(1)                      = coor(2,is1)
         z(1)                      = coor(3,is1)
c
         x(2)                      = coor(1,is2)
         y(2)                      = coor(2,is2)
         z(2)                      = coor(3,is2)
c
         x(3)                      = coor(1,is3)
         y(3)                      = coor(2,is3)
         z(3)                      = coor(3,is3)
c
         x(4)                      = coor(1,is4)
         y(4)                      = coor(2,is4)
         z(4)                      = coor(3,is4)
c     
c
c        Computing the basis function gradients
c
         CALL GRADFB(x, y, z, b, c, d, vol6)
c
         ait                       = us6/volt(jt)
c
         dbx(1)                    = b(1)*ait
         dbx(2)                    = b(2)*ait
         dbx(3)                    = b(3)*ait
         dbx(4)                    = b(4)*ait
c
         dby(1)                    = c(1)*ait
         dby(2)                    = c(2)*ait
         dby(3)                    = c(3)*ait
         dby(4)                    = c(4)*ait
c
         dbz(1)                    = d(1)*ait
         dbz(2)                    = d(2)*ait
         dbz(3)                    = d(3)*ait
         dbz(4)                    = d(4)*ait
c
c
c        Computing the P1-gradients on each tetraedra
c
c
         DO 60 ivar=1,5
c        --------------
            dxt(ivar)              = uph(ivar,1)*dbx(1) +
     &                               uph(ivar,2)*dbx(2) +
     &                               uph(ivar,3)*dbx(3) +
     &                               uph(ivar,4)*dbx(4)
c
            dyt(ivar)              = uph(ivar,1)*dby(1) +
     &                               uph(ivar,2)*dby(2) +
     &                               uph(ivar,3)*dby(3) +
     &                               uph(ivar,4)*dby(4)
c
            dzt(ivar)              = uph(ivar,1)*dbz(1) +
     &                               uph(ivar,2)*dbz(2) +
     &                               uph(ivar,3)*dbz(3) +
     &                               uph(ivar,4)*dbz(4)
c
60       CONTINUE
c        --------
c
         IF (nordre .EQ. 1) GOTO 1110
         IF (nordre .EQ. 4) GOTO 1105
c
c        Computing the hermitian nodal gradients
c
         ait                       = us4*volt(jt)
c
         DO 95 ivar=1,5
c        --------------
            dx(ivar,is1)           = dx(ivar,is1) + dxt(ivar)*ait
            dx(ivar,is2)           = dx(ivar,is2) + dxt(ivar)*ait
            dx(ivar,is3)           = dx(ivar,is3) + dxt(ivar)*ait
            dx(ivar,is4)           = dx(ivar,is4) + dxt(ivar)*ait
c
            dy(ivar,is1)           = dy(ivar,is1) + dyt(ivar)*ait
            dy(ivar,is2)           = dy(ivar,is2) + dyt(ivar)*ait
            dy(ivar,is3)           = dy(ivar,is3) + dyt(ivar)*ait
            dy(ivar,is4)           = dy(ivar,is4) + dyt(ivar)*ait
c     
            dz(ivar,is1)           = dz(ivar,is1) + dzt(ivar)*ait
            dz(ivar,is2)           = dz(ivar,is2) + dzt(ivar)*ait
            dz(ivar,is3)           = dz(ivar,is3) + dzt(ivar)*ait
            dz(ivar,is4)           = dz(ivar,is4) + dzt(ivar)*ait
c
95       CONTINUE
c        --------
c
1105     CONTINUE    
c
         IF (nordre .EQ. 4) THEN
c
            DO 1090 k=1,4
c           -------------
               is                  = nu(k,jt)
c
               DO 1095 ivar=1,5
c
                  dx(ivar,is)      = (1.0d0 - pentel(is))*dxt(ivar) + 
     &            pentel(is)*0.5d0*(SIGN(1.0d0, dxt(ivar)) +
     &                            SIGN(1.0d0, dx(ivar,is)))*
     &            MIN(ABS(dxt(ivar)), ABS(dx(ivar,is)))
c
                  dy(ivar,is)      = (1.0d0 - pentel(is))*dyt(ivar) + 
     &            pentel(is)*0.5d0*(SIGN(1.0d0, dyt(ivar)) +
     &                            SIGN(1.0d0, dy(ivar,is)))*
     &            MIN(ABS(dyt(ivar)), ABS(dy(ivar,is)))
c
                  dz(ivar,is)      = (1.0d0 - pentel(is))*dzt(ivar) + 
     &            pentel(is)*0.5d0*(SIGN(1.0d0, dzt(ivar)) +
     &                            SIGN(1.0d0, dz(ivar,is)))*
     &            MIN(ABS(dzt(ivar)), ABS(dz(ivar,is)))
c
                  pentel(is)       = 1.0d0
c
1095           CONTINUE
c
1090        CONTINUE
c           --------
c
         ENDIF    !          IF (nordre .EQ. 4) 
c
1110     CONTINUE
c
1000  CONTINUE    ! FIN DE LA BOUCLE SUR LES ELEMENTS
c     ========
c
c
c
c
c
c     Completing the non-limited nodal gradients
c
      IF ((nordre .EQ. 2) .OR. (nordre .EQ. 3)) THEN 
c
        is = 0
C$AD II-LOOP
        DO 1040 ivar=1,5
c        ================
            DO 1050 is=1,ns
c           ---------------
               ais                 = 1.0d0/vols(is)
c
               dx(ivar,is)         = dx(ivar,is)*ais
               dy(ivar,is)         = dy(ivar,is)*ais
               dz(ivar,is)         = dz(ivar,is)*ais
c
1050        CONTINUE
c           --------
1040     CONTINUE
c        ========
c
      ENDIF !       IF ((nordre .EQ. 2) .OR. (nordre .EQ. 3)) THEN 
c
c
c
      RETURN
      END
