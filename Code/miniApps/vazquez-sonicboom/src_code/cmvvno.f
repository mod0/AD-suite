      SUBROUTINE CMVVNO
c---------------------------------------------------------------------   
c Updates the normals to control/volumes boundaries
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------
c     Local variables definition  
      INTEGER nvrm
      PARAMETER (nvrm    = 6)
      INTEGER iseg  , jt, is1, is2, is3, is4, is
      INTEGER ind   , jc, jc3
      INTEGER nub1  , nub2   , k  , kv , in, i1
      INTEGER nor(6), nex(6) , nusg(6,2)
      REAL*8    eps3  , eta1   , eta2    , eta3
      REAL*8    et1   , et2    , et3     , esign
      REAL*8    us36dt
      REAL*8    x(4,2), y(4,2) , z(4,2)
      REAL*8    w(6,3), x12(3) , x13(3)
      REAL*8    x14(3), x23(3) , x24(3)
      REAL*8    x34(3), y12(3) , y24(3)  
      REAL*8    y13(3), y14(3) , y23(3)
      REAL*8    y34(3), z12(3) , z13(3)
      REAL*8    z14(3), z23(3) , z24(3), z34(3)
      REAL*8    dbijx(6,4)
      REAL*8    dbijy(6,4)
      REAL*8    dbijz(6,4)
c
c     Initializations
c
      nor(1)                       = 2 
      nex(1)                       = 1
      nusg(1,1)                    = 3
      nusg(1,2)                    = 4
c
      nor(2)                       = 3
      nex(2)                       = 1
      nusg(2,1)                    = 4
      nusg(2,2)                    = 2
c
      nor(3)                       = 4
      nex(3)                       = 1
      nusg(3,1)                    = 2
      nusg(3,2)                    = 3
c
      nor(4)                       = 3
      nex(4)                       = 2
      nusg(4,1)                    = 1
      nusg(4,2)                    = 4
c
      nor(5)                       = 4
      nex(5)                       = 2
      nusg(5,1)                    = 3
      nusg(5,2)                    = 1
c
      nor(6)                       = 4
      nex(6)                       = 3
      nusg(6,1)                    = 1
      nusg(6,2)                    = 2
c
      DO 1 is=1,ns
c
         coco(1,is)                = coor(1,is) 
         coco(2,is)                = coor(2,is) 
         coco(3,is)                = coor(3,is)
c
1     CONTINUE
c
      dt                           = 0.1
c
      DO 26 i1=1,4
         DO 25 jt=1,nt
            tvno(jt,i1,1)          = 0.0
            tvno(jt,i1,2)          = 0.0
            tvno(jt,i1,3)          = 0.0
25       CONTINUE
26    CONTINUE
c
      DO 10 iseg=1,nseg
         vnocl(1,iseg)             = 0.0
         vnocl(2,iseg)             = 0.0
         vnocl(3,iseg)             = 0.0
         sigma(iseg)               = 0.0
10    CONTINUE
c
      us36dt                       = 1.0/(dt*36.0)
      eps3                         = 1.0/(24.0*3.0)
      esign                        = 1.0
c
      DO 1000 jt=1,nt
c
         is1                       = nu(1,jt)
         is2                       = nu(2,jt)
         is3                       = nu(3,jt)
         is4                       = nu(4,jt)
c
         x(1,1)                    = coor(1,is1)
         x(2,1)                    = coor(1,is2)
         x(3,1)                    = coor(1,is3)
         x(4,1)                    = coor(1,is4)
c
         y(1,1)                    = coor(2,is1)
         y(2,1)                    = coor(2,is2)
         y(3,1)                    = coor(2,is3)
         y(4,1)                    = coor(2,is4)
c
         z(1,1)                    = coor(3,is1)
         z(2,1)                    = coor(3,is2)
         z(3,1)                    = coor(3,is3)
         z(4,1)                    = coor(3,is4)
c
         x(1,2)                    = coco(1,is1)
         x(2,2)                    = coco(1,is2)
         x(3,2)                    = coco(1,is3)
         x(4,2)                    = coco(1,is4)
c
         y(1,2)                    = coco(2,is1)
         y(2,2)                    = coco(2,is2)
         y(3,2)                    = coco(2,is3)
         y(4,2)                    = coco(2,is4)
c
         z(1,2)                    = coco(3,is1)
         z(2,2)                    = coco(3,is2)
         z(3,2)                    = coco(3,is3)
         z(4,2)                    = coco(3,is4)
c     
         DO 20 jc=1,2
c
            x12(jc)                = x(2,jc) - x(1,jc) 
            x13(jc)                = x(3,jc) - x(1,jc) 
            x14(jc)                = x(4,jc) - x(1,jc) 
            x23(jc)                = x(3,jc) - x(2,jc) 
            x24(jc)                = x(4,jc) - x(2,jc) 
            x34(jc)                = x(4,jc) - x(3,jc) 
c
            y12(jc)                = y(2,jc) - y(1,jc) 
            y13(jc)                = y(3,jc) - y(1,jc) 
            y14(jc)                = y(4,jc) - y(1,jc) 
            y23(jc)                = y(3,jc) - y(2,jc) 
            y24(jc)                = y(4,jc) - y(2,jc) 
            y34(jc)                = y(4,jc) - y(3,jc) 
c
            z12(jc)                = z(2,jc) - z(1,jc) 
            z13(jc)                = z(3,jc) - z(1,jc) 
            z14(jc)                = z(4,jc) - z(1,jc) 
            z23(jc)                = z(3,jc) - z(2,jc) 
            z24(jc)                = z(4,jc) - z(2,jc) 
            z34(jc)                = z(4,jc) - z(3,jc) 
c      
20       CONTINUE
c
         DO 30 in=1,4
c
            jc                     = in
            jc3                    = 3 - jc
c
            IF (in .EQ. 3) THEN
               jc                  = 2
               jc3                 = 2
            ENDIF
c
            IF (in .EQ. 4) THEN
               jc                  = 1 
               jc3                 = 1
            ENDIF
c
            dbijx(1,in)            = y34(jc)*(z14(jc3)  + z24(jc3)) -
     &                               z34(jc)*(y14(jc3)  + y24(jc3))
            dbijx(2,in)            = y24(jc)*(-z12(jc3) + z23(jc3)) - 
     &                               z24(jc)*(-y12(jc3) + y23(jc3))
            dbijx(3,in)            =-y23(jc)*(-z13(jc3) + z34(jc3)) + 
     &                               z23(jc)*(-y13(jc3) + y34(jc3))
            dbijx(4,in)            = y14(jc)*(z24(jc3)  + z34(jc3)) - 
     &                               z14(jc)*(y24(jc3)  + y34(jc3)) 
            dbijx(5,in)            = y13(jc)*(z12(jc3)  + z14(jc3)) - 
     &                               z13(jc)*(y12(jc3)  + y14(jc3)) 
            dbijx(6,in)            =-y12(jc)*(z23(jc3)  + z24(jc3)) +
     &                               z12(jc)*(y23(jc3)  + y24(jc3))
c
            dbijy(1,in)            = z34(jc)*(x14(jc3)  + x24(jc3)) - 
     &                               x34(jc)*(z14(jc3)  + z24(jc3))
            dbijy(2,in)            = z24(jc)*(-x12(jc3) + x23(jc3)) - 
     &                               x24(jc)*(-z12(jc3) + z23(jc3))
            dbijy(3,in)            =-z23(jc)*(-x13(jc3) + x34(jc3)) + 
     &                               x23(jc)*(-z13(jc3) + z34(jc3)) 
            dbijy(4,in)            = z14(jc)*(x24(jc3)  + x34(jc3)) - 
     &                               x14(jc)*(z24(jc3)  + z34(jc3))
            dbijy(5,in)            = z13(jc)*(x12(jc3)  + x14(jc3)) - 
     &                               x13(jc)*(z12(jc3)  + z14(jc3))
            dbijy(6,in)            =-z12(jc)*(x23(jc3)  + x24(jc3)) + 
     &                               x12(jc)*(z23(jc3)  + z24(jc3)) 
c
            dbijz(1,in)            = x34(jc)*(y14(jc3)  + y24(jc3)) - 
     &                               y34(jc)*(x14(jc3)  + x24(jc3))
            dbijz(2,in)            = x24(jc)*(-y12(jc3) + y23(jc3)) - 
     &                               y24(jc)*(-x12(jc3) + x23(jc3))
            dbijz(3,in)            =-x23(jc)*(-y13(jc3) + y34(jc3)) + 
     &                               y23(jc)*(-x13(jc3) + x34(jc3)) 
            dbijz(4,in)            = x14(jc)*(y24(jc3)  + y34(jc3)) - 
     &                               y14(jc)*(x24(jc3)  + x34(jc3))
            dbijz(5,in)            = x13(jc)*(y12(jc3)  + y14(jc3)) - 
     &                               y13(jc)*(x12(jc3)  + x14(jc3)) 
            dbijz(6,in)            =-x12(jc)*(y23(jc3)  + y24(jc3)) + 
     &                               y12(jc)*(x23(jc3)  + x24(jc3))
     &                               
c
30       CONTINUE
c
         DO 40 k=1,6
c
            is1                    = nu(nor(k),jt)
            is2                    = nu(nex(k),jt)
            is3                    = nu(nusg(k,1),jt)
            is4                    = nu(nusg(k,2),jt)
c
            w(k,1)                 = us36dt*
     &                               ((coco(1,is1) - coor(1,is1) + 
     &                                 coco(1,is2) - coor(1,is2))*13.0 +
     &                                (coco(1,is3) - coor(1,is3) + 
     &                                 coco(1,is4) - coor(1,is4))*5.0)
c
            w(k,2)                 = us36dt*
     &                               ((coco(2,is1) - coor(2,is1) + 
     &                                 coco(2,is2) - coor(2,is2))*13.0 +
     &                                (coco(2,is3) - coor(2,is3) + 
     &                                 coco(2,is4) - coor(2,is4))*5.0)
c
            w(k,3)                 = us36dt*
     &                               ((coco(3,is1) - coor(3,is1) +
     &                                 coco(3,is2) - coor(3,is2))*13.0 +
     &                                (coco(3,is3) - coor(3,is3) + 
     &                                 coco(3,is4) - coor(3,is4))*5.0)
c
40       CONTINUE
c
         DO 500 k=1,6
c
            is1                    = nu(nor(k),jt)
            is2                    = nu(nex(k),jt)
c
            et1                    = eps3*
     &                              (dbijx(k,3) + dbijx(k,4) +
     &                               0.5*(dbijx(k,1) + dbijx(k,2)))
c
            et2                    = eps3*
     &                              (dbijy(k,3) + dbijy(k,4) +
     &                               0.5*(dbijy(k,1) + dbijy(k,2)))
c
            et3                    = eps3*
     &                              (dbijz(k,3) + dbijz(k,4) + 
     &                               0.5*(dbijz(k,1) + dbijz(k,2)))
c
            tvno(jt,nor(k),1)      = tvno(jt,nor(k),1) + et1 
            tvno(jt,nor(k),2)      = tvno(jt,nor(k),2) + et2 
            tvno(jt,nor(k),3)      = tvno(jt,nor(k),3) + et3 
c
            tvno(jt,nex(k),1)      = tvno(jt,nex(k),1) - et1 
            tvno(jt,nex(k),2)      = tvno(jt,nex(k),2) - et2 
            tvno(jt,nex(k),3)      = tvno(jt,nex(k),3) - et3 
c
            ind                    = 0
c
            DO 100 kv=1,ndmax
c
               IF (jaret(is1,kv) .EQ. 0) GOTO 110
c
               iseg                = jaret(is1,kv)
c
               nub1                = nubo(1,iseg)
               nub2                = nubo(2,iseg)
c
               IF ((is1 .EQ. nub1) .AND. (is2 .EQ. nub2)) THEN
                  ind              = 1
                  esign            = ABS(esign)
                  GOTO 200
               ELSE
                  IF ((is1 .EQ. nub2) .AND. (is2 .EQ. nub1)) THEN
                     ind           = 1
                     esign         =-ABS(esign)
                     GOTO 200
                  ENDIF
               ENDIF
c
100         CONTINUE
c
110         IF (ind .EQ. 0) THEN 
c
               WRITE(6, *) 'Error in CMVVNO'
               WRITE(6, *) is1, is2, nub1, nub2
c
               CALL TILT
c
            ENDIF
c
200         CONTINUE
c
            eta1                   = esign*et1
            eta2                   = esign*et2
            eta3                   = esign*et3
c
            vnocl(1,iseg)          = vnocl(1,iseg) - eta1 
            vnocl(2,iseg)          = vnocl(2,iseg) - eta2 
            vnocl(3,iseg)          = vnocl(3,iseg) - eta3 
c
            sigma(iseg)            = sigma(iseg) - 
     &                               (eta1*w(k,1) + 
     &                                eta2*w(k,2) + eta3*w(k,3))
c
500      CONTINUE
c
1000  CONTINUE
c
      DO 1025 iseg=1,nseg
         eta1                      = vnocl(1,iseg)*vnocl(1,iseg) +
     &                               vnocl(2,iseg)*vnocl(2,iseg) + 
     &                               vnocl(3,iseg)*vnocl(3,iseg)
         eta1                      = SQRT(eta1)
         IF (eta1 .LT. 1.0e-16) THEN 
            vnocl(1,iseg)          = 1.0e-16
            vnocl(2,iseg)          = 1.0e-16
            vnocl(3,iseg)          = 1.0e-16
         ENDIF
1025  CONTINUE  
c
      RETURN
      END
