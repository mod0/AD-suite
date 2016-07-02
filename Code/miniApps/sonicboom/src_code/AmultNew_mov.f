      SUBROUTINE AmultNew_mov(v1x,v1y,v1z,Av1x,Av1y,Av1z,
     .           si,coor,nubo,nu,
     &           logfr,mvlgfr,kspring,ns,nt,nseg)
c
c
c
      implicit none
c
      include 'param.h'
c
      real*8 v1x(*), v1y(*), v1z(*), Av1x(*), Av1y(*),Av1z(*)
      real*8 si
c
      integer ns, nt, nseg, nu(4,ntmax)
      integer logfr(nsmax), mvlgfr(nsmax), kspring, nubo(2,nsegmax)
      real*8 coor(3,nsmax)

c
      real*8 newsi
c
      INTEGER iseg,isglow,isgup, nubo1, nubo2
      INTEGER is
c
      real*8 epsilon
c
      real*8 tb11, tb12, tb13,
     $     tb21, tb22, tb23,
     $     tb31, tb32, tb33
c
      real*8 xlocbar1,
     $     xlocbar2,
     $     xlocbar3
c
      real*8 zlocbar1,
     $     zlocbar2,
     $     zlocbar3
c
      real*8 xmag, ymag, zmag
c
      real*8 kbar11, kbar12, kbar13, kbar14, kbar15, kbar16,
     $     kbar22, kbar23, kbar24, kbar25, kbar26,
     $     kbar33, kbar34, kbar35, kbar36,
     $     kbar44, kbar45, kbar46,
     $     kbar55, kbar56,
     $     kbar66
c
      real*8 a, b, c
c
      real*8 mult1, mult2, mult3
c
      integer ele
c
c
      real*8 u1, u2, u3, u4
      real*8 v1, v2, v3, v4
      real*8 w1, w2, w3, w4
c
      integer toflag, barflag
c
c
c     ***********************************************************
c
      epsilon = 1.0d-08
c
      barflag = 0
      toflag  = 0
c
      if(kspring.eq.1) then

         barflag = 1

      endif
c
      if(kspring.eq.2) then

         toflag = 1

      endif
c
c     ------------------------------------------------------------
c
      isglow                    = 1
      isgup                     = nseg
c
c
      DO is = 1, ns

        Av1x(is)          = 0.0d0
        Av1y(is)          = 0.0d0
        Av1z(is)          = 0.0d0

      END DO
c
c
      if (barflag.EQ.1) then
c
c
c
c     Loop on local list of edges
c
      DO iseg = isglow,isgup
c
c        Local indexing of the vertices of the current edge
c
         nubo1                     = nubo(1,iseg)
         nubo2                     = nubo(2,iseg)
c
         newsi = si
c
c        Computing the stiffness of the current edge
c
         u1     =     v1x(nubo1)
         v1     =     v1y(nubo1)
         w1     =     v1z(nubo1)
c
         u2     =     v1x(nubo2)
         v2     =     v1y(nubo2)
         w2     =     v1z(nubo2)
c
c
c     compute local x direction, along edge
c
      xlocbar1 = coor(1,nubo2) - coor(1,nubo1)
      xlocbar2 = coor(2,nubo2) - coor(2,nubo1)
      xlocbar3 = coor(3,nubo2) - coor(3,nubo1)
c
      xmag = sqrt(xlocbar1*xlocbar1+xlocbar2*xlocbar2+xlocbar3*xlocbar3)
c
      tb11 = xlocbar1/xmag
      tb21 = xlocbar2/xmag
      tb31 = xlocbar3/xmag

c
c     pick local y to be orthogonal to local x
c
      if(abs(xlocbar1).le.epsilon.and.abs(xlocbar2).le.epsilon) then
c
c     this means that the current xlocal is in the z direction
c  
        tb12     =  1.0d0
        tb22     =  0.0d0
        tb32     =  0.0d0
c
*        ymag = 1.0d0
c
        zlocbar1 =   0.0d0
        zlocbar2 =   tb31*tb12
        zlocbar3 = - tb21*tb12      
c
        zmag = sqrt(zlocbar2*zlocbar2+zlocbar3*zlocbar3)
c
        tb13 =  0.0d0
        tb23 =  zlocbar2/zmag
        tb33 =  zlocbar3/zmag

c
      else
c
c
        ymag = sqrt(xlocbar2*xlocbar2+xlocbar1*xlocbar1)
c
        tb12 =  xlocbar2/ymag 
        tb22 = -xlocbar1/ymag 
        tb32 =  0.0d0

c
        zlocbar1 = -tb31*tb22
        zlocbar2 =  tb31*tb12
        zlocbar3 =  tb11*tb22 - tb21*tb12      
c
        zmag = sqrt(zlocbar1*zlocbar1+zlocbar2*zlocbar2+
     $              zlocbar3*zlocbar3)
c

*      write(6,*)'---> the magnitude of z is', zmag
c
        tb13 =  zlocbar1/zmag
        tb23 =  zlocbar2/zmag
        tb33 =  zlocbar3/zmag 
c
      endif
c
      a =  1.0d0/xmag
      b =  a
      c = -a
c
c
      mult1  = tb11*tb11
      kbar11 = a*mult1
      kbar12 = a*tb11*tb21
      kbar13 = a*tb11*tb31
      kbar14 = c*mult1
      kbar15 = c*tb11*tb21
      kbar16 = c*tb11*tb31
c
c
      mult2  = tb21*tb21
   
      kbar22 = a*mult2
      kbar23 = a*tb21*tb31
      kbar24 = c*tb11*tb21  
      kbar25 = c*mult2
      kbar26 = c*tb21*tb31
c
c
      mult3  = tb31*tb31
    
      kbar33 = a*mult3
      kbar34 = c*tb11*tb31
      kbar35 = c*tb21*tb31
      kbar36 = c*mult3
c
c
    
      kbar44 = b*tb11*tb11
      kbar45 = b*tb11*tb21
      kbar46 = b*tb11*tb31
c
c
c
      kbar55 = b*tb21*tb21
      kbar56 = b*tb21*tb31
c
c  

      kbar66 = b*tb31*tb31
c
c
c        superblub
c
         IF((logfr(nubo1).EQ.0) .OR.   (mvlgfr(nubo1).GT.0)) THEN
c
c
             Av1x(nubo1) = Av1x(nubo1) +
     $                 newsi*(kbar11*u1+kbar12*v1+kbar13*w1+
     $                     kbar14*u2+kbar15*v2+kbar16*w2)

             Av1y(nubo1) = Av1y(nubo1) +
     $                 newsi*(kbar12*u1+kbar22*v1+kbar23*w1+
     $                     kbar24*u2+kbar25*v2+kbar26*w2)

             Av1z(nubo1) = Av1z(nubo1) +
     $                 newsi*(kbar13*u1+kbar23*v1+kbar33*w1+
     $                     kbar34*u2+kbar35*v2+kbar36*w2)
c
c
c
         ENDIF
c
         IF ((logfr(nubo2).EQ.0) .OR. (mvlgfr(nubo2).GT.0)) THEN
c
c
c
             Av1x(nubo2) = Av1x(nubo2) +
     $                 newsi*(kbar14*u1+kbar24*v1+kbar34*w1+
     $                     kbar44*u2+kbar45*v2+kbar46*w2)

             Av1y(nubo2) = Av1y(nubo2) +
     $                 newsi*(kbar15*u1+kbar25*v1+kbar35*w1+
     $                     kbar45*u2+kbar55*v2+kbar56*w2)

             Av1z(nubo2) = Av1z(nubo2) +
     $                 newsi*(kbar16*u1+kbar26*v1+kbar36*w1+
     $                     kbar46*u2+kbar56*v2+kbar66*w2)
c

         ENDIF
c
      ENDDO
c
c
c
c     end of barflag
c     |
      endif

c
c
c     tetra computation
c
c
      if(toflag.EQ.1) then
c
c
c
c
c     now we need to add the tetra elements stifnesses
c
c     LOOP over local list of TETRAS
c
      do ele = 1, nt
c
c  
         u1  =  v1x(nu(1,ele))
         v1  =  v1y(nu(1,ele))
         w1  =  v1z(nu(1,ele))
c
         u2  =  v1x(nu(2,ele))
         v2  =  v1y(nu(2,ele))
         w2  =  v1z(nu(2,ele))
c
         u3  =  v1x(nu(3,ele))
         v3  =  v1y(nu(3,ele))
         w3  =  v1z(nu(3,ele))
c
         u4  =  v1x(nu(4,ele))
         v4  =  v1y(nu(4,ele))
         w4  =  v1z(nu(4,ele))
c
c
         IF ((logfr(nu(1,ele)).EQ.0) .OR. (mvlgfr(nu(1,ele)).GT.0))THEN
c
           Av1x(nu(1,ele)) = Av1x(nu(1,ele))+
     $             si*(ktet1_1(ele)*u1+ktet1_2(ele)*v1+ktet1_3(ele)*w1+
     $             ktet1_4(ele)*u2+ktet1_5(ele)*v2+ktet1_6(ele)*w2+
     $             ktet1_7(ele)*u3+ktet1_8(ele)*v3+ktet1_9(ele)*w3+
     $             ktet1_10(ele)*u4+ktet1_11(ele)*v4+ktet1_12(ele)*w4)
c
          Av1y(nu(1,ele)) = Av1y(nu(1,ele))+
     $             si*(ktet1_2(ele)*u1+ktet2_2(ele)*v1+ktet2_3(ele)*w1+
     $             ktet2_4(ele)*u2+ktet2_5(ele)*v2+ktet2_6(ele)*w2+
     $             ktet2_7(ele)*u3+ktet2_8(ele)*v3+ktet2_9(ele)*w3+
     $             ktet2_10(ele)*u4+ktet2_11(ele)*v4+ktet2_12(ele)*w4)
c
          Av1z(nu(1,ele)) = Av1z(nu(1,ele))+
     $             si*(ktet1_3(ele)*u1+ktet2_3(ele)*v1+ktet3_3(ele)*w1+
     $             ktet3_4(ele)*u2+ktet3_5(ele)*v2+ktet3_6(ele)*w2+
     $             ktet3_7(ele)*u3+ktet3_8(ele)*v3+ktet3_9(ele)*w3+
     $             ktet3_10(ele)*u4+ktet3_11(ele)*v4+ktet3_12(ele)*w4)
c
c
         endif
c
c
c
         IF ((logfr(nu(2,ele)).EQ.0) .OR. (mvlgfr(nu(2,ele)).GT.0)) THEN
c
c
          Av1x(nu(2,ele)) = Av1x(nu(2,ele))+
     $              si*(ktet1_4(ele)*u1+ktet2_4(ele)*v1+ktet3_4(ele)*w1+
     $              ktet4_4(ele)*u2+ktet4_5(ele)*v2+ktet4_6(ele)*w2+
     $              ktet4_7(ele)*u3+ktet4_8(ele)*v3+ktet4_9(ele)*w3+
     $              ktet4_10(ele)*u4+ktet4_11(ele)*v4+ktet4_12(ele)*w4)
c
          Av1y(nu(2,ele)) = Av1y(nu(2,ele))+
     $              si*(ktet1_5(ele)*u1+ktet2_5(ele)*v1+ktet3_5(ele)*w1+
     $              ktet4_5(ele)*u2+ktet5_5(ele)*v2+ktet5_6(ele)*w2+
     $              ktet5_7(ele)*u3+ktet5_8(ele)*v3+ktet5_9(ele)*w3+
     $              ktet5_10(ele)*u4+ktet5_11(ele)*v4+ktet5_12(ele)*w4)
c
          Av1z(nu(2,ele)) = Av1z(nu(2,ele))+
     $              si*(ktet1_6(ele)*u1+ktet2_6(ele)*v1+ktet3_6(ele)*w1+
     $              ktet4_6(ele)*u2+ktet5_6(ele)*v2+ktet6_6(ele)*w2+
     $              ktet6_7(ele)*u3+ktet6_8(ele)*v3+ktet6_9(ele)*w3+
     $              ktet6_10(ele)*u4+ktet6_11(ele)*v4+ktet6_12(ele)*w4)
c
c
c
         endif
c
c
         IF ((logfr(nu(3,ele)).EQ.0) .OR. (mvlgfr(nu(3,ele)).GT.0))THEN
c
c
          Av1x(nu(3,ele)) = Av1x(nu(3,ele))+
     $             si*(ktet1_7(ele)*u1+ktet2_7(ele)*v1+ktet3_7(ele)*w1+
     $             ktet4_7(ele)*u2+ktet5_7(ele)*v2+ktet6_7(ele)*w2+
     $             ktet7_7(ele)*u3+ktet7_8(ele)*v3+ktet7_9(ele)*w3+
     $             ktet7_10(ele)*u4+ktet7_11(ele)*v4+ktet7_12(ele)*w4)
c
          Av1y(nu(3,ele)) = Av1y(nu(3,ele))+
     $             si*(ktet1_8(ele)*u1+ktet2_8(ele)*v1+ktet3_8(ele)*w1+
     $             ktet4_8(ele)*u2+ktet5_8(ele)*v2+ktet6_8(ele)*w2+
     $             ktet7_8(ele)*u3+ktet8_8(ele)*v3+ktet8_9(ele)*w3+
     $             ktet8_10(ele)*u4+ktet8_11(ele)*v4+ktet8_12(ele)*w4)
c
          Av1z(nu(3,ele)) = Av1z(nu(3,ele))+
     $             si*(ktet1_9(ele)*u1+ktet2_9(ele)*v1+ktet3_9(ele)*w1+
     $             ktet4_9(ele)*u2+ktet5_9(ele)*v2+ktet6_9(ele)*w2+
     $             ktet7_9(ele)*u3+ktet8_9(ele)*v3+ktet9_9(ele)*w3+
     $             ktet9_10(ele)*u4+ktet9_11(ele)*v4+ktet9_12(ele)*w4)
c
c
         endif
c
c
c
         IF ((logfr(nu(4,ele)).EQ.0) .OR. (mvlgfr(nu(4,ele)).GT.0))THEN
c
          Av1x(nu(4,ele)) = Av1x(nu(4,ele))+
     $          si*(ktet1_10(ele)*u1+ktet2_10(ele)*v1+ktet3_10(ele)*w1+
     $          ktet4_10(ele)*u2+ktet5_10(ele)*v2+ktet6_10(ele)*w2+
     $          ktet7_10(ele)*u3+ktet8_10(ele)*v3+ktet9_10(ele)*w3+
     $          ktet10_10(ele)*u4+ktet10_11(ele)*v4+ktet10_12(ele)*w4)
c
          Av1y(nu(4,ele)) = Av1y(nu(4,ele))+
     $          si*(ktet1_11(ele)*u1+ktet2_11(ele)*v1+ktet3_11(ele)*w1+
     $          ktet4_11(ele)*u2+ktet5_11(ele)*v2+ktet6_11(ele)*w2+
     $          ktet7_11(ele)*u3+ktet8_11(ele)*v3+ktet9_11(ele)*w3+
     $          ktet10_11(ele)*u4+ktet11_11(ele)*v4+ktet11_12(ele)*w4)
c
          Av1z(nu(4,ele)) = Av1z(nu(4,ele))+
     $          si*(ktet1_12(ele)*u1+ktet2_12(ele)*v1+ktet3_12(ele)*w1+
     $          ktet4_12(ele)*u2+ktet5_12(ele)*v2+ktet6_12(ele)*w2+
     $          ktet7_12(ele)*u3+ktet8_12(ele)*v3+ktet9_12(ele)*w3+
     $          ktet10_12(ele)*u4+ktet11_12(ele)*v4+ktet12_12(ele)*w4)
c
c
c
         ENDIF
c
c     end of loop over tetras
c     |
      enddo                
c
c
c
c     end of toflag
c     |
      endif
c
c
c


c      DO is=1,ns

c        IF (mvlgfr(is).GT.0) then
     
c           CALL PRJDEP(is, Av1x(is), Av1y(is), Av1z(is))

c        ENDIF

c      ENDDO

      RETURN

      END









































