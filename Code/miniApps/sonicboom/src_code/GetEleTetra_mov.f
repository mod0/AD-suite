c     *******************************************************************
c 
      subroutine GetEleTetra_mov(nufaces,nucycles,lowlambda,highlambda,
     &           map_pointer,coor,nu,nt)
c
c
      implicit none
c
      include 'param.h'
c
      integer ele
c
      real*8 t11, t12, t13, 
     $     t21, t22, t23,
     $     t31, t32, t33
c
      real*8 xi, yi, zi,
     $     xj, yj, zj,
     $     xk, yk, zk
c
      real*8 ktetra(12,12)
c
      real*8 xlocj
      real*8 xlock, ylock, yfac     
c
      real*8 lambda
c
      real*8 ktetglob(12,12)
c
      real*8 lowlambda,highlambda
      integer nufaces,nucycles,map_pointer(4,3,4)
      integer nt,nu(4,ntmax)
      real*8 coor(3,nsmax)
c
c    --------------------------variables for serial version--------------
c
      real*8 rx12, ry12, rz12
      real*8 rx13, ry13, rz13
      real*8 rx14, ry14, rz14
      real*8 rx23, ry23, rz23
      real*8 rx24, ry24, rz24
      real*8 rx34, ry34, rz34
c
      real*8 avx(4), avy(4), avz(4)
c
      real*8 bvx(4,3), bvy(4,3), bvz(4,3)
      real*8 cvx(4,3), cvy(4,3), cvz(4,3)
      real*8 rhx(4,3), rhy(4,3), rhz(4,3)
c
      integer face, cycle
c
      real*8 num, den
c
      real*8  r11, r12, r14, r15, r16,
     $      r22, r23, r24, r25, r26,
     $      r31, r32, r33, r34, r35, r36
c 
       real*8 c1, c2, c3, cfac
c
       real*8 b13, b23
       real*8 a12, a13, a23
c
       real*8 x23
c
       real*8 lsq12, lsq13, lsq23
c
      real*8 xqj, yqj, zqj,
     $     magrqj
c
      real*8 magcross
c
      real*8 xik, yik, zik
c
      real*8 area, areasq
      real*8 lsq1, ls1, lasq, la1la
c
      real*8         k11, k12, k13, k14, k15, k16,
     $             k22, k23, k24, k25, k26,
     $             k33, k34, k35, k36,
     $             k44, k45, k46,
     $             k55, k56,
     $             k66
c
      real*8               kg11,kg12,kg13,kg14,kg15,kg16,kg17,kg18,kg19,
     $                   kg22,kg23,kg24,kg25,kg26,kg27,kg28,kg29,
     $                   kg33,kg34,kg35,kg36,kg37,kg38,kg39,
     $                   kg44,kg45,kg46,kg47,kg48,kg49,
     $                   kg55,kg56,kg57,kg58,kg59,
     $                   kg66,kg67,kg68,kg69,
     $                   kg77,kg78,kg79,
     $                   kg88,kg89,  
     $                   kg99   
c
      real*8 short1, short3, short4, short5, short6, 
     $     short8, short9, short10, short11, short12, short13, short14,
     $     short15
c
      real*8 fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9, fac10
      real*8 fac11, fac12, fac13, fac14, fac15, fac16, fac17, fac18, 
     $     fac19, fac20
      real*8 fac21, fac22, fac23, fac24, fac25, fac26, fac27, fac28, 
     $     fac29, fac30
c
      integer p, q, r, s
c
      integer the_row, the_col
c 
      integer ro, col
c
      integer row_1, row_2, row_3
      integer col_1, col_2, col_3
c
c
c
c    --------------------------end of variables for serial version--------------
c
c     loop over all tetra elements and compute the respective tetra stiffness matrix
c
      do ele = 1, nt
c
c
c            initialize current tetra stiffness matrix
c
             do ro = 1, 12

                do col = 1, 12
          
                   ktetglob(ro,col) = 0.0

                enddo

             enddo
c
c
          rx12 =  coor(1,nu(2,ele)) - coor(1,nu(1,ele))
          ry12 =  coor(2,nu(2,ele)) - coor(2,nu(1,ele))
          rz12 =  coor(3,nu(2,ele)) - coor(3,nu(1,ele))
c
c
c
          rx13 =  coor(1,nu(3,ele)) - coor(1,nu(1,ele))
          ry13 =  coor(2,nu(3,ele)) - coor(2,nu(1,ele))
          rz13 =  coor(3,nu(3,ele)) - coor(3,nu(1,ele))
c
c
c
          rx14 =  coor(1,nu(4,ele)) - coor(1,nu(1,ele))
          ry14 =  coor(2,nu(4,ele)) - coor(2,nu(1,ele))
          rz14 =  coor(3,nu(4,ele)) - coor(3,nu(1,ele))
c
c
c
          rx23 =  coor(1,nu(3,ele)) - coor(1,nu(2,ele))
          ry23 =  coor(2,nu(3,ele)) - coor(2,nu(2,ele))
          rz23 =  coor(3,nu(3,ele)) - coor(3,nu(2,ele))
c
c
c
          rx24 =  coor(1,nu(4,ele)) - coor(1,nu(2,ele)) 
          ry24 =  coor(2,nu(4,ele)) - coor(2,nu(2,ele))
          rz24 =  coor(3,nu(4,ele)) - coor(3,nu(2,ele))
c
c
c
          rx34 =  coor(1,nu(4,ele)) - coor(1,nu(3,ele))
          ry34 =  coor(2,nu(4,ele)) - coor(2,nu(3,ele))
          rz34 =  coor(3,nu(4,ele)) - coor(3,nu(3,ele))
c
c
c
c         now compute the normal vectors of the 4 tetra faces
c
c     face 1-2-3
c
      avx(1)  =   ry12*rz13 - rz12*ry13
      avy(1)  = -(rx12*rz13 - rz12*rx13)
      avz(1)  =   rx12*ry13 - ry12*rx13     
c
c     face 2-1-4
c
      avx(2)  =   ry14*rz12 - rz14*ry12
      avy(2)  = -(rx14*rz12 - rz14*rx12)
      avz(2)  =   rx14*ry12 - ry14*rx12     
c
c     face 4-1-3
c
      avx(3)  =   ry13*rz14 - rz13*ry14
      avy(3)  = -(rx13*rz14 - rz13*rx14)
      avz(3)  =   rx13*ry14 - ry13*rx14     
c
c     face 2-4-3
c
      avx(4)  =   ry24*rz23 - rz24*ry23
      avy(4)  = -(rx24*rz23 - rz24*rx23)
      avz(4)  =   rx24*ry23 - ry24*rx23     
c
c     --------------------------------------------------------------
c     map the computed vectors into the ones needed
c     each loop over the faces and their corresponding cycles
c     the first index is the face number and the second one is the 
c     cycle number
c     --------------------------------------------------------------
c
c     face 1, cycle 1
c
c
      bvx(1,1)    =  rx24
      bvy(1,1)    =  ry24
      bvz(1,1)    =  rz24
c
      cvx(1,1)    =  rx13
      cvy(1,1)    =  ry13
      cvz(1,1)    =  rz13
c
      rhx(1,1)    = -rx12
      rhy(1,1)    = -ry12
      rhz(1,1)    = -rz12
c
c     face 1, cycle 2
c
      bvx(1,2)    =  rx14
      bvy(1,2)    =  ry14
      bvz(1,2)    =  rz14
c
      cvx(1,2)    = -rx23
      cvy(1,2)    = -ry23
      cvz(1,2)    = -rz23
c
      rhx(1,2)    =  rx13
      rhy(1,2)    =  ry13
      rhz(1,2)    =  rz13 
c
c     face 1, cycle 3
c
      bvx(1,3)    =  rx34
      bvy(1,3)    =  ry34
      bvz(1,3)    =  rz34
c
      cvx(1,3)    = -rx12
      cvy(1,3)    = -ry12
      cvz(1,3)    = -rz12
c
      rhx(1,3)    = -rx23
      rhy(1,3)    = -ry23
      rhz(1,3)    = -rz23
c
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
c
c     face 2, cycle 1
c
c
c
      bvx(2,1)    =  rx13
      bvy(2,1)    =  ry13
      bvz(2,1)    =  rz13
c
      cvx(2,1)    =  rx24
      cvy(2,1)    =  ry24
      cvz(2,1)    =  rz24
c
      rhx(2,1)    =  rx12
      rhy(2,1)    =  ry12
      rhz(2,1)    =  rz12 
c
c     face 2, cycle 2
c
      bvx(2,2)    =  rx23
      bvy(2,2)    =  ry23
      bvz(2,2)    =  rz23
c
      cvx(2,2)    = -rx14
      cvy(2,2)    = -ry14
      cvz(2,2)    = -rz14
c
      rhx(2,2)    =  rx24
      rhy(2,2)    =  ry24
      rhz(2,2)    =  rz24 
c
c     face 2, cycle 3
c
      bvx(2,3)    = -rx34
      bvy(2,3)    = -ry34
      bvz(2,3)    = -rz34
c
      cvx(2,3)    =  rx12
      cvy(2,3)    =  ry12
      cvz(2,3)    =  rz12
c
      rhx(2,3)    = -rx14
      rhy(2,3)    = -ry14
      rhz(2,3)    = -rz14
c
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c     face 3, cycle 1
c
c
c
      bvx(3,1)    =  rx12
      bvy(3,1)    =  ry12
      bvz(3,1)    =  rz12
c
      cvx(3,1)    = -rx34
      cvy(3,1)    = -ry34
      cvz(3,1)    = -rz34
c
      rhx(3,1)    =  rx14
      rhy(3,1)    =  ry14
      rhz(3,1)    =  rz14 
c
c     face 3, cycle 2
c
      bvx(3,2)    = -rx24
      bvy(3,2)    = -ry24
      bvz(3,2)    = -rz24
c
      cvx(3,2)    = -rx13
      cvy(3,2)    = -ry13
      cvz(3,2)    = -rz13
c
      rhx(3,2)    = -rx34
      rhy(3,2)    = -ry34
      rhz(3,2)    = -rz34
c
c     face 3 , cycle 3
c
      bvx(3,3)    = -rx23
      bvy(3,3)    = -ry23
      bvz(3,3)    = -rz23
c
      cvx(3,3)    =  rx14
      cvy(3,3)    =  ry14
      cvz(3,3)    =  rz14
c
      rhx(3,3)    = -rx13
      rhy(3,3)    = -ry13
      rhz(3,3)    = -rz13
c 
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c     face 4, cycle 1
c
c
c
      bvx(4,1)    = -rx14
      bvy(4,1)    = -ry14
      bvz(4,1)    = -rz14
c
      cvx(4,1)    =  rx23
      cvy(4,1)    =  ry23
      cvz(4,1)    =  rz23
c
      rhx(4,1)    = -rx24
      rhy(4,1)    = -ry24
      rhz(4,1)    = -rz24        
c
c     face 4, cycle 2
c
      bvx(4,2)    = -rx13
      bvy(4,2)    = -ry13
      bvz(4,2)    = -rz13
c
      cvx(4,2)    = -rx24
      cvy(4,2)    = -ry24
      cvz(4,2)    = -rz24
c
      rhx(4,2)    =  rx34
      rhy(4,2)    =  ry34
      rhz(4,2)    =  rz34 
c
c     face 4, cycle 3
c
      bvx(4,3)    = -rx12
      bvy(4,3)    = -ry12
      bvz(4,3)    = -rz12
c
      cvx(4,3)    =  rx34
      cvy(4,3)    =  ry34
      cvz(4,3)    =  rz34
c
      rhx(4,3)    =  rx23
      rhy(4,3)    =  ry23
      rhz(4,3)    =  rz23 
c
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c
c   ---------------------------------------------------------------------
c
c
c
      do face = 1, nufaces
c
c     2nd COMPUTE ylocal
c            
      magcross = sqrt(avx(face)*avx(face) + 
     $                avy(face)*avy(face) + 
     $                avz(face)*avz(face))
c
      t12 = avx(face)/magcross
      t22 = avy(face)/magcross
      t32 = avz(face)/magcross
c         
         do cycle = 1, nucycles
c
            xi = coor(1,nu(map_pointer(face,cycle,2),ele))   
            yi = coor(2,nu(map_pointer(face,cycle,2),ele))
            zi = coor(3,nu(map_pointer(face,cycle,2),ele))
c
c
            xk = coor(1,nu(map_pointer(face,cycle,4),ele)) 
            yk = coor(2,nu(map_pointer(face,cycle,4),ele))
            zk = coor(3,nu(map_pointer(face,cycle,4),ele))
c
c
c
c
            num = (avz(face)*bvy(face,cycle)*rhx(face,cycle) - 
     $             avy(face)*bvz(face,cycle)*rhx(face,cycle) -
     $             avz(face)*bvx(face,cycle)*rhy(face,cycle) + 
     $             avx(face)*bvz(face,cycle)*rhy(face,cycle) + 
     $             avy(face)*bvx(face,cycle)*rhz(face,cycle) - 
     $             avx(face)*bvy(face,cycle)*rhz(face,cycle))
c
c
            den = (-(avz(face)*bvy(face,cycle)*cvx(face,cycle))  + 
     $               avy(face)*bvz(face,cycle)*cvx(face,cycle)   + 
     $               avz(face)*bvx(face,cycle)*cvy(face,cycle)   - 
     $               avx(face)*bvz(face,cycle)*cvy(face,cycle)   - 
     $               avy(face)*bvx(face,cycle)*cvz(face,cycle)   + 
     $               avx(face)*bvy(face,cycle)*cvz(face,cycle))
c
c
c
c     check if local y direction concides with the base face
c     normal vector, in that case we pick the kinematic triangle with
c     the minimal area, see notes from 03/27/1998
c
      if (abs(den).lt.1.0e-20) then
c
c

          num = -rhx(face,cycle)*cvx(face,cycle)
     $          -rhy(face,cycle)*cvy(face,cycle)
     $          -rhz(face,cycle)*cvz(face,cycle)

          den = cvx(face,cycle)*cvx(face,cycle)+ 
     $          cvy(face,cycle)*cvy(face,cycle)+
     $          cvz(face,cycle)*cvz(face,cycle)
c
c
      endif


      lambda = num/den       

c
c
c    lambda is equal to the lumping factor
c
c
c
c
c    now, we need to compute the actual intersection point
c
      xj = coor(1,nu(map_pointer(face,cycle,1),ele)) + 
     $     lambda*cvx(face,cycle)

      yj = coor(2,nu(map_pointer(face,cycle,1),ele)) + 
     $     lambda*cvy(face,cycle)

      zj = coor(3,nu(map_pointer(face,cycle,1),ele)) + 
     $     lambda*cvz(face,cycle)
c
c
c
c     change lambda if necessary, that means if the intersection lies
c     outside the actual tetrahedron
c     the change of lambda will be accomplished according to the
c     following rule
c     lambda > 1.0 --> lambda = 1.0
c     lambda < 0.0 --> lambda = 0.0
c     this is only done to get the right force lumping
c     the actual intersection point will be computed using the actual
c     lambda !!!!
c
c
c
      if(lambda.lt.lowlambda) then

         lambda = 0.0

      endif
c
      if(lambda.gt.highlambda) then

         lambda = 1.0

      endif
c
c

       ls1          =      (1.0 - lambda)
       lsq1         =      ls1*ls1
       lasq         =      lambda*lambda
       la1la        =      ls1*lambda
c
c
c     **************************************************************
c
c     this routine calculates the transformation matrix
c     necessary to get the global triangle stiffness matrix 
c
c     1st COMPUTE xlocal
c
      xqj  = xj - coor(1,nu(map_pointer(face,cycle,2),ele))
      yqj  = yj - coor(2,nu(map_pointer(face,cycle,2),ele))
      zqj  = zj - coor(3,nu(map_pointer(face,cycle,2),ele))
c 
      magrqj = sqrt(xqj*xqj + yqj*yqj + zqj*zqj)
c
      t11 = xqj/magrqj
      t21 = yqj/magrqj
      t31 = zqj/magrqj
c
c
c
c     3rd COMPUTE zlocal
c
      t13 =   t21*t32-t31*t22
      t23 = -(t11*t32-t31*t11)
      t33 =   t11*t22-t21*t11
c 
c     choose local coordinates of node i to be 0,0,0      
c     transform rij and rik into the local system to obtain
c     the local triangle co-ordinates
c
c
      xik = xk - xi
      yik = yk - yi
      zik = zk - zi
c
      xlocj = t11*(xj - xi) + t21*(yj - yi) + t31*(zj - zi)
c
      xlock = t11*xik + t21*yik + t31*zik
      ylock = t12*xik + t22*yik + t32*zik
c
c
c
c
c     ************************************************************************
c
c    
c
c
c      epsilon = 1.0e-06
c
c            
c      compute geometric properties for the matrices R and C
c
c
c
c      i = 1, j = 2, k = 3
c
c-----------> NOTE: in the local triangle co-ordinate system we have
c
c      xloci = yloci = zloci = 0
c      and ylocj = 0
c
c
       x23  = xlock - xlocj
c
c
       yfac = ylock*ylock
c      
       lsq12 = xlocj*xlocj
       lsq13 = xlock*xlock+yfac
       lsq23 = x23*x23+yfac
c
       a12 =  xlocj/lsq12

       a13 =  xlock/lsq13

       a23 =  x23/lsq23
c
c      
       b13 =  ylock/lsq13
   
       b23 =  ylock/lsq23
c
       r11 =  b13
       r12 =  a12 - a13
       
       r14 = -a12
       r15 = -b13
       r16 =  a13
c
       
       r22 = -a12
       r23 = -b23
       r24 =  a12 + a23
       r25 =  b23
       r26 = -a23
c
       r31 = -b13
       r32 =  a13
       r33 =  b23
       r34 = -a23
       r35 =  b13 - b23
       r36 = -a13 + a23
c
c
c
       area = 0.5*xlocj*ylock
c
       areasq =  area * area 
c    
       cfac = lsq12/areasq

       c1 = cfac*lsq13
       c2 = cfac*lsq23
       c3 = lsq13*lsq23/areasq
c
c
c
c      generate the local triangle stiffness matrix 
c
c
c
c     optimize this routine by computing commonly used factors beforehand 
c
c
      short1   =  c1*r11
   
      short3   =  c3*r31

      k11 = short1*r11              + short3*r31
      k12 = short1*r12              + short3*r32
      k13 =                           short3*r33
      k14 = short1*r14              + short3*r34
      k15 = short1*r15              + short3*r35
      k16 = short1*r16              + short3*r36
c
      short4   =  c1*r12
      short5   =  c2*r22
      short6   =  c3*r32

      k22 = short4*r12 + short5*r22 + short6*r32
      k23 =              short5*r23 + short6*r33
      k24 = short4*r14 + short5*r24 + short6*r34
      k25 = short4*r15 + short5*r25 + short6*r35
      k26 = short4*r16 + short5*r26 + short6*r36
c
      
      short8   =   c2*r23
      short9   =   c3*r33

      k33 =              short8*r23 + short9*r33
      k34 =              short8*r24 + short9*r34
      k35 =              short8*r25 + short9*r35
      k36 =              short8*r26 + short9*r36
c
      short10   =  c1*r14
      short11   =  c2*r24
      short12   =  c3*r34

      k44 = short10*r14 + short11*r24 + short12*r34
      k45 = short10*r15 + short11*r25 + short12*r35
      k46 = short10*r16 + short11*r26 + short12*r36
c
      short13    =   c1*r15
      short14    =   c2*r25
      short15    =   c3*r35

      k55 = short13*r15 + short14*r25 + short15*r35
      k56 = short13*r16 + short14*r26 + short15*r36
c
c
c 
      k66 = c1*r16*r16 + c2*r26*r26 + c3*r36*r36
c
c
c
c  **************************************************************
c
c  **************************************************************
c
c     1
c

      fac1    =   k11*t11 + k12*t12
      fac2    =   k12*t11 + k22*t12

      kg11    =   t11*fac1 + t12*fac2
      kg12    =   t21*fac1 + t22*fac2
      kg13    =   t31*fac1 + t32*fac2

      fac3    =   k13*t11 + k23*t12
      fac4    =   k14*t11 + k24*t12

      kg14    =   t11*fac3 + t12*fac4
      kg15    =   t21*fac3 + t22*fac4
      kg16    =   t31*fac3 + t32*fac4

      fac5    =   k15*t11 + k25*t12
      fac6    =   k16*t11 + k26*t12

      kg17    =   t11*fac5 + t12*fac6
      kg18    =   t21*fac5 + t22*fac6
      kg19    =   t31*fac5 + t32*fac6
c
c     2
c

      fac7    =  k11*t21 + k12*t22
      fac8    =  k12*t21 + k22*t22

      kg22    =  t21*fac7 + t22*fac8
      kg23    =  t31*fac7 + t32*fac8

      fac9    =  k13*t21 + k23*t22
      fac10   =  k14*t21 + k24*t22

      kg24    =  t11*fac9 + t12*fac10
      kg25    =  t21*fac9 + t22*fac10
      kg26    =  t31*fac9 + t32*fac10

      fac11   =  k15*t21 + k25*t22
      fac12   =  k16*t21 + k26*t22

      kg27    =  t11*fac11 + t12*fac12
      kg28    =  t21*fac11 + t22*fac12
      kg29    =  t31*fac11 + t32*fac12
c
c     3
c
      kg33     =  t31*(k11*t31 + k12*t32) + t32*(k12*t31 + k22*t32)

      fac13    =  k13*t31 + k23*t32
      fac14    =  k14*t31 + k24*t32

      kg34     =  t11*fac13 + t12*fac14
      kg35     =  t21*fac13 + t22*fac14
      kg36     =  t31*fac13 + t32*fac14

      fac15    =  k15*t31 + k25*t32
      fac16    =  k16*t31 + k26*t32

      kg37     =  t11*fac15 + t12*fac16
      kg38     =  t21*fac15 + t22*fac16
      kg39     =  t31*fac15 + t32*fac16
c
c     4
c
      fac17    =  k33*t11 + k34*t12
      fac18    =  k34*t11 + k44*t12

      kg44     =  t11*fac17 + t12*fac18
      kg45     =  t21*fac17 + t22*fac18
      kg46     =  t31*fac17 + t32*fac18

      fac19    =  k35*t11 + k45*t12
      fac20    =  k36*t11 + k46*t12 

      kg47     =  t11*fac19 + t12*fac20
      kg48     =  t21*fac19 + t22*fac20
      kg49     =  t31*fac19 + t32*fac20
c
c     5
c
      fac21    =  k33*t21 + k34*t22
      fac22    =  k34*t21 + k44*t22

      kg55     =  t21*fac21 + t22*fac22
      kg56     =  t31*fac21 + t32*fac22

      fac23    =  k35*t21 + k45*t22
      fac24    =  k36*t21 + k46*t22

      kg57     =  t11*fac23 + t12*fac24
      kg58     =  t21*fac23 + t22*fac24
      kg59     =  t31*fac23 + t32*fac24
c
c     6
c
       kg66     =  t31*(k33*t31 + k34*t32) + t32*(k34*t31 + k44*t32)

       fac25    =  k35*t31 + k45*t32
       fac26    =  k36*t31 + k46*t32


       kg67     =  t11*fac25 + t12*fac26
       kg68     =  t21*fac25 + t22*fac26
       kg69     =  t31*fac25 + t32*fac26
c
c      7
c
        fac27    =  k55*t11 + k56*t12
        fac28    =  k56*t11 + k66*t12

        kg77     =  t11*fac27 + t12*fac28
        kg78     =  t21*fac27 + t22*fac28
        kg79     =  t31*fac27 + t32*fac28
c
c       8
c
         fac29     =  k55*t21 + k56*t22
         fac30     =  k56*t21 + k66*t22

         kg88      =  t21*fac29 + t22*fac30
         kg89      =  t31*fac29 + t32*fac30
c
c       9
c
         kg99      =  t31*(k55*t31 + k56*t32) + t32*(k56*t31 + k66*t32)
c
c
c      computing the global tetrahedron stiffness matrix in
c      the generic p-q-r-s system using the lumping factor lambda
c

c
c
c
       ktetra(1,1)  = kg44*lsq1
       ktetra(1,2)  = kg45*lsq1
       ktetra(1,3)  = kg46*lsq1
       ktetra(1,4)  = kg14*ls1
       ktetra(1,5)  = kg24*ls1
       ktetra(1,6)  = kg34*ls1
       ktetra(1,7)  = kg44*la1la
       ktetra(1,8)  = kg45*la1la
       ktetra(1,9)  = kg46*la1la
       ktetra(1,10) = kg47*ls1
       ktetra(1,11) = kg48*ls1
       ktetra(1,12) = kg49*ls1
c
c
       ktetra(2,2)  = kg55*lsq1
       ktetra(2,3)  = kg56*lsq1
       ktetra(2,4)  = kg15*ls1
       ktetra(2,5)  = kg25*ls1
       ktetra(2,6)  = kg35*ls1
       ktetra(2,7)  = kg45*la1la
       ktetra(2,8)  = kg55*la1la
       ktetra(2,9)  = kg56*la1la
       ktetra(2,10) = kg57*ls1
       ktetra(2,11) = kg58*ls1
       ktetra(2,12) = kg59*ls1
c
c
       ktetra(3,3)  = kg66*lsq1
       ktetra(3,4)  = kg16*ls1
       ktetra(3,5)  = kg26*ls1
       ktetra(3,6)  = kg36*ls1
       ktetra(3,7)  = kg46*la1la
       ktetra(3,8)  = kg56*la1la
       ktetra(3,9)  = kg66*la1la
       ktetra(3,10) = kg67*ls1
       ktetra(3,11) = kg68*ls1
       ktetra(3,12) = kg69*ls1
c
c

       ktetra(4,7)  = kg14*lambda
       ktetra(4,8)  = kg15*lambda
       ktetra(4,9)  = kg16*lambda

c
c

       ktetra(5,7)  = kg24*lambda
       ktetra(5,8)  = kg25*lambda
       ktetra(5,9)  = kg26*lambda

c
c

       ktetra(6,7)  = kg34*lambda
       ktetra(6,8)  = kg35*lambda
       ktetra(6,9)  = kg36*lambda

c
c

       ktetra(7,7)  = kg44*lasq
       ktetra(7,8)  = kg45*lasq
       ktetra(7,9)  = kg46*lasq
       ktetra(7,10) = kg47*lambda
       ktetra(7,11) = kg48*lambda
       ktetra(7,12) = kg49*lambda
c
c
       ktetra(8,8)  = kg55*lasq
       ktetra(8,9)  = kg56*lasq
       ktetra(8,10) = kg57*lambda
       ktetra(8,11) = kg58*lambda
       ktetra(8,12) = kg59*lambda
c
c
       ktetra(9,9)  = kg66*lasq
       ktetra(9,10) = kg67*lambda
       ktetra(9,11) = kg68*lambda
       ktetra(9,12) = kg69*lambda
c
c
 
c
c      *********************************************************************
c
      p = map_pointer(face,cycle,1) 
      q = map_pointer(face,cycle,2) 
      r = map_pointer(face,cycle,3)
      s = map_pointer(face,cycle,4)
c
c
c     p-p block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(p-1)
      the_col = 3*(p-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(1,1)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(1,2)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(1,3)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(1,2)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(2,2)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(2,3)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(1,3)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(2,3)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(3,3)
c
c     p-q block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(p-1)
      the_col = 3*(q-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(1,4)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(1,5)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(1,6)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(2,4)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(2,5)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(2,6)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(3,4)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(3,5)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(3,6)
c
c     p-r block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(p-1)
      the_col = 3*(r-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(1,7)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(1,8)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(1,9)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(2,7)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(2,8)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(2,9)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(3,7)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(3,8)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(3,9)
c
c     p-s block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(p-1)
      the_col = 3*(s-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(1,10)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(1,11)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(1,12)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(2,10)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(2,11)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(2,12)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(3,10)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(3,11)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(3,12)
c
c     q-p block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(q-1)
      the_col = 3*(p-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(1,4)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(2,4)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(3,4)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(1,5)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(2,5)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(3,5)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(1,6)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(2,6)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(3,6)
c
c     q-q block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(q-1)
      the_col = 3*(q-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + kg11 
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + kg12
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + kg13
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + kg12  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + kg22
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + kg23
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + kg13  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + kg23
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + kg33
c
c     q-r block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(q-1)
      the_col = 3*(r-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(4,7)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(4,8)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(4,9)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(5,7)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(5,8)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(5,9)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(6,7)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(6,8)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(6,9)
c
c     q-s block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(q-1)
      the_col = 3*(s-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + kg17 
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + kg18
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + kg19
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + kg27  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + kg28
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + kg29
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + kg37 
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + kg38
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + kg39
c
c     r-p block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(r-1)
      the_col = 3*(p-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(1,7)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(2,7)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(3,7)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(1,8)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(2,8)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(3,8)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(1,9)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(2,9)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(3,9)
c
c     r-q block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(r-1)
      the_col = 3*(q-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(4,7)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(5,7)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(6,7)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(4,8)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(5,8)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(6,8)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(4,9)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(5,9)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(6,9)
c
c     r-r block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(r-1)
      the_col = 3*(r-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(7,7)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(7,8)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(7,9)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(7,8)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(8,8)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(8,9)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(7,9)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(8,9)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(9,9)
c
c     r-s block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(r-1)
      the_col = 3*(s-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(7,10)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(7,11)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(7,12)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(8,10)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(8,11)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(8,12)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(9,10)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(9,11)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(9,12)
c
c     s-p block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(s-1)
      the_col = 3*(p-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(1,10)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(2,10)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(3,10)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(1,11)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(2,11)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(3,11)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(1,12)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(2,12)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(3,12)
c
c     s-q block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(s-1)
      the_col = 3*(q-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + kg17 
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + kg27
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + kg37
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + kg18  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + kg28
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + kg38
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + kg19 
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + kg29
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + kg39
c
c     s-r block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(s-1)
      the_col = 3*(r-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + ktetra(7,10)  
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + ktetra(8,10)
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + ktetra(9,10)
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + ktetra(7,11)  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + ktetra(8,11)
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + ktetra(9,11)
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + ktetra(7,12)  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + ktetra(8,12)
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + ktetra(9,12)
c
c
c     s-s block
c
c     set matrix pointer in upper left corner of current block
c
c
      the_row = 3*(s-1)
      the_col = 3*(s-1)
c
      row_1   = the_row + 1
      row_2   = the_row + 2
      row_3   = the_row + 3
c
      col_1   = the_col + 1
      col_2   = the_col + 2
      col_3   = the_col + 3
c
      ktetglob(row_1,col_1) = ktetglob(row_1,col_1) + kg77 
      ktetglob(row_1,col_2) = ktetglob(row_1,col_2) + kg78
      ktetglob(row_1,col_3) = ktetglob(row_1,col_3) + kg79
c
      ktetglob(row_2,col_1) = ktetglob(row_2,col_1) + kg78  
      ktetglob(row_2,col_2) = ktetglob(row_2,col_2) + kg88
      ktetglob(row_2,col_3) = ktetglob(row_2,col_3) + kg89
c
      ktetglob(row_3,col_1) = ktetglob(row_3,col_1) + kg79  
      ktetglob(row_3,col_2) = ktetglob(row_3,col_2) + kg89
      ktetglob(row_3,col_3) = ktetglob(row_3,col_3) + kg99
c
c     end of loop over cycles
c     |
      enddo
c
c     end of loop over faces
c     |
      enddo
c
c     assign current element stiffness matrix to kteti_j(ele)
c
c
       ktet1_1(ele)  = ktetglob(1,1)
       ktet1_2(ele)  = ktetglob(1,2)
       ktet1_3(ele)  = ktetglob(1,3)
       ktet1_4(ele)  = ktetglob(1,4)
       ktet1_5(ele)  = ktetglob(1,5)
       ktet1_6(ele)  = ktetglob(1,6)
       ktet1_7(ele)  = ktetglob(1,7)
       ktet1_8(ele)  = ktetglob(1,8)
       ktet1_9(ele)  = ktetglob(1,9)
       ktet1_10(ele) = ktetglob(1,10)
       ktet1_11(ele) = ktetglob(1,11)
       ktet1_12(ele) = ktetglob(1,12)
c
       ktet2_2(ele)  = ktetglob(2,2)
       ktet2_3(ele)  = ktetglob(2,3)
       ktet2_4(ele)  = ktetglob(2,4)
       ktet2_5(ele)  = ktetglob(2,5)
       ktet2_6(ele)  = ktetglob(2,6)
       ktet2_7(ele)  = ktetglob(2,7)
       ktet2_8(ele)  = ktetglob(2,8)
       ktet2_9(ele)  = ktetglob(2,9)
       ktet2_10(ele)  = ktetglob(2,10)
       ktet2_11(ele)  = ktetglob(2,11)
       ktet2_12(ele)  = ktetglob(2,12)
c
       ktet3_3(ele)  = ktetglob(3,3)
       ktet3_4(ele)  = ktetglob(3,4)
       ktet3_5(ele)  = ktetglob(3,5)
       ktet3_6(ele)  = ktetglob(3,6)
       ktet3_7(ele)  = ktetglob(3,7)
       ktet3_8(ele)  = ktetglob(3,8)
       ktet3_9(ele)  = ktetglob(3,9)
       ktet3_10(ele)  = ktetglob(3,10)
       ktet3_11(ele)  = ktetglob(3,11)
       ktet3_12(ele)  = ktetglob(3,12)
c
       ktet4_4(ele)  = ktetglob(4,4)
       ktet4_5(ele)  = ktetglob(4,5)
       ktet4_6(ele)  = ktetglob(4,6)
       ktet4_7(ele)  = ktetglob(4,7)
       ktet4_8(ele)  = ktetglob(4,8)
       ktet4_9(ele)  = ktetglob(4,9)
       ktet4_10(ele)  = ktetglob(4,10)
       ktet4_11(ele)  = ktetglob(4,11)
       ktet4_12(ele)  = ktetglob(4,12)
c
       ktet5_5(ele)  = ktetglob(5,5)
       ktet5_6(ele)  = ktetglob(5,6)
       ktet5_7(ele)  = ktetglob(5,7)
       ktet5_8(ele)  = ktetglob(5,8)
       ktet5_9(ele)  = ktetglob(5,9)
       ktet5_10(ele)  = ktetglob(5,10)
       ktet5_11(ele)  = ktetglob(5,11)
       ktet5_12(ele)  = ktetglob(5,12)
c
       ktet6_6(ele)  = ktetglob(6,6)
       ktet6_7(ele)  = ktetglob(6,7)
       ktet6_8(ele)  = ktetglob(6,8)
       ktet6_9(ele)  = ktetglob(6,9)
       ktet6_10(ele)  = ktetglob(6,10)
       ktet6_11(ele)  = ktetglob(6,11)
       ktet6_12(ele)  = ktetglob(6,12)
c
       ktet7_7(ele)  = ktetglob(7,7)
       ktet7_8(ele)  = ktetglob(7,8)
       ktet7_9(ele)  = ktetglob(7,9)
       ktet7_10(ele)  = ktetglob(7,10)
       ktet7_11(ele)  = ktetglob(7,11)
       ktet7_12(ele)  = ktetglob(7,12)
c
       ktet8_8(ele)  = ktetglob(8,8)
       ktet8_9(ele)  = ktetglob(8,9)
       ktet8_10(ele)  = ktetglob(8,10)
       ktet8_11(ele)  = ktetglob(8,11)
       ktet8_12(ele)  = ktetglob(8,12)
c
       ktet9_9(ele)  = ktetglob(9,9)
       ktet9_10(ele)  = ktetglob(9,10)
       ktet9_11(ele)  = ktetglob(9,11)
       ktet9_12(ele)  = ktetglob(9,12)
c
       ktet10_10(ele)  = ktetglob(10,10)
       ktet10_11(ele)  = ktetglob(10,11)
       ktet10_12(ele)  = ktetglob(10,12)
c
       ktet11_11(ele)  = ktetglob(11,11)
       ktet11_12(ele)  = ktetglob(11,12)
c
       ktet12_12(ele)  = ktetglob(12,12)
c
   
c
c     end of loop over local tetra elements
c     |
      enddo
c
c
      return

      end





