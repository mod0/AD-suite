      
                  
      SUBROUTINE CMVFAC
c---------------------------------------------------------------------
c Computes/updates the cell and tetraedra volumes  
c Computes/updates the mesh tetraedra and cell volumes 
c Computes/updates the normals at the physical boundary faces
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c      INCLUDE 'Paramopt3D.h'
c---------------------------------------------------------------------
c     Local variables definition  
      INTEGER k , is , jt, ks  , ifac
      INTEGER nsf1(3),  nsf2(3), nsf3(3), nsf4(3), nut(4) , nsf(4,3)
      REAL*8    x(4)   , y(4)    , z(4)
      REAL*8    xx(3,2), yy(3,2) , zz(3,2), xx12(4), yy12(4), zz12(4)
      REAL*8    xx13(4), yy13(4) , zz13(4)
      REAL*8    xet1   , xet2    , xet3
      REAL*8    z12    , z13     , z14    , z23    , z24    , z34
      REAL*8    b1, b2 , b3, b4     
      REAL*8    dbij11 , dbij12  , dbij13 , dbij21 , dbij22 
      REAL*8    dbij23 , dbij31  , dbij32 , dbij33 , dbij41
      REAL*8    dbij42 , dbij43
      REAL*8    us3    , us6     , us3dt  , eps
      REAL*8    vol6   , xww1    , xww2   , xww3
c
c     Initialisations
c
      nsf1(1)                      = 1
      nsf1(2)                      = 2
      nsf1(3)                      = 3
c
      nsf2(1)                      = 1
      nsf2(2)                      = 4
      nsf2(3)                      = 2
c
      nsf3(1)                      = 1
      nsf3(2)                      = 3
      nsf3(3)                      = 4
c
      nsf4(1)                      = 2
      nsf4(2)                      = 4
      nsf4(3)                      = 3
c
      dt                           = 0.1
c
      us3                          = 1.0/3.0
      us3dt                        = 1.0/(dt*3.0)
      us6                          = 1.0/6.0
      eps                          = 1.0e-6
c
      DO 1 k=1,3
         nsf(1,k)                  = nsf1(k)
         nsf(2,k)                  = nsf2(k)
         nsf(3,k)                  = nsf3(k)
         nsf(4,k)                  = nsf4(k)
1     CONTINUE
c
      DO 5 is=1,ns
         vols(is)                  = 0.0
         coco(1,is)                = coor(1,is) 
         coco(2,is)                = coor(2,is) 
         coco(3,is)                = coor(3,is)
5     CONTINUE
c
      DO 100 jt=1,nt
c
         DO 10 ks=1,4
c
            nut(ks)                = nu(ks,jt)
c
            is                     = nut(ks)
c
            x(ks)                  = coco(1,is)
            y(ks)                  = coco(2,is)
            z(ks)                  = coco(3,is)
c
10       CONTINUE
c
         z12                       = z(2) - z(1)
         z13                       = z(3) - z(1)
         z14                       = z(4) - z(1)
         z23                       = z(3) - z(2)
         z24                       = z(4) - z(2)
         z34                       = z(4) - z(3)
c
         b1                        = y(2)*z34 - y(3)*z24 + y(4)*z23
         b2                        =-y(1)*z34 + y(3)*z14 - y(4)*z13
         b3                        = y(1)*z24 - y(2)*z14 + y(4)*z12
         b4                        =-y(1)*z23 + y(2)*z13 - y(3)*z12
c
         vol6                      = x(1)*b1 + x(2)*b2 + 
     &                               x(3)*b3 + x(4)*b4
c
c        Updates the volumes of the cells attached to the  
c        current tetraedra vertices
c
         DO 20 ks=1,4
            vols(nut(ks))         = vols(nut(ks)) + vol6/24.0
20       CONTINUE
c
c        Updates the current tetraedra volume 
c
         volt(jt)                  = us6*vol6
c
100   CONTINUE
c
      DO 200 ifac=1,nfac
c
         xx(1,1)                   = coor(1,nsfac(1,ifac))
         yy(1,1)                   = coor(2,nsfac(1,ifac))
         zz(1,1)                   = coor(3,nsfac(1,ifac))
c
         xx(1,2)                   = coco(1,nsfac(1,ifac))
         yy(1,2)                   = coco(2,nsfac(1,ifac))
         zz(1,2)                   = coco(3,nsfac(1,ifac))
c
         xx(2,1)                   = coor(1,nsfac(2,ifac))
         yy(2,1)                   = coor(2,nsfac(2,ifac))
         zz(2,1)                   = coor(3,nsfac(2,ifac))
c
         xx(2,2)                   = coco(1,nsfac(2,ifac))
         yy(2,2)                   = coco(2,nsfac(2,ifac))
         zz(2,2)                   = coco(3,nsfac(2,ifac))
c
         xx(3,1)                   = coor(1,nsfac(3,ifac))
         yy(3,1)                   = coor(2,nsfac(3,ifac))
         zz(3,1)                   = coor(3,nsfac(3,ifac))
c
         xx(3,2)                   = coco(1,nsfac(3,ifac))
         yy(3,2)                   = coco(2,nsfac(3,ifac))
         zz(3,2)                   = coco(3,nsfac(3,ifac))
c
         xww1                      = us3dt*(xx(1,2) - xx(1,1) + 
     &                                      xx(2,2) - xx(2,1) + 
     &                                      xx(3,2) - xx(3,1))
c
         xww2                      = us3dt*(yy(1,2) - yy(1,1) + 
     &                                      yy(2,2) - yy(2,1) + 
     &                                      yy(3,2) - yy(3,1))
c
         xww3                      = us3dt*(zz(1,2) - zz(1,1) + 
     &                                      zz(2,2) - zz(2,1) + 
     &                                      zz(3,2) - zz(3,1))
c
         xx12(1)                   = xx(2,1) - xx(1,1) 
         xx13(1)                   = xx(3,1) - xx(1,1) 
c
         yy12(1)                   = yy(2,1) - yy(1,1)
         yy13(1)                   = yy(3,1) - yy(1,1) 
c
         zz12(1)                   = zz(2,1) - zz(1,1) 
         zz13(1)                   = zz(3,1) - zz(1,1)
c 
         xx12(2)                   = xx(2,2) - xx(1,2) 
         xx13(2)                   = xx(3,2) - xx(1,2) 
c
         yy12(2)                   = yy(2,2) - yy(1,2) 
         yy13(2)                   = yy(3,2) - yy(1,2) 
c
         zz12(2)                   = zz(2,2) - zz(1,2) 
         zz13(2)                   = zz(3,2) - zz(1,2) 
c
         dbij11                    = yy13(2)*zz12(1) - yy12(1)*zz13(2)
         dbij21                    = yy13(1)*zz12(2) - yy12(2)*zz13(1)
         dbij31                    = yy13(2)*zz12(2) - yy12(2)*zz13(2)
         dbij41                    = yy13(1)*zz12(1) - yy12(1)*zz13(1)
c
         dbij12                    =-xx13(2)*zz12(1) + xx12(1)*zz13(2)
         dbij22                    =-xx13(1)*zz12(2) + xx12(2)*zz13(1)
         dbij32                    =-xx13(2)*zz12(2) + xx12(2)*zz13(2)
         dbij42                    =-xx13(1)*zz12(1) + xx12(1)*zz13(1)
c
         dbij13                    = xx13(2)*yy12(1) - xx12(1)*yy13(2)
         dbij23                    = xx13(1)*yy12(2) - xx12(2)*yy13(1)
         dbij33                    = xx13(2)*yy12(2) - xx12(2)*yy13(2)
         dbij43                    = xx13(1)*yy12(1) - xx12(1)*yy13(1)
c
c        
         xet1                      = us3*(dbij31 + dbij41 + 
     &                                    0.5*(dbij11 + dbij21))
         xet2                      = us3*(dbij32 + dbij42 + 
     &                                    0.5*(dbij12 + dbij22))
         xet3                      = us3*(dbij33 + dbij43 + 
     &                                    0.5*(dbij13 + dbij23))
c
         vnfac(1,ifac)             = us6*xet1
         vnfac(2,ifac)             = us6*xet2
         vnfac(3,ifac)             = us6*xet3
c
         vnsig(ifac)               = us6*(xet1*xww1 + xet2*xww2 + 
     &                                    xet3*xww3)
c
200   CONTINUE
C
      RETURN
      END
