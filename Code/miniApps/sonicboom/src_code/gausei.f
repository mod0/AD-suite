
      SUBROUTINE GAUSEI
c---------------------------------------------------------------------   
c Solves the linear implict system using the Gauss-Seidel method
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
c     Local variables definition
      INTEGER i, j, icc , ia   , ks  , is, id1
      INTEGER iseg, isgs, incgs, nsg , nsg2 
      REAL*8    x1  , x2  , x3   , x4  , x5
      REAL*8    res , res0, alpha
c
      alpha                        = 1.0
c
      DO 15 ia=1,5
c
         DO 20 is=1,ns
c
            dy(ia,is)              = diag(is,ia,1)*ce(1,is) + 
     &                               diag(is,ia,2)*ce(2,is) +
     &                               diag(is,ia,3)*ce(3,is) + 
     &                               diag(is,ia,4)*ce(4,is) +
     &                               diag(is,ia,5)*ce(5,is)
c
20       CONTINUE
c
15    CONTINUE
c
      IF (nbrel .GT. 0) THEN
c
         DO 300 nit=1,nbrel
c
            res                    = 0.0
c
            DO 25 ia=1,5
c
               DO 30 is=1,ns 
c
                  dx(ia,is)        = ce(ia,is)
c
30             CONTINUE
c
25          CONTINUE
c
            DO 65 is=1,ns
               dx(5,is)            = adh1(is)*dx(5,is)
65          CONTINUE  
c
            IF (MOD(nit, 2) .EQ. 0) THEN
               isgs                = ns + 1
               incgs               =-1
            ELSE
               isgs                = 0
               incgs               = 1
            ENDIF
c
            DO 35 ks=1,ns
c
               isgs                = isgs + incgs
               nsg2                = ndeg(isgs)
c
               DO 40 nsg=1,nsg2
c
                  iseg             = jaret(isgs,nsg)
c
                  icc              = 0
c
                  IF (isgs .EQ. nubo(1,iseg)) THEN
                     icc           = 2
                     id1           = 50*(iseg - 1) + 25
                  ELSE
                     icc           = 1
                     id1           = 50*(iseg - 1)
                  ENDIF
c
                  i                = nubo(3-icc,iseg)
                  j                = nubo(icc,iseg)
c
                  dx(1,i)          = dx(1,i) + stmat(id1+1)*dy(1,j)    +
     &                                         stmat(id1+2)*dy(2,j)    +
     &                                         stmat(id1+3)*dy(3,j)    +
     &                                         stmat(id1+4)*dy(4,j)    +
     &                                         stmat(id1+5)*dy(5,j)
c
                  dx(2,i)          = dx(2,i) + adh1(i)*(
     &                                         stmat(id1+5+1)*dy(1,j)  +
     &                                         stmat(id1+5+2)*dy(2,j)  +
     &                                         stmat(id1+5+3)*dy(3,j)  +
     &                                         stmat(id1+5+4)*dy(4,j)  +
     &                                         stmat(id1+5+5)*dy(5,j))
c     
                  dx(3,i)          = dx(3,i) + adh1(i)*(
     &                                         stmat(id1+10+1)*dy(1,j) +
     &                                         stmat(id1+10+2)*dy(2,j) +
     &                                         stmat(id1+10+3)*dy(3,j) +
     &                                         stmat(id1+10+4)*dy(4,j) +
     &                                         stmat(id1+10+5)*dy(5,j)) 
c
                  dx(4,i)          = dx(4,i) + adh1(i)*(
     &                                         stmat(id1+15+1)*dy(1,j) +
     &                                         stmat(id1+15+2)*dy(2,j) +
     &                                         stmat(id1+15+3)*dy(3,j) +
     &                                         stmat(id1+15+4)*dy(4,j) +
     &                                         stmat(id1+15+5)*dy(5,j))
c
                  dx(5,i)          = dx(5,i) + adh1(i)*(
     &                                         stmat(id1+20+1)*dy(1,j) + 
     &                                         stmat(id1+20+2)*dy(2,j) +
     &                                         stmat(id1+20+3)*dy(3,j) +
     &                                         stmat(id1+20+4)*dy(4,j) +
     &                                         stmat(id1+20+5)*dy(5,j))
c
40             CONTINUE
c
               is                  = isgs
c
               x1                  = diag(is,1,1)*dx(1,is) + 
     &       diag(is,1,2)*dx(2,is) + diag(is,1,3)*dx(3,is) + 
     &       diag(is,1,4)*dx(4,is) + diag(is,1,5)*dx(5,is)
               x2                  = diag(is,2,1)*dx(1,is) + 
     &       diag(is,2,2)*dx(2,is) + diag(is,2,3)*dx(3,is) + 
     &       diag(is,2,4)*dx(4,is) + diag(is,2,5)*dx(5,is)
               x3                  = diag(is,3,1)*dx(1,is) + 
     &       diag(is,3,2)*dx(2,is) + diag(is,3,3)*dx(3,is) + 
     &       diag(is,3,4)*dx(4,is) + diag(is,3,5)*dx(5,is)
               x4                  = diag(is,4,1)*dx(1,is) + 
     &       diag(is,4,2)*dx(2,is) + diag(is,4,3)*dx(3,is) + 
     &       diag(is,4,4)*dx(4,is) + diag(is,4,5)*dx(5,is)
               x5                  = diag(is,5,1)*dx(1,is) + 
     &       diag(is,5,2)*dx(2,is) + diag(is,5,3)*dx(3,is) + 
     &       diag(is,5,4)*dx(4,is) + diag(is,5,5)*dx(5,is)
c
               dx(1,is)            = x1
               dx(2,is)            = x2
               dx(3,is)            = x3
               dx(4,is)            = x4
               dx(5,is)            = x5
c 
c              Computing the residual 
c
               res                 = res + 
     &                               (dx(1,is) - dy(1,is))**2 + 
     &                               (dx(2,is) - dy(2,is))**2 +
     &                               (dx(3,is) - dy(3,is))**2 + 
     &                               (dx(4,is) - dy(4,is))**2 +
     &                               (dx(5,is) - dy(5,is))**2
c
               dy(1,is)            = (1.0 - alpha)*dy(1,is) + 
     &                               alpha*dx(1,is)
               dy(2,is)            = (1.0 - alpha)*dy(2,is) + 
     &                               alpha*dx(2,is)
               dy(3,is)            = (1.0 - alpha)*dy(3,is) + 
     &                               alpha*dx(3,is)
               dy(4,is)            = (1.0 - alpha)*dy(4,is) + 
     &                               alpha*dx(4,is)
               dy(5,is)            = (1.0 - alpha)*dy(5,is) + 
     &                               alpha*dx(5,is)
c
35          CONTINUE
c
            IF (nit .EQ. 1) res0   = res
c
            res                    = SQRT(res/res0)
c
            IF ((nit .EQ. nbrel) .OR. (res .LT. err)) GOTO 310
c
300      CONTINUE
c
      ENDIF
c
310   CONTINUE
c
c      WRITE(17, 1200) cfl, nit, res
1200  FORMAT(f12.2,2x,i6,2x,e12.5)
c
      RETURN
      END
