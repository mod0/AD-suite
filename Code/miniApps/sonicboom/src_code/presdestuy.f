





      SUBROUTINE PRESDESTUY(ctrl)
C
c   Calcul de la pression desiree si DEFORTUY=0 (i.e que l'on calcule la 
c   pression sur une deformation simulee par transpiration : ITRANS=1, contr=1)
c
c   Comme on a mis contr=1, CTRL=f(x)-5. On va alors le remettre a 0.
c
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c
      INTEGER ISP, IS,ii 
      REAL*8 ctrl(nnsp)
      REAL*8 sol(nsmax,8),coefr(5),vcarre,denom
C
      do isp = 1,nsp
         is = node2D3D(isp)
         pdesp(is) = gam1*(ua(5,is) - 0.5*(ua(2,is)**2 +
     $        ua(3,is)**2 + ua(4,is)**2)/ua(1,is))
         ctrl(isp) = 0.
      end do
C

      coefr(1)                     = rhoref
      coefr(2)                     = rhoref*vref
      coefr(3)                     = rhoref*vref
      coefr(4)                     = rhoref*vref
      coefr(5)                     = pref
c
      DO is = 1,ns
         sol(is,1) = ua(1,is)
         sol(is,2) = ua(2,is)/ua(1,is)
         sol(is,3) = ua(3,is)/ua(1,is)
         sol(is,4) = ua(4,is)/ua(1,is)
         sol(is,5) = ua(5,is)
         sol(is,6) = gam1*(sol(is,5)*coefr(5) - 0.5*sol(is,1)*(sol(is,2)
     &        *sol(is,2) + sol(is,3)*sol(is,3) + sol(is,4)*sol(is,4)))
         sol(is,7) = SQRT( (sol(is,2)*sol(is,2)+sol(is,3)*sol(is,3)+
     &               sol(is,4)*sol(is,4))/(gam*sol(is,6)/sol(is,1)))
         vcarre = sol(is,2)**2+sol(is,3)**2+sol(is,4)
         denom = roin*0.5*vcarre
         sol(is,8) = (pin - sol(is,6))/denom
      END DO
c
c     sol(is,1) ---> Rho
c     sol(is,2) ---> u
c     sol(is,3) ---> v
c     sol(is,4) ---> w
c     sol(is,5) ---> E
c     sol(is,6) ---> P
c     sol(is,7) ---> Mach
c     sol(is,8) ---> Cp
c
      Do ii = 1,8
         Do is = 1,ns
           WRITE(21,113) sol(is,ii)
         End Do
      End Do
c
113   FORMAT(e15.8)

      RETURN
      END
