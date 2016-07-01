
      SUBROUTINE ECRITVIGIE(numvigie)
C
      INCLUDE 'Param3D.h'
C
      INTEGER is, ii, numvigie
      REAL*8 sol(nsmax,9),coefr(5),vcarre,denom
C
c      Le numero de fichier numvigie debute a fort.22. Il est ensuite
c      incremente de 1.
c      Ecriture dans numvigie des solutions au format hdb
c      pour la visualisation avec VIGIE :
c      rho, u, v, w, E, Pression, Mach, Cp, Entropie.
c

      
      write(6,*)
     .  'Subroutine ECRITVIGIE returned.. writting only in GID files'
      
      return
      
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
         sol(is,9) = (sol(is,6)/pin)*(roin/sol(is,1))**gam - 1.
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
c     sol(is,9) ---> Entropie
c
      Do ii = 1,9
         Do is = 1,ns
           WRITE(numvigie,113) sol(is,ii)
         End Do
      End Do
c
113   FORMAT(e15.8)
c
c      DO 1000 is=1,ns
c         WRITE(58, 115) is, ua(1,is), ua(2,is), ua(3,is),
c     &                      ua(4,is), ua(5,is)
c1000  CONTINUE
c  
c      WRITE(58, 112) kt, t
c
c112   FORMAT(i6,2x,e12.5)
c115   FORMAT(i6,5(e15.8,1x))
c
      return
      end
