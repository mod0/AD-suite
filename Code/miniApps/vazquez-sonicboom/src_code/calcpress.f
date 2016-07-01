

      SUBROUTINE calcpress
c
c *** Calcul des pressions 

      include 'Param3D.h'
c
      INTEGER is
      REAL*8 usro,pres(nsmax)

      DO 1 is = 1,ns

         usro=1.0/ua(1,is)
         pres(is)=gam1*(ua(5,is)-0.5*(ua(2,is)*ua(2,is)
     1              +ua(3,is)*ua(3,is)+ua(4,is)*ua(4,is))*usro)

1     CONTINUE

c      pression          : p= pres(is)
c      densite           : rho= ua(1,is)
c      vitesse suivant x : u = ua(2,is)/ro
c      vitesse suivant y : v = ua(3,is)/ro
c      vitesse suivant z : w = ua(4,is)/ro
c      Energie           : e = ua(5,is)

      RETURN
      END
