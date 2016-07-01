      SUBROUTINE AEROF1(cl)
c     -----------------------------------------------------------------
c     Computation of the lift coefficient in the case of a 
c     non-viscous flow
c     -----------------------------------------------------------------
      INCLUDE 'Param3D.h'
c     -----------------------------------------------------------------
      REAL*8     cl
c     Local variables definition
      INTEGER  is1   , is2 , is3 , if1 , ifac
      REAL*8     pm    , ds2
      REAL*8     rhom  , rhum, rhvm, rhwm, rhem  
c
      ds2                          = 2.0
c
c     Looping over the submesh faces placed ont the body surface
c
      cl                           = 0.0
c
      DO 100 if1=1,nf1
c
         ifac                      = noe1(nf11+if1)
c
         is1                       = nsfac(1,ifac)
         is2                       = nsfac(2,ifac)
         is3                       = nsfac(3,ifac)
c
         rhom                      = ua(1,is1) 
         rhum                      = ua(2,is1)
         rhvm                      = ua(3,is1) 
         rhwm                      = ua(4,is1)
         rhem                      = ua(5,is1) 
c   
         pm                        = gam1*(rhem - 
     &                                     0.5*(rhum*rhum + rhvm*rhvm +
     &                                         rhwm*rhwm)/rhom)
c
         cl                        = cl + ds2*vnfac(2,ifac)*pm
c
         rhom                      = ua(1,is2) 
         rhum                      = ua(2,is2)
         rhvm                      = ua(3,is2) 
         rhwm                      = ua(4,is2)
         rhem                      = ua(5,is2) 
c
         pm                        = gam1*(rhem - 
     &                                     0.5*(rhum*rhum + rhvm*rhvm +
     &                                         rhwm*rhwm)/rhom)
c
         cl                        = cl + ds2*vnfac(2,ifac)*pm
c
         rhom                      = ua(1,is3) 
         rhum                      = ua(2,is3)
         rhvm                      = ua(3,is3) 
         rhwm                      = ua(4,is3)
         rhem                      = ua(5,is3) 
c
         pm                        = gam1*(rhem - 
     &                                     0.5*(rhum*rhum + rhvm*rhvm +
     &                                         rhwm*rhwm)/rhom)
c
         cl                        = cl + ds2*vnfac(2,ifac)*pm
c
100   CONTINUE
c
      RETURN
      END
