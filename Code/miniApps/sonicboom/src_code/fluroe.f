      SUBROUTINE FLUROE
c---------------------------------------------------------------------   
c Computes the convective fluxes using the approximate Riemann
c solver of Roe adapted to dynamic meshes
c A Van Albada limiter is used in the M.U.S.C.L. procedure
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
c     Local variables definition
      INTEGER is     , iseg   , nubo1  , nubo2
      REAL*8    gamo   , usg0   , pow    , coeff  , dtpred 
      REAL*8    ro     , usro   , u ,v, w, p, rnorm
      REAL*8    aix    , aiy    , aiz    , delta  , switch
      REAL*8    squsr1 , squsr2
      REAL*8    beta2  , beta3  , beta
      REAL*8    dpm    , dpor   , dpex   , aux1   , aux2 , e2 
      REAL*8    xsigm  , vp1    , vp4    , vp5
      REAL*8    uas1(2), uas2(2), uas3(2), uas4(2), uas5(2)
      REAL*8    vno(3)
      REAL*8    prod1  , prod2  , bsign
      REAL*8    tet1   , tet2   , tet3
      REAL*8    cr     , cr2    , ener1  , ener2  , pror , qir
      REAL*8    uar1   , uar2   , uar3   , uar4   , uar5
      REAL*8    dif1   , dif2   , dif3   , dif4   , dif5 
      REAL*8    flur1  , flur2  , flur3  , flur4  , flur5
      REAL*8    fltr1  , fltr2  , fltr3  , fltr4  , fltr5
      REAL*8    flum1  , flum2  , flum3  , flum4  , flum5
c
      bsign                        = 1.0d0
c
      IF (nordre .EQ. 2) bsign     =-1.0d0   
c
      beta                         = 0.5d0
c
      beta2                        = beta
      beta3                        = 0.5d0*(1.0d0 - 2.0d0*beta)
c
      e2                           = 1.0d-16
c
      gamo                         = 0.5d0*(gam - 1.0d0)
      usg0                         = 1.0d0/gamo
      pow                          = 1.0d0/(2.0d0*gam)
      coeff                        = gam1/(gam + 1.0d0)
c
C$AD II-LOOP
      DO 5 is=1,ns
c
         dtpred                    = 0.5d0*dtl(is)*ipred
c
         ro                        = ua(1,is)
         usro                      = 1.0d0/ro
         u                         = ua(2,is)*usro
         v                         = ua(3,is)*usro
         w                         = ua(4,is)*usro
         p                         = gam1*(ua(5,is) - 
     &                                     0.5d0*ro*(u*u + v*v + w*w))
c
         un(1,is)                  = ro - dtpred*(u*dx(1,is) + 
     &                                            v*dy(1,is) + 
     &                                            w*dz(1,is) + 
     &                             ro*(dx(2,is) + dy(3,is) + dz(4,is)))
         un(2,is)                  = u  - dtpred*(u*dx(2,is) + 
     &                                            v*dy(2,is) + 
     &                                            w*dz(2,is) + 
     &                                            dx(5,is)*usro)
         un(3,is)                  = v  - dtpred*(u*dx(3,is) + 
     &                                            v*dy(3,is) + 
     &                                            w*dz(3,is) +
     &                                            dy(5,is)*usro)
         un(4,is)                  = w  - dtpred*(u*dx(4,is) + 
     &                                            v*dy(4,is) + 
     &                                            w*dz(4,is) +
     &                                            dz(5,is)*usro)
         un(5,is)                  = p  - dtpred*(u*dx(5,is) + 
     &                                            v*dy(5,is) + 
     &                                            w*dz(5,is) +
     &                                            gam*p*(dx(2,is) + 
     &                                                   dy(3,is) + 
     &                                                   dz(4,is)))
c
5     CONTINUE
c
C$AD II-LOOP
      DO 10 iseg=1,nseg
c
         nubo1                     = nubo(1,iseg)
         nubo2                     = nubo(2,iseg)
c
c        Indirect addressing on vertices physical states
c
         uas1(1)                   = un(1,nubo1)
         uas1(2)                   = un(1,nubo2)
         uas2(1)                   = un(2,nubo1)
         uas2(2)                   = un(2,nubo2)
         uas3(1)                   = un(3,nubo1)
         uas3(2)                   = un(3,nubo2)
         uas4(1)                   = un(4,nubo1)
         uas4(2)                   = un(4,nubo2)
         uas5(1)                   = un(5,nubo1)
         uas5(2)                   = un(5,nubo2)
c
         flur1                     = 0.0d0
         flur2                     = 0.0d0
         flur3                     = 0.0d0
         flur4                     = 0.0d0
         flur5                     = 0.0d0
c         
         fltr1                     = 0.0d0
         fltr2                     = 0.0d0
         fltr3                     = 0.0d0
         fltr4                     = 0.0d0
         fltr5                     = 0.0d0
c
         IF (nordre .EQ. 1) GOTO 275
c
         aix                       = coor(1,nubo2) - coor(1,nubo1)
         aiy                       = coor(2,nubo2) - coor(2,nubo1)
         aiz                       = coor(3,nubo2) - coor(3,nubo1)
c
         flur1                     = beta2*
     &         (aix*dx(1,nubo1) + aiy*dy(1,nubo1) + aiz*dz(1,nubo1)) +
     &                               beta3*(uas1(2) - uas1(1))
         flur2                     = beta2*
     &         (aix*dx(2,nubo1) + aiy*dy(2,nubo1) + aiz*dz(2,nubo1)) +
     &                               beta3*(uas2(2) - uas2(1))
         flur3                     = beta2*
     &         (aix*dx(3,nubo1) + aiy*dy(3,nubo1) + aiz*dz(3,nubo1)) + 
     &                               beta3*(uas3(2) - uas3(1))
         flur4                     = beta2*
     &         (aix*dx(4,nubo1) + aiy*dy(4,nubo1) + aiz*dz(4,nubo1)) +
     &                               beta3*(uas4(2) - uas4(1))
         flur5                     = beta2*
     &         (aix*dx(5,nubo1) + aiy*dy(5,nubo1) + aiz*dz(5,nubo1)) +
     &                               beta3*(uas5(2) - uas5(1))
c
         fltr1                     = beta2*
     &         (aix*dx(1,nubo2) + aiy*dy(1,nubo2) + aiz*dz(1,nubo2)) +
     &                               beta3*(uas1(2) - uas1(1))
         fltr2                     = beta2*
     &         (aix*dx(2,nubo2) + aiy*dy(2,nubo2) + aiz*dz(2,nubo2)) +
     &                               beta3*(uas2(2) - uas2(1))
         fltr3                     = beta2*
     &         (aix*dx(3,nubo2) + aiy*dy(3,nubo2) + aiz*dz(3,nubo2)) +
     &                               beta3*(uas3(2) - uas3(1))
         fltr4                     = beta2*
     &         (aix*dx(4,nubo2) + aiy*dy(4,nubo2) + aiz*dz(4,nubo2)) +
     &                               beta3*(uas4(2) - uas4(1))
         fltr5                     = beta2*
     &         (aix*dx(5,nubo2) + aiy*dy(5,nubo2) + aiz*dz(5,nubo2)) +
     &                               beta3*(uas5(2) - uas5(1))
c
         IF ((nordre .EQ. 3) .OR. (nordre .EQ. 4)) GOTO 275
c
c        Auxiliary Values for the Van Albada Procedure
c
         dpm                       =-(uas1(2) - uas1(1))
         dpex                      =-4.0d0*flur1 - dpm
         aux1                 = 0.25*(1.0d0 + SIGN(1.0d0, dpex*dpm))
         dpor                      =-4.0d0*fltr1 - dpm
         aux2                 = 0.25*(1.0d0 + SIGN(1.0d0, dpor*dpm))
c
         flur1                     = aux1*
     &                               ((dpex*dpex + e2)*dpm +
     &                                (dpm*dpm   + e2)*dpex)/
     &                               (dpex*dpex  + dpm*dpm + 2.0d0*e2)
         fltr1                     = aux2*
     &                               ((dpor*dpor + e2)*dpm +
     &                                (dpm*dpm   + e2)*dpor)/
     &                               (dpor*dpor  + dpm*dpm + 2.0d0*e2)
c
         dpm                       =-(uas2(2) - uas2(1))
         dpex                      =-4.0d0*flur2 - dpm
         aux1                 = 0.25*(1.0d0 + SIGN(1.0d0, dpex*dpm))
         dpor                      =-4.0d0*fltr2 - dpm
         aux2                 = 0.25*(1.0d0 + SIGN(1.0d0, dpor*dpm))
c
         flur2                     = aux1*
     &                               ((dpex*dpex + e2)*dpm +
     &                                (dpm*dpm   + e2)*dpex)/
     &                               (dpex*dpex  + dpm*dpm + 2.0d0*e2)
         fltr2                     = aux2*
     &                               ((dpor*dpor + e2)*dpm +
     &                                (dpm*dpm   + e2)*dpor)/
     &                               (dpor*dpor  + dpm*dpm + 2.0d0*e2)
c
         dpm                       =-(uas3(2) - uas3(1))
         dpex                      =-4.0d0*flur3 - dpm
         aux1                 = 0.25*(1.0d0 + SIGN(1.0d0, dpex*dpm))
         dpor                      =-4.0d0*fltr3 - dpm
         aux2                 = 0.25*(1.0d0 + SIGN(1.0d0, dpor*dpm))
c
         flur3                     = aux1*
     &                               ((dpex*dpex + e2)*dpm +
     &                                (dpm*dpm   + e2)*dpex)/
     &                               (dpex*dpex  + dpm*dpm + 2.0d0*e2)
         fltr3                     = aux2*
     &                               ((dpor*dpor + e2)*dpm +
     &                                (dpm*dpm   + e2)*dpor)/
     &                               (dpor*dpor  + dpm*dpm + 2.0d0*e2)
c
         dpm                       =-(uas4(2) - uas4(1))
         dpex                      =-4.0d0*flur4 - dpm
         aux1                 = 0.25*(1.0d0 + SIGN(1.0d0, dpex*dpm))
         dpor                      =-4.0d0*fltr4 - dpm
         aux2                 = 0.25*(1.0d0 + SIGN(1.0d0, dpor*dpm))
c      
         flur4                     = aux1*
     &                               ((dpex*dpex + e2)*dpm +
     &                                (dpm*dpm   + e2)*dpex)/
     &                               (dpex*dpex  + dpm*dpm + 2.0d0*e2)
         fltr4                     = aux2*
     &                               ((dpor*dpor + e2)*dpm +
     &                                (dpm*dpm   + e2)*dpor)/
     &                               (dpor*dpor  + dpm*dpm + 2.0d0*e2)
c
         dpm                       =-(uas5(2) - uas5(1))
         dpex                      =-4.0d0*flur5 - dpm
         aux1                 = 0.25*(1.0d0 + SIGN(1.0d0, dpex*dpm))
         dpor                      =-4.0d0*fltr5 - dpm
         aux2                 = 0.25*(1.0d0 + SIGN(1.0d0, dpor*dpm))
c
         flur5                     = aux1*
     &                               ((dpex*dpex + e2)*dpm +
     &                                (dpm*dpm   + e2)*dpex)/
     &                               (dpex*dpex  + dpm*dpm + 2.0d0*e2)
         fltr5                     = aux2*
     &                               ((dpor*dpor + e2)*dpm +
     &                                (dpm*dpm   + e2)*dpor)/
     &                               (dpor*dpor  + dpm*dpm + 2.0d0*e2)
c
275      CONTINUE
c
         uas1(1)                   = uas1(1) + bsign*flur1
         uas2(1)                   = uas2(1) + bsign*flur2
         uas3(1)                   = uas3(1) + bsign*flur3
         uas4(1)                   = uas4(1) + bsign*flur4
         uas5(1)                   = uas5(1) + bsign*flur5
c
         uas1(2)                   = uas1(2) - bsign*fltr1
         uas2(2)                   = uas2(2) - bsign*fltr2
         uas3(2)                   = uas3(2) - bsign*fltr3
         uas4(2)                   = uas4(2) - bsign*fltr4
         uas5(2)                   = uas5(2) - bsign*fltr5
c
         rnorm                     = 1.0d0/SQRT(
     &                               vnocl(1,iseg)*vnocl(1,iseg) +
     &                               vnocl(2,iseg)*vnocl(2,iseg) + 
     &                               vnocl(3,iseg)*vnocl(3,iseg))
c
         vno(1)                    =-vnocl(1,iseg)*rnorm
         vno(2)                    =-vnocl(2,iseg)*rnorm
         vno(3)                    =-vnocl(3,iseg)*rnorm
c
         xsigm                     =-sigma(iseg)*rnorm
c
         prod1                     = uas2(1)*vno(1) + uas3(1)*vno(2) + 
     &                               uas4(1)*vno(3)
c
         ener1                     = uas5(1)/gam1 + 
     &                               0.5d0*uas1(1)*(uas2(1)*uas2(1) + 
     &                                            uas3(1)*uas3(1) + 
     &                                            uas4(1)*uas4(1))
c 
         prod2                     = uas2(2)*vno(1) + uas3(2)*vno(2) + 
     &                               uas4(2)*vno(3)
c
         ener2                     = uas5(2)/gam1 + 
     &                               0.5d0*uas1(2)*(uas2(2)*uas2(2) +
     &                                            uas3(2)*uas3(2) + 
     &                                            uas4(2)*uas4(2))
c
         flum1                     = uas1(1)*(prod1 - xsigm) + 
     &                               uas1(2)*(prod2 - xsigm)
c
         flum2                     = uas1(1)*uas2(1)*(prod1 - xsigm) + 
     &                               uas1(2)*uas2(2)*(prod2 - xsigm) + 
     &                               (uas5(1) + uas5(2))*vno(1)
     &                               
c
         flum3                     = uas1(1)*uas3(1)*(prod1 - xsigm) +
     &                               uas1(2)*uas3(2)*(prod2 - xsigm) + 
     &                               (uas5(1) + uas5(2))*vno(2)
c
         flum4                     = uas1(1)*uas4(1)*(prod1 - xsigm) + 
     &                               uas1(2)*uas4(2)*(prod2 - xsigm) + 
     &                               (uas5(1) + uas5(2))*vno(3)
c
         flum5                     = ener1*(prod1 - xsigm) + 
     &                               ener2*(prod2 - xsigm) + 
     &                               uas5(1)*prod1 + uas5(2)*prod2 
c
         squsr1                    = SQRT(uas1(1))
         squsr2                    = SQRT(uas1(2))
c
         usro                      = 1.0d0/(squsr1 + squsr2)
c
         uar1                      = (squsr1*uas1(1) + 
     &                                squsr2*uas1(2))*usro
c
         uar2                      = (squsr1*uas2(1) + 
     &                                squsr2*uas2(2))*usro
c
         uar3                      = (squsr1*uas3(1) + 
     &                                squsr2*uas3(2))*usro
c
         uar4                      = (squsr1*uas4(1) + 
     &                                squsr2*uas4(2))*usro
c
         uar5                      = ((ener1 + uas5(1))/
     &                                squsr1 + 
     &                                (ener2 + uas5(2))/
     &                                squsr2)*usro
c
         pror                      = vno(1)*uar2 + vno(2)*uar3 +
     &                               vno(3)*uar4
c
         qir                       = 0.5d0*(uar2*uar2 + uar3*uar3 +
     &                                    uar4*uar4)
c
         tet1                      = vno(3)*uar3 - vno(2)*uar4
         tet2                      = vno(1)*uar4 - vno(3)*uar2
         tet3                      = vno(2)*uar2 - vno(1)*uar3
c
         cr2                       = gam1*(uar5 - qir)
         cr                        = SQRT(cr2)
         cr2                       = 1.0d0/cr2
c
         dif1                      = uas1(1) - uas1(2)
         dif2                      = uas1(1)*uas2(1) - uas1(2)*uas2(2)
         dif3                      = uas1(1)*uas3(1) - uas1(2)*uas3(2)
         dif4                      = uas1(1)*uas4(1) - uas1(2)*uas4(2)
         dif5                      = ener1 - ener2
c
         vp1                       = pror - xsigm
         vp4                       = pror + cr - xsigm
         vp5                       = pror - cr - xsigm
c

         IF (ient .EQ. 1) THEN
c
           delta                  = dabs(vp4)/100.d0
c
           vp1                    = dsqrt(vp1*vp1+delta)
           vp4                    = dsqrt(vp4*vp4+delta)
           vp5                    = dsqrt(vp5*vp5+delta)
c
c
         ENDIF
c
         flur1                     = ABS(vp1)*
     &              ((vno(1)*(1.0d0 - gam1*qir*cr2) - tet1)*dif1 + 
     &               (vno(1)*gam1*uar2*cr2)*dif2  +  
     &               (vno(3)  + (vno(1)*gam1*uar3*cr2))*dif3   +
     &               (-vno(2) + (vno(1)*gam1*uar4*cr2))*dif4   -
     &               (vno(1)*gam1*cr2)*dif5)
c
         flur2                     = ABS(vp1)*
     &              ((vno(2)*(1.0d0 - gam1*qir*cr2) - tet2)*dif1 + 
     &               (-vno(3) + (vno(2)*gam1*uar2*cr2))*dif2   + 
     &               (vno(2)*gam1*uar3*cr2)*dif3  + 
     &               (vno(1)  + (vno(2)*gam1*uar4*cr2))*dif4   -
     &               (vno(2)*gam1*cr2)*dif5)
c
         flur3                     = ABS(vp1)*
     &              ((vno(3)*(1.0d0 - gam1*qir*cr2) - tet3)*dif1 + 
     &               (vno(2)  + (vno(3)*gam1*uar2*cr2))*dif2   + 
     &               (-vno(1) + (vno(3)*gam1*uar3*cr2))*dif3   + 
     &               (vno(3)*gam1*uar4*cr2)*dif4  - 
     &               (vno(3)*gam1*cr2)*dif5)
c
         flur4                     = ABS(vp4)*
     &              ((-cr*pror   + gam1*qir)*dif1  + 
     &               ( cr*vno(1) - gam1*uar2)*dif2 + 
     &               ( cr*vno(2) - gam1*uar3)*dif3 + 
     &               ( cr*vno(3) - gam1*uar4)*dif4 + 
     &               gam1*dif5)
c
         flur5                     = ABS(vp5)*
     &              (( cr* pror  + gam1* qir)*dif1 + 
     &               (-cr*vno(1) - gam1*uar2)*dif2 + 
     &               (-cr*vno(2) - gam1*uar3)*dif3 +
     &               (-cr*vno(3) - gam1*uar4)*dif4 + 
     &               gam1*dif5)
c
         fltr1                     = vno(1)*flur1 + vno(2)*flur2 +
     &                               vno(3)*flur3 + 
     &                               0.5d0*(flur4 + flur5)*cr2 
c
         fltr2                     = 
     &              (uar2*vno(1))*flur1 +
     &              (uar2*vno(2) - vno(3))*flur2 + 
     &              (uar2*vno(3) + vno(2))*flur3 + 
     &              0.5d0*vno(1)*(flur4 - flur5)/cr + 
     &              0.5d0*uar2*(flur4 + flur5)*cr2 
c
         fltr3                     = 
     &              (uar3*vno(1) + vno(3))*flur1 +
     &              (uar3*vno(2))*flur2 + 
     &              (uar3*vno(3) - vno(1))*flur3 +
     &              0.5d0*vno(2)*(flur4 - flur5)/cr + 
     &              0.5d0*uar3*(flur4 + flur5)*cr2
c
         fltr4                     = 
     &              (uar4*vno(1) - vno(2))*flur1 +
     &              (uar4*vno(2) + vno(1))*flur2 +
     &              (uar4*vno(3))*flur3 + 
     &              0.5d0*vno(3)*(flur4 - flur5)/cr + 
     &              0.5d0*uar4*(flur4 + flur5)*cr2
c
         fltr5                     = 
     &              (qir*vno(1) + tet1)*flur1 +
     &              (qir*vno(2) + tet2)*flur2 +
     &              (qir*vno(3) + tet3)*flur3 +
     &              0.5d0*pror*(flur4 - flur5)/cr + 
     &              0.5d0*uar5*(flur4 + flur5)*cr2
c
         rnorm                     = 0.5d0/rnorm
c
         ce(1,nubo1)               = ce(1,nubo1) - (flum1 + fltr1)*rnorm
         ce(2,nubo1)               = ce(2,nubo1) - (flum2 + fltr2)*rnorm
         ce(3,nubo1)               = ce(3,nubo1) - (flum3 + fltr3)*rnorm  
         ce(4,nubo1)               = ce(4,nubo1) - (flum4 + fltr4)*rnorm
         ce(5,nubo1)               = ce(5,nubo1) - (flum5 + fltr5)*rnorm
c
         ce(1,nubo2)               = ce(1,nubo2) + (flum1 + fltr1)*rnorm
         ce(2,nubo2)               = ce(2,nubo2) + (flum2 + fltr2)*rnorm
         ce(3,nubo2)               = ce(3,nubo2) + (flum3 + fltr3)*rnorm
         ce(4,nubo2)               = ce(4,nubo2) + (flum4 + fltr4)*rnorm
         ce(5,nubo2)               = ce(5,nubo2) + (flum5 + fltr5)*rnorm
c
10    CONTINUE
c
      RETURN
      END
