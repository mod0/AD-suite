

      SUBROUTINE IMPROE
c---------------------------------------------------------------------   
c Computes the Jacobian matrix of the approximate Riemann solver 
c of Roe
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
      INTEGER ia  , icc ,j , iseg, id1
      INTEGER is  , is1   , is2 
      REAL*8    p, p1, p2, q1, q2
      REAL*8    delta       , switch
      REAL*8    ro  , usro
      REAL*8    uqro, usro1 , usro2
      REAL*8    u, v, w, q  , c   , enth
      REAL*8    vp1 , vp4   , vp5
      REAL*8    vno1  , vno2, vno3
      REAL*8    gc 
      REAL*8    rnorm       , usc2, unx
      REAL*8    xsigm       , xsig
      REAL*8    tmw(5,5)
      REAL*8    dmat(5,5)
      REAL*8    dmatr(5,5)
      REAL*8    dmatf(5,5)
c
c     Initializations
c
      DO 1000 iseg=1,nseg
c
         is1                       = nubo(1,iseg)
         is2                       = nubo(2,iseg)
c
         rnorm                     = SQRT(vnocl(1,iseg)*vnocl(1,iseg) +
     &                                    vnocl(2,iseg)*vnocl(2,iseg) + 
     &                                    vnocl(3,iseg)*vnocl(3,iseg))
c
         vno1                      =-vnocl(1,iseg)/rnorm
         vno2                      =-vnocl(2,iseg)/rnorm
         vno3                      =-vnocl(3,iseg)/rnorm
c
         xsigm                     =-sigma(iseg)/rnorm
c
         usro1                     = SQRT(ua(1,is1))
         usro2                     = SQRT(ua(1,is2))
         uqro                      = 1.0/(usro1 + usro2)
c
         ro                        = (usro1*ua(1,is1) + 
     &                                usro2*ua(1,is2))*uqro
         u                         = (usro1*ua(2,is1)/ua(1,is1) +
     &                                usro2*ua(2,is2)/ua(1,is2))*uqro
         v                         = (usro1*ua(3,is1)/ua(1,is1) +
     &                                usro2*ua(3,is2)/ua(1,is2))*uqro
         w                         = (usro1*ua(4,is1)/ua(1,is1) +
     &                                usro2*ua(4,is2)/ua(1,is2))*uqro
c        
         q1                        = (ua(2,is1)*ua(2,is1) + 
     &                                ua(3,is1)*ua(3,is1) + 
     &                                ua(4,is1)*ua(4,is1))/
     &                               (ua(1,is1)*ua(1,is1))
c        
         q2                        = (ua(2,is2)*ua(2,is2) + 
     &                                ua(3,is2)*ua(3,is2) + 
     &                                ua(4,is2)*ua(4,is2))/
     &                               (ua(1,is2)*ua(1,is2))
c        
         p1                        = gam1*(ua(5,is1) - 0.5*ua(1,is1)*q1)
         p2                        = gam1*(ua(5,is2) - 0.5*ua(1,is2)*q2)
c
         enth                      = ((ua(5,is1) + p1)/usro1 +
     &                                (ua(5,is2) + p2)/usro2)*uqro
c
         q                         = 0.5*(u*u + v*v + w*w)
         c                         = SQRT(gam1*(enth - q))
c
         unx                       = vno1*u + vno2*v + vno3*w
         gc                        = gam1/(c*c)
         usc2                      = 0.5/(c*c)
c
         vp1                       = ABS(rnorm*(unx - xsigm))
         vp4                       = ABS(rnorm*(unx + c - xsigm))
         vp5                       = ABS(rnorm*(unx - c - xsigm))
c
         IF (ient .EQ. 1) THEN
c
            delta                  = epsiim*vp4  
            switch                 = 0.5 + SIGN(0.5, ABS(vp1) - delta) 
     &                               
            vp1                    = switch*ABS(vp1) + 
     &                               0.5*(1.0 - switch)*
     &                               (delta + vp1*vp1)/delta
c
            switch                 = 0.5 + SIGN(0.5, ABS(vp4) - delta) 
            vp4                    = switch*ABS(vp4) + 
     &                               0.5*(1.0 - switch)*
     &                               (delta + vp4*vp4)/delta
c
            switch                 = 0.5 + SIGN(0.5, ABS(vp5) - delta)
            vp5                    = switch*ABS(vp5) + 
     &                               0.5*(1.0 - switch)*
     &                               (delta + vp5*vp5)/delta
c
         ENDIF
c
         tmw(1,1)                  = vno1 + vno2*w - v*vno3 - gc*vno1*q  
         tmw(1,2)                  = gc*vno1*u  
         tmw(1,3)                  = gc*vno1*v + vno3
         tmw(1,4)                  = gc*vno1*w - vno2
         tmw(1,5)                  =-gc*vno1
c
         tmw(2,1)                  = vno2 + vno3*u - w*vno1 - gc*vno2*q
         tmw(2,2)                  = gc*vno2*u - vno3 
         tmw(2,3)                  = gc*vno2*v 
         tmw(2,4)                  = gc*vno2*w + vno1 
         tmw(2,5)                  =-gc*vno2 
c
         tmw(3,1)                  = vno3 + vno1*v - u*vno2 - gc*vno3*q 
         tmw(3,2)                  = gc*vno3*u + vno2 
         tmw(3,3)                  = gc*vno3*v - vno1 
         tmw(3,4)                  = gc*vno3*w 
         tmw(3,5)                  =-gc*vno3 
c       
         tmw(4,1)                  = gam1*q - unx*c 
         tmw(4,2)                  =-gam1*u + c*vno1 
         tmw(4,3)                  =-gam1*v + c*vno2 
         tmw(4,4)                  =-gam1*w + c*vno3 
         tmw(4,5)                  = gam1
c
         tmw(5,1)                  = gam1*q + unx*c 
         tmw(5,2)                  =-gam1*u - c*vno1 
         tmw(5,3)                  =-gam1*v - c*vno2 
         tmw(5,4)                  =-gam1*w - c*vno3
         tmw(5,5)                  = gam1
c
         DO 70 j=1,5
c
            dmat(1,j)              = vp1*(vno1*tmw(1,j) +
     &                                    vno2*tmw(2,j) +
     &                                    vno3*tmw(3,j))
     &                             + vp4*usc2*tmw(4,j)
     &                             + vp5*usc2*tmw(5,j)
c
            dmat(2,j)              = vp1*(vno1*u*tmw(1,j) +
     &                                    (vno2*u - vno3)*tmw(2,j) +
     &                                    (vno3*u + vno2)*tmw(3,j))
     &                             + vp4*(usc2*(u + c*vno1))*tmw(4,j)
     &                             + vp5*(usc2*(u - c*vno1))*tmw(5,j)
c
            dmat(3,j)              = vp1*(vno2*v*tmw(2,j)   +
     &                                    (vno1*v + vno3)*tmw(1,j) +
     &                                    (vno3*v - vno1)*tmw(3,j))
     &                             + vp4*(usc2*(v + c*vno2))*tmw(4,j)
     &                             + vp5*(usc2*(v - c*vno2))*tmw(5,j)
c
            dmat(4,j)              = vp1*(vno3*w*tmw(3,j) +
     &                                    (vno1*w - vno2)*tmw(1,j) +
     &                                    (vno2*w + vno1)*tmw(2,j))
     &                             + vp4*(usc2*(w + c*vno3))*tmw(4,j)
     &                             + vp5*(usc2*(w - c*vno3))*tmw(5,j)
c
            dmat(5,j)              = vp1*((vno1*q + 
     &                                     (vno3*v - vno2*w))*tmw(1,j) +
     &                                    (vno2*q + 
     &                                     (vno1*w - vno3*u))*tmw(2,j) +
     &                                    (vno3*q +
     &                                     (vno2*u - vno1*v))*tmw(3,j))
     &                             + vp4*(usc2*(enth + unx*c))*tmw(4,j)
     &                             + vp5*(usc2*(enth - unx*c))*tmw(5,j)
c
70       CONTINUE
c
         DO 80 j=1,5
c
            dmatr(1,j)             = 0.5*dmat(1,j)
            dmatr(2,j)             = 0.5*dmat(2,j)
            dmatr(3,j)             = 0.5*dmat(3,j)
            dmatr(4,j)             = 0.5*dmat(4,j)
            dmatr(5,j)             = 0.5*dmat(5,j)
c
80       CONTINUE
c
         DO 90 icc=1,2
c
            is                     = nubo(icc,iseg)
c
            ro                     = ua(1,is)
            usro                   = 1.0/ro
c
            u                      = ua(2,is)*usro
            v                      = ua(3,is)*usro
            w                      = ua(4,is)*usro
            q                      = 0.5*(u*u + v*v + w*w)
c
            p                      = gam1*(ua(5,is) - ro*q)
c   
            c                      = SQRT(gam*p*usro)
c
            enth                   = (ua(5,is) + p)*usro
c
            unx                    = vno1*u + vno2*v + vno3*w
c
            gc                     = gam1/(c*c)
            xsig                   = enth
            usc2                   = 0.5/(c*c)
c
            tmw(1,1)               =-xsigm
            tmw(1,2)               = vno1
            tmw(1,3)               = vno2
            tmw(1,4)               = vno3
            tmw(1,5)               = 0.0
c
            tmw(2,1)               = gam1*vno1*q - u*unx
            tmw(2,2)               = (1.0 - gam1)*vno1*u + unx - xsigm
            tmw(2,3)               =-gam1*vno1*v + u*vno2
            tmw(2,4)               =-gam1*vno1*w + u*vno3
            tmw(2,5)               = gam1*vno1
c
            tmw(3,1)               = gam1*vno2*q - v*unx
            tmw(3,2)               =-gam1*vno2*u + v*vno1
            tmw(3,3)               = (1.0 - gam1)*vno2*v + unx - xsigm
            tmw(3,4)               =-gam1*vno2*w + v*vno3
            tmw(3,5)               = gam1*vno2
c
            tmw(4,1)               = gam1*vno3*q - w*unx
            tmw(4,2)               =-gam1*vno3*u + w*vno1
            tmw(4,3)               =-gam1*vno3*v + w*vno2
            tmw(4,4)               = (1.0 - gam1)*vno3*w + unx - xsigm
            tmw(4,5)               = gam1*vno3
c
            tmw(5,1)               = (gam1*q - xsig)*unx
            tmw(5,2)               = xsig*vno1 - gam1*u*unx
            tmw(5,3)               = xsig*vno2 - gam1*v*unx
            tmw(5,4)               = xsig*vno3 - gam1*w*unx
            tmw(5,5)               = (1.0 + gam1)*unx - xsigm
c
            IF (icc .EQ. 1) THEN
c
               DO 110 j=1,5
c
                  dmatf(1,j)       = 0.5*rnorm*tmw(1,j) + dmatr(1,j)
                  dmatf(2,j)       = 0.5*rnorm*tmw(2,j) + dmatr(2,j)
                  dmatf(3,j)       = 0.5*rnorm*tmw(3,j) + dmatr(3,j)
                  dmatf(4,j)       = 0.5*rnorm*tmw(4,j) + dmatr(4,j)
                  dmatf(5,j)       = 0.5*rnorm*tmw(5,j) + dmatr(5,j)
c
110            CONTINUE
c
c               IF (istok .EQ. 1) THEN
c
                  DO 120 j=1,5
c
                     id1           = 50*(iseg - 1) + j
c
                     stmat(id1)    = stmat(id1)    + dmatf(1,j)
                     stmat(id1+5)  = stmat(id1+5)  + dmatf(2,j)
                     stmat(id1+10) = stmat(id1+10) + dmatf(3,j)
                     stmat(id1+15) = stmat(id1+15) + dmatf(4,j)
                     stmat(id1+20) = stmat(id1+20) + dmatf(5,j)
c
120               CONTINUE
c
c               ENDIF
c
            ENDIF
c
            IF (icc .EQ. 2) THEN
c
               DO 130 j=1,5
c
                  dmatf(1,j)       =-0.5*rnorm*tmw(1,j) + dmatr(1,j)
                  dmatf(2,j)       =-0.5*rnorm*tmw(2,j) + dmatr(2,j)
                  dmatf(3,j)       =-0.5*rnorm*tmw(3,j) + dmatr(3,j)
                  dmatf(4,j)       =-0.5*rnorm*tmw(4,j) + dmatr(4,j)
                  dmatf(5,j)       =-0.5*rnorm*tmw(5,j) + dmatr(5,j)
c
130            CONTINUE
c
c               IF (istok .EQ. 1) THEN
c
                  DO 140 j=1,5
c
                     id1           = 50*(iseg - 1) + 25 + j
c
                     stmat(id1)    = stmat(id1)    + dmatf(1,j)
                     stmat(id1+5)  = stmat(id1+5)  + dmatf(2,j)
                     stmat(id1+10) = stmat(id1+10) + dmatf(3,j)
                     stmat(id1+15) = stmat(id1+15) + dmatf(4,j)
                     stmat(id1+20) = stmat(id1+20) + dmatf(5,j)
c
140               CONTINUE
c
c               ENDIF
c
            ENDIF
c
            DO 145 ia=1,5
c
               is                  = nubo(icc,iseg) 
c
               diag(is,ia,1)       = diag(is,ia,1) + dmatf(ia,1)
               diag(is,ia,2)       = diag(is,ia,2) + dmatf(ia,2)
               diag(is,ia,3)       = diag(is,ia,3) + dmatf(ia,3)
               diag(is,ia,4)       = diag(is,ia,4) + dmatf(ia,4)
               diag(is,ia,5)       = diag(is,ia,5) + dmatf(ia,5)
c
145         CONTINUE
c
90       CONTINUE
c
1000  CONTINUE
c
c
      RETURN
      END
