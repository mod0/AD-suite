      SUBROUTINE CONDBORDS(ctrlno)
C
C*** Traitement des bords de la matrice implicite
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
C     variables d'appel
C
       LOGICAL logique
C
C     Variables locales
C
      INTEGER ia, ib, j, i, is, ifac, if11, if1, if3,ilogf
      REAL*8 dvu(5),ctrlno
      REAL*8 rho, usrho, u, v, w, q, p, enth 
      REAL*8 vnx, vny, vnz, rnorm, vno1, vno2, vno3
      REAL*8 xsig, unx, cs2, gc, usc2
      REAL*8 vp1 ,vp4, vp5
      REAL*8 tmw(5,5), rmat(5,5)
      REAL*8 c
C
      IF (nf1 .GT. 0) THEN
c
         DO 500 if1=1,nf1
c
            ifac                   = noe1(nf11+if1)
c
c*** Inclusion des conditions de transpiration : ITRANS=1
c    On ne les inclut que sur la carlingue (cas du falcon) ou bien l'aile
c    ou bien la tuyere.
c    La carlingue, l'aile et la tuyere sont reperees si logfac = 200 ou 
c    si logfac = -2 ou si logfac = 7.
c    Test : Si transpiration, on traite le plan de symetrie comme avant
c                             et on applique la transpiration sur la coque.
c           Si non transpiration, on traite la coque et l'axe de symetrie
c                                 de la meme facon.
c    
            logique = .false.
c
            ilogf=logfac(ifac)
            if (ilogf.eq.-2 .and. ctrlno.lt.1.d-10) ilogf=2            
            
            if (itrans.eq.0) logique = .true.
            if (ITRANS.EQ.1) then
              if (ilogf.ne.-2 )then
                  if (ilogf.ne.200) then
                     if (ilogf.ne.7) then
                        logique = .true.
                     endif
                  endif
               endif
            endif
c
c             if (((ITRANS.EQ.1).AND.(ilogf.ne.-2).and.
c     $           (ilogf.ne.200).and.(ilogf.ne.7)
c     $             .OR. (ITRANS.EQ.0)) then
c
            if (logique) then
c
            dvu(2)                 = gam1*vnfac(1,ifac)
            dvu(3)                 = gam1*vnfac(2,ifac)
            dvu(4)                 = gam1*vnfac(3,ifac)
            dvu(5)                 = gam1*vnsig(ifac)
c     
            DO 510 j=1,3
c
               is                  = nsfac(j,ifac)
c
               rho                  = ua(1,is) 
               usrho                = 1.0/rho
               u                   = ua(2,is)*usrho
               v                   = ua(3,is)*usrho
               w                   = ua(4,is)*usrho
               q                   = 0.5*(u*u + v*v + w*w)
c     
               DO 525 i=2,5
c
                  diag(is,i,1)     = diag(is,i,1) + q*dvu(i)
                  diag(is,i,2)     = diag(is,i,2) - u*dvu(i)
                  diag(is,i,3)     = diag(is,i,3) - v*dvu(i)
                  diag(is,i,4)     = diag(is,i,4) - w*dvu(i)
                  diag(is,i,5)     = diag(is,i,5) + dvu(i)
c
525            CONTINUE
c
510         CONTINUE
c
         endif
c
500      CONTINUE
C
      ENDIF
c
      IF ((nf2 + nf3) .GT. 0) THEN
c
         DO 600 if3=1,nf2+nf3
c
            ifac                   = noe1(nf11+nf1+if3)
c
            vnx                    = vnfac(1,ifac)
            vny                    = vnfac(2,ifac)
            vnz                    = vnfac(3,ifac)
c
            rnorm                  = SQRT(vnx*vnx+ vny*vny + vnz*vnz)
c
            vno1                   = vnx/rnorm
            vno2                   = vny/rnorm
            vno3                   = vnz/rnorm
c
            xsig                   = vnsig(ifac)
c
            DO 610 i=1,3
c
               is                  = nsfac(i,ifac)
c
               rho                  = ua(1,is)
               usrho                = 1.0/rho
c
               u                   = ua(2,is)*usrho
               v                   = ua(3,is)*usrho
               w                   = ua(4,is)*usrho
c
               q                   = 0.5*(u*u + v*v + w*w)
c
               p                   = gam1*(ua(5,is) - rho*q)
c
               c                   = SQRT(gam*p*usrho)
c            
               enth                = (ua(5,is) + p)*usrho
c
               unx                 = vno1*u + vno2*v + vno3*w
c
               cs2                 = 1.0/(c*c)
               gc                  = gam1*cs2 
               usc2                = 0.5*cs2
c 
               vp1                 = MAX(rnorm*unx - xsig, 0.0d0)
               vp4                 = MAX(rnorm*(unx + c) - xsig, 0.0d0)
               vp5                 = MAX(rnorm*(unx - c) - xsig, 0.0d0)
c
               tmw(1,1)            = vno1 + usrho*(vno2*w - v*vno3) - 
     &                               gc*vno1*q
               tmw(1,2)            = gc*vno1*u
               tmw(1,3)            = gc*vno1*v + usrho*vno3
               tmw(1,4)            = gc*vno1*w - usrho*vno2
               tmw(1,5)            =-gc*vno1
c
               tmw(2,1)            = vno2 + usrho*(vno3*u - w*vno1) -
     &                               gc*vno2*q
               tmw(2,2)            = gc*vno2*u - usrho*vno3
               tmw(2,3)            = gc*vno2*v
               tmw(2,4)            = gc*vno2*w + usrho*vno1
               tmw(2,5)            =-gc*vno2
c   
               tmw(3,1)            = vno3 + usrho*(vno1*v  - u*vno2) -
     &                               gc*vno3*q
               tmw(3,2)            = gc*vno3*u + usrho*vno2
               tmw(3,3)            = gc*vno3*v - usrho*vno1
               tmw(3,4)            = gc*vno3*w
               tmw(3,5)            =-gc*vno3
c   
               tmw(4,1)            = gam1*q - unx*c
               tmw(4,2)            =-gam1*u + c*vno1
               tmw(4,3)            =-gam1*v + c*vno2
               tmw(4,4)            =-gam1*w + c*vno3
               tmw(4,5)            = gam1
c
               tmw(5,1)            = gam1*q + unx*c
               tmw(5,2)            =-gam1*u - c*vno1
               tmw(5,3)            =-gam1*v - c*vno2
               tmw(5,4)            =-gam1*w - c*vno3
               tmw(5,5)            = gam1
c
               DO 700 j=1,5
c
                  rmat(1,j)        = vp1*(tmw(1,j)*vno1  + 
     &                                    tmw(2,j)*vno2  +
     &                                    tmw(3,j)*vno3) +
     &                               vp4*tmw(4,j)*usc2 +
     &                               vp5*tmw(5,j)*usc2
c
                  rmat(2,j)        = vp1*(tmw(1,j)*vno1*u  +
     &                                    tmw(2,j)*(vno2*u - rho*vno3) + 
     &                                    tmw(3,j)*(vno3*u + rho*vno2))  
     &                             + vp4*tmw(4,j)*(usc2*(u + c*vno1))  
     &                             + vp5*tmw(5,j)*(usc2*(u - c*vno1))
c
                  rmat(3,j)       = vp1*(tmw(1,j)*(vno1*v  + rho*vno3) + 
     &                                    tmw(2,j)*(vno2*v) + 
     &                                    tmw(3,j)*(vno3*v - rho*vno1))
     &                             + vp4*tmw(4,j)*(usc2*(v + c*vno2))
     &                             + vp5*tmw(5,j)*(usc2*(v - c*vno2))
c
                  rmat(4,j)        = vp1*(tmw(1,j)*(vno1*w - rho*vno2) +
     &                                    tmw(2,j)*(vno2*w + rho*vno1) +
     &                                    tmw(3,j)*(vno3*w))
     &                             + vp4*tmw(4,j)*(usc2*(w + c*vno3))
     &                             + vp5*tmw(5,j)*(usc2*(w - c*vno3))
c
                  rmat(5,j)        = vp1*(tmw(1,j)*(vno1*q  + 
     &                                      rho*(vno3*v - vno2*w)) +
     &                                    tmw(2,j)*(vno2*q  + 
     &                                       rho*(vno1*w - vno3*u)) + 
     &                                    tmw(3,j)*(vno3*q  +
     &                                           rho*(vno2*u - vno1*v)))
     &                             + vp4*tmw(4,j)*(usc2*(enth + unx*c))
     &                             + vp5*tmw(5,j)*(usc2*(enth - unx*c))
c
700            CONTINUE
c
               DO 720 ia=1,5
c
                  Do 721 ib=1,5
c
                  diag(is,ia,ib)    = diag(is,ia,ib) + rmat(ia,ib)
c
721               Continue
c
720            CONTINUE
c
610         CONTINUE
c
600      CONTINUE  
c
      ENDIF
c
      IF (ivis .EQ. 1) THEN 
c
         DO 785 if11=1,nf11
c
            ifac                   = noe1(if11)
c
            DO 790 ia=1,3
c
               is                  = nsfac(ia,ifac)
c
               diag(is,2,1)        = 0.0
               diag(is,2,2)        = 1.0
               diag(is,2,3)        = 0.0
               diag(is,2,4)        = 0.0
               diag(is,2,5)        = 0.0
c
               diag(is,3,1)        = 0.0
               diag(is,3,2)        = 0.0
               diag(is,3,3)        = 1.0
               diag(is,3,4)        = 0.0
               diag(is,3,5)        = 0.0
c
               diag(is,4,1)        = 0.0
               diag(is,4,2)        = 0.0
               diag(is,4,3)        = 0.0
               diag(is,4,4)        = 1.0
               diag(is,4,5)        = 0.0
c
               diag(is,5,1)        =-tbrd
               diag(is,5,2)        = 0.0
               diag(is,5,3)        = 0.0
               diag(is,5,4)        = 0.0
               diag(is,5,5)        = 1.0
c
790         CONTINUE  
c
785      CONTINUE  
c
      ENDIF
c
c
      RETURN
      END
