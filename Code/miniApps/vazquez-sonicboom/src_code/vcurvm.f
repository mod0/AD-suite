      SUBROUTINE VCURVM(ctrlno)
c---------------------------------------------------------------------   
c        Traitement des conditions aux bords
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c---------------------------------------------------------------------   
c     Local variables definition
c
      INTEGER if1 , if2  , if3  , is1 , is2 , is3, ifac, is , j,ilogf
      REAL*8    rhom, rhum , rhvm , rhwm, rhem, pm , vnx , vny, vnz
      REAL*8    sig , pm1  , pm2  , pm3
      REAL*8    cap , capd1, capd2, capd3
      REAL*8    ro  , usro , u , v, w, p, c, e, usc, unx
      REAL*8    vp1 , vp4  , vp5  , vpp , vpm
      REAL*8    x1  , x2   , x3   , x4  , x5  , vit, c2  , cc
      REAL*8    xn1 , xn2  , ff, gg
      REAL*8    fgp(5) , fgm(5),ctrlno
C
c*** Traitement des bords glissants
c
      
      IF (nf1 .GT. 0) THEN
c         
C$AD II-LOOP
         DO 1000 if1=1,nf1
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
c*** Non transpiration :
c      Psi_bords = (0,P.nx,P.ny,P.nz,0)
c

            ilogf=logfac(ifac)
            if (ilogf.eq.-2 .and. ctrlno.lt.1.d-10) ilogf=2
            
            if (((ITRANS.EQ.1).AND.(ilogf.ne.-2).and.
     $        (ilogf.ne.200).and.(ilogf.ne.7))
     $        .OR. (ITRANS.EQ.0)) then
c
              is1                    = nsfac(1,ifac)
            is2                    = nsfac(2,ifac)
            is3                    = nsfac(3,ifac)
c
            rhom                   = ua(1,is1) 
            rhum                   = ua(2,is1)
            rhvm                   = ua(3,is1) 
            rhwm                   = ua(4,is1)
            rhem                   = ua(5,is1) 
c
            pm1                    = gam1*(rhem - 
     &        0.5d0*(rhum*rhum + rhvm*rhvm +
     &        rhwm*rhwm)/rhom)
c
            rhom                   = ua(1,is2) 
            rhum                   = ua(2,is2) 
            rhvm                   = ua(3,is2) 
            rhwm                   = ua(4,is2) 
            rhem                   = ua(5,is2) 
c
            pm2                    = gam1*(rhem - 
     &        0.5d0*(rhum*rhum + rhvm*rhvm +
     &        rhwm*rhwm)/rhom)
c
            rhom                   = ua(1,is3)
            rhum                   = ua(2,is3)
            rhvm                   = ua(3,is3)
            rhwm                   = ua(4,is3)
            rhem                   = ua(5,is3)
c
            pm3                    = gam1*(rhem - 
     &        0.5d0*(rhum*rhum + rhvm*rhvm +
     &        rhwm*rhwm)/rhom)
c
            CE(2,is1)              = CE(2,is1) - vnfac(1,ifac)*pm1
            CE(3,is1)              = CE(3,is1) - vnfac(2,ifac)*pm1
            CE(4,is1)              = CE(4,is1) - vnfac(3,ifac)*pm1
            CE(5,is1)              = CE(5,is1) - vnsig(ifac)*pm1
c
            CE(2,is2)              = CE(2,is2) - vnfac(1,ifac)*pm2
            CE(3,is2)              = CE(3,is2) - vnfac(2,ifac)*pm2
            CE(4,is2)              = CE(4,is2) - vnfac(3,ifac)*pm2
            CE(5,is2)              = CE(5,is2) - vnsig(ifac)*pm2
c
            CE(2,is3)              = CE(2,is3) - vnfac(1,ifac)*pm3
            CE(3,is3)              = CE(3,is3) - vnfac(2,ifac)*pm3
            CE(4,is3)              = CE(4,is3) - vnfac(3,ifac)*pm3
            CE(5,is3)              = CE(5,is3) - vnsig(ifac)*pm3
c
            endif
c
1000     CONTINUE
c
      ENDIF
c
c*** Traitement des conditions a l'infini (Steger-Warming)
c
      IF (nf3 .GT. 0) THEN
c
C$AD II-LOOP
         DO 2000 if3=1,nf3
c
            ifac                   = noe1(nf11+nf1+nf2+if3)
c     
            vnx                    = vnfac(1,ifac)
            vny                    = vnfac(2,ifac)
            vnz                    = vnfac(3,ifac)
c  
c            sig                    = vnsig(ifac)
c
            cap                    = SQRT(vnx*vnx+ vny*vny + vnz*vnz)
            capd1                  = vnx/cap
            capd2                  = vny/cap
            capd3                  = vnz/cap
c              
            DO 2100 j=1,3    
c
c   Partie +
c
               is                  = nsfac(j,ifac)
c
               ro                  = ua(1,is)
               usro                = 1.0d0/ro
c
               u                   = ua(2,is)*usro
               v                   = ua(3,is)*usro
               w                   = ua(4,is)*usro
c
               vit                 = 0.5d0*(u*u + v*v + w*w)
               pm                  = gam1*(ua(5,is) - ro*vit)
c      
               c2                  = gam*pm*usro
               c                   = SQRT(c2)
               usc                 = 1.0d0/c
               cc                  = vit + c2/gam1               
               unx                 = capd1*u + capd2*v + capd3*w

c
c
               x1                  = ro
               x2                  = ro*u
               x3                  = ro*v
               x4                  = ro*w
               x5                  = pm/gam1 + ro*vit  
c                   
               vp1                 = MAX(cap*unx, 0.0d0)
               vp4                 = MAX(cap*(unx + c), 0.0d0)
               vp5                 = MAX(cap*(unx - c), 0.0d0)
c
               vpp                 = vp1 - 0.5d0*(vp4 + vp5)
               vpm                 = 0.5d0*(vp5 - vp4)
c
               xn1                 = gam1*usc*(-vit*x1 + 
     &                                         u*x2 + v*x3 + w*x4 - x5)
               xn2                 = unx*x1 - 
     &                               capd1*x2 - capd2*x3 - capd3*x4
c
               ff                  = usc*(vpp*xn1 + vpm*xn2)
               gg                  = vpm*xn1 + vpp*xn2
c
               fgp(1)              = vp1*x1 + ff
               fgp(2)              = vp1*x2 + u*ff  + capd1*gg
               fgp(3)              = vp1*x3 + v*ff  + capd2*gg
               fgp(4)              = vp1*x4 + w*ff  + capd3*gg
               fgp(5)              = vp1*x5 + cc*ff + unx*gg

c
c   Partie -
c
               is                  = nsfac(j,ifac)
c
               ro                  = ub3(1,if3)               
               usro                = 1.0d0/ro
               u                   = ub3(2,if3)
               v                   = ub3(3,if3)
               w                   = ub3(4,if3)
               vit                 = 0.5d0*(u*u + v*v + w*w)               
               pm                   = ub3(5,if3)
c
               c2                  = gam*pm*usro
               c                   = SQRT(c2)
               usc                 = 1.0d0/c
               cc                  = vit + c2/gam1               
               unx                 = capd1*u + capd2*v + capd3*w
c
               x1                  = ro
               x2                  = ro*u
               x3                  = ro*v
               x4                  = ro*w
               x5                  = pm/gam1 + ro*vit  
c                   
               vp1                 = MIN(cap*unx, 0.0d0)
               vp4                 = MIN(cap*(unx + c) , 0.0d0)
               vp5                 = MIN(cap*(unx - c) , 0.0d0)
               
c
               vpp                 = vp1 - 0.5d0*(vp4 + vp5)
               vpm                 = 0.5d0*(vp5 - vp4)
c
               xn1                 = gam1*usc*(-vit*x1 + 
     &                                         u*x2 + v*x3 + w*x4 - x5)
               xn2                 = unx*x1 - 
     &                               capd1*x2 - capd2*x3 - capd3*x4
c
               ff                  = usc*(vpp*xn1 + vpm*xn2)
               gg                  = vpm*xn1 + vpp*xn2
c
c
c
               fgm(1)              = vp1*x1 + ff
               fgm(2)              = vp1*x2 + u*ff  + capd1*gg
               fgm(3)              = vp1*x3 + v*ff  + capd2*gg
               fgm(4)              = vp1*x4 + w*ff  + capd3*gg
               fgm(5)              = vp1*x5 + cc*ff + unx*gg
c
               CE(1,is)            = CE(1,is) - fgp(1) - fgm(1)
               CE(2,is)            = CE(2,is) - fgp(2) - fgm(2)
               CE(3,is)            = CE(3,is) - fgp(3) - fgm(3)
               CE(4,is)            = CE(4,is) - fgp(4) - fgm(4)
               CE(5,is)            = CE(5,is) - fgp(5) - fgm(5)

               
2100        CONTINUE  
c
2000     CONTINUE
c
      ENDIF
c

      
c      stop
      
      END
