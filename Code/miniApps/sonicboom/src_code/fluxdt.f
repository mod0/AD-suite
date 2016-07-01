      SUBROUTINE FLUXDT
c---------------------------------------------------------------------   
c Computes the local time steps
c Computes and gathers the viscous fluxes
c The hermitian nodal gradients are computed in GRADNOD, not here.
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------
c     Local variables definition
      INTEGER ivar  , k     , jt   , ia  , ib
      INTEGER is    , is1   , is2  , is3 , is4  , id1 , id2
      INTEGER nex(6), nor(6), iseg , nub1, nub2 , js1 , js2
      REAL*8    usro  , usgam1
      REAL*8    gam53 , usrey , gampr
      REAL*8    vol6  , ait   , ais  , clhaut
      REAL*8    romi  , roma  , ema  , ecmi, vma
      Real*8    qxmi  , qymi  , qzmi , qxma, qyma , qzma
      REAL*8    xivis , ds3   , us4  , us6
      REAL*8    x(4)  , y(4)  , z(4) , b(4), c(4) , d(4)
      REAL*8    dbx(5), dby(5), dbz(5)
      REAL*8    dxt(5), dyt(5), dzt(5)
      REAL*8    pma   , dtvis , dtmax(4)
      REAL*8    uph(5,4)      , ei(4), um(3)
      REAL*8    eix  , eiy    , eiz
      REAL*8    txx  , txy    , txz   
      REAL*8    tyz  , tzz    , tyy
      REAL*8    r5   , s5     , q5
      REAL*8    a11  , a22    , a33  , a44
      REAL*8    dbxr(4)       , dbyr(4)    , dbzr(4)
      REAL*8    rmat(5,5,4)   , smat(5,5,4), qmat(5,5,4)
      REAL*8    amat(5,5,6)   , bmat(5,5,6)
c
c     Initialisations
c
      nex(1)                       = 1
      nex(2)                       = 1
      nex(3)                       = 1
      nex(4)                       = 2
      nex(5)                       = 2
      nex(6)                       = 3
c
      nor(1)                       = 2
      nor(2)                       = 3
      nor(3)                       = 4
      nor(4)                       = 3
      nor(5)                       = 4
      nor(6)                       = 4
c
      IF (ivis .EQ. 0) THEN 
         xivis                     = 0.0
      ELSE
         xivis                     = 1.0
      ENDIF
c
      dt                           = 1.0e+20
c
      ds3                          = 2.0/3.0
      us4                          = 1.0/4.0
      us6                          = 1.0/6.0
c
      usgam1                       = 1.0/gam1
      gam53                        = (gam - 5.0/3.0)/gam1
      usrey                        = 1.0/rey
      gampr                        = gam/pr
c
      DO 10 is=1,ns
         dtl(is)                = 1.0e+20
10    CONTINUE
c
      DO 1000 jt=1,nt
c
         DO 30 k=1,4
c
            is                     = nu(k,jt)
c
            uph(1,k)               = ua(1,is)
            usro                   = 1.0/uph(1,k)
c
            uph(2,k)               = ua(2,is)*usro
            uph(3,k)               = ua(3,is)*usro
            uph(4,k)               = ua(4,is)*usro
c
            uph(5,k)               = gam1*(ua(5,is) -
     &                                     0.5*uph(1,k)*(
     &                                     uph(2,k)*uph(2,k) +
     &                                     uph(3,k)*uph(3,k) +
     &                                     uph(4,k)*uph(4,k)))
c
c           Specific internal energy
c
            ei(k)                  = usgam1*usro*uph(5,k)
c
30       CONTINUE
c
c        Computing the local time steps
c
         is1                       = nu(1,jt)
         is2                       = nu(2,jt)
         is3                       = nu(3,jt)
         is4                       = nu(4,jt)
c
         roma                      = MAX(ABS(ua(1,is1)), ABS(ua(1,is2)))
         roma                      = MAX(roma, ABS(ua(1,is3)))
         roma                      = MAX(roma, ABS(ua(1,is4)))
c
         qxmi                      = MIN(ABS(ua(2,is1)), ABS(ua(2,is2)))
         qxmi                      = MIN(qxmi, ABS(ua(2,is3)))
         qxmi                      = MIN(qxmi, ABS(ua(2,is4)))
c
         qymi                      = MIN(ABS(ua(3,is1)), ABS(ua(3,is2)))
         qymi                      = MIN(qymi, ABS(ua(3,is3)))
         qymi                      = MIN(qymi, ABS(ua(3,is4)))
c
         qzmi                      = MIN(ABS(ua(4,is1)), ABS(ua(4,is2)))
         qzmi                      = MIN(qzmi, ABS(ua(4,is3)))
         qzmi                      = MIN(qzmi, ABS(ua(4,is4)))
c
         ema                       = MAX(ABS(ua(5,is1)), ABS(ua(5,is2)))
         ema                       = MAX(ema, ABS(ua(5,is3)))
         ema                       = MAX(ema, ABS(ua(5,is4)))
c
         ecmi                      = 0.5*(qxmi*qxmi + qymi*qymi + 
     &                                    qzmi*qzmi)/roma
c
         pma                       = gam1*ABS(ema - ecmi)
c
         romi                      = MIN(ABS(ua(1,is1)), ABS(ua(1,is2)))
         romi                      = MIN(romi, ABS(ua(1,is3)))
         romi                      = MIN(romi, ABS(ua(1,is4)))
c
         qxma                      = MAX(ABS(ua(2,is1)), ABS(ua(2,is2)))
         qxma                      = MAX(qxma, ABS(ua(2,is3)))
         qxma                      = MAX(qxma, ABS(ua(2,is4)))
c
         qyma                      = MAX(ABS(ua(3,is1)), ABS(ua(3,is2)))
         qyma                      = MAX(qyma, ABS(ua(3,is3)))
         qyma                      = MAX(qyma, ABS(ua(3,is4)))
c
         qzma                      = MAX(ABS(ua(4,is1)), ABS(ua(4,is2)))
         qzma                      = MAX(qzma, ABS(ua(4,is3)))
         qzma                      = MAX(qzma, ABS(ua(4,is4)))
c
         vma                       = SQRT(qxma*qxma + qyma*qyma + 
     &                                    qzma*qzma)/romi
c
c        Minimal wave velocity for the current tetraedra
c
         pma                       = 1.0/(vma + SQRT(gam*pma/romi))
c
c        Viscous contribution to the time step
c
         dtvis                     = 2.0*xivis*usrey*
     &                               MAX(gampr, 1.0d0)/romi
c
c        Computing the basis function gradients
c
         x(1)                      = coor(1,is1)
         y(1)                      = coor(2,is1)
         z(1)                      = coor(3,is1)
c
         x(2)                      = coor(1,is2)
         y(2)                      = coor(2,is2)
         z(2)                      = coor(3,is2)
c
         x(3)                      = coor(1,is3)
         y(3)                      = coor(2,is3)
         z(3)                      = coor(3,is3)
c
         x(4)                      = coor(1,is4)
         y(4)                      = coor(2,is4)
         z(4)                      = coor(3,is4)
c     
         CALL GRADFB(x, y, z, b, c, d, vol6)
c
         ait                       = us6/volt(jt)
c
         dbx(1)                    = b(1)*ait
         dbx(2)                    = b(2)*ait
         dbx(3)                    = b(3)*ait
         dbx(4)                    = b(4)*ait
c
         dby(1)                    = c(1)*ait
         dby(2)                    = c(2)*ait
         dby(3)                    = c(3)*ait
         dby(4)                    = c(4)*ait
c
         dbz(1)                    = d(1)*ait
         dbz(2)                    = d(2)*ait
         dbz(3)                    = d(3)*ait
         dbz(4)                    = d(4)*ait
c
c        Computing the minimal time steps for each 
c        vertex of a tetraedra
c
         clhaut                    = dbx(1)*dbx(1) + dby(1)*dby(1) +
     &                               dbz(1)*dbz(1)
c
         dtmax(1)                  = pma/
     &                               (SQRT(clhaut) + pma*dtvis*clhaut)
     &                                
c
         clhaut                    = dbx(2)*dbx(2) + dby(2)*dby(2) +
     &                               dbz(2)*dbz(2)
c
         dtmax(2)                  = pma/
     &                               (SQRT(clhaut) + pma*dtvis*clhaut)
c
         clhaut                    = dbx(3)*dbx(3) + dby(3)*dby(3) +
     &                               dbz(3)*dbz(3)
c
         dtmax(3)                  = pma/
     &                               (SQRT(clhaut) + pma*dtvis*clhaut)
c
         clhaut                    = dbx(4)*dbx(4) + dby(4)*dby(4) +
     &                               dbz(4)*dbz(4)
c
         dtmax(4)                  = pma/
     &                               (SQRT(clhaut) + pma*dtvis*clhaut)
c
c        Computing the local time step as the minimum of 
c        the contributions from each of the attached tetraedras
c
         dtl(is1)                  = MIN(dtl(is1), dtmax(1))
         dtl(is2)                  = MIN(dtl(is2), dtmax(2))
         dtl(is3)                  = MIN(dtl(is3), dtmax(3))
         dtl(is4)                  = MIN(dtl(is4), dtmax(4))
cc
cc        Computing the P1-gradients on each tetraedra
cc
cc         DO 60 ivar=1,5
cc
cc            dxt(ivar)              = uph(ivar,1)*dbx(1) +
cc     &                               uph(ivar,2)*dbx(2) +
cc     &                               uph(ivar,3)*dbx(3) +
cc     &                               uph(ivar,4)*dbx(4)
cc
cc            dyt(ivar)              = uph(ivar,1)*dby(1) +
cc     &                               uph(ivar,2)*dby(2) +
cc     &                               uph(ivar,3)*dby(3) +
cc     &                               uph(ivar,4)*dby(4)
cc
cc            dzt(ivar)              = uph(ivar,1)*dbz(1) +
cc     &                               uph(ivar,2)*dbz(2) +
cc     &                               uph(ivar,3)*dbz(3) +
cc     &                               uph(ivar,4)*dbz(4)
cc
cc60       CONTINUE
cc
         IF (ivis .EQ. 1) THEN 
c
            um(1)                  = us4*(uph(2,1) + uph(2,2) + 
     &                                    uph(2,3) + uph(2,4))
            um(2)                  = us4*(uph(3,1) + uph(3,2) + 
     &                                    uph(3,3) + uph(3,4))
            um(3)                  = us4*(uph(4,1) + uph(4,2) + 
     &                                    uph(4,3) + uph(4,4))
c
            eix                    = dbx(1)*ei(1) + dbx(2)*ei(2) +
     &                               dbx(3)*ei(3) + dbx(4)*ei(4)
c
            eiy                    = dby(1)*ei(1) + dby(2)*ei(2) +
     &                               dby(3)*ei(3) + dby(4)*ei(4)
c
            eiz                    = dbz(1)*ei(1) + dbz(2)*ei(2) +
     &                               dbz(3)*ei(3) + dbz(4)*ei(4)
c
            txx                    = ds3*(2.0*dxt(2) - dyt(3) - dzt(4))
            tyy                    = ds3*(2.0*dyt(3) - dxt(2) - dzt(4))
            tzz                    = ds3*(2.0*dzt(4) - dxt(2) - dyt(3))
c       
            txy                    = dyt(2) + dxt(3)
            txz                    = dxt(4) + dzt(2)
            tyz                    = dyt(4) + dzt(3)
c
            r5                     = um(1)*txx + um(2)*txy +
     &                               um(3)*txz + gampr*eix
c
            s5                     = um(1)*txy + um(2)*tyy +
     &                               um(3)*tyz + gampr*eiy
c
            q5                     = um(1)*txz + um(2)*tyz +
     &                               um(3)*tzz + gampr*eiz
c
            ait                    = usrey*volt(jt)
c
            DO 240 k=1,4
c
               is                  = nu(k,jt)
c
               ce(2,is)            = ce(2,is) - ait*(dbx(k)*txx +
     &                                               dby(k)*txy +
     &                                               dbz(k)*txz)
               ce(3,is)            = ce(3,is) - ait*(dbx(k)*txy +
     &                                               dby(k)*tyy +
     &                                               dbz(k)*tyz)
               ce(4,is)            = ce(4,is) - ait*(dbx(k)*txz +
     &                                               dby(k)*tyz +
     &                                               dbz(k)*tzz)
               ce(5,is)            = ce(5,is) - ait*(dbx(k)*r5  +
     &                                               dby(k)*s5  +
     &                                               dbz(k)*q5)
c
240         CONTINUE
c
c           Contribution of the viscous fluxes to the implicit matrix
c
            IF (nexp .EQ. 0) THEN 
c
               dbxr(1)             = dbx(1)/uph(1,1)
               dbxr(2)             = dbx(2)/uph(1,2)
               dbxr(3)             = dbx(3)/uph(1,3)
               dbxr(4)             = dbx(4)/uph(1,4)
c
               dbyr(1)             = dby(1)/uph(1,1)
               dbyr(2)             = dby(2)/uph(1,2)
               dbyr(3)             = dby(3)/uph(1,3)
               dbyr(4)             = dby(4)/uph(1,4)
c
               dbzr(1)             = dbz(1)/uph(1,1)
               dbzr(2)             = dbz(2)/uph(1,2)
               dbzr(3)             = dbz(3)/uph(1,3)
               dbzr(4)             = dbz(4)/uph(1,4)
c
               DO 400 k=1,4
c
                  DO 410 ivar=1,5
                     rmat(1,ivar,k)= 0.0
                     smat(1,ivar,k)= 0.0
                     qmat(1,ivar,k)= 0.0
410               CONTINUE
c
                  rmat(2,1,k)      =-ds3*(
     &                               2.0*uph(2,k)*dbxr(k) -
     &                               uph(3,k)*dbyr(k) - 
     &                               uph(4,k)*dbzr(k))
                  rmat(2,2,k)      = 2.0*ds3*dbxr(k)
                  rmat(2,3,k)      =-ds3*dbyr(k)
                  rmat(2,4,k)      =-ds3*dbzr(k)
                  rmat(2,5,k)      = 0.0
c
                  rmat(3,1,k)      =-(   
     &                               uph(3,k)*dbxr(k) +
     &                               uph(2,k)*dbyr(k))
                  rmat(3,2,k)      = dbyr(k)
                  rmat(3,3,k)      = dbxr(k)
                  rmat(3,4,k)      = 0.0
                  rmat(3,5,k)      = 0.0
c
                  rmat(4,1,k)      =-(
     &                               uph(4,k)*dbxr(k) +
     &                               uph(2,k)*dbzr(k))
                  rmat(4,2,k)      = dbzr(k)
                  rmat(4,3,k)      = 0.0
                  rmat(4,4,k)      = dbxr(k)
                  rmat(4,5,k)      = 0.0
c
c
                  rmat(5,1,k)      =-us4*(
     &                               uph(2,k)*txx + 
     &                               uph(3,k)*txy +
     &                               uph(4,k)*txz)/
     &                               uph(1,k) +
     &                               um(1)*rmat(2,1,k) + 
     &                               um(2)*rmat(3,1,k) +
     &                               um(3)*rmat(4,1,k) - 
     &                               gampr*dbxr(k)*ei(k)
c
                  rmat(5,2,k)      = us4*txx/uph(1,k) +
     &                               um(1)*rmat(2,2,k) +
     &                               um(2)*rmat(3,2,k) + 
     &                               um(3)*rmat(4,2,k) -
     &                               gampr*dbxr(k)*uph(2,k)
c
                  rmat(5,3,k)      = us4*txy/uph(1,k) +
     &                               um(1)*rmat(2,3,k) +
     &                               um(2)*rmat(3,3,k) +
     &                               um(3)*rmat(4,3,k) - 
     &                               gampr*dbxr(k)*uph(3,k)
c
                  rmat(5,4,k)      = us4*txz/uph(1,k) +
     &                               um(1)*rmat(2,4,k) +
     &                               um(2)*rmat(3,4,k) + 
     &                               um(3)*rmat(4,4,k) -
     &                               gampr*dbxr(k)*uph(4,k)
c
                  rmat(5,5,k)      = gampr*dbxr(k)
c
                  smat(3,1,k)      =-ds3*(
     &                               2.0*uph(3,k)*dbyr(k) -
     &                               uph(2,k)*dbxr(k) -
     &                               uph(4,k)*dbzr(k))
c
                  smat(3,2,k)      =-ds3*dbxr(k)
                  smat(3,3,k)      = 2.0*ds3*dbyr(k)
                  smat(3,4,k)      =-ds3*dbzr(k)
                  smat(3,5,k)      = 0.0
c
                  smat(4,1,k)      =-(
     &                               uph(4,k)*dbyr(k) +
     &                               uph(3,k)*dbzr(k))
                  smat(4,2,k)      = 0.0
                  smat(4,3,k)      = dbzr(k)
                  smat(4,4,k)      = dbyr(k)
                  smat(4,5,k)      = 0.0
c
                  smat(5,1,k)      =-us4*(
     &                               uph(2,k)*txy +
     &                               uph(3,k)*tyy +
     &                               uph(4,k)*tyz)/
     &                               uph(1,k) +
     &                               um(1)*rmat(3,1,k) +
     &                               um(2)*smat(3,1,k) +
     &                               um(3)*smat(4,1,k) -
     &                               gampr*dbyr(k)*ei(k)
c
                  smat(5,2,k)      = us4*txy/uph(1,k) +
     &                               um(1)*rmat(3,2,k) +
     &                               um(2)*smat(3,2,k) +
     &                               um(3)*smat(4,2,k) -
     &                               gampr*dbyr(k)*uph(2,k)
c
                  smat(5,3,k)      = us4*tyy/uph(1,k) +
     &                               um(1)*rmat(3,3,k) +
     &                               um(2)*smat(3,3,k) +
     &                               um(3)*smat(4,3,k) -
     &                               gampr*dbyr(k)*uph(3,k)
c
                  smat(5,4,k)      = us4*tyz/uph(1,k) +
     &                               um(1)*rmat(3,4,k) +
     &                               um(2)*smat(3,4,k) +
     &                               um(3)*smat(4,4,k) -
     &                               gampr*dbyr(k)*uph(4,k)
c
                  smat(5,5,k)      = gampr*dbyr(k)
c
                  qmat(4,1,k)      =-ds3*(
     &                               2.0*uph(4,k)*dbzr(k) -
     &                               uph(2,k)*dbxr(k) -
     &                               uph(3,k)*dbyr(k))
                  qmat(4,2,k)      = smat(3,2,k)
                  qmat(4,3,k)      = rmat(2,3,k)
                  qmat(4,4,k)      = 2.0*ds3*dbzr(k)
                  qmat(4,5,k)      = 0.0
c
                  qmat(5,1,k)      =-us4*(
     &                               uph(2,k)*txz +
     &                               uph(3,k)*tyz +
     &                               uph(4,k)*tzz)/
     &                               uph(1,k) +
     &                               um(1)*rmat(4,1,k) +
     &                               um(2)*smat(4,1,k) +
     &                               um(3)*qmat(4,1,k) -
     &                               gampr*dbzr(k)*ei(k)
c
                  qmat(5,2,k)      = us4*txz/uph(1,k) +
     &                               um(1)*rmat(4,2,k) +
     &                               um(2)*smat(4,2,k) +
     &                               um(3)*qmat(4,2,k) -
     &                               gampr*dbzr(k)*uph(2,k)
c
                  qmat(5,3,k)      = us4*tyz/uph(1,k) +
     &                               um(1)*rmat(4,3,k) +
     &                               um(2)*smat(4,3,k) +
     &                               um(3)*qmat(4,3,k) -
     &                               gampr*dbzr(k)*uph(3,k)
c
                  qmat(5,4,k)      = us4*tzz/uph(1,k) +
     &                               um(1)*rmat(4,4,k) +
     &                               um(2)*smat(4,4,k) +
     &                               um(3)*qmat(4,4,k) -
     &                               gampr*dbzr(k)*uph(4,k)
c
                  qmat(5,5,k)      = gampr*dbzr(k)
c    
                  DO 220 ivar=1,5
                     smat(2,ivar,k)= rmat(3,ivar,k)
                     qmat(2,ivar,k)= rmat(4,ivar,k)
                     qmat(3,ivar,k)= smat(4,ivar,k)
220               CONTINUE
c
400            CONTINUE
c
c              Treating the diagonal terms
c
               DO 230 ia=2,5
c
                  DO 235 ib=1,5
c
                     a11           = rmat(ia,ib,1)*dbx(1) + 
     &                               smat(ia,ib,1)*dby(1) + 
     &                               qmat(ia,ib,1)*dbz(1)
                     a22           = rmat(ia,ib,2)*dbx(2) + 
     &                               smat(ia,ib,2)*dby(2) + 
     &                               qmat(ia,ib,2)*dbz(2)
                     a33           = rmat(ia,ib,3)*dbx(3) + 
     &                               smat(ia,ib,3)*dby(3) + 
     &                               qmat(ia,ib,3)*dbz(3)
                     a44           = rmat(ia,ib,4)*dbx(4) +
     &                               smat(ia,ib,4)*dby(4) + 
     &                               qmat(ia,ib,4)*dbz(4)
c
                     diag(is1,ia,ib)    = diag(is1,ia,ib) + ait*a11
                     diag(is2,ia,ib)    = diag(is2,ia,ib) + ait*a22
                     diag(is3,ia,ib)    = diag(is3,ia,ib) + ait*a33
                     diag(is4,ia,ib)    = diag(is4,ia,ib) + ait*a44
c
235               CONTINUE
c
230            CONTINUE
c
c              Treating the extra-diagonal terms
c
               DO 250 ia=2,5
c
                  DO 260 ib=1,5
c
                    amat(ia,ib,1)  = ait*(
     &                               rmat(ia,ib,2)*dbx(1) + 
     &                               smat(ia,ib,2)*dby(1) +
     &                               qmat(ia,ib,2)*dbz(1))
                    bmat(ia,ib,1)  = ait*(
     &                               rmat(ia,ib,1)*dbx(2) + 
     &                               smat(ia,ib,1)*dby(2) +
     &                               qmat(ia,ib,1)*dbz(2))
                    amat(ia,ib,2)  = ait*(
     &                               rmat(ia,ib,3)*dbx(1) +
     &                               smat(ia,ib,3)*dby(1) +
     &                               qmat(ia,ib,3)*dbz(1))
                    bmat(ia,ib,2)  = ait*(
     &                               rmat(ia,ib,1)*dbx(3) +
     &                               smat(ia,ib,1)*dby(3) +
     &                               qmat(ia,ib,1)*dbz(3))
                    amat(ia,ib,3)  = ait*(
     &                               rmat(ia,ib,4)*dbx(1) +
     &                               smat(ia,ib,4)*dby(1) +
     &                               qmat(ia,ib,4)*dbz(1))
                    bmat(ia,ib,3)  = ait*(
     &                               rmat(ia,ib,1)*dbx(4) +
     &                               smat(ia,ib,1)*dby(4) +
     &                               qmat(ia,ib,1)*dbz(4))
                    amat(ia,ib,4)  = ait*(
     &                               rmat(ia,ib,3)*dbx(2) +
     &                               smat(ia,ib,3)*dby(2) +
     &                               qmat(ia,ib,3)*dbz(2))
                    bmat(ia,ib,4)  = ait*(
     &                               rmat(ia,ib,2)*dbx(3) + 
     &                               smat(ia,ib,2)*dby(3) +
     &                               qmat(ia,ib,2)*dbz(3))
                    amat(ia,ib,5)  = ait*(
     &                               rmat(ia,ib,4)*dbx(2) +
     &                               smat(ia,ib,4)*dby(2) +
     &                               qmat(ia,ib,4)*dbz(2))
                    bmat(ia,ib,5)  = ait*(
     &                               rmat(ia,ib,2)*dbx(4) +
     &                               smat(ia,ib,2)*dby(4) +
     &                               qmat(ia,ib,2)*dbz(4))
                    amat(ia,ib,6)  = ait*(
     &                               rmat(ia,ib,4)*dbx(3) +
     &                               smat(ia,ib,4)*dby(3) +
     &                               qmat(ia,ib,4)*dbz(3))
                    bmat(ia,ib,6)  = ait*(
     &                               rmat(ia,ib,3)*dbx(4) +
     &                               smat(ia,ib,3)*dby(4) +
     &                               qmat(ia,ib,3)*dbz(4))
c
260               CONTINUE
c
250            CONTINUE
c
               DO 270 k=1,6
c
                  js1              = nu(nor(k),jt)
                  js2              = nu(nex(k),jt)
c
                  DO 275 ivar=1,ndeg(js1)
c
                     iseg          = jaret(js1,ivar)
c
                     nub1          = nubo(1,iseg)
                     nub2          = nubo(2,iseg)
c
                     id1           = 50*(iseg - 1)
                     id2           = 50*(iseg - 1) + 25
c
                     IF ((js1 .EQ. nub1) .AND. (js2 .EQ. nub2)) THEN 
c
                        DO 280 ia=2,5
c
                           DO 285 ib=1,5
c
                              stmat(id2+5*(ia-1)+ib)   = 
     &                        stmat(id2+5*(ia-1)+ib) - bmat(ia,ib,k)
                              stmat(id1+5*(ia-1)+ib)   = 
     &                        stmat(id1+5*(ia-1)+ib) - amat(ia,ib,k)
c
285                        CONTINUE
c
280                     CONTINUE
c
                        GOTO 270
c
                     ENDIF
c
                     IF ((js1 .EQ. nub2) .AND. (js2 .EQ. nub1)) THEN 
c
                        DO 290 ia=2,5
c
                           DO 295 ib=1,5
c
                              stmat(id1+5*(ia-1)+ib)   = 
     &                        stmat(id1+5*(ia-1)+ib) - bmat(ia,ib,k)
                              stmat(id2+5*(ia-1)+ib)   = 
     &                        stmat(id2+5*(ia-1)+ib) - amat(ia,ib,k)
c
295                        CONTINUE
c
290                     CONTINUE
c
                        GOTO 270
c
                     ENDIF
c
275               CONTINUE
c
270            CONTINUE
c
             ENDIF
c
           ENDIF
c
cc         IF (nordre .EQ. 1) GOTO 1110
cc         IF (nordre .EQ. 4) GOTO 1105
cc
cc        Computing the hermitian nodal gradients
cc
cc         ait                       = us4*volt(jt)
cc
cc        DO 95 ivar=1,5
cc
cc            dx(ivar,is1)           = dx(ivar,is1) + dxt(ivar)*ait
cc            dx(ivar,is2)           = dx(ivar,is2) + dxt(ivar)*ait
cc            dx(ivar,is3)           = dx(ivar,is3) + dxt(ivar)*ait
cc            dx(ivar,is4)           = dx(ivar,is4) + dxt(ivar)*ait
cc
cc            dy(ivar,is1)           = dy(ivar,is1) + dyt(ivar)*ait
cc            dy(ivar,is2)           = dy(ivar,is2) + dyt(ivar)*ait
cc            dy(ivar,is3)           = dy(ivar,is3) + dyt(ivar)*ait
cc            dy(ivar,is4)           = dy(ivar,is4) + dyt(ivar)*ait
cc     
cc            dz(ivar,is1)           = dz(ivar,is1) + dzt(ivar)*ait
cc            dz(ivar,is2)           = dz(ivar,is2) + dzt(ivar)*ait
cc            dz(ivar,is3)           = dz(ivar,is3) + dzt(ivar)*ait
cc            dz(ivar,is4)           = dz(ivar,is4) + dzt(ivar)*ait
cc
cc95       CONTINUE
cc
cc1105     CONTINUE    
cc
cc         IF (nordre .EQ. 4) THEN
cc
cc            DO 1090 k=1,4
ccc
cc               is                  = nu(k,jt)
ccc
cc               DO 1095 ivar=1,5
cc
cc                  dx(ivar,is)      = (1.0 - pentel(is))*dxt(ivar) + 
cc     &            pentel(is)*0.5*(SIGN(1.0, dxt(ivar)) +
cc     &                            SIGN(1.0, dx(ivar,is)))*
cc     &            MIN(ABS(dxt(ivar)), ABS(dx(ivar,is)))
cc
cc                  dy(ivar,is)      = (1.0 - pentel(is))*dyt(ivar) + 
cc     &            pentel(is)*0.5*(SIGN(1.0, dyt(ivar)) +
cc     &                            SIGN(1.0, dy(ivar,is)))*
cc     &            MIN(ABS(dyt(ivar)), ABS(dy(ivar,is)))
cc
cc                  dz(ivar,is)      = (1.0 - pentel(is))*dzt(ivar) + 
cc     &            pentel(is)*0.5*(SIGN(1.0, dzt(ivar)) +
cc     &                            SIGN(1.0, dz(ivar,is)))*
cc     &            MIN(ABS(dzt(ivar)), ABS(dz(ivar,is)))
cc
cc                  pentel(is)       = 1.0
cc
cc1095           CONTINUE
cc
c1090        CONTINUE
cc
cc         ENDIF
cc
cc1110     CONTINUE
cc
           
 1000    CONTINUE
         
      DO 1010 is=1,ns
         dt                        = MIN(dt, dtl(is))  
1010  CONTINUE  
c
      dt                           = cfl*dt
c
      dt                           = MIN(dt, tmax - t)
c
      IF (iloc .EQ. 1) THEN
         DO 1020 is=1,ns
            dtl(is)                = cfl*dtl(is)
1020     CONTINUE
      ELSE
         DO 1030 is=1,ns
            dtl(is)                = dt
1030     CONTINUE
      ENDIF
cc
cc     Completing the non-limited nodal gradients
cc
cc      IF ((nordre .EQ. 2) .OR. (nordre .EQ. 3)) THEN 
cc
cc         DO 1040 ivar=1,5
cc
cc            DO 1050 is=1,ns
cc
cc               ais                 = 1.0/vols(is)
cc
cc               dx(ivar,is)         = dx(ivar,is)*ais
cc               dy(ivar,is)         = dy(ivar,is)*ais
cc               dz(ivar,is)         = dz(ivar,is)*ais
cc
cc1050        CONTINUE
cc
cc1040     CONTINUE
cc
cc      ENDIF
cc

      
      END
