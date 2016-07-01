
      SUBROUTINE RESU3D
c---------------------------------------------------------------------
c Output extremum values of the physical solution   
c Saves the physical solution of file
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
      INTEGER nps   , is, ii, ifac, jj
      REAL*8    rmax  , rmin   , pmax  , pmin
      REAL*8    uxmax , uxmin  , vxmax , vxmin,  wxmax , wxmin
      REAL*8    xcmax , xcmin, emin, emax
      REAL*8    rm, ux, vx, wx , pm, xc
      REAL*8    coefr(5), pression
c
      coefr(1)                     = rhoref
      coefr(2)                     = rhoref*vref
      coefr(3)                     = rhoref*vref
      coefr(4)                     = rhoref*vref
      coefr(5)                     = pref
c
      rmax                         =-1.0e+30
      rmin                         = 1.0e+30
      uxmax                        =-1.0e+30
      uxmin                        = 1.0e+30
      vxmax                        =-1.0e+30
      vxmin                        = 1.0e+30
      wxmax                        =-1.0e+30
      wxmin                        = 1.0e+30
      pmax                         =-1.0e+30
      pmin                         = 1.0e+30
      xcmax                        =-1.0e+30
      xcmin                        = 1.0e+30
      emax                         =-1.0e+30
      emin                         = 1.0e+30
c
      nps                          = 0
c
      DO 10 is=1,ns
c
         rm                        = ua(1,is)*coefr(1)
         rmax                      = MAX(rmax, rm)
         rmin                      = MIN(rmin, rm)
c
         ux                        = ua(2,is)*coefr(2)/rm
         vx                        = ua(3,is)*coefr(2)/rm
         wx                        = ua(4,is)*coefr(2)/rm
         uxmax                     = MAX(uxmax, ux)
         uxmin                     = MIN(uxmin, ux)
         vxmax                     = MAX(vxmax, vx)
         vxmin                     = MIN(vxmin, vx)
         wxmax                     = MAX(wxmax, wx)
         wxmin                     = MIN(wxmin, wx)
c
         emax                      = MAX(ua(5,is),emax)
         emin                      = MIN(ua(5,is),emin)
c
         pm                        = gam1*(ua(5,is)*coefr(5) - 
     &                               0.5*rm*(ux*ux + vx*vx + wx*wx))
         pmax                      = MAX(pmax, pm)
         pmin                      = MIN(pmin, pm)
c
         xc                        = SQRT((ux*ux + vx*vx + wx*wx)/
     &                                    (gam*pm/rm))
         IF (xc .GE. 1.0) nps      = nps + 1
         xcmax                     = MAX(xcmax, xc)
         xcmin                     = MIN(xcmin, xc)
c
10    CONTINUE
c
c
      pmnew= pmax
      
      WRITE(6, *) 
     &   '          --------------------------------------------------'
      WRITE(6, 110) ' Iteration number             : ', kt
      WRITE(6, 109) ' Physical time                : ', t
      WRITE(6, 110) ' Number of supersonic points  : ', nps
      WRITE(6, 111) ' Min and Max Density          : ', rmin,  rmax
      WRITE(6, 111) ' Min and Max Pressure         : ', pmin,  pmax
      WRITE(6, 111) ' Min and Max Mach number      : ', xcmin, xcmax
      WRITE(6, 111) ' Min and Max X-Velocity       : ', uxmin, uxmax
      WRITE(6, 111) ' Min and Max Y-Velocity       : ', vxmin, vxmax
      WRITE(6, 111) ' Min and Max Z-Velocity       : ', wxmin, wxmax
      WRITE(6, 111) ' Min and Max Energy           : ', emin, emax
      WRITE(6, *) ' '
c      WRITE(6, *) ' '
c
109   FORMAT(10x,a32,e14.7)
110   FORMAT(10x,a32,i8)
111   FORMAT(10x,a32,e14.7,2x,e14.7)
      
c
      RETURN
      END
