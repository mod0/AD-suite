    

      SUBROUTINE MOVMSH
c     ----------------------------------------------------------------- 
c     Moves the mesh and distributes the new mesh vertices 
c     coordinates to neighboring vertices
c     -----------------------------------------------------------------
      include 'Param3D.h'
c     -----------------------------------------------------------------
c     Local variables definition
      INTEGER ns0, kjac0
      INTEGER is, iseg    , nubo1   , nubo2
      REAL*8    rsjac0
      REAL*8    stfseg, aix , aiy     , aiz
      REAL*8    cosT  , sinT
      REAL*8    stifsm(nsmax)
c
c     Initializations
c
      maxjac                       = 4
      rsjacf                       = 1.0e-6
      ns0                          = 0
c
c     Initialization of the elastic center 
c
      xce                          = 0.2
      yce                          = 0.0 
      zce                          = 0.0
c
      DO 10 is=1,ns
         stifsm(is)                = 0.0
c
         xt0(is)                   = coor(1,is) 
         yt0(is)                   = coor(2,is) 
         zt0(is)                   = coor(3,is)
c 
         deltx(1,is)               = 0.0
         delty(1,is)               = 0.0
         deltz(1,is)               = 0.0
         deltx(2,is)               = 0.0
         delty(2,is)               = 0.0
         deltz(2,is)               = 0.0
c
10    CONTINUE
c  
      cosT                         = cos(tetaT)
      sinT                         = sin(tetaT)
c
c     Updating the position of the mesh vertices placed on 
c     the body surface using the rotation of angle tetaT and 
c     center of rotation xce, yce in the x-y plane
c
      DO 20 is=1,ns
c
c        Selecting the mesh vertices placed on the body
c
         IF (logfr(is) .EQ. -2) THEN 
c
            coco(1,is)             = (xt0(is) - xce)*cosT - 
     &                               (yt0(is) - yce)*sinT + xce
            coco(2,is)             = (xt0(is) - xce)*sinT + 
     &                               (yt0(is) - yce)*cosT + yce
c
         ENDIF   
c          
20    CONTINUE
c
c     Computing the predicted displacements
c     They must be equal to zero for the mesh vertices 
c     placed on the farfield boundaries
c     They are known for the mesh vertices placed on the body suface
c     
      DO 30 is=1,ns
c
         deltxp(is)                = 0.0
         deltyp(is)                = 0.0
         deltzp(is)                = 0.0
c
c        Selecting the submesh vertices placed on the body
c
         IF (logfr(is) .EQ. -2) THEN 
c
            deltxp(is)             = coco(1,is) - coor(1,is)
            deltyp(is)             = coco(2,is) - coor(2,is)
            deltzp(is)             = coco(3,is) - coor(3,is)
c
            deltx(1,is)            = deltxp(is)
            delty(1,is)            = deltyp(is)
            deltz(1,is)            = deltzp(is)
c
         ENDIF
c
c        Computing the predicted displacements for the interior 
c        mesh vertices 
c
         IF (logfr(is) .EQ. 0) THEN 
c
            ns0                    = ns0 + 1
c
            deltxp(is)             = 2.0*deltx(1,is) - deltx(2,is)
            deltyp(is)             = 2.0*delty(1,is) - delty(2,is)
            deltzp(is)             = 2.0*deltz(1,is) - deltz(2,is)
c
c           Updating the interior mesh vertices displacements 
c           at physical time iteration kt-1
c
            deltx(2,is)            = deltx(1,is)
            delty(2,is)            = delty(1,is)
            deltz(2,is)            = deltz(1,is)
c
c           Initializations for the Jacobi iterations
c
            deltx(1,is)            = 0.0
            delty(1,is)            = 0.0
            deltz(1,is)            = 0.0
c
         ENDIF
c          
30    CONTINUE
c
c     Loop on local list of edges
c
      DO 50 iseg=1,nseg
c
c        Local indexing of the vertices of the current edge
c
         nubo1                     = nubo(1,iseg)
         nubo2                     = nubo(2,iseg)
c        
c        Computing the stiffness of the current edge
c
         aix                       = coor(1,nubo2) - coor(1,nubo1)
         aiy                       = coor(2,nubo2) - coor(2,nubo1)
         aiz                       = coor(3,nubo2) - coor(3,nubo1)
c
         stfseg                    = 1.0/SQRT(aix*aix + 
     &                                        aiy*aiy + aiz*aiz)
c
c        Gathering the stiffness of the current edge into
c        the corresponding edge end-points stiffnesses
c
         stifsm(nubo1)             = stifsm(nubo1) + stfseg 
c
         stifsm(nubo2)             = stifsm(nubo2) + stfseg 
c
50    CONTINUE  
c
c     Updating the coordinates of the interior mesh vertices 
c     Jacobi iterations
c
      kjac                         = 0
      kjac0                        = 1
c
1000  kjac                         = kjac + 1
c
c     Loop on local list of edges
c
      DO 100 iseg=1,nseg
c
c        Local indexing of the vertices of the current edge
c
         nubo1                     = nubo(1,iseg)
         nubo2                     = nubo(2,iseg)
c        
c        Computing the stiffness of the current edge
c
         aix                       = coor(1,nubo2) - coor(1,nubo1)
         aiy                       = coor(2,nubo2) - coor(2,nubo1)
         aiz                       = coor(3,nubo2) - coor(3,nubo1)
c
         stfseg                    = 1.0/SQRT(aix*aix + 
     &                                        aiy*aiy + aiz*aiz)
c
c        Computing the displacements for the interior    
c        active submesh vertices         
c
         IF (logfr(nubo1) .EQ. 0) THEN 
c
             deltx(1,nubo1)        = deltx(1,nubo1) + 
     &                               stfseg*deltxp(nubo2)
             delty(1,nubo1)        = delty(1,nubo1) + 
     &                               stfseg*deltyp(nubo2)
             deltz(1,nubo1)        = deltz(1,nubo1) + 
     &                               stfseg*deltzp(nubo2)
c
         ENDIF
c
         IF (logfr(nubo2) .EQ. 0) THEN 
c
             deltx(1,nubo2)        = deltx(1,nubo2) + 
     &                               stfseg*deltxp(nubo1)
             delty(1,nubo2)        = delty(1,nubo2) + 
     &                               stfseg*deltyp(nubo1)
            deltz(1,nubo2)         = deltz(1,nubo2) + 
     &                               stfseg*deltzp(nubo1)
c
          ENDIF
c
100   CONTINUE  
c
c     Updating the displacements for the interior mesh vertices
c
      DO 150 is=1,ns
c
         IF (logfr(is) .EQ. 0) THEN 
c
            deltx(1,is)            = deltx(1,is)/stifsm(is)        
            delty(1,is)            = delty(1,is)/stifsm(is)        
            deltz(1,is)            = deltz(1,is)/stifsm(is)        
c
         ENDIF         
c
150   CONTINUE
c
c     Computing the residual 
c
      rsjac                        = 0.0
c
      DO 160 is=1,ns
c
         IF (logfr(is) .EQ. 0) THEN
c
            rsjac                  = rsjac + 
     &                               (deltxp(is) - deltx(1,is))*
     &                               (deltxp(is) - deltx(1,is)) +
     &                               (deltyp(is) - delty(1,is))*
     &                               (deltyp(is) - delty(1,is)) +
     &                               (deltzp(is) - deltz(1,is))*
     &                               (deltzp(is) - deltz(1,is))
c
         ENDIF
c
160   CONTINUE
c
      rsjac                        = SQRT(rsjac)/FLOAT(ns0)
c
      IF (kjac .EQ. kjac0) rsjac0  = rsjac
c
      rsjac                        = rsjac/rsjac0

      rsjac                        = SQRT(rsjac)/FLOAT(ns0)
c
      IF (kjac .EQ. kjac0) rsjac0  = rsjac
c
      rsjac                        = rsjac/rsjac0
c
c     Testing for convergence or maximum number of Jacobi iterations
c
      IF (kjac  .EQ. maxjac) GOTO 1100
      IF (rsjac .LT. rsjacf) GOTO 1100
c
c     Initializations for the next Jacobi iteration
c
      DO 170 is=1,ns
c
         IF (logfr(is) .EQ. 0) THEN 
c
            deltxp(is)             = deltx(1,is)
            deltyp(is)             = delty(1,is)
            deltzp(is)             = deltz(1,is)
            deltx(1,is)            = 0.0
            delty(1,is)            = 0.0
            deltz(1,is)            = 0.0
c
         ENDIF
c
170   CONTINUE
c
      GOTO 1000
c
1100  CONTINUE
c
c     Updating the positions of the interior submesh vertices
c
      DO 2000 is=1,ns
c
         IF (logfr(is) .EQ. 0) THEN 
c
            coco(1,is)             = coor(1,is) + deltx(1,is) 
            coco(2,is)             = coor(2,is) + delty(1,is) 
            coco(3,is)             = coor(3,is) + deltz(1,is) 
c
         ENDIF
c
2000  CONTINUE
c
      RETURN
      END
