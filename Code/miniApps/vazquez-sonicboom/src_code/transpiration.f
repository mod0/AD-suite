      SUBROUTINE TRANSPIRATION(PSI,CTRL)
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER ISP,IS,I
      REAL*8 QT, PRESSION, U, V, W
      REAL*8 PI, AALPHA
      REAL*8 PSI(5,NSMAX), CTRL(NNSP), VNCOQ(3,nnsp)
      REAL*8 ray, normalis, x, ddx
c
c*** Inclusion des conditions de transpiration : ITRANS=1
c    On ne les inclut que sur la carlingue (cas du falcon) ou bien l'aile.
c    La carlingue et l'aile sont reperes si lo1fac = 200 ou si logfac = -2.
c    Transpiration : on traite le plan de symetrie comme avant
c                       et on applique la transpiration sur la coque.
C
C**** cas ou l'on se trouve soit sur la carlingue, soit sur l'aile
C     ==> TRANSPIRATION
C
C     PSI_bords = Q(T)*(RHO,RHOU,RHOV,RHOW,E) + (0,P.nx,P.ny,P.nz,P.Q(T))
c                -->   -->  --->    -->
C     avec Q(T) = V . (NO - N(T))   NO = (nx,ny,nz) est stoque dans VNO
C                                   --->
c                                   N(T) est stoque dans VNCOQ
C
c      IF ((COEFM1.NE.10).AND.(NPRES.EQ.0)) THEN 
c   
c  Cas ou l'on ne se trouve pas sur la tuyere et ou l'on veut calculer
c  la pression desiree par incidence simulee par transpiration.
C
C*** INCIDENCE DE 3.06 DEGRES. ROTATION 2D, DANS LE PLAN (X,Y).
C
c         PI = 4.0 * ATAN(1.0)
c         AALPHA = 3.06*PI/180
C
c         DO ISP = 1,NSP
c          VNCOQ(1,ISP) = VNO(1,ISP)*COS(AALPHA)+VNO(2,ISP)*SIN(AALPHA)
c          VNCOQ(2,ISP) = -VNO(1,ISP)*SIN(AALPHA)+VNO(2,ISP)*COS(AALPHA)
c          VNCOQ(3,ISP) = VNO(3,ISP)
c         END DO
c
c      ENDIF
C
c      IF ((COEFM1.EQ.10).OR.((COEFM1.NE.10).AND.(NPRES.NE.0))) THEN
c
         CALL NORMCOQ(CTRL,VNCOQ)
c
c$$$         if (kt.eq.1) then
c$$$c
c$$$            print*,'Transpiration.'
c$$$         do isp = 1,nsp
c$$$            is = node2d3d(isp)
c$$$            write(6,222) isp,vncoq(1,isp),vncoq(2,isp),coor(1,is)
c$$$     $                   ,coor(2,is)
c$$$c
c$$$ 222        format(i2,1x,4(e10.3,1x))
c$$$         end do
c$$$c
c$$$         write(6,*)'  '
c$$$         write(6,*) 'Normales specifiees'
c$$$c
c$$$      endif
c
c$$$
c$$$c
c$$$               ray=2.99719
c$$$c
c$$$               DO ISP = 1,NSP
c$$$c
c$$$                  NORMALIS = 0.
c$$$                  DO I = 1,3
c$$$                     NORMALIS = NORMALIS + VNO(I,ISP)*VNO(I,ISP)
c$$$                  END DO
c$$$                  NORMALIS = SQRT(NORMALIS)
c$$$c
c$$$                  VNCOQ(1,ISP) = 0.
c$$$                  VNCOQ(2,ISP) = - NORMALIS
c$$$                  VNCOQ(3,ISP) = 0.
c$$$c
c$$$                  IS = NODE2D3D(ISP) 
c$$$                  IF (((COOR(1,IS).GE.0.).AND.(COOR(1,IS).LE.2.)).AND.
c$$$     $                 (COOR(2,IS).EQ.5.)) THEN
c$$$                     x = 0.5*coor(1,is) + 0.5
c$$$                     if (abs(x-1.0).le.0.5) then
c$$$                     ddx=x-1.0
c$$$                     vncoq(1,isp) = ddx/sqrt(ray*ray-ddx*ddx)
c$$$                     ddx = sqrt(vncoq(1,isp)**2 + 1.)
c$$$                     vncoq(1,isp) = -vncoq(1,isp)*normalis/ddx
c$$$                     vncoq(2,isp) = normalis/ddx
c$$$c
c$$$
c$$$                     endif
c$$$                  ENDIF
c$$$c
c$$$         if (kt.eq.1) then
c$$$            write(6,222) isp,vncoq(1,isp),vncoq(2,isp),coor(1,is)
c$$$     $                   ,coor(2,is)
c$$$         endif
c$$$
c$$$               END DO
      
c
c      ENDIF
C



         isp=0
C$AD II-LOOP
         DO ISP = 1,NSP
C
            IS = NODE2D3D(isp)
c
            U = UA(2,IS)/UA(1,IS)
            V = UA(3,IS)/UA(1,IS)
            W = UA(4,IS)/UA(1,IS)
c
            QT = U*(VNO(1,ISP)-VNCOQ(1,ISP)) +
     $           V*(VNO(2,ISP)-VNCOQ(2,ISP)) +
     $           W*(VNO(3,ISP)-VNCOQ(3,ISP))
c
            PRESSION = GAM1*(UA(5,IS)-(UA(2,IS)**2 +
     $           UA(3,IS)**2 + UA(4,IS)**2)/(2.d0*UA(1,IS)))
C
            PSI(1,IS) = PSI(1,IS) - UA(1,IS)*QT
            PSI(2,IS) = PSI(2,IS) - UA(2,IS)*QT - PRESSION*VNO(1,ISP)
            PSI(3,IS) = PSI(3,IS) - UA(3,IS)*QT - PRESSION*VNO(2,ISP)
            PSI(4,IS) = PSI(4,IS) - UA(4,IS)*QT - PRESSION*VNO(3,ISP)
            PSI(5,IS) = PSI(5,IS) - (UA(5,IS) + PRESSION)*QT
         END DO            
C 
      RETURN
      END
