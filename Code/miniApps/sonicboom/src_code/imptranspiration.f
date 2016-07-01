      SUBROUTINE IMPTRANSPIRATION(CTRL)
c
c*** Conditions de transpiration en implicite.
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c
      INTEGER ISP,IS,I
      REAL*8 UG,VG,WG,QG,QT,DV2,DV3,DV4,DV5
      REAL*8 PI, AALPHA
      REAL*8 VNCOQ(3,NNSP), PRESSION, CTRL(NNSP)
      REAL*8 DQDRHO(NSMAX), DQDRHOU(NSMAX), DQDRHOV(NSMAX),
     $     DQDRHOW(NSMAX), DQDE(NSMAX)
      REAL*8 ray, normalis, x, ddx
c
C**** cas ou l'on se trouve soit sur la carlingue, soit sur l'aile, soit
C     sur la tuyere,
C     ==> TRANSPIRATION
C
C     PSI_bords = Q(T)*(RHO,RHOU,RHOV,RHOW,E) + (0,P.nx,P.ny,P.nz,P.Q(T))
c                -->   -->  --->    -->
C     avec Q(T) = V . (NO - N(T))   NO = (nx,ny,nz) est stoque dans VNO
C                                   --->
c                                   N(T) est stoque dans VNCOQ
c
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
c$$$
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
c$$$c
c$$$               END DO

c
c      ENDIF
c
         DO ISP = 1,NSP
c
            IS = NODE2D3D(ISP)
C
            UG = UA(2,IS)/UA(1,IS)
            VG = UA(3,IS)/UA(1,IS)
            WG = UA(4,IS)/UA(1,IS)
c
            QT = UG*(VNO(1,ISP)-VNCOQ(1,ISP)) +
     $           VG*(VNO(2,ISP)-VNCOQ(2,ISP)) +
     $           WG*(VNO(3,ISP)-VNCOQ(3,ISP))

            QG = 0.5 * (UG**2 + VG**2 + WG**2)
C
            PRESSION = GAM1*(UA(5,IS)-(UA(2,IS)**2 +
     $           UA(3,IS)**2 + UA(4,IS)**2)/(2.*UA(1,IS)))
c
            DV2 = GAM1*VNO(1,ISP)
            DV3 = GAM1*VNO(2,ISP)
            DV4 = GAM1*VNO(3,ISP)
            DV5 = GAM1 * QT
c
            DQDRHO(IS) = - QT/UA(1,IS)
            DQDRHOU(IS) = (VNO(1,ISP)-VNCOQ(1,ISP))/UA(1,IS)
            DQDRHOV(IS) = (VNO(2,ISP)-VNCOQ(2,ISP))/UA(1,IS)
            DQDRHOW(IS) = (VNO(3,ISP)-VNCOQ(3,ISP))/UA(1,IS)
            DQDE(IS) = 0.
C
            DIAG(IS,1,1) = DIAG(IS,1,1) + QT + UA(1,IS)*DQDRHO(IS)
            DIAG(IS,1,2) = DIAG(IS,1,2) + UA(1,IS)*DQDRHOU(IS)
            DIAG(IS,1,3) = DIAG(IS,1,3) + UA(1,IS)*DQDRHOV(IS)
            DIAG(IS,1,4) = DIAG(IS,1,4) + UA(1,IS)*DQDRHOW(IS)
            DIAG(IS,1,5) = DIAG(IS,1,5) + UA(1,IS)*DQDE(IS)

C
            DIAG(IS,2,1) = DIAG(IS,2,1) + QG*DV2 + UA(2,IS)*DQDRHO(IS)
            DIAG(IS,2,2) = DIAG(IS,2,2) + QT - UG*DV2 +
     $                                       UA(2,IS)*DQDRHOU(IS)
            DIAG(IS,2,3) = DIAG(IS,2,3) - VG*DV2 +
     $                                       UA(2,IS)*DQDRHOV(IS)
            DIAG(IS,2,4) = DIAG(IS,2,4) - WG*DV2 +
     $                                       UA(2,IS)*DQDRHOW(IS)
            DIAG(IS,2,5) = DIAG(IS,2,5) + DV2 + UA(2,IS)*DQDE(IS)
C
            DIAG(IS,3,1) = DIAG(IS,3,1) + QG*DV3 + UA(3,IS)*DQDRHO(IS)
            DIAG(IS,3,2) = DIAG(IS,3,2) - UG*DV3 + UA(3,IS)*DQDRHOU(IS)
            DIAG(IS,3,3) = DIAG(IS,3,3) + QT - VG*DV3 +
     $                                          UA(3,IS)*DQDRHOV(IS)
            DIAG(IS,3,4) = DIAG(IS,3,4) - WG*DV3 +
     $                                          UA(3,IS)*DQDRHOW(IS)
            DIAG(IS,3,5) = DIAG(IS,3,5) + DV3 + UA(3,IS)*DQDE(IS)
C
            DIAG(IS,4,1) = DIAG(IS,4,1) + QG*DV4 + UA(4,IS)*DQDRHO(IS)
            DIAG(IS,4,2) = DIAG(IS,4,2) - UG*DV4 + UA(4,IS)*DQDRHOU(IS)
            DIAG(IS,4,3) = DIAG(IS,4,3) - VG*DV4 + UA(4,IS)*DQDRHOV(IS)
            DIAG(IS,4,4) = DIAG(IS,4,4) + QT - WG*DV4 +
     $                                       UA(4,IS)*DQDRHOW(IS)
            DIAG(IS,4,5) = DIAG(IS,4,5) + DV4 + UA(4,IS)*DQDE(IS)
C
            DIAG(IS,5,1) = DIAG(IS,5,1) + QG*DV5 +
     $                   (UA(5,IS) + PRESSION)*DQDRHO(IS)
            DIAG(IS,5,2) = DIAG(IS,5,2) - UG*DV5 +
     $                   (UA(5,IS) + PRESSION)*DQDRHOU(IS)
            DIAG(IS,5,3) = DIAG(IS,5,3) - VG*DV5 +
     $                   (UA(5,IS) + PRESSION)*DQDRHOV(IS)
            DIAG(IS,5,4) = DIAG(IS,5,4) - WG*DV5 +
     $                   (UA(5,IS) + PRESSION)*DQDRHOW(IS)
            DIAG(IS,5,5) = DIAG(IS,5,5) + QT + DV5 +
     $                   (UA(5,IS) + PRESSION)*DQDE(IS)
C
         END DO
C
         RETURN
         END
