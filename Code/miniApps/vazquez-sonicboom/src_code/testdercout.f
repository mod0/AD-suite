
      SUBROUTINE TESTDERCOUT(CTRL)
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER ISP,KVAR,IS,ISP0,ISV0,I
      REAL*8 EPS,PRESSION,DIFF,COUTO,DMAX,Cout1,CoutPres
      REAL*8 CTRL(NNSP), VNCOQ(3,NNSP),CX,CY,CCD,CCL,DENOMINATEUR
C
      EPS = 1.E-03
C
      CX = 0.
      CY = 0.
C
c      CALL NORMCOQ(CTRL,VNCOQ)
C
      COUT1 = 0.
c
      DO 10 ISP = 1,NSP
         IS = NODE2D3D(ISP)
         PRESSION = GAM1 * (UA(5,IS) - 0.5*(UA(2,IS)**2 + UA(3,IS)**2
     $              + UA(4,IS)**2)/UA(1,IS))
c         cx = cx + (vncoq(1,isp)*pression)
c         cy = cy + (vncoq(2,isp)*pression)
      COUT1 = COUT1 + AIRESP0(ISP)/3. *(PRESSION - PDESP(IS))**2
 10   CONTINUE
c
c      Denominateur = roin*(uxin**2+uyin**2+uzin**2)/2.
c
c      CX = CX/Denominateur
c      CY = CY/Denominateur
C
c      CCL= -sin(tetacdcl)*CX + cos(tetacdcl)*CY
c      CCD= cos(tetacdcl)*CX + sin(tetacdcl)*CY
C
c      COUTO = coeftrainee*(CCD-CDTARGET)**2 + (CCL - CLTARGET)**2
c     $        + COEFPRES*Cout1
C
c      PRINT*,'LA FONCTIONNELLE COUT A POUR VALEUR : ',COUTO
      PRINT*,'LA FONCTIONNELLE COUT A POUR VALEUR : ',COUT1
c
      DO I = 1,8
c
      DMAX = 0.
C
      DO ISP0 =1,NSP
         ISV0 = NODE2D3D(ISP0)
         DO KVAR = 1,5
            UA(KVAR,ISV0) = UA(KVAR,ISV0) + EPS
c            CX = 0.
c            Cy = 0.
            CoutPres = 0.
            DO  ISP = 1,NSP
               IS = NODE2D3D(ISP)
               PRESSION = GAM1 * (UA(5,IS) - 0.5*(UA(2,IS)**2 +
     $              UA(3,IS)**2+ UA(4,IS)**2)/UA(1,IS))
c
c               cx = cx + (vncoq(1,isp)*pression)
c               cy = cy + (vncoq(2,isp)*pression)
c
       COUTPRES = COUTPRES + AIRESP0(ISP)/3. *(PRESSION - PDESP(IS))**2
c
            END DO
c
c            CX = CX/Denominateur
c            CY = CY/Denominateur
C
c            CCL= -sin(tetacdcl)*CX + cos(tetacdcl)*CY
c            CCD= cos(tetacdcl)*CX + sin(tetacdcl)*CY
c
c            COUT = coeftrainee*(CCD-CDTARGET)**2 + (CCL - CLTARGET)**2 +
c     $             coefpres*coutpres
c
c            DIFF = (COUT-COUTO)/EPS
            DIFF = (COUTPRES-COUT1)/EPS
            DMAX = MAX(DMAX, ABS(DIFF-DJDW(KVAR,ISV0)))
C
                  UA(KVAR,ISV0) = UA(KVAR,ISV0) - EPS
         END DO
      END DO
C
      PRINT*,'DMAX =',DMAX,' EPSILON = ',EPS
C
      EPS = EPS/10.
c
      END DO
C
      RETURN
      END
