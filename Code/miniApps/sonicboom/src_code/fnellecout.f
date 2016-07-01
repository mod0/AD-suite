      SUBROUTINE FNELLECOUT(CCD,CCL)
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      REAL*8 CCD, CCL, COUT1
      INTEGER ISP, IS
      REAL*8 PRESSION
c
C**** Calcul de la fonctionnelle cout J(W). 
C            nsp
C            --  1                                                   2
C     J(W) = \   -  Aire(isp) * (Pression(is) - Pression_desiree(is)) 
C            /   3
C            --
C          isp=1
C     is=node2d3d(isp)
C                                         2                  2
C            + coeftrainee*(CD - CDTARGET)  + (CL - CLTARGET)
C
C     La fonctionnelle est stockee dans le scalaire Cout
C
      COUT1 = 0.
C
      DO 10 ISP = 1,NSP
         IS = NODE2D3D(ISP)
         PRESSION = GAM1 * (UA(5,IS) - 0.5*(UA(2,IS)**2 + UA(3,IS)**2
     $              + UA(4,IS)**2)/UA(1,IS))
      COUT1 = COUT1 + AIRESP0(ISP)/3. *(PRESSION - PDESP(IS))**2
 10   CONTINUE
C
      COUT = COEFTRAINEE*(CCD-CDTARGET)**2 + (CCL - CLTARGET)**2 +
     $       COEFPRES*COUT1
c
      RETURN
      END
