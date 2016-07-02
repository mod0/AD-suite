

      SUBROUTINE DCOUTDW(CTRL,CCL,CCD)
C
C*** Cette procedure calcule la derivee de la fonctionnelle cout J(W,gamma) par
C    rapport a W. 
C
C        nsp       
C         __
C    dJ   \   2              dPression
C    -- = /  ---*AIRES(isp)* --------- (is)*(Pression(is)-Pres_des(is))
C    dW   --  3                  dW
C       isp=1
C   is=node2d3d(isp)
C
c                                 __                   
c                                 \          dPression   
c      + 2*omega*( CD - CDTARGET)*/   A(isp)*--------- 
c                                 --             dW    
c                                isp        
c                          __
c                          \         dPression
c      + 2*(CL - CLTARGET)*/  B(isp)*---------
c                          --           dW
c                          isp
c
c     A(isp) = cos(teta)*N(1,isp) + sin(teta)*N(2,isp)
c     B(isp) = cos(teta)*N(2,isp) - sin(teta)*N(1,isp)
c
C   La derivee est stockee dans DJDW(1:5,1:ns)
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER ISP, IS, K
      REAL*8 U, V, W, PRESSION, DPDW(5,NSMAX), DIFF,DIFF1, CCL, CCD
      REAL*8 CTRL(NNSP), VNCOQ(3,NNSP), COEF, DIFF2
c
      write(6,*) 'Entree dans DCoutDW.'
c
      COEF = ROIN*(UXIN**2+UYIN**2+UZIN**2)/2.
C
      CALL NORMCOQ(CTRL,VNCOQ)
C
      DO 10 IS = 1,NS
         DO 20 K =1,5
            DJDW(K,IS) = 0.
 20      CONTINUE
 10   CONTINUE
C
      DO 30 ISP = 1,NSP
c
         IS = NODE2D3D(ISP)
c
         U = UA(2,IS)/UA(1,IS)
         V = UA(3,IS)/UA(1,IS)
         W = UA(4,IS)/UA(1,IS)
c
         PRESSION = GAM1*(UA(5,IS) - 0.5*UA(1,IS)*(U**2 + V**2 + W**2))
c
c
c
         DPDW(1,IS) =  GAM1*0.5*(U**2 + V**2 + W**2)
         DPDW(2,IS) = -GAM1*U
         DPDW(3,IS) = -GAM1*V
         DPDW(4,IS) = -GAM1*W
         DPDW(5,IS) = GAM1
c
         DIFF = (COS(TETACDCL)*VNCOQ(1,ISP)+SIN(TETACDCL)
     $            *VNCOQ(2,ISP))/COEF
         DIFF1 = (COS(TETACDCL)*VNCOQ(2,ISP)-SIN(TETACDCL)
     $            *VNCOQ(1,ISP))/COEF
         DIFF2 = PRESSION - PDESP(IS)
c
         DO 40 K = 1,5
            DJDW(K,IS) = DJDW(K,IS)
     $           + 2*COEFTRAINEE*(CCD-CDTARGET)*DIFF
     $             *DPDW(K,IS) + 2*(CCL-CLTARGET)*DIFF1*DPDW(K,IS)
     $             + COEFPRES*2.*AIRESP0(ISP)*DIFF2*DPDW(K,IS)/3.
 40      CONTINUE
C
 30   CONTINUE
C
      write(6,*) 'Sortie de DCoutdW.'
c
      RETURN
      END
