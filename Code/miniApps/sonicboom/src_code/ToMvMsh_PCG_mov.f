      SUBROUTINE MOVMSH_PCG_mov(coor,nu,logfr,nubo,kspring,ns,nt,
     &                      nseg,xw)
c     ----------------------------------------------------------------- 
c
c     Moves the mesh and distributes the new mesh vertices 
c     coordinates to neighboring vertices. Use of PCG.
c
c     ----------------------------------------------------------------- 
c
      implicit none
c
      include 'param.h'
c
      integer nu(4,ntmax), logfr(nsmax), nubo(2,nsegmax)
      integer kspring, ns, nt, nseg
      real*8 coor(3,nsmax), xw(3,nsmax)
      real*8 deltx(nsmax,2), delty(nsmax,2), deltz(nsmax,2)
      real*8 rx(nsmax), ry(nsmax), rz(nsmax)
      real*8 rsjacf
      integer imvmsh, kpcg, maxjac
      real*8 Apx(nsmax), Apy(nsmax), Apz(nsmax), respcg
      
c
      INTEGER isglow , isgup
      INTEGER is  , iseg, nubo1 , nubo2
      REAL*8    res0PCG, pAp, zr
      REAL*8    betaPCG, alphaPCG
      REAL*8    zx(nsmax), zy(nsmax), zz(nsmax)
      REAL*8    px(nsmax), py(nsmax), pz(nsmax)
      REAL*8    rp
c
      real*8 epsilon

c
c
c
      real*8 tb11, tb21, tb31
c
c
      real*8 xlocbar1,
     $     xlocbar2,
     $     xlocbar3
c
      real*8 xmag
c
      real*8 kbar11, kbar22, kbar33, kbar44, kbar55, kbar66
c
c
      real*8 a
c
c
      real*8 fac1, fac2, fac3
c
      real*8 lowlambda,highlambda
c
      integer ele
      integer toflag, barflag
      integer printres
      integer nufaces,nucycles,map_pointer(4,3,4),mvlgfr(nsmax)
c
      real*8 bdot_mov
      external bdot_mov
c
c     ***********************************************************
c
      open(20,file='moco',status='unknown')
      imvmsh=3
      maxjac=50
      rsjacf=1.e-9
c
      do is=1,ns
        mvlgfr(is)=0
      enddo
c
      call MovMap_mov(nufaces,nucycles,lowlambda,highlambda,map_pointer)

c
      barflag   = 0
      toflag    = 0
      printres  = 1
c
c
c
      if(kspring.eq.1) then

         barflag = 1

      endif
c
      if(kspring.eq.2) then

         toflag = 1

      endif
c
      epsilon = 1.0e-08
c
      isglow                    = 1
      isgup                     = nseg
c
c     Initializations
c
      DO is=1,ns

         stifsmx(is)                  = 0.0d0
         stifsmy(is)                  = 0.0d0
         stifsmz(is)                  = 0.0d0
c
      ENDDO
c  
c     Computing the predicted displacements
c     Taking into account the logical boundary 
c     They must be equal to zero for the mesh vertices 
c     placed on the farfield boundaries
c     They are known for the mesh vertices placed on the body suface
c     
      DO is=1,ns
c
CML Commented this out. It may help in case of motion. Perhaps a second order predictor
C would be even better...

         deltx(is,1)                = 0.0d0
         delty(is,1)                = 0.0d0
         deltz(is,1)                = 0.0d0

c
C Second order:
c

c        tmp = deltx(is,2)
c        deltx(is,2) = deltx(is,1)
c        deltx(is,1) = 2*deltx(is,1)-tmp
c        tmp = delty(is,2)
c        delty(is,2) = delty(is,1)
c        delty(is,1) = 2*delty(is,1)-tmp
c        tmp = deltz(is,2)
c        deltz(is,2) = deltz(is,1)
c        deltz(is,1) = 2*deltz(is,1)-tmp
c
c        Selecting the submesh vertices placed on the structure 
c        logfr =-2 : slipping vertex (Euler flow)
c        logfr =-3 : non-slipping vertex (Navier-Stokes flow) 
c
         IF (logfr(is) .LT. 0) THEN 
c
            deltx(is,1)     = xw(1,is) - coor(1,is)
            delty(is,1)     = xw(2,is) - coor(2,is)
            deltz(is,1)     = xw(3,is) - coor(3,is)
c
         ENDIF
c
         
      ENDDO
c
45    CONTINUE
c
c     Computation of D (system to solve Au = b, D = diag(A))
c     D = stifsm()
c
c     compute geometric properties for the computation of the tetra
c     stiffness matrix
c

      if (toflag.EQ.1) then

          call GetEleTetra_mov(nufaces,nucycles,lowlambda,highlambda,
     &           map_pointer,coor,nu,nt)
      endif

c
      if(barflag.EQ.1) then
c
c     Loop on local list of edges
c
      DO iseg=isglow,isgup
c
c        Local indexing of the vertices of the current edge
c
         nubo1                     = nubo(1,iseg)
         nubo2                     = nubo(2,iseg)
c 
c  -------------------->       
c        Computing the stiffness of the current edge
c
c

         xlocbar1                       = coor(1,nubo2) - coor(1,nubo1)
         xlocbar2                       = coor(2,nubo2) - coor(2,nubo1)
         xlocbar3                       = coor(3,nubo2) - coor(3,nubo1)
c
c
c
      xmag = sqrt(xlocbar1*xlocbar1+xlocbar2*xlocbar2+xlocbar3*xlocbar3)
c
      tb11 = xlocbar1/xmag
      tb21 = xlocbar2/xmag
      tb31 = xlocbar3/xmag
c
      a =  1.0d0/xmag
c
      fac1 = tb11*tb11
      fac2 = tb21*tb21
      fac3 = tb31*tb31
c
      kbar11 = a*fac1
c
      kbar22 = a*fac2
c
      kbar33 = a*fac3
c
      kbar44 = a*fac1
c
      kbar55 = a*fac2
c  
      kbar66 = a*fac3
c
c
c
c        Gathering the stiffness of the current edge into
c        the corresponding edge end-points stiffnesses
c
         stifsmx(nubo1)             = stifsmx(nubo1) + kbar11
         stifsmy(nubo1)             = stifsmy(nubo1) + kbar22
         stifsmz(nubo1)             = stifsmz(nubo1) + kbar33 
c
         stifsmx(nubo2)             = stifsmx(nubo2) + kbar44
         stifsmy(nubo2)             = stifsmy(nubo2) + kbar55
         stifsmz(nubo2)             = stifsmz(nubo2) + kbar66 
c

c
      ENDDO
c
c     end of barflag
c     |
      endif
c
      if(toflag.EQ.1) then
c
c     LOOP ON LOCAL LIST OF TETRAS
c
      do ele = 1, nt
c
c
c        Gathering the stiffness of the current tetra into
c        the corresponding tetra nodal stiffnesses
c
         stifsmx(nu(1,ele))    = stifsmx(nu(1,ele)) + 
     $                           ktet1_1(ele)
         stifsmy(nu(1,ele))    = stifsmy(nu(1,ele)) + 
     $                           ktet2_2(ele)
         stifsmz(nu(1,ele))    = stifsmz(nu(1,ele)) + 
     $                           ktet3_3(ele)
c
         stifsmx(nu(2,ele))    = stifsmx(nu(2,ele)) + 
     $                           ktet4_4(ele)
         stifsmy(nu(2,ele))    = stifsmy(nu(2,ele)) + 
     $                           ktet5_5(ele)
         stifsmz(nu(2,ele))    = stifsmz(nu(2,ele)) + 
     $                           ktet6_6(ele)
c
         stifsmx(nu(3,ele))    = stifsmx(nu(3,ele)) + 
     $                           ktet7_7(ele)
         stifsmy(nu(3,ele))    = stifsmy(nu(3,ele)) + 
     $                           ktet8_8(ele)
         stifsmz(nu(3,ele))    = stifsmz(nu(3,ele)) + 
     $                           ktet9_9(ele)
c
         stifsmx(nu(4,ele))    = stifsmx(nu(4,ele)) + 
     $                           ktet10_10(ele)
         stifsmy(nu(4,ele))    = stifsmy(nu(4,ele)) + 
     $                           ktet11_11(ele)
         stifsmz(nu(4,ele))    = stifsmz(nu(4,ele)) + 
     $                           ktet12_12(ele)
c
c     end of loop over tetras
c     |
      enddo
c
c     end of toflag
c     |
      endif
c
c
c     ---------------------------------------------------------------
c     ---------------------------------------------------------------
c
c
c     Computation of b (system to solve Au = b)
c     b = rx(),ry(),rz()
c
c
      call AmultNew_mov(deltx,delty,deltz,rx,ry,rz,1.0d0,coor,nubo,nu,
     &           logfr,mvlgfr,kspring,ns,nt,nseg)
c
c
c     Computation of the l2_norm of rx(), ry(), rz()
c
      res0PCG     = BDOT_mov(rx,ry,rz,rx,ry,rz,logfr,mvlgfr,ns)
c
      
      res0PCG     = SQRT(res0PCG)


c
c
      IF (res0PCG.LT.rsjacf) GOTO 1100
c
c     Computation of D^{-1}r, stored in zx(),zy(),zz()
c
      IF (imvmsh.EQ.2) THEN
c
        DO is=1,ns

          IF ((logfr(is).EQ.0) .OR. (mvlgfr(is).GT.0)) THEN

            zx(is) = rx(is)
            zy(is) = ry(is)
            zz(is) = rz(is)

          ENDIF

        ENDDO
c
      ENDIF
c
      IF (imvmsh.EQ.3) THEN
c
        DO is=1,ns

          IF ((logfr(is).EQ.0) .OR. (mvlgfr(is).GT.0)) THEN

            zx(is) = rx(is)/stifsmx(is)
            zy(is) = ry(is)/stifsmy(is)
            zz(is) = rz(is)/stifsmz(is)

          ENDIF

        ENDDO
c
      ENDIF
c
      
      DO is=1,ns
c
        IF ((logfr(is).EQ.0) .OR. (mvlgfr(is).GT.0)) THEN

          px(is) = zx(is)
          py(is) = zy(is)
          pz(is) = zz(is)

        ELSE

          px(is) = 0.0d0
          py(is) = 0.0d0
          pz(is) = 0.0d0

        ENDIF
c
      ENDDO
c
c     Computation of the scalar product  z().r() = zx().rx() + zy().ry()
c                                                            + zz().rz()
c
      zr = BDOT_mov(zx,zy,zz,rx,ry,rz,logfr,mvlgfr,ns)

c     Computation of Ap, stored in Apx(), Apy(), Apz()
c
      call AmultNew_mov(px,py,pz,Apx,Apy,Apz,-1.0d0,coor,nubo,nu,
     &           logfr,mvlgfr,kspring,ns,nt,nseg)

c
c     Computation of the scalar product pAp = px().Apx() + py().Apy()
c                                                        + pz().Apz()
c
      pAp = BDot_mov(px,py,pz,Apx,Apy,Apz,logfr,mvlgfr,ns)
c
c     Computation of alphaPCG
c
      alphaPCG = zr/pAp

c
c     Computation of deltx, delty, deltz
c
      DO is=1,ns
c
         IF ((logfr(is) .EQ. 0) .OR. (mvlgfr(is) .GT. 0)) THEN
c
           deltx(is,1) = deltx(is,1) + alphaPCG*px(is)
           delty(is,1) = delty(is,1) + alphaPCG*py(is)
           deltz(is,1) = deltz(is,1) + alphaPCG*pz(is)
c
c
         ENDIF         
c
      ENDDO
c
c     Computation of rx, ry, rz
c
      DO is=1,ns
c
         IF ((logfr(is) .EQ. 0) .OR. (mvlgfr(is) .GT. 0)) THEN
c
           rx(is) = rx(is) - alphaPCG*Apx(is)
           ry(is) = ry(is) - alphaPCG*Apy(is)
           rz(is) = rz(is) - alphaPCG*Apz(is)
c
         ENDIF
c
      ENDDO
c
c     Beginning of the PCG iterations (#2, #3, ...)
c     ---------------------------------------------
c
      kPCG = 1
c

1000  kPCG = kPCG + 1

c

c
c     Computation of the norm of the residual #(kPCG-1)
c
      resPCG = BDot_mov(rx,ry,rz,rx,ry,rz,logfr,mvlgfr,ns)
c
      resPCG                       = SQRT(resPCG)
c
      resPCG = resPCG/res0PCG

      if(printres.EQ.1) then

        if (kpcg.eq.2)
     .    write(6,*) 'First kPCG = ', kPCG,'  resPCG = ', resPCG
        call flush(6)
 
      endif

      WRITE(42,*) kPCG, resPCG
      CALL FLUSH(42)

      IF (resPCG.LE.rsjacf .or. kPCG.EQ.(maxjac+1)) GOTO 1100
c
c     Computation of D^{-1}r
c
      IF (imvmsh.EQ.2) THEN
c
        DO is=1,ns
c
          IF ((logfr(is).EQ.0) .OR. (mvlgfr(is).GT.0)) THEN

            zx(is) = rx(is)
            zy(is) = ry(is)
            zz(is) = rz(is)

          ENDIF
c
        ENDDO
c
      ENDIF
c
      IF (imvmsh.EQ.3) THEN
c
        DO is=1,ns
c
          IF ((logfr(is).EQ.0) .OR. (mvlgfr(is).GT.0)) THEN

            zx(is) = rx(is)/stifsmx(is)
            zy(is) = ry(is)/stifsmy(is)
            zz(is) = rz(is)/stifsmz(is)

          ENDIF
c
        ENDDO
c
      ENDIF
c
c     Computation of betaPCG
c
c
      rp = 1.0d0/zr

      zr = bdot_mov(zx,zy,zz,rx,ry,rz,logfr,mvlgfr,ns)

      betaPCG = -bdot_mov(Apx,Apy,Apz,zx,zy,zz,logfr,mvlgfr,ns)/pAp
c
c     Computation of px(), py(), pz()
c
      DO is=1,ns
c
         IF ((logfr(is) .EQ. 0) .OR. (mvlgfr(is) .GT. 0)) THEN
c
           px(is) = zx(is) + betaPCG*px(is)
           py(is) = zy(is) + betaPCG*py(is)
           pz(is) = zz(is) + betaPCG*pz(is)
c
         ENDIF        
c
      ENDDO
c
c     Computation of Ap, stored in Apx(), Apy(), Apz()
c
      call AmultNew_mov(px,py,pz,Apx,Apy,Apz,-1.0d0,coor,nubo,nu,
     &           logfr,mvlgfr,kspring,ns,nt,nseg)
c
c     Computation of the scalar product p().Ap() = px().Apx() 
c                                     + py().Apy() + pz().Apz()
      pAp = bdot_mov(px,py,pz,Apx,Apy,Apz,logfr,mvlgfr,ns)
c
c     Computation of alphaPCG
c
      IF (ABS(pAp).GT.1.e-20) THEN

        alphaPCG = zr/pAp

      ELSE

        alphaPCG = 0.0d0

      ENDIF
c
c     Computation of deltx, delty, deltz
c
      DO is=1,ns
c
         IF ((logfr(is) .EQ. 0) .OR. (mvlgfr(is) .GT. 0)) THEN
c
           deltx(is,1) = deltx(is,1) + alphaPCG*px(is)
           delty(is,1) = delty(is,1) + alphaPCG*py(is)
           deltz(is,1) = deltz(is,1) + alphaPCG*pz(is)
c
         ENDIF         
c
      ENDDO
c
c     Computation of rx, ry, rz
c
      DO is=1,ns
c
         IF ((logfr(is) .EQ. 0) .OR. (mvlgfr(is) .GT. 0)) THEN
c
           rx(is) = rx(is) - alphaPCG*Apx(is)
           ry(is) = ry(is) - alphaPCG*Apy(is)
           rz(is) = rz(is) - alphaPCG*Apz(is)
c
         ENDIF
c
      ENDDO

      rp = BDot_mov(rx,ry,rz,rx,ry,rz,logfr,mvlgfr,ns)

c
      GOTO 1000
c

1100  CONTINUE

c
c
c
c
c     Updating the positions of the interior submesh vertices
c
      DO 2000 is=1,ns
c
c         IF ((logfr(is) .EQ. 0) .OR. (mvlgfr(is) .GT. 0)) THEN 
         IF ((logfr(is) .GE. 0) .OR. (mvlgfr(is) .GT. 0)) THEN 
c
c           Projection of the updated displacements for those 
c           vertices that are constrained
c
c            IF (mvlgfr(is) .GT. 0) 
c     &        CALL PRJDEP(is, deltx(is,1), delty(is,1), deltz(is,1))
c
            xw(1,is)               = coor(1,is) + deltx(is,1) 
            xw(2,is)               = coor(2,is) + delty(is,1) 
            xw(3,is)               = coor(3,is) + deltz(is,1) 
c
         ENDIF
c
2000  CONTINUE
c

      if(printres.EQ.1) then
        if (kpcg.eq.2)
     .    write(6,*) 'Last kPCG = ', kPCG,'  resPCG = ', resPCG
        call flush(6)
      endif

      
      RETURN
      END







