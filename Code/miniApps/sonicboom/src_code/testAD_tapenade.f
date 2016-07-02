    

      SUBROUTINE ETAT(CTRL,tempsmax)
c---------------------------------------------------------------------   
c Tridimensional Euler and Navier-Stokes equations solver adapted 
c to dynamic meshes
c Explicit multi-steps or implicit time integration process
c Upwind schemes and linear interpolation method for the computation 
c of the convective fluxes using a finite volume formulation
c Classical central Galerkin P1-finite element method for the 
c computation of the diffusive fluxes
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
      INCLUDE 'DiffCall_flow_d.h'
c      INCLUDE 'DiffCall_gamma_b.h'
      INCLUDE 'DiffCall_flow_b.h'
c---------------------------------------------------------------------
c     Local variables definition
c
      INTEGER ivar , i, is , ia, ib,isp, tempsmax, iprint
      REAL*8    aux
      REAL*8    pi
      REAL*8    tetaA, domega, cl
      REAL*8    tetaTd,ctrlno
      REAL*8    ceeps(5, nsmax), diffdiv
c
      INTEGER j,iseg, niter, IETAT, jj
      REAL*8 CTRL(NNSP),ctrld(nnsp),ctrlb(nnsp)
      REAL*8 residu,vnorme(3), vnormeprec(3), vnormestock(3)
C
      INTEGER NCOMPT,kvar
      REAL*8 PRESSION,coorpbar(3,nnsp),vncoqbar(3,nnsp)

      real*8 ceSave(5,12966),uaSave(5,nsmax),cedSave(5,nsmax)
      real*8 uadSave(5,nsmax), ctrlSave(nnsp)
      real*8 scalTest1,scalTest2,scalTest3,tmpTest
c llh: needed for derivatives by DD:
      real*8 ddeps, epszero
      integer phase, ddfile
      common /comphase/ phase, ddfile, ddeps, epszero
      real*8 tt0, cpuOrig, cpuTgt, cpuAdj
      REAL*4 ctjac
C
c     Initialisations
c

      ctrlno=0.d0
      do isp=1,nsp
        ctrlno= ctrlno+ctrl(isp)*ctrl(isp)
      end do

      ifre                         = 1000
c
      pi                           = 4.0*ATAN(1.0)
c
      tetaA                        = 2.0
      domega                        = 4.0
c
      dt                           = 0.1
      kt                           = kt0
      t                            = t0
c
      som                          = 1.0
c
c     Initialisation of the timers
c
      DO 20 i=1,icpmax
         time(i)                   = 0.0
20    CONTINUE
c
c      REWIND(17)
      REWIND(19)
c
c     Boucle non lineaire en temps
c     ============================
c
1000  CONTINUE
c
      kt                           = kt + 1
c
c     Computing the cfl
c
      aux                          = MAX(ycfl, zcfl*(kt - kt0))
c
      cfl                          = MIN(MAX(aux, xcfl/som), cflmax) 
c
c     Initialisation des termes de la matrice implicite.
c     Van Leer : zm(5,5,2,nseg) = termes extradiagonaux
c     Roe : stmat(5*5*2*nseg) = termes extradiagonaux
c     Van Leer et Roe : diag(5,5,ns) = termes diagonaux
c
C$DOACROSS LOCAL(ISEG,I,J),SHARE(ZM)
      DO iseg = 1,nseg
         Do i = 1,5
            Do j =1,5
               ZM(i,j,1,iseg) = 0.
               ZM(i,j,2,iseg) = 0.
            End Do
         End Do
      End Do
c

C$DOACROSS LOCAL(IS,IA,IB),SHARE(DIAG)
      DO  is=1,ns
         DO  ia=1,5
            DO  ib=1,5
               diag(is,ia,ib)      = 0.0
            END DO
         END DO
      END DO
c
      if (iflux.gt.1) then
C$DOACROSS LOCAL(ISEG)
        DO iseg=1,nsegs
          stmat(iseg)               = 0.0
        END DO
      end if

c
c     Computing psi
c      
      IF (iflux .EQ. 1) THEN                                !van leer
        
C$DOACROSS LOCAL(IS,IA),SHARE(DX,DY,DZ)
        DO  is=1,ns
          DO  ia=1,5
            dx(ia,is)            = 0.0
            dy(ia,is)            = 0.0
            dz(ia,is)            = 0.0
          END DO
        END DO
c
c     Initialisation of the nodal fluxes 
c
C$DOACROSS LOCAL(IS),SHARE(CE)
        DO is=1,ns
          ce(1,is)                  = 0.0
          ce(2,is)                  = 0.0
          ce(3,is)                  = 0.0
          ce(4,is)                  = 0.0
          ce(5,is)                  = 0.0
        END DO

        CALL GRADNOD

        CALL FLUX( NS, NT, NSEG, NORDRE, NU, NUBO,
     .    DX, DY, DZ, VNOCL, COOR, UA, GAM, CE)
        
        CALL VCURVM(ctrlno)
        
        IF (itrans.eq.1 .and. ctrlno.gt.1.0d-10)
     .    CALL TRANSPIRATION(CE,CTRL)        
        CALL CONDDIRFLUX        
        
        
      ELSE if (iflux.eq.2) then                             !roe

C**************************************************************
C This is the call of the function you want to differentiate,
C Suggestion: differentiate it for iteration kt==5.
C**************************************************************
C      VALIDATION BETWEEN DIVIDED DIFFERENCES, FORWARD MODE AD, AND REVERSE MODE AD
C      To activate at time step N, put (kt.eq.N) :
        if (kt.eq.5) then
          
          write (6,*) 'VALIDATION TESTS FOR vazquez-sonicboom:'
C Snapshot the initial (current) state:
          do isp=1,nsp
            ctrlSave(isp) = ctrl(isp)
          end do
          do is=1,ns
             do jj=1,5
                uaSave(jj,is) = ua(jj,is)
                ceSave(jj,is) = ce(jj,is)
             enddo
          enddo
C Snapshot the input forward derivatives:
          do is=1,ns
             do jj=1,5
                uadSave(jj,is)= 1.0
                cedSave(jj,is)= 1.0
             enddo
          enddo
C Initialise the input forward derivatives :
          do is=1,ns
             do jj=1,5
                uad(jj,is)=uadSave(jj,is)
                ced(jj,is)=cedSave(jj,is)
             enddo
          enddo
C Divided differences: build initial state, with epsilon modification, from the snapshot:
          ddeps = 1.d-6
          do isp=1,nsp
            ctrl(isp) = ctrlSave(isp)
          enddo
          do is=1,ns
             do jj=1,5
                ua(jj,is) = uaSave(jj,is) + ddeps*uadSave(jj,is)
             enddo
          enddo
C Divided differences: run the original code at point x+epsilon
          call cputim(tt0)
          call psiroe(ctrl,ctrlno)
          call cputim(cpuOrig)
          cpuOrig = cpuOrig-tt0
          write (6,*) 'Time of  original  function: ',cpuOrig
C Divided differences: save the result, i.e. f(x+eps)
          do is=1,ns
             do jj=1,5
                ceeps(jj,is) = ce(jj,is)
             enddo
          enddo
C Divided differences: rebuild the initial state,
C  without the epsilon modification, from the snapshot:
          do isp=1,nsp
            ctrl(isp) = ctrlSave(isp)
          enddo
          do is=1,ns
             do jj=1,5
                ua(jj,is) = uaSave(jj,is)
                ce(jj,is) = ceSave(jj,is)
             enddo
          enddo
C Divided differences: Re-initialise the input forward derivatives :
          do is=1,ns
             do jj=1,5
                uad(jj,is)=uadSave(jj,is)
                ced(jj,is)=cedSave(jj,is)
             enddo
          enddo
C Divided differences: run the fwd AD code at point x
          call cputim(tt0)
          call psiroe_d(ctrl,ctrlno)
          call cputim(cpuTgt)
          cpuTgt = cpuTgt-tt0
          write (6,*) 'Time of tangent AD function: ',cpuTgt
C Divided differences: compute the square norm of f(x+eps)-f(x)/eps:
          scalTest1=0.0
          do is=1,ns
             do jj=1,5
                diffdiv = (ceeps(jj,is) - ce(jj,is))/ddeps
                scalTest1=scalTest1+diffdiv*diffdiv
             enddo
          enddo
C Forward AD mode: compute the square norm of \dot{f}(x), which should be the same:
          scalTest2=0.0
          do is=1,ns
            do jj=1,5
              scalTest2=scalTest2+ced(jj,is)*ced(jj,is)
            end do
          end do
C Reverse AD mode: rebuild the initial state:
          do isp=1,nsp
            ctrl(isp) = ctrlSave(isp)
          enddo
          do is=1,ns
             do jj=1,5
                ua(jj,is) = uaSave(jj,is)
                ce(jj,is) = ceSave(jj,is)
             enddo
          enddo
C Reverse AD mode used for the dot-product test:
C  initialize the input reverse derivatives with the output fwd derivatives:
          do is=1,ns
             do jj=1,5
                uab(jj,is) = 0.0
                ceb(jj,is) = ced(jj,is)
             enddo
          enddo
C Reverse AD mode: run the reverse AD code at point x and print time and stack:
          call cputim(tt0)
          call psiroe_b(ctrl,ctrlno)
          call cputim(cpuAdj)
          cpuAdj = cpuAdj-tt0
          write (6,*) 'Time of adjoint AD function: ',cpuAdj

          if (cpuOrig.ne.0.0)
     +        write (6,1003) 
     +            ' Slowdown factors: tangent',cpuTgt/cpuOrig,
     +            '    adjoint',cpuAdj/cpuOrig
          call printstackmax()
          call PRINTTRAFFIC()

C Reverse AD mode: dot-product test: display <xb|xd>, should be the same
C  as the above <yd|yd>.
          scalTest3=0.0
          do is=1,ns
            do jj=1,5
              scalTest3=scalTest3+uab(jj,is)*uadSave(jj,is)
            end do
          end do
          write (6,*) ''
          write (6,1001) 'Divided differences   = ',scalTest1,
     +         ' (epsilon=',ddeps,')'
          write (6,1002) 'AD Forward derivative = ',scalTest2 
          write (6,1002) 'AD Reverse <xb|xd>    = ',scalTest3
C The resolution is now broken by the tests above: we can't go on...
          stop
       endif

       call psiroe(ctrl,ctrlno)

C**************************************************************

1001  format(a,e26.20,a,e8.2,a)
1002  format(a,e26.20)
1003  format(a,f5.2,a,f5.2)
        
      ENDIF

c     Computation of the local time steps
c     Computation and gathering of the viscous fluxes
c
      CALL FLUXDT


      t                            = t + dt

c     Computation of the instantaneous angle of attack
c
      tetaTd                       = tetaA*SIN(domega*t)
      tetaT                        = tetaTd*pi/180.0

c
      IF (ivis .EQ. 1) CALL CDLNS
c
c     Computing the residual
c
      som                          = 0.0
      dro                          = 0.0
c
      DO 55 is=1,ns
         som                       = som + ce(5,is)*ce(5,is)
         dro                       = MAX(dro, ABS(ce(5,is)))
55    CONTINUE
c
      som                          = SQRT(som)
c
c     Normalisation du residu :
c     -------------------------
c
      IF (kt .EQ. 1 .or. faczz(19).lt.-1.0d0) THEN 
	if (faczz(19).lt.-1.0d0) then
           faczz(19)=som 
           faczz(18)=dro 
        end if
c         som1                      = som
c         dro1                      = dro

         som1                      = faczz(19)
         dro1                      = faczz(18)
c
c      write(57,110) som1, dro1  ! En cas de continuation du code
c
c110   FORMAT(2e14.7)
      ENDIF
c
      som                          = som/som1
      dro                          = dro/dro1
c
 1918 FORMAT(10x,'Iteration en temps:',I4,
     $' residu non lineaire:',e13.6)
c
      call flunow(6)                                        !dummy in interactive but very useful
                                                            !running in batch...
c      open(16,access='append')
c      WRITE(16, 90) kt, log10(som), log10(dro), itopt
c      close(16)
c
90    FORMAT(i5,2x,2(e12.5,2x),i5)
c
      IF (nexp .EQ. 0) GOTO 1250
c
c     Explicit time integration
c
c     Updating the physical solution
c

C$DOACROSS LOCAL(IS,IVAR),SHARE(UA,CE,DTL,VOLS,UN)
      DO is=1,ns
         DO ivar=1,5
            un(ivar,is) = ua(ivar,is) + ce(ivar,is)*dtl(is)
     $                                 /vols(is)
         END DO
      END DO
c
      GOTO 1500
c
1250  CONTINUE
c
c     Implicit time integration
c
c     Computes the Jacobian matrix 
c

      IF (iflux .EQ. 1) THEN 
         CALL MATVL(NEQUATION, NSMAX, NSGMAX, NUBO, VNOCL, GAM, 
     &                        UA,DIAG,ZM,ns,nseg)
      ENDIF
c
      IF (iflux .EQ. 2) THEN 
         CALL IMPROE
      ENDIF
c

      CALL CONDBORDS(ctrlno)
      IF (itrans.eq.1 .and. ctrlno.gt.1.0d-10)
     .  CALL IMPTRANSPIRATION(CTRL)
      if (iflux.eq.1) then
        CALL CDMAT(NEQUATION,NSMAX,NSGMAX,NUBO,LOGFR,ZM,DIAG,ns,nseg)
      else
        call CDMAT(NEQUATION,NSMAX,NSGMAX,NUBO,LOGFR,stmat,DIAG,ns,nseg)
      endif

c      if (bound.eq.1) then
c         CALL CONDBORDS(ctrlno)
c         IF (itrans.eq.1 .and. ctrlno.gt.1.0d-10)
c     .     CALL IMPTRANSPIRATION(CTRL)
c      endif
c      if (diric.eq.1) then
c         CALL CONDBORDS(ctrlno)
c         IF (itrans.eq.1 .and. ctrlno.gt.1.0d-10)
c     .     CALL IMPTRANSPIRATION(CTRL)
c         if (iflux.eq.1) then
c       CALL CDMAT(NEQUATION,NSMAX,NSGMAX,NUBO,LOGFR,ZM,DIAG,ns,nseg)
c         else
c       call CDMAT(NEQUATION,NSMAX,NSGMAX,NUBO,LOGFR,stmat,DIAG,ns,nseg)
c         endif
c      endif

c
      if (Newton.eq.0) CALL DIAGAJOUT
C
c     Inversion des termes diagonaux
c
      CALL INVERSION
c
c     Solving the implicit linear system using the Jacobi method
c
C$DOACROSS LOCAL(IS,IA)
      DO  is=1,ns
         DO  ia=1,5
            dx(ia,is)            = 0.0
         END DO
      END DO
c
      IETAT = 1
c
      IF (iflux.eq.1) THEN
          IF (irlax .EQ. 1) THEN 
             IF (ONESHOT.EQ.1) THEN
                IF (ITOPT.EQ.1) THEN
                   CALL JACOBI(NEQUATION,nsmax,nsgmax,nubo,nbrel,
     $     zm,diag,ce,err,dx,residu, niter, IETAT,ctjac,ns,nseg)
                ELSE
                   CALL JACOBI(NEQUATION,nsmax,nsgmax,nubo,nbrelos,
     $     zm,diag,ce,err,dx,residu, niter, IETAT,ctjac,ns,nseg)
                ENDIF
             ELSE
                CALL JACOBI(NEQUATION,nsmax,nsgmax,nubo,nbrel,
     $     zm,diag,ce,err,dx,residu, niter, IETAT,ctjac,ns,nseg)
             ENDIF
c
c       WRITE(6,1492)residu,niter
c
c1492   FORMAT('          Jacobi: residu=',
c     $        e10.4,' nbre de relax.=',i5)
c
c              open(18,access='append')
c              write(18,180)residu, niter, kt, itopt
c              close(18)
c
180   FORMAT(e14.7,1x,3(i5,1x))
c
          ENDIF
      ENDIF
c
      IF (iflux.eq.2) THEN
         IF (irlax .EQ. 1) THEN
            CALL JACOBIROE
         ENDIF
      ENDIF

c
c     Updating the physical solution
c
      DO is=1,ns
         DO ivar=1,5
            if (Newton.eq.1) then 
               un(ivar,is) = ua(ivar,is) + 0.5*dx(ivar,is)
            else
               un(ivar,is) = ua(ivar,is) + dx(ivar,is)
            endif
         END DO
         call clippin(is,un(1,is),ua(1,is))
      END DO
c
1500  CONTINUE

c
c     Swapping the new and old physical states
c
C$DOACROSS LOCAL(IS,IVAR),SHARE(UA,UN)
      DO is=1,ns
         DO ivar=1,5
            ua(ivar,is) = un(ivar,is)
         END DO
      END DO
c
      IF (MOD(kt, ifre)  .EQ. 0) THEN
c
c        Computing the lift coefficient
c
         CALL AEROF1(cl)
c
         WRITE(19, '(i8,3f14.6)') kt, t, tetaTd, cl
c
      ENDIF
c
      if (Newton.eq.1) then
c
      Do i = 1,3
         vnormeprec(i) = 1.
         vnorme(i) = 0.
      End Do
c
C$DOACROSS LOCAL(IS), REDUCTION(VNORME(1),VNORME(2),VNORME(3)),
C$&        SHARE(CE)
      Do is = 1,ns
        if (logfr(is).ne.1 .and. logfr(is).ne.5) then
          vnorme(1) = vnorme(1) + ce(1,is)*ce(1,is)
          vnorme(2) = vnorme(2) + ce(2,is)*ce(2,is)
          vnorme(3) = vnorme(3) + ce(5,is)*ce(5,is)
        endif
      End Do
c
            Do i = 1,3
               vnorme(i) = sqrt(vnorme(i))
               vnormestock(i) = vnorme(i)/vnormeprec(i)
               vnormeprec(i) = vnorme(i)
            End Do
c
      open(34,access='append')
      write(34,210) vnormestock(1),vnormestock(2),vnormestock(3),kt,cfl
      close(34)
c
210   FORMAT(3(e14.7,1x),1x,i5,1x,e14.7)
c
      endif
c
      IF ((kt .LT. tempsmax) .AND. 
     &    (ABS(t - tmax) .GT. 1.0e-06) .AND.
     &    (som .GT. resf)) GOTO 1000
c
      IPRINT=1
      IF(IPRINT.EQ.1)CALL RESU3D
c
      IF (coefm1.eq.1) THEN  ! Cas du tube a choc
c
      open(31)
      open(32)
      open(33)
c
      Do is=1,ns
         if ((abs(coor(2,is)-0.5).lt.1e-04).and.(abs(coor(3,is)-0.5)
     &                   .lt.1e-04)) then
         WRITE(31,116) coor(1,is), ua(1,is)
         WRITE(32,116) coor(1,is), ua(2,is)/ua(1,is)
         WRITE(33,116) coor(1,is), ua(3,is)/ua(1,is)
         endif
      End do 
c
c   Remise a 0 de kt0
c 
      kt0 =0
c
      CLOSE(31)
      CLOSE(32)
      CLOSE(33)
c
      ENDIF
c
116   FORMAT(2e15.8)
c
c
3000  FORMAT(a38,f12.6)
c
c      CLOSE(17)
      CLOSE(19)

      RETURN
      END
