
      SUBROUTINE ETATADJOINT(CTRL)
C                                           ___
C*** Cette procedure calcule l'etat-adjoint ||  stocke dans Piadj(1:5,1:ns)
C    par resolution de l'equation lineaire :
C                                         
C      ( dPsi )*         ___    dJ
C      ( ---- )          ||  =  --  avec la methode de Jacobi.
C      (  dW  )vanleer1         dW
C
C*** Appel de la procedure ADJOINTE qui calcule la matrice adjointe stockee
C    dans ZM et DIAG. 
C    Appel de la procedure Jacobi qui prend en entrees ZM, DIAG, dJdW et
C    qui ressort l'etat-adjoint Piadj.
C
C
C
C    DO n=1,...,p
C                                         
C      ( dPsi )*              ___ n          dJ       ( dPsi )*     ___ n-1
C      ( ---- )       ( Delta || )        =  --   -   ( ---- )      || 
C      (  dW  )roe1                          dW       (  dW  )roe2
C
C                                                    |                      |
C                                                     ----------------------
C                                                              ceb  <<<---  mode inverse
C
C                                           |                               |
C                                            -------------------------------
C                                                         rhsadj
C
C
C      Pi^n  =  Pi^n-1 + Delta Pi^n
c
c
C    END DO
C      
            
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
      INCLUDE 'DiffCall_flow_d.h'
c      INCLUDE 'DiffCall_gamma_b.h'
      INCLUDE 'DiffCall_flow_b.h'

      REAL*8 RES, err1,erritor,ctrlno,sumrhs,sumceb,sumdel
      REAL*8 CTRL(NNSP),ctrld(nnsp),ctrlb(nnsp)
      REAL*8 PRESSION,coorpbar(3,nnsp),vncoqbar(3,nnsp)
      INTEGER ITER, IETAT, nitor, iitor,is,j,isp

      real*8 delpi(5,nsmax),rhsadj(5,nsmax),piaux(5,nsmax)
      real*4 ctjac
      
      
C [llh Nan]      call second(ct1)

      write(6,*) 'Entree dans EtatAdjoint.'      

      nitor  = 1
      erritor= 1.0d-5

      do j=1,5
        do is=1,ns
          delpi(j,is)  =  piadj(j,is)
          rhsadj(j,is) =  djdw(j,is)
        end do          
      end do
c
c...  For roe flux, assign
c
c           zm(5,5,2,nsgmax) <--- stmat(nsegs)
c
c     remember nsegs=50*nsgmax!
c      
      if (iflux.eq.2) then
        ctrlno=0.d0
        do isp=1,nsp
          ctrlno= ctrlno+ctrl(isp)*ctrl(isp)
        end do
        do j=1,5
          do is=1,ns
            delpi(j,is)  =  0.d0
            piaux(j,is)  =  piadj(j,is)
          end do          
        end do        
        nitor= 10
c        call stmat2zm(stmat,zm,5,5,2,nsgmax,nseg,nsegs)        
      end if

c
c...  Evaluate the adjoint matrix
c
      
      call adjointe
      
C [llh Nan]      call second(ct3)

c
c...  Loop to compute Delta Pi (roe) or directly compute Pi (van leer) 
c

      do iitor=1,nitor

        
c...  Jacobi cycle 

        IETAT=0
        if (iflux.eq.1) then
          CALL JACOBI(
     .      NEQUATION,nsmax,nsgmax,nubo,nbrelPi,   zm,diag,rhsadj,
     .      errJACPI,delpi,res,iter,IETAT,ctjac,ns,nseg)

        else if (iflux.eq.2) then
          CALL JACOBI(
     .      NEQUATION,nsmax,nsgmax,nubo,nbrelPi,stmat,diag,rhsadj,
     .      errJACPI,delpi,res,iter,IETAT,ctjac,ns,nseg)

          sumdel=0.d0
          sumrhs=0.d0
          sumceb=0.d0
          do j=1,5
            do is=1,ns
              sumrhs= sumrhs+rhsadj(j,is)*rhsadj(j,is)
              sumdel= sumdel+delpi(j,is)*delpi(j,is)
              sumceb= sumceb+ceb(j,is)*ceb(j,is)
              piadj(j,is) = piadj(j,is)-delpi(j,is)
              ceb(j,is)   = piadj(j,is)
              uab(j,is)   = 0.d0
              delpi(j,is) = 0.d0
            end do          
          end do                    
          
          write(6,*) 'ctrlno=  ',ctrlno,'  sumceb=  ',sumceb,
     .      '  sumrhs=  ',sumrhs
          call psiroe_b(ctrl,ctrlno)
        
          sumceb=0.d0
          sumrhs=0.d0
          do j=1,5
            do is=1,ns
              rhsadj(j,is) = djdw(j,is) - ceb(j,is)
              sumceb= sumceb+ceb(j,is)*ceb(j,is)
              sumrhs= sumrhs+rhsadj(j,is)*rhsadj(j,is)
            end do          
          end do          

          write(6,*) 'ctrlno=  ',ctrlno,'  sumceb=  ',sumceb
          if (iitor.eq.1)
     .      write(6,*) 'Solving 2nd. order adjoint state...'
          write(6,*)
     .      '  Jacobi Residual=  ',res,
     .      '  Delta Pi norm  =  ',sumdel,
     .      '  Outer iteration= ',iitor
          
        end if


      end do

        
C [llh Nan]      call second(ct4)
C [llh Nan]      ctjad= ctjad+ct4-ct3
      
C
      WRITE(6,*) '  '
      WRITE(6,*) 'Dans la resolution de l"etat-adjoint :'
      WRITE(6,*) '--------------------------------------'
      WRITE(6,*) '  '
      WRITE(6,*) 'Le residu vaut ',res,'  pour ',iter,' relaxations.'
      WRITE(6,*) '  '
c
              open(19,access='append')
              write(19,180)res,iter,itopt
              close(19)
c
180   FORMAT(e14.7,1x,2i5)
C

C [llh Nan]      call second(ct2)
C [llh Nan]      ctato= ctato+ct2-ct1
      
      write(6,*) 'Sortie de EtatAdjoint.'

      END
