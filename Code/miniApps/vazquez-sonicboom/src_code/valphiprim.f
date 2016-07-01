




      SUBROUTINE VALPHIPRIM(CTRL)
C
c  Cette procedure verifie le calcul de la matrice implicite 
c  de la facon suivante :
c
c     Phi(W + eps dW) - Phi(W)
c     ------------------------ - Phi'(W).dW = O(eps)
c               eps
c
c                Phi(W) . dt
c   avec : dw = -------------
c                  volume
c
c   On se place en calcul explicite. Cette procedure est appelee dans Etat.f
c
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      REAL*8 PSI1(NEQUATION,nsmax)
      REAL*8 DELTAW(NEQUATION,nsmax), DIFF(NEQUATION,nsmax)
      REAL*8 DIAGO(nsmax,NEQUATION,NEQUATION)
      REAL*8 TMP(NEQUATION,nsmax), CTRL(NNSP)
      REAL*8 DMAX, NORME, epsilon, ZM1(NEQUATION,NEQUATION,2,nsgmax)
C
      INTEGER ISEG, I, J, IS, K, isa, isb, ii
C
c      COMMON/TABLOC/TMP, PSI1
c

      write (6,*) 'STOPPED IN VALPHIPRIM:'
      write (6,*) 'CORRECT ctrlno AS DONE IN ETAT SUBROUTINE!!'

      stop
      
      epsilon = 1.d-05
c
c   Calcul de Phi(W) -->  stocke dans PSI1
c   Calcul de dW     -->  stocke dans DELTAW
c
               CALL FLUX(NSMAX, NTMAX, NSGMAX, NORDRE, NU, NUBO,
     $               DX, DY, DZ, VNOCL, COOR, UA, GAM, CE)
c
               if (bound.eq.1) then
                  CALL VCURVM
         IF (itrans.eq.1) CALL TRANSPIRATION(CE,CTRL)
               endif
c
               if (diric.eq.1) then
                  CALL VCURVM
         IF (itrans.eq.1) CALL TRANSPIRATION(CE,CTRL)
                  CALL CONDDIRFLUX
               endif
c
               DO IS = 1,NSMAX
                  DO K = 1,NEQUATION
                     PSI1(K,IS) = CE(K,IS)
c                     DELTAW(K,IS) = DTL(IS)*PSI1(K,IS)/VOLS(IS)
                    DELTAW(K,IS) = 1.
                     DO J=1,NEQUATION
                        DIAG(IS,K,J) = 0.
                     END DO
                  END DO
               END DO
C
c   Calcul de Phi'(W)  -->  stocke dans ZM et DIAG
c
               DO iseg = 1,nsgmax
                  Do i = 1,NEQUATION
                     Do j =1,NEQUATION
                        ZM1(i,j,1,iseg) = 0.
                        ZM1(i,j,2,iseg) = 0.
                     End Do
                  End Do
               End Do
c            
             CALL MATVL(NEQUATION, NSMAX, NSGMAX,
     &                             NUBO, VNOCL, GAM, UA,
     &                                DIAG,ZM1,ns,nseg) 
c
         if (bound.eq.1) then 
         CALL CONDBORDS
         IF (itrans.eq.1) CALL IMPTRANSPIRATION(CTRL)
c
         endif
c
         if (diric.eq.1) then
         CALL CONDBORDS
         IF (itrans.eq.1) CALL IMPTRANSPIRATION(CTRL)
         CALL CDMAT(NEQUATION,NSMAX,NSGMAX,NUBO,LOGFR,ZM1,DIAG,ns,nseg)
         endif
c
c   Calcul de Phi'(W).dW  -->  stocke dans tmp
c
               do is = 1,nsmax
                  do i =1,5
                     tmp(i,is) = 0.
                  end do
               end do
c
c     Termes extra-diagonaux
c
         do iseg = 1,nsgmax
            isa = nubo(1,iseg)
            isb = nubo(2,iseg)
               do k = 1,5
                  do j = 1,5
                     tmp(k,isa) = tmp(k,isa) +
     $                    zm1(k,j,2,iseg)*deltaW(j,isb)
                     tmp(k,isb) = tmp(k,isb) +
     $                    zm1(k,j,1,iseg)*deltaW(j,isa)
                  end do
               end do
         end do
c
c     Termes diagonaux
c
         do is = 1,nsmax
               do i = 1,5
                  do j = 1,5
               tmp(i,is) = tmp(i,is) + diag(is,i,j)*deltaW(j,is)
                  end do
            end do
         end do
c
c   Calcul de W + eps dW  -->  stocke dans ua
c
      DO 60 i=1,5
         DO 65 is=1,nsmax
c            ua(i,is) = ua(i,is) + epsilon*ce(i,is)*dtl(is)
c     $                    /vols(is)
            ua(i,is) = ua(i,is) + epsilon*deltaw(i,is)
65       CONTINUE
60    CONTINUE
c
c   Calcul de Phi(W + eps dW)  -->  stocke dans ce
c
               CALL FLUX(NSMAX, NTMAX, NSGMAX, NORDRE, NU, NUBO,
     $               DX, DY, DZ, VNOCL, COOR, UA, GAM, CE)
c
               if (bound.eq.1) then
                  CALL VCURVM
         IF (itrans.eq.1) CALL TRANSPIRATION(CE,CTRL)
               endif
c
               if (diric.eq.1) then
                  CALL VCURVM
         IF (itrans.eq.1) CALL TRANSPIRATION(CE,CTRL)
                  CALL CONDDIRFLUX
               endif
c
c   Calcul de (Phi(W+epsdW)-Phi(W))/eps
c
c   Attention : dans PSI1 est stocke en fait -Phi(W)
c               dans CE   "    "       "     -Phi(W+epsdW)
c
               DO IS = 1,Nsmax
                  DO K = 1,5
                     DIFF(K,IS) = (PSI1(K,IS) - CE(K,IS))/epsilon
                   END DO
               END DO
C
c
         dmax = -1.
         norme = 0.
         do is = 1,nsmax
            if (logfr(is).ne.1) then
               do i = 1,5
                  dmax = max(dmax,abs(tmp(i,is)-diff(i,is)))
                  norme = norme + (tmp(i,is)-diff(i,is))**2
               write(100+i,222) is, logfr(is), tmp(i,is), diff(i,is),
     $       coor(1,is),coor(2,is),coor(3,is),abs(tmp(i,is)-diff(i,is))
               end do
            endif
         end do
         print*,'epsilon = ',epsilon
         print*,'MAX = ',dmax
         print*,'norme =',sqrt(norme)
c
222      FORMAT(2i5,1x,6(1x,e14.7))
c

               do is = 1,ns
                  do i = 1,5
                     do j = 1,5
                        diag(is,i,j) = 0.
                     end do
                  end do
               end do
c
               DO iseg = 1,nsgmax
                  Do i = 1,NEQUATION
                     Do j =1,NEQUATION
                        ZM1(i,j,1,iseg) = 0.
                        ZM1(i,j,2,iseg) = 0.
                     End Do
                  End Do
               End Do
c
         RETURN
         END
