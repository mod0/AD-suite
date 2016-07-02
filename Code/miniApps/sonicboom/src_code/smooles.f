      subroutine smooles
c
c     This subroutine performs a least square smoothing of the normals
c
c     W --> element defined normals, discontinued function. It is the "wrinkled" 
c           normal to be smoothed
c
c     S --> smoothed normal, defined on the skin points 
c
c     M --> mass matrix: Int_Omega ( N^i N^j dOmega ), N^i is the i-th node shape function
c
c
      include 'Param3D.h'
      include 'Paramopt3D.h'

      real*8 sh,vnox,vnono,coorpbar(3,nnsp)
      integer jtp, isp, idim,inod,izztask,ite,jn1,jn2,jn3


C    I HAVE TO CHECK THIS SUB...!!!!
      

c...  Initialization 

      do idim=1,3
        do isp = 1,nsp
          vnon(idim,isp)     = 0.d0
          vno(idim,isp)      = 0.d0
          coorpbar(idim,isp) = 0.d0
        end do
      end do
      
c...  Evaluate R (stored in coorpbar)

      sh  = 1.d0/3.d0

      do idim=1,3
        do inod=1,3
          do jtp = 1,ntp            
            isp=nup(inod,jtp)
            jn1=nup(1,jtp)
            jn2=nup(2,jtp)
            jn3=nup(3,jtp)
            vnox= sh*(vno(idim,jn1)+vno(idim,jn2)+vno(idim,jn3))

            coorpbar(idim,isp) =
     .        coorpbar(idim,isp)+sh*vnox*airtp(jtp)

          end do
        end do
      end do

c...  Evaluate S^0 (stored in vno)

      do isp=1,nsp
        do idim=1,3
          vno(idim,isp) = coorpbar(idim,isp) / airesp(isp) 
        end do
      end do


c...  Iterative Jacobi process
      
c      do ite=1,10
c        
cc...  Evaluate M S^i (stored in vnon)
c
c        do idim=1,3
c          do inod=1,3
c            do jtp = 1,ntp            
c              isp=nup(inod,jtp)
c              jn1=nup(1,jtp)
c              jn2=nup(2,jtp)
c              jn3=nup(3,jtp)
c              vnox= sh*(vno(idim,jn1)+vno(idim,jn2)+vno(idim,jn3))
c              vnon(idim,isp) = vnon(idim,isp)+sh*vnox*airtp(jtp)
c            end do
c          end do
c        end do
c
cc...  Evaluate S^i+1 = S^i + L^{-1} ( R - M S^i )    (stored in vno)
c
c        do isp=1,nsp
c          do idim=1,3
c            vno(idim,isp) = vno(idim,isp)
c     .        + (coorpbar(idim,isp) - vnon(idim,isp)) / airesp(isp) 
c            vnon(idim,isp)=0.d0
c          end do
c        end do
c        
c        
c      end do


      
c...  Normalize
      
      do isp=1,nsp
        vnono= 0.d0
        do idim=1,3
          vnono         = vnono + vno(idim,isp)*vno(idim,isp)
        end do
        vnono= sqrt(vnono)
        do idim=1,3
          vnon(idim,isp) = vno(idim,isp) / vnono
        end do        
      end do

      

      end
