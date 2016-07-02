      subroutine smlesgr
c
c     This subroutine performs a least square smoothing of the gradient
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

      real*8 sh,grax
      integer jtp,isp,inod,jn1,jn2,jn3

      

c...  Initialization 

      do isp = 1,nsp
        delta(isp)      = 0.d0
      end do
      
c...  Evaluate R = M * W (stored in vno)

      sh  = 1.d0/3.d0
      
      do inod=1,3
        do jtp = 1,ntp            
          isp=nup(inod,jtp)
          jn1=nup(1,jtp)
          jn2=nup(2,jtp)
          jn3=nup(3,jtp)
          
          grax= sh*(grad(jn1)+grad(jn2)+grad(jn3))
          
          delta(isp) =
     .      delta(isp)+sh*grax*airtp(jtp)
          
        end do
      end do
      
c...  Evaluate S = R * L^-1

      do isp=1,nsp
        grad(isp) = delta(isp) / airesp(isp) 
      end do


      end
