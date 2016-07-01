      subroutine smlesdi(ctrl)
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

      real*8 sh,coorx,vnono,coorpbar(3,nnsp),ctrl(nnsp)
      integer jtp, isp, is,idim,inod,izztask,ite,jn1,jn2,jn3

      

c...  Initialization 

      do idim=1,3
        do isp = 1,nsp
          vno(idim,isp)      = 0.d0
          coorpbar(idim,isp) = ctrl(isp)*vnon(idim,isp)
        end do
      end do
      
c...  Evaluate R (stored in vno)

      sh  = 1.d0/3.d0

      do idim=1,3
        do inod=1,3
          do jtp = 1,ntp            
            isp=nup(inod,jtp)
            jn1=nup(1,jtp)
            jn2=nup(2,jtp)
            jn3=nup(3,jtp)

            coorx= sh*
     .        (coorpbar(idim,jn1)+coorpbar(idim,jn2)
     .        +coorpbar(idim,jn3))

            vno(idim,isp) =
     .        vno(idim,isp)+sh*coorx*airtp(jtp)

          end do
        end do
      end do

c...  Evaluate S 

      do isp=1,nsp
        do idim=1,3
          coorpbar(idim,isp) = vno(idim,isp) / airesp(isp) 
        end do
      end do

c...  Add coorpbar to coorp

      do isp=1,nsp
        is= node2d3d(isp)
        do idim=1,3
          coorp(idim,isp) = coorp(idim,isp) + coorpbar(idim,isp)
        end do
        ctrl(isp)=0.0d0
      end do
      

      end
