      subroutine intnor
c
c     
c     This sub evaluates the normals using the divergence theorem
c     in each node.
c     
c
c     
c
c     
c
c
      include 'Param3D.h'
      include 'Paramopt3D.h'

      real*8 sh,vnox,vnono,sol(nsmax,9),
     .  x(4),y(4),z(4),b(4),c(4),d(4),vol6,us6,u,v,w,volttra
      integer jtp, isp, is, idim,inod,izztask,jt,kn

      
      do is=1,ns
        sol(is,1)=0.d0
        sol(is,2)=0.d0
        sol(is,3)=0.d0
      end do
      
      do jt=1,nt
        
        do kn=1,4
            x(kn)= coor(1,nu(kn,jt))
            y(kn)= coor(2,nu(kn,jt))
            z(kn)= coor(3,nu(kn,jt))                      
        end do

        call gradfb(x,y,z,b,c,d,vol6)

        do kn=1,4
          is= nu(kn,jt)
          sol(is,1)=sol(is,1)+b(kn)/6.d0
          sol(is,2)=sol(is,2)+c(kn)/6.d0
          sol(is,3)=sol(is,3)+d(kn)/6.d0
        end do
                
        
      end do

      do is=1,ns
        isp= node3d2d(is)
        if (isp.gt.0) then
          write (6,*) isp,vnon(1,isp),vnon(2,isp),vnon(3,isp)
          vno(1,isp)= sol(is,1)
          vno(2,isp)= sol(is,2)
          vno(3,isp)= sol(is,3)
          vnono=
     .      vno(1,isp)*vno(1,isp)+vno(2,isp)*vno(2,isp)
     .      + vno(3,isp)*vno(3,isp)
          vnono= sqrt(vnono)
          vnon(1,isp)= sol(is,1)/vnono
          vnon(2,isp)= sol(is,2)/vnono
          vnon(3,isp)= sol(is,3)/vnono
          write (6,*) isp,vnon(1,isp),vnon(2,isp),vnon(3,isp)
        end if        
      end do



      
      end
