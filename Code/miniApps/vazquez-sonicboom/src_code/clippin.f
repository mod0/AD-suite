      subroutine clippin(ipoin,unew,uold)
      
      INCLUDE 'Param3D.h'

      real*8 uold(5),unew(5),rm,pm,ux,vx,wx,pmo
      integer ipoin,nclip

      rm  = unew(1)
      ux  = unew(2)/rm
      vx  = unew(3)/rm
      wx  = unew(4)/rm
c     
      pm  = gam1*(unew(5) -
     &     0.5*rm*(ux*ux + vx*vx + wx*wx))
      
      if (pm.lt.1.0e-6) then
         
         unew(1)=uold(1)
         unew(2)=uold(2)
         unew(3)=uold(3)
         unew(4)=uold(4)
         unew(5)=uold(5)
         rm  = unew(1)
         ux  = unew(2)/rm
         vx  = unew(3)/rm
         wx  = unew(4)/rm
c     
         pmo  = gam1*(unew(5) -
     &        0.5*rm*(ux*ux + vx*vx + wx*wx))

c         write (6,*) 'Clipping p node: ',ipoin,pm,pmo

      end if
      
      
      end
