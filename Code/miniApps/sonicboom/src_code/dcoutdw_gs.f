      subroutine dcoutdw(CTRL,CCL,CCD)
C
C*** Cette procedure calcule la derivee de la fonctionnelle cout J(W,gamma) par
C    rapport a W (programmed by M. Vazquez, feb.2002). 
c    It takes into account the pressure gradient below the object.
c
c    There are two loops. The first one is done over the boundary nodes
c    and evaluates the drag and lift variations. The second loop is done
c    over the elements, keeping those under a given y-coordinate value.
c    While the first integral is done over the surface elements, the
c    second is over the volume elements.
C
C
      
      include 'Param3D.h'
      include 'Paramopt3D.h'

      integer*4 isp,is,kn,kw,jt,ia,ib,iwhico
      real*8
     .  u2,pression(4),entropy(4),
     .  adjoi(4),dpdw(5,4),grap(3),xgrapn,xdpdw(5),
     .  diff,diff1,diff2,ccl,ccd,usro,yx,xx,xshape,ybelow,
     .  edrag,elift,epres,frali,dicd,dicl,dixma,xfin,edixm,
     .  ctrl(nnsp),vncoq(3,nnsp),coef,uph(5,4),um(3),
     .  x(4),y(4),z(4),b(4),c(4),d(4),vol6,us6,u,v,w,volttra


      edrag = faczz(1)
      elift = faczz(2)
      epres = faczz(3)
      frali = faczz(4)
      ybelow= faczz(5)

      iwhico= int(faczz(6))                                 
      
      write(6,*)
     .  'Entree dans DCoutDW - GP'
      call flunow(6)
c

      coef = roin*(uxin*uxin + uyin*uyin + uzin*uzin)/2.0d0

      call normcoq(ctrl,vncoq)
      
      us6=1.0d0/6.0d0
      xshape=0.25d0
      
      do is=1,ns
        do kw=1,5
          djdw(kw,is)=0.0d0
        end do
      end do

c...  1. Loop over the coque nodes to calculate drag-lift variation contribution

      
      do isp=1,nsp

        is   = node2d3d(isp)
        
        u = ua(2,is)/ua(1,is)
        v = ua(3,is)/ua(1,is)
        w = ua(4,is)/ua(1,is)
        u2= (u*u+v*v+w*w)
                
        dpdw(1,1)=  gam1*0.5d0*u2
        dpdw(2,1)= -gam1*u
        dpdw(3,1)= -gam1*v
        dpdw(4,1)= -gam1*w
        dpdw(5,1)=  gam1

        diff =
     .    (cos(tetacdcl)*vncoq(1,isp)
     .    +sin(tetacdcl)*vncoq(iwhico,isp))/coef        
        diff1=
     .    (cos(tetacdcl)*vncoq(iwhico,isp)
     .    -sin(tetacdcl)*vncoq(1,isp))/coef        

        dicd= (ccd-cdtarget)
        dicl= (ccl-cltarget/frali)
        
        do kw=1,5
          djdw(kw,is)=djdw(kw,is)
     .      + edrag * 2.0d0 
     .      *     dicd * diff  * dpdw(kw,1)
     .      + elift * 2.0d0 
     .      *     dicl * diff1 * dpdw(kw,1)
        end do
        
      end do

      
c...  2. Loop over the elements below the wing (i.e. whose central point y-component
c     is lower than ybelow) to calculate the volume integral which corresponds
c     to the gradient norm.
      

      do jt=1,nt

        xx =
     .    xshape*(coor(1,nu(1,jt))+coor(1,nu(2,jt))
     .    +coor(1,nu(3,jt))+coor(1,nu(4,jt)))
        yx =
     .    xshape*(coor(iwhico,nu(1,jt))+coor(iwhico,nu(2,jt))
     .    +coor(iwhico,nu(3,jt))+coor(iwhico,nu(4,jt)))

        if (yx.lt.ybelow) then

          grap(1) = 0.0d0
          grap(2) = 0.0d0
          grap(3) = 0.0d0

          do kn=1,4
            
            is = nu(kn,jt)
            
            u = ua(2,is)/ua(1,is)
            v = ua(3,is)/ua(1,is)
            w = ua(4,is)/ua(1,is)
            u2= (u*u+v*v+w*w)
            
            
            pression(kn)=
     .        gam1*(ua(5,is)-0.5d0*ua(1,is)*u2)
            entropy(kn) =
     .        (pression(kn)/pin)*(roin/ ua(1,is))**gam - 1.d0
            
            dpdw(1,kn)=
     .        - gam * (roin)**gam * pression(kn)
     .        / ua(1,is)**(gam+1) / pin

            dpdw(2,kn)= 0.d0
            dpdw(3,kn)= 0.d0
            dpdw(4,kn)= 0.d0
            dpdw(5,kn)= 0.d0
            
            x(kn)= coor(1,nu(kn,jt))
            y(kn)= coor(2,nu(kn,jt))
            z(kn)= coor(3,nu(kn,jt))
            
            
          end do
          
          call gradfb(x,y,z,b,c,d,vol6)
          
          
          do kn=1,4
            
            grap(1)= grap(1) + entropy(kn) * ( b(kn)/vol6 )
            grap(2)= grap(2) + entropy(kn) * ( c(kn)/vol6 )
            grap(3)= grap(3) + entropy(kn) * ( d(kn)/vol6 )

          end do
          
          xgrapn = grap(1)*grap(1)+grap(2)*grap(2)+grap(3)*grap(3)
          xgrapn = sqrt(xgrapn)

          volttra= us6*vol6

          do kn=1,4
c...  Assembly: 
            is = nu(kn,jt)
            do kw=1,5
              djdw(kw,is)=djdw(kw,is)
     .          + epres*coefpres*2.0d0*volttra
     .          *      xgrapn*xshape*dpdw(kw,kn)
            end do
          end do
          
        end if

      end do
      
      end
      
