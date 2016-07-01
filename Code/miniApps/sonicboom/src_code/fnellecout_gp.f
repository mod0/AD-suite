      subroutine fnellecout(CCD,CCL)
c
C**** Calcul de la fonctionnelle cout J(W). 
C            nt
C            --  1                    2
C     J(W) = \   -  Volu(k) * |Grad P| 
C            /   4
C            --
C            k
C     
C                                         2                  2
C            + coeftrainee*(CD - CDTARGET)  + (CL - CLTARGET)
C
C     La fonctionnelle est stockee dans le scalaire Cout
C


      include 'Param3D.h'
      include 'Paramopt3D.h'

      integer*4 isp,is,kn,jt,ia,ib,iwhico,iwhiex
      real*8
     .  u2,pression(4),adjoi(4),grap(3),graj(3),xgrapn2,
     .  dicd,dicl,ccl,ccd,usro,yx,xe,xshape,dixma,xfin,edixm,
     .  coef,uph(5,4),um(3),u,v,w,voltet,
     .  edrag,elift,epres,frali,ybelow,yextra,
     .  x(4),y(4),z(4),b(4),c(4),d(4),vol6,us6,cout1

c      elift = 1.0d0
c      epres = 1.0d1
c      elift = 1.0d-1
c      epres = 1.0d0
c      frali = 1.0d1
c      ybelow= -0.06d0

      edrag = faczz(1)
      elift = faczz(2)
      epres = faczz(3)
      frali = faczz(4)
      ybelow= faczz(5)
      yextra= faczz(7)
      iwhico= int(faczz(6))
      iwhiex= int(faczz(8))
            

      edixm = 0.0d0
      xfin  = 0.0d0
      
      
      cout  = 0.0d0
      cout1 = 0.0d0
      us6   = 1.0d0/6.0d0
      xshape= 0.25d0

      do jt=1,nt                                            !loop over the elements

        xe =
     .    xshape*(coor(iwhiex,nu(1,jt))+coor(iwhiex,nu(2,jt))
     .    +coor(iwhiex,nu(3,jt))+coor(iwhiex,nu(4,jt)))
        yx =
     .    xshape*(coor(iwhico,nu(1,jt))+coor(iwhico,nu(2,jt))
     .    +coor(iwhico,nu(3,jt))+coor(iwhico,nu(4,jt)))

c        if (yx.lt.-0.5d0) then
c        if (yx.lt.-1.2d0) then
c        if (yx.lt.ybelow
c     .    .and. xe.lt.0.8d0) then
c        if (yx.lt.ybelow .and. xe.gt.yextra) then
        if (yx.lt.ybelow) then
          grap(1) = 0.0d0
          grap(2) = 0.0d0
          grap(3) = 0.0d0
          graj(1) = 0.0d0
          graj(2) = 0.0d0
          graj(3) = 0.0d0

          do kn=1,4
            
            is = nu(kn,jt)

            u = ua(2,is)/ua(1,is)
            v = ua(3,is)/ua(1,is)
            w = ua(4,is)/ua(1,is)
            u2= (u*u+v*v+w*w)
            
            pression(kn)= gam1*(ua(5,is)-0.5d0*ua(1,is)*u2)

            adjoi(kn)   = piadj(1,is)/2.d0
            
            x(kn)= coor(1,nu(kn,jt))
            y(kn)= coor(2,nu(kn,jt))
            z(kn)= coor(3,nu(kn,jt))
            
          end do
          
          call gradfb(x,y,z,b,c,d,vol6)
          
          do kn=1,4
            grap(1)= grap(1) + pression(kn) * ( b(kn)/vol6 )
            grap(2)= grap(2) + pression(kn) * ( c(kn)/vol6 )
            grap(3)= grap(3) + pression(kn) * ( d(kn)/vol6 )
            graj(1)= graj(1) +    adjoi(kn) * ( b(kn)/vol6 )
            graj(2)= graj(2) +    adjoi(kn) * ( c(kn)/vol6 )
            graj(3)= graj(3) +    adjoi(kn) * ( d(kn)/vol6 )
          end do

          
          xgrapn2 = grap(1)*grap(1)+grap(2)*grap(2)+grap(3)*grap(3)

c          xgrapn2 = grap(1)*grap(1)                         !ojo, prueba...
          
          
c... hay que pensar en esto para la adjunta: ni idea que es dadj/dw..!
c     .      + graj(1)*graj(1)+graj(2)*graj(2)+graj(3)*graj(3)

          voltet = us6*vol6

          cout1= cout1 + xgrapn2*voltet
          
        end if

      end do

      dicd  = (ccd-cdtarget)
      dicl  = (ccl-cltarget/frali)
      dixma = (xfin-xman)
      
      cout =
     .  edrag*dicd*dicd
     .  + elift*dicl*dicl
     .  + epres*coefpres*cout1
     .  + edixm*dixma*dixma

c... This is to fix a maximum movement ahead..      
c     .  + edixm*dixma*dixma
      
 
      
      end
