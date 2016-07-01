      subroutine zzoutsur(ctrl,izztask)
C
C*** Write the trace of the results at the optimized boundary
C

      include 'Param3D.h'
      include 'Paramopt3D.h'

      integer  isp,is,i,ii,izztask,
     .  imate,nnodb,ndime,idime,ifac,inodb,ifaw,ielem
      real*8   ctrl(nnsp),sol(nsmax,9),coefr(5),coorpbar(3,nnsp),
     .  vcarre,denom,reval,dimv(3),ctrans,cmovem,cvcoe,cpcoe,tinf
      real*4   tact4
      character*4 wopos
      character*8 sumsh,vomsh,sures,vores
      integer*4 lunit

      vomsh='volu.msh'
      sumsh='surf.msh'
      vores='volu.res'
      sures='surf.res'

      if (kmov.gt.0) then
        vomsh='vmm'//char(kmov+65)//'.msh'
        vores='vmm'//char(kmov+65)//'.res'
c        sumsh='smm'//char(kmov+65)//'.msh'
c        sures='smm'//char(kmov+65)//'.res'
      end if
      
      open (11,file=vores,status='unknown',access='append')
      open (13,file=sures,status='unknown',access='append')

      open (12,file='conv.data',status='unknown',access='append')
      
      tact4=dble(itopt)

      coefr(1)                     = rhoref
      coefr(2)                     = rhoref*vref
      coefr(3)                     = rhoref*vref
      coefr(4)                     = rhoref*vref
      coefr(5)                     = pref
c

      
      if (itopt.eq.0 .or. izztask.eq.2 .or. imvmsh.gt.0) then

        if (itopt.eq.0 .or. izztask.eq.2) then
          write(12,*)
     .      '# CONVERGENCE FILE'
          write(12,*)
     .      '# itopt  gradj  cout  cl  cd  cltarget  cdtarget '
        end if
          
        if (ncontopt.eq.0 .or. imvmsh.gt.0) then

c...     Surface mesh for gid:

                    
          open (15,file=sumsh,status='unknown')

          imate=2
          nnodb=3
          ndime=3
          
          write(15,2000) ndime,nnodb

c     .      'MESH ',' dimension ',ndime,
c     .      ' Elemtype ','Triangle', ' Nnode ',nnodb
          
          
          write(15,*) 'coordinates'
          do isp=1,nsp
            is=node2d3d(isp)
            write(15,100) isp,(coin(idime,is),idime=1,3)
          end do
          write(15,*) 'end coordinates'
          write(15,*) 'elements'
          ifaw=0
          do ifac=1,nfac
            if (logfac(ifac).eq.-2) then 
              ifaw=ifaw+1
              write(15,*)
     .          ifaw,(node3d2d(nsfac(inodb,ifac)),inodb=1,3),imate
            end if
          end do
          write(15,*) 'end elements'
          
          
          close(15)

c...     Volume mesh for gid:

          open (15,file=vomsh,status='unknown')

          imate=2
          nnodb=4
          ndime=3
          
          write(15,2100) ndime,nnodb

c     .      'MESH ',' dimension ',ndime,
c     .      ' Elemtype ','Tetrahedra', ' Nnode ',nnodb
          
          
          write(15,*) 'coordinates'
          do is=1,ns
            write(15,100) is,(coor(idime,is),idime=1,3)
          end do
          write(15,*) 'end coordinates'
          write(15,*) 'elements'
          ifaw=0
          do ielem=1,nt
            write(15,*)
     .        ielem,(nu(inodb,ielem),inodb=1,nnodb),imate
          end do
          write(15,*) 'end elements'
          
          
          close(15)

          
          write (6,*)
     .      '--> zzoutsur: writing GiD files ',vomsh,' and ',sumsh
        
          
        end if

      end if
      if (izztask.eq.1 .or. izztask.eq.2) then

        write(12,*)
     .    itopt,reazz(1),reazz(2),reazz(3),reazz(4),reazz(5),reazz(6)

        do is = 1,ns
          sol(is,1) = ua(1,is)
          sol(is,2) = ua(2,is)/ua(1,is)
          sol(is,3) = ua(3,is)/ua(1,is)
          sol(is,4) = ua(4,is)/ua(1,is)
          sol(is,5) = ua(5,is)
          sol(is,6) = gam1*(sol(is,5)*coefr(5)
     .      - 0.5*sol(is,1)*(sol(is,2)
     .      *sol(is,2) + sol(is,3)*sol(is,3) + sol(is,4)*sol(is,4)))
          sol(is,7) = SQRT( (sol(is,2)*sol(is,2)+sol(is,3)*sol(is,3)
     .      + sol(is,4)*sol(is,4))/(gam*sol(is,6)/sol(is,1)))
          vcarre = sol(is,2)**2+sol(is,3)**2+sol(is,4)
          denom = roin*0.5*vcarre
          sol(is,8) = (pin - sol(is,6))/denom
          sol(is,9) = (sol(is,6)/pin)*(roin/sol(is,1))**gam - 1.
        end do
c
c     sol(is,1) ---> Rho
c     sol(is,2) ---> u
c     sol(is,3) ---> v
c     sol(is,4) ---> w
c     sol(is,5) ---> E
c     sol(is,6) ---> P
c     sol(is,7) ---> Mach
c     sol(is,8) ---> Cp
c     sol(is,9) ---> Entropie
c
        
         
        write (6,*)
     .    '--> zzoutsur: writing GiD files ',vores,' and ',sures
        
c...   volu.res, volume:
        
        wopos='PRES'
        write(11,1000) wopos,1,tact4,1,1,0                        
        reval=1.0d0
        do is=1,ns
          write(11,1200) is,sol(is,6)/reval
        end do
        
        wopos='MACH'
        write(11,1000) wopos,1,tact4,1,1,0                        
        reval=1.0d0
        do is=1,ns
          write(11,1200) is,sol(is,7)/reval
        end do

        wopos='DENS'
        write(11,1000) wopos,1,tact4,1,1,0                        
        reval=1.0d0
        do is=1,ns
          write(11,1200) is,sol(is,1)/reval
        end do

        wopos='ENER'
        write(11,1000) wopos,1,tact4,1,1,0
        do is=1,ns
          write(11,1200) is,sol(is,5)
        end do

        wopos='VELO'
        write(11,1000) wopos,1,tact4,2,1,0                        
        do is=1,ns
          write(11,1200) is,(sol(is,i),i=2,4)
        end do

        
c        wopos='ENTR'
c        write(11,1000) wopos,1,tact4,1,1,0                        
c        reval=1.0d0
c        do is=1,ns
c          write(11,1200) is,sol(is,9)/reval
c        end do
c
c        wopos='ADCO'
c        write(11,1000) wopos,1,tact4,1,1,0                        
c        do is=1,ns
c          write(11,1200) is,piadj(1,is)
c        end do
c
c        wopos='ADEN'
c        write(11,1000) wopos,1,tact4,1,1,0                        
c        do is=1,ns
c          write(11,1200) is,piadj(5,is)
c        end do
c
c        wopos='ADMO'
c        write(11,1000) wopos,1,tact4,2,1,0                        
c        do is=1,ns
c          write(11,1200) is,(piadj(i,is),i=2,4)
c        end do

        if (imvmsh.eq.1) then
          wopos='DIMV'
          write(11,1000) wopos,1,tact4,2,1,0                        
          reval=1.0d0
          do is=1,ns
            do i = 1,3
c              dimv(i) = (coor(i,is) - coin(i,is))
              dimv(i) = 0.0d0
            end do
            write(11,1200) is,(dimv(i)/reval,i=1,3)
          end do
        end if
        
c...   surf.res, surface: 

        wopos='DISU'
        write(13,1000) wopos,1,tact4,2,1,0                        
        reval=1.0d0
        do isp=1,nsp
          is= node2d3d(isp)
          do i = 1,3
            ctrans=  ctrl(isp)  * vnon(i,isp)
            cmovem=  coor(i,is) - coin(i,is)
c            cmovem= 0.0d0
            coorpbar(i,isp) = ctrans + cmovem
          end do
          write(13,1200) isp,(coorpbar(i,isp)/reval,i=1,3)
        end do

        wopos='CTRL'
        write(13,1000) wopos,1,tact4,1,1,0                        
        reval=1.0d0
        do isp=1,nsp
          is= node2d3d(isp)
          write(13,1200) isp,ctrl(isp)
        end do
        
        wopos='PRES'
        write(13,1000) wopos,1,tact4,1,1,0                        
        reval=1.0d0
        do isp=1,nsp
          is= node2d3d(isp)
          write(13,1200) isp,sol(is,6)/reval
        end do
        
        wopos='MACH'
        write(13,1000) wopos,1,tact4,1,1,0                        
        reval=1.0d0
        do isp=1,nsp
          is= node2d3d(isp)
          write(13,1200) isp,sol(is,7)/reval
        end do

        wopos='GRAD'
        write(13,1000) wopos,1,tact4,1,1,0                        
        reval=1.0d0
        do isp=1,nsp
          is= node2d3d(isp)
          write(13,1200) isp,grad(isp)
        end do

      end if
      
      close(11)
      close(12)
      close(13)
      
      lunit=11
      call flunow(lunit)
      lunit=12
      call flunow(lunit)
      lunit=13
      call flunow(lunit)

  100 format(i8,3f20.7)
      
 1000 format(a15,i5,e12.5,3i5)
 1500 format(a5,10x,i5,e12.8,3i5)
 1200 format(1x,i8,3(2x,e12.5))
 1400 format(1x,i8,2x,f12.5)
 2000 format(17h MESH   dimension,
     .  i3,28h Elemtype  Triangle    Nnode,i3)
 2100 format(17h MESH   dimension,
     .  i3,28h Elemtype  Tetrahedra  Nnode,i3)
     

      end
      
