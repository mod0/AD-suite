      SUBROUTINE MAIL3D(coextr,rwhiex)
c---------------------------------------------------------------------   
c Reads the mesh definition and constructs the metrics
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------
c     Global variables definition
      INTEGER itrs, itrt
      COMMON/comtri/itrs, itrt 
c
c     Local variables definition 
      INTEGER ndle,ifus,iala,iotr,isli,iwhiex,idir,kdifi,idifi
      PARAMETER (ndle = 4)
      INTEGER numf
      INTEGER is, is1, is2, is3, i, k, ifac, jt, it, iseg,c,j
      INTEGER klog,npore,irea
      INTEGER idrap(ntmax)
      real*4 tact4
      REAL*8    xmin, ymin  , zmin, zmiri,zmile,abs1,abs2,abs3
      REAL*8    xmax, ymax  , zmax, zmari,zmale,dzm,dxm,dym,
     .  coextr,rwhiex,xtract(3),a1,a2,a3,pian,
     .  aenorm(20),opnorm(20),aecor(3),opcor(3),flexi(3,20)
c
c     Initializations
c
      iwhiex=int(rwhiex)                                    !if zero, do not extract anything
      
      DO 1 is=1,ns
         ndeg(is)                  = 0
         DO 2 k=1,ndmax
            jaret(is,k)            = 0
2        CONTINUE
1     CONTINUE

      numf                         = 60

      do i=1,20
        opnorm(i)  =0.d0
        aenorm(i)  =0.d0
        flexi(1,i)=0.d0
        flexi(2,i)=0.d0
        flexi(3,i)=0.d0
      end do
      
      REWIND(numf)

      READ(numf, *) ns, nt, nfac

      IF (ns .GT. nsmax) THEN
         write (6,*) ns,nsmax
         WRITE(6, *) 'NSMAX too small'
         CALL TILT
      ENDIF

      READ(numf, *) ((coor(i,is)   , i=1,3), is=1,ns)
      READ(numf, *) ((nu(k,jt)     , k=1,4), jt=1,nt)
      READ(numf, *) (logfac(ifac)  , ifac=1,nfac)
      READ(numf, *) ((nsfac(k,ifac), k=1,3), ifac=1,nfac)     

      CLOSE(numf)

C [llh]      write (6,*) 
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) '                   2- Lecture du maillage :'
C [llh]      WRITE(6,*) '                   ************************ '
C [llh]      WRITE(6,*) ' '
C [llh]      IF (ns.eq.909) WRITE(6,*)'Travail sur le tube a choc a 909 noeuds'
C [llh]      IF (ns.eq.2203) WRITE(6,*) 'Travail sur l"aile M6 '
C [llh]      IF (ns.eq.10188) WRITE(6,*) 'Travail sur le demi avion Falcon' 
C [llh]      IF (ns.eq.225) WRITE(6,*)'Travail sur la tuyere'
C [llh]      if (ns.eq.30514) WRITE(6,*) 'Travail sur le complete avion Falcon'
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) 'Nombre de noeuds :',ns
C [llh]      WRITE(6,*) 'Nombre de tetraedres :',nt
C [llh]      WRITE(6,*) 'Nombre de facettes :',nfac
C [llh]      WRITE(6,*) ' '     
c
c     Keep original coordinates
c

      do i=1,3
        do is=1,ns
          coin(i,is)=coor(i,is)
        end do
      end do
      
c      
c     Correct the coordinates for initial global grid displacements
c
c     They are contained in these files: 
c 
c          disp.aero.(CHAR)    ---->   containing aeroelastic DISPLACEMENTS (ONLY AEROELASTIC!)
c                                      produced by AERO, relatives to its own
c                                      cor0.data ( = corf.aero.(CHAR-1))
c                                      They will give a (NON-OPTIMIZED) CRUISE GRID.
c
c          corf.aero.(CHAR-1)  ---->   containing PHYSICAL COORDINATES
c                                      produced by MARCO'S in the preceding 
c                                      iteration, i.e. the PREVIOUD OPTIMIZED GRID. It corresponds to 
c                                      the last OPTIMIZED LANDED GRID.
c
c          crui.aero.(CHAR-1)  ---->   containing PHYSICAL COORDINATES
c                                      produced by MARCO'S in the preceding 
c                                      iteration. It corresponds to 
c                                      the last OPTIMIZED CRUISE GRID.
c
c
c
c     Then, the new coor will be:
c
c          coor = coorf(CHAR-1) + disp(CHAR) 
c
c
c
c
c
c     CHAR is A, B, C, ...  CHAR-1 means the precedent letter. 
c     For the initial state, CHAR = A. Therefore, it won't read any corf.aero.(CHAR-1)
c     Recall that char(65) = A
c
c

      kdifi= intzz(3)                             

      if (kdifi.lt.0) then

        open (numf,file='mesh.init',status='old')          

        do is=1,ns
          read(numf, 115) i, coco(1,is),coco(2,is),coco(3,is)
          ua(1,is)=coor(1,is)
          ua(2,is)=coor(2,is)
          ua(3,is)=coor(3,is)
        end do          
                
        close(numf)

        opnorm(1)= 0.d0
        do i=1,3
          do is=1,ns
            opcor(i)      = coco(i,is)-ua(i,is)
            opnorm(1) = opnorm(1)+opcor(i)*opcor(i)
            ua(i,is)      = 0.d0
            coor(i,is)    = coco(i,is)                      !initial mesh 
          end do            
        end do          
        
        opnorm(1)= sqrt(opnorm(1))/ns
        
        
      end if

      
      if (kdifi.gt.1) then
c
c...  Read LANDED previously optimized meshes (in corf.aero.*-1)        
c        
        do idifi=2,kdifi
          open (numf,file='corf.aero.'//char(idifi-2+65),status='old')          
          do is=1,ns
            read(numf, 115) i, coco(1,is),coco(2,is),coco(3,is)
          end do          

          close(numf)
          
          opnorm(idifi)= 0.d0
          do i=1,3
            do is=1,ns
              opcor(i)      = coco(i,is)-coin(i,is)
              opnorm(idifi) = opnorm(idifi)+opcor(i)*opcor(i)

              if (idifi.eq.kdifi) coor(i,is)= coco(i,is)    !initial mesh is the previous optimized (landed)

            end do            
          end do          
          
c...  mean optimization modifications norm: ||(Omega_opt^i+1 - Omega_opt^i)|| / ns
          
          opnorm(idifi)= sqrt(opnorm(idifi))/ns
          
        end do
        
  115   format(i6,3(e15.8,1x))
        
      end if
      
      if (kdifi.gt.0) then
c
c...  Read new CRUISE aeroelastic deformations (in disp.aero.*)
c        
        do idifi=1,kdifi
          open (numf,file='disp.aero.'//char(idifi-1+65),status='old')
          
          read(numf,'(a)')
          read(numf,*) npore
          if (npore.eq.ns) then
            irea=0              
            do while(irea.eq.0) 
              read(numf,*,end=555) tact4
              do i=1,npore
                read(numf,*) (coco(j,i),j=1,3)
              end do
            end do
          else
            write (6,*) 
            write (6,*) 
            write (6,*) '================================'
            write (6,*) 'ERROR WHEN READING disp.aero!!!!'
            write (6,*) '================================'
            write (6,*) 'Different number of points...'
            write (6,*) 'Skipped, zero initial displacements!'
            write (6,*) 
            write (6,*) 
          end if
          
  555     close(numf)

          aenorm(idifi)= 0.d0
          do i=1,3
            do is=1,ns
              aecor(i)      = coco(i,is)
              aenorm(idifi) = aenorm(idifi)+aecor(i)*aecor(i)
              ua(i,is)      = coco(i,is)
            end do            
          end do          
          
c...  mean aeroelastic modifications norm ||(Delta Omega_aer^i+1 - Delta Omega_aer^i)|| / ns
          
          aenorm(idifi)= sqrt(aenorm(idifi))/ns
          

          if (ns.eq.22014) then
c...  only the agard problem

            pian= 3.1415927d0
            
            a1= coin(1,4881)-coin(1,4915)                   ! initial Dx = x1-x2
            a2= coin(3,4881)+coco(3,4881)                   ! deformed z1
            a3= coin(3,4915)+coco(3,4915)                   ! deformed z2
            flexi(1,idifi)=180.0d0*asin((a3-a2)/a1)/pian
            
            a1= coin(2,4915)-coin(2,1)                      ! initial Dy = y2-y3
            a2= coin(3,4915)+coco(3,4915)                   ! deformed z2
            a3= coin(3,1)  +coco(3,1)                       ! deformed z3
            flexi(2,idifi)=180.0d0*asin((a2-a3)/a1)/pian
            
            a1= coin(2,4881)-coin(2,101)                    ! initial Dy = y1-y4
            a2= coin(3,4881)+coco(3,4881)                   ! deformed z1
            a3= coin(3,101)   +coco(3,101)                  ! deformed z4
            flexi(3,idifi)=180.0d0*asin((a2-a3)/a1)/pian


          else if (ns.eq.12963) then
c...  only the sbj's delta wings            

            pian= 3.1415927d0
            
            a1= coin(1,10134)-coin(1,8762)                  ! initial Dx = x1-x2
            a2= coin(3,10134)+coco(3,10134)                 ! deformed z1
            a3= coin(3,8762)+coco(3,8762)                   ! deformed z2
            flexi(1,idifi)=180.0d0*asin((a3-a2)/a1)/pian
            
            a1= coin(2,8762)-coin(2,2839)                   ! initial Dy = y2-y3
            a2= coin(3,8762)+coco(3,8762)                   ! deformed z2
            a3= coin(3,2839)  +coco(3,2839)                 ! deformed z3
            flexi(2,idifi)=180.0d0*asin((a2-a3)/a1)/pian
            
            a1= coin(2,10134)-coin(2,11589)                 ! initial Dy = y1-y4
            a2= coin(3,10134)+coco(3,10134)                 ! deformed z1
            a3= coin(3,11589)   +coco(3,11589)              ! deformed z4
            flexi(3,idifi)=180.0d0*asin((a2-a3)/a1)/pian

            
          end if
          
        end do

        open (numf,file='aeop.conv',status='unknown')
        write(numf,*)
     .    '# CONVERGENCE FILE FOR AEROELASTIC COUPLED OPTIMIZATION'
        write(numf,*) '# 1. global iteration'
        write(numf,*) '# 2. optimization norm'
        write(numf,*) '# 3. aeroelastic norm'
        write(numf,*) '# 4. aeroelastic torsion'
        write(numf,*) '# 5. aeroelastic flexion in the leading edge'
        write(numf,*) '# 6. aeroelastic flexion in the trailing edge'
        
        do idifi=1,kdifi
          write(numf,*) idifi,opnorm(idifi),aenorm(idifi),
     .      flexi(1,idifi),flexi(2,idifi),flexi(3,idifi)
        end do
        close(numf)

      end if

      if (kdifi.gt.0) then
        
        open (numf,file='crui.aero.'//char(kdifi-1+65),status='unknown')
        
        if (npore.eq.ns) then
          do i=1,npore
            coor(1,i)= coor(1,i)+ua(1,i)
            coor(2,i)= coor(2,i)+ua(2,i)
            coor(3,i)= coor(3,i)+ua(3,i)
            write(numf,*) i,coor(1,i),coor(2,i),coor(3,i)
          end do
        end if
        
        close(numf)
        
C [llh]        write (6,*) 
C [llh]        write (6,*) ' --->>>> INITIAL MESH MOVEMENT GIVEN:'
C [llh]        write (6,*) 
C [llh]        write (6,*) 'General iteration number  :   ',kdifi
C [llh]        write (6,*) 
C [llh]        if (kdifi.gt.1)
C [llh]     .    write (6,*) 'Initial coordinates file:     ',
C [llh]     .    'corf.aero.'//char(kdifi-2+65)
C [llh]        write (6,*) 'Initial displacements file:   ',
C [llh]     .    'disp.aero.'//char(kdifi-1+65)
C [llh]        write (6,*) 
C [llh]        write (6,*) 'Initial mesh movement files read ok...'
C [llh]        write (6,*) 
      else if (kdifi.lt.0) then
C [llh]        write (6,*) 
C [llh]        write (6,*) ' --->>>> INITIAL MESH MOVEMENT GIVEN:'
C [llh]        write (6,*) 
C [llh]        if (kdifi.gt.1)
C [llh]     .    write (6,*) 'Initial coordinates file:     ',
C [llh]     .    'mesh.init'
C [llh]        write (6,*) 
C [llh]        write (6,*) 'Initial mesh movement files read ok...'
C [llh]        write (6,*) 
      end if

c
c     Computing the coordinates extrema
c
      xmin                         = 1.0e+16
      ymin                         = 1.0e+16
      zmin                         = 1.0e+16
      xmax                         =-1.0e+16
      ymax                         =-1.0e+16
      zmax                         =-1.0e+16

      if (ns.eq.30514) then

        zmiri= 1.4d0
        zmari= 9.1d0
        zmile=-1.4d0
        zmale=-9.1d0

        zmin = 1.4d0
        zmax = 9.1d0

        ymax = -0.27                                        !forget engines and tail..
        
        xmax = -4.0                                         !but include the nose
        
      else
        
        DO 100 is=1,ns
          xmin                      = MIN(xmin, coor(1,is)) 
          ymin                      = MIN(ymin, coor(2,is)) 
          zmin                      = MIN(zmin, coor(3,is)) 
          xmax                      = MAX(xmax, coor(1,is)) 
          ymax                      = MAX(ymax, coor(2,is)) 
          zmax                      = MAX(zmax, coor(3,is)) 
  100   CONTINUE  
c

        dzm = (zmax-zmin)/1.e5
        zmin=  zmin+dzm
        
      end if
        
c
c*** Verification :
c
c      IF (((NS.EQ.2203).AND.(COEFM1.NE.5)).OR.((NS.EQ.10188).AND.
c     $     (COEFM1.NE.5))) THEN
c         PRINT*,'INCOMPATIBILITE POUR COEFM1'
c         STOP
c      ENDIF
C
c      IF ((NS.EQ.225).AND.(COEFM1.NE.10)) THEN
c         PRINT*,'INCOMPATIBILITE POUR COEFM1'
c         STOP
c      ENDIF
c
C [llh]      write (6,*) 
C [llh]      WRITE(6,*) 'MIN et MAX des coordonnees du maillage : '
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6, 1000) xmin, xmax
C [llh]      WRITE(6, 1100) ymin, ymax
C [llh]      WRITE(6, 1200) zmin, zmax
1000  FORMAT(' Xmin = ', f15.6, ' Xmax = ', f15.6)
1100  FORMAT(' Ymin = ', f15.6, ' Ymax = ', f15.6)
1200  FORMAT(' Zmin = ', f15.6, ' Zmax = ', f15.6)
c      
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) '          Fin de lecture des donnees du maillage '
C [llh]      WRITE(6,*) '          ************************************** '
C [llh]      WRITE(6,*) ' '
C [llh]      WRITE(6,*) ' '
c

c      open (7,file='moco',status='unknown')
      
      xtract(1)=xmin+coextr
      xtract(2)=ymin+coextr
      xtract(3)=zmin+coextr       
         
      IF (coefm1 .EQ. 5) THEN 
c
c*** Reperage des facettes de la carlingue :
c          logfac=200 : carlingue
c          logfac=53  : axe de symetrie
c          logfac=11  : infini
c
         DO 85 is=1,ns
            logfr(is)              = 0            
85       CONTINUE 
c

c
c*** modifying logfac for the aero-falcon problem
c
         ifus=0
         iala=0
         iotr=0
         isli=0
         idir=0
         
         do ifac=1,nfac
           if (ivis.eq.0  .and. (logfac(ifac).eq. 3)) logfac(ifac)=2
           if (ivis.eq.0  .and. (logfac(ifac).eq.-3)) logfac(ifac)=2
           if (diric.eq.1 .and. (logfac(ifac).eq. 4)) logfac(ifac)=5 !force st-w to diri
           if (iwhiex.gt.0.and. (logfac(ifac).eq.-2)) logfac(ifac)=2 !force to extract
         end do
         

         DO 90 ifac=1,nfac
c
c**** On recupere dans log1fac les logfac egaux a 200
c
c
           if ((logfac(ifac).eq.200).or.(logfac(ifac).eq.53).or.
     $       (logfac(ifac).eq.2)) then
             logfr(nsfac(1,ifac))   = 2
             logfr(nsfac(2,ifac))   = 2
             logfr(nsfac(3,ifac))   = 2
             isli=isli+1
           else if ((logfac(ifac).eq.11).or.(logfac(ifac).eq.4)) then
             logfr(nsfac(1,ifac))   = 4                     !steger - warming
             logfr(nsfac(2,ifac))   = 4
             logfr(nsfac(3,ifac))   = 4
             iotr=iotr+1             
           else if (logfac(ifac).eq.5) then
             logfr(nsfac(1,ifac))   = 5                     !fixed (dirichlet)
             logfr(nsfac(2,ifac))   = 5
             logfr(nsfac(3,ifac))   = 5
             idir=idir+1
           endif


c
c**** Cas de l'aile : on cherche a extraire l'aile du maillage 3D
c     Elle est reperee pour logfac=-2
c     This is done ONLY when iwhiex <> 0.
c     
c
c

             
           abs1=coor(iwhiex,nsfac(1,ifac))
           abs2=coor(iwhiex,nsfac(2,ifac))
           abs3=coor(iwhiex,nsfac(3,ifac))
           if (abs1.lt.0.d0) abs1=-abs1
           if (abs2.lt.0.d0) abs2=-abs2
           if (abs3.lt.0.d0) abs3=-abs3

           IF (logfac(ifac) .EQ. 2 .and. iwhiex.gt.0) THEN

c
c...  Check if there are faces 2 to extract and put them to -2
c             
                          
             IF (
     $         (abs1 .GT. xtract(iwhiex)) .and. 
     $         (abs2 .GT. xtract(iwhiex)) .and. 
     $         (abs3 .GT. xtract(iwhiex))) THEN 
c     
               
               if (ns.eq.30514) then
                 
                 if ((coor(2,nsfac(1,ifac)).lt.ymax) .or.   !falcon's wings
     .             (coor(2,nsfac(2,ifac)).lt.ymax) .or.
     .             (coor(2,nsfac(3,ifac)).lt.ymax)) then
                   
                   logfac(ifac)        =-2
                 end if
                 
               else
                 
c                 if (((coor(1,nsfac(1,ifac)).lt.26800.d0) .or.   !OJO!!!!! INCLUYE EL FUSELAJE ENTRE ELLAS PARA EL FALSUP
c     .                (coor(1,nsfac(2,ifac)).lt.26800.d0) .or.
c     .                (coor(1,nsfac(3,ifac)).lt.26800.d0)) 
c     .                           .and. 
c     .               ((coor(1,nsfac(1,ifac)).gt. 8900.d0) .or.
c     .                (coor(1,nsfac(2,ifac)).gt. 8900.d0) .or.
c     .                (coor(1,nsfac(3,ifac)).gt. 8900.d0))) then
c                       
c                       logfac(ifac) = -2                 
c
c                 end if

                 logfac(ifac)        =-2

               end if    

             end if

c             if ((ns.eq.30514) .and.                        !falcon's nose
c     .         ((coor(1,nsfac(1,ifac)).lt.xmax) .or.
c     .         (coor(1,nsfac(2,ifac)).lt.xmax) .or.
c     .         (coor(1,nsfac(3,ifac)).lt.xmax))) then
c               
c               logfac(ifac)        =-2
c               ifus=ifus+1
c
c             end if
             
             
           end if

           if (logfac(ifac) .EQ. -2 .and. iwhiex.gt.0) then
c
c...  Check if there are nodes -2 that should be 2 (out of the "extract" range)
c                          
c             if (abs1 .lt. xtract(iwhiex)) then
c               logfr(nsfac(1,ifac)) = 2
c             end if
c             if (abs2 .lt. xtract(iwhiex)) then
c               logfr(nsfac(2,ifac)) = 2
c             end if
c             if (abs3 .lt. xtract(iwhiex)) then
c               logfr(nsfac(3,ifac)) = 2
c             end if
             
           end if

c

c           if (logfac(ifac).eq.-2) then
c             if ((coor(1,nsfac(1,ifac)).gt.0.25d0) .or. 
c     .         (coor(1,nsfac(2,ifac)).gt.0.25d0) .or.
c     .         (coor(1,nsfac(3,ifac)).gt.0.25d0)) then
c               
c               logfac(ifac)        = 2
c               ifus=ifus+1
c             end if                          
c           end if


           if (logfac(ifac).eq.-2) iala=iala+1
           
   90    CONTINUE  
c
       ENDIF
c

       iala= iala-ifus
C [llh]       write(6,*) 'The surface to optimize comprises: '
C [llh]       write(6,*) '   Wings facets     :', iala
C [llh]       write(6,*) '   Fuselage facets  :', ifus
C [llh]       write(6,*) '   Total            :', iala+ifus
C [llh]       write(6,*) 


c
      IF (itrt .EQ. 1) CALL CALDEG
c
      IF (coefm1 .NE. 2) GOTO 500
c
c***  Si coefm1=2, on travaille sur l'aile.
c     Defining the table of logical references for the 
c     boundary mesh vertices logical 
c
      DO 3 is=1,ns
         logfr(is)                 = 0            
3     CONTINUE 
c
      DO 4 ifac=1,nfac
c
C [llh]        write (6,*) 
        
         is1                       = nsfac(1,ifac)
         is2                       = nsfac(2,ifac)
         is3                       = nsfac(3,ifac)
c
c        If the current boundary face is a slipping/non-slipping 
c        face then test if it is placed on the wing.
c
         IF (logfac(ifac) .EQ. 2) THEN
c
           abs1=coor(iwhiex,nsfac(3,is1))
           abs2=coor(iwhiex,nsfac(3,is2))
           abs3=coor(iwhiex,nsfac(3,is3))
           if (abs1.lt.0.d0) abs1=-abs1
           if (abs2.lt.0.d0) abs2=-abs2
           if (abs3.lt.0.d0) abs3=-abs3
           
c           IF ((abs(coor(iwhiex,is1) .GT. xtract(iwhiex))) .OR. 
c     &          (abs(coor(iwhiex,is2) .GT. xtract(iwhiex))) .OR. 
c     &          (abs(coor(iwhiex,is3) .GT. xtract(iwhiex)))) THEN 
c
           IF ((abs1 .GT. xtract(iwhiex)) .OR. 
     &       (abs2 .GT. xtract(iwhiex)) .OR. 
     &       (abs3 .GT. xtract(iwhiex))) THEN 
             
             logfac(ifac)        =-2
             logfr(is1)          =-2
             logfr(is2)          =-2
             logfr(is3)          =-2
c
           ENDIF
c
         ENDIF
c
c
c        If the current boundary face is a freestream 
c        face then test if it needs to become a slipping face
c
         IF (logfac(ifac) .EQ. 4 .or. logfac(ifac).eq.5) THEN 

           abs1=coor(iwhiex,nsfac(3,is1))- xtract(iwhiex)
           abs2=coor(iwhiex,nsfac(3,is2))- xtract(iwhiex)
           abs3=coor(iwhiex,nsfac(3,is3))- xtract(iwhiex)
           if (abs1.lt.0.d0) abs1=-abs1
           if (abs2.lt.0.d0) abs2=-abs2
           if (abs3.lt.0.d0) abs3=-abs3


c            IF ((ABS(coor(3,is1) - zmin) .LT. 1.0e-016) .AND. 
c     &          (ABS(coor(3,is2) - zmin) .LT. 1.0e-016) .AND.
c     &          (ABS(coor(3,is3) - zmin) .LT. 1.0e-016)) 
            IF ((ABS1 .LT. 1.0e-016) .AND. 
     &          (ABS2 .LT. 1.0e-016) .AND.
     &          (ABS3 .LT. 1.0e-016)) 
     &         logfac(ifac)        = 2 
c
c            IF ((ABS(coor(3,is1) - zmax) .LT. 1.0e-016) .AND. 
c     &          (ABS(coor(3,is2) - zmax) .LT. 1.0e-016) .AND.
c     &          (ABS(coor(3,is3) - zmax) .LT. 1.0e-016)) 
c     &         logfac(ifac)        = 2 
c
            logfr(is1)             = logfac(ifac)
            logfr(is2)             = logfac(ifac)
            logfr(is3)             = logfac(ifac)
c
         ENDIF
c         
4     CONTINUE 
c
500   CONTINUE
c

      if (coefm1.eq.5) then                                 !for me: ALWAYS 5!!!!!

        do ifac=1,nfac
          if (logfac(ifac).eq.4) then
             logfr(nsfac(1,ifac))   = 4                     !correct: s-w second best
             logfr(nsfac(2,ifac))   = 4
             logfr(nsfac(3,ifac))   = 4            
          end if            
        end do

        do ifac=1,nfac
          if (logfac(ifac).eq.5) then
             logfr(nsfac(1,ifac))   = 5                     !correct: dir takes all
             logfr(nsfac(2,ifac))   = 5
             logfr(nsfac(3,ifac))   = 5            
          end if            
        end do

        
        do ifac=1,nfac
          if (logfac(ifac).eq.-2) then                      !final check for logfac=-2
            if (logfr(nsfac(1,ifac)).eq.0)
     .        logfr(nsfac(1,ifac))  = -2    
            if (logfr(nsfac(2,ifac)).eq.0)
     .        logfr(nsfac(2,ifac))  = -2    
            if (logfr(nsfac(3,ifac)).eq.0)
     .        logfr(nsfac(3,ifac))  = -2    
          end if            
        end do

      end if
      

      IF (coefm1.eq.10) THEN
         DO IFAC = 1,NFAC
            IF (((COOR(1,NSFAC(1,IFAC)).EQ.XMIN).AND.
     $           (COOR(1,NSFAC(2,IFAC)).EQ.XMIN).AND.
     $           (COOR(1,NSFAC(3,IFAC)).EQ.XMIN)).OR.
     $           (COOR(1,NSFAC(1,IFAC)).EQ.XMAX).AND.
     $           (COOR(1,NSFAC(2,IFAC)).EQ.XMAX).AND.
     $           (COOR(1,NSFAC(3,IFAC)).EQ.XMAX)) THEN
               LOGFAC(IFAC)=4
            ELSE
               IF (((COOR(2,NSFAC(1,IFAC)).EQ.YMAX).AND.
     $             (COOR(2,NSFAC(2,IFAC)).EQ.YMAX).AND.
     $             (COOR(2,NSFAC(3,IFAC)).EQ.YMAX)).AND.
     $          (((COOR(1,NSFAC(1,IFAC)).GE.0.).AND.
     $              (COOR(1,NSFAC(1,IFAC)).LE.2.))
     $              .AND.
     $          ((COOR(1,NSFAC(2,IFAC)).GE.0.).AND.
     $              (COOR(1,NSFAC(2,IFAC)).LE.2.))
     $              .AND.
     $          ((COOR(1,NSFAC(3,IFAC)).GE.0.).AND.
     $              (COOR(1,NSFAC(3,IFAC)).LE.2.)))) THEN
                  LOGFAC(IFAC) = 7
               ELSE
                  LOGFAC(IFAC) = 2
               ENDIF
            ENDIF
         END DO
c
         DO is=1,ns
            logfr(is)              = 0            
         end do
c
         do ifac = 1,nfac
            if ((logfac(ifac).eq.2).or.(logfac(ifac).eq.7)) then
               logfr(nsfac(1,ifac))   = 2
               logfr(nsfac(2,ifac))   = 2
               logfr(nsfac(3,ifac))   = 2
            endif
         end do
c
         do ifac = 1,nfac
                  if (logfac(ifac).eq.4) then
                     logfr(nsfac(1,ifac))   = 4
                     logfr(nsfac(2,ifac))   = 4
                     logfr(nsfac(3,ifac))   = 4
                  endif
         end do
c
      ENDIF
C

      

      CALL CALFRO

      CALL SEG3D
c
      DO 105 is=1,ns
         ndeg(is)                  = 0
105   CONTINUE  
c
c     Forms the list of edges attached to each vertex 
c

      DO 5 iseg=1,nseg
c
         is1                       = nubo(1,iseg)
         is2                       = nubo(2,iseg)
         ndeg(is1)                 = ndeg(is1) + 1
         IF (ndeg(is1) .GT. ndmax) GOTO 15
         ndeg(is2)                 = ndeg(is2) + 1
         IF (ndeg(is2) .GT. ndmax) GOTO 15
         jaret(is1,ndeg(is1))      = iseg
         jaret(is2,ndeg(is2))      = iseg

5     CONTINUE
c
      
      GOTO 1010
c
15    CONTINUE
c
      WRITE(6, *) 'Error in MAIL3D'
      WRITE(6, *) 'Increase NDMAX : ', ndmax
c
      CALL TILT
c
1010  CONTINUE
c
      CALL REPT3D
c
      DO 10 it=1,nt
         idrap(it)                 = 0
10    CONTINUE  
c
      DO 20 ifac=1,nfac
c
         is                        = nsfac(1,ifac)
c
         IF ((fv(is) .EQ. 3) .OR. (fv(is) .EQ. 5) .OR.  
     &       (fv(is) .EQ. 7)) THEN 
            DO 30 it=1,nbvoi(is)
               idrap(ivoi(is,it))  = 1
30          CONTINUE  
         ENDIF
c
         is                        = nsfac(2,ifac)
c
         IF ((fv(is) .EQ. 3) .OR. (fv(is) .EQ. 5) .OR. 
     &       (fv(is) .EQ. 7)) THEN 
            DO 40 it=1,nbvoi(is)
               idrap(ivoi(is,it))  = 1
40          CONTINUE  
         ENDIF
c
         IF ((fv(is) .EQ. 3) .OR. (fv(is) .EQ. 5) .OR.  
     &       (fv(is) .EQ. 7)) THEN 
            is                     = nsfac(3,ifac)
            DO 50 it=1,nbvoi(is)
               idrap(ivoi(is,it))  = 1
50          CONTINUE  
         ENDIF

20    CONTINUE
  
c      open(10,file='uee')
c      write(10,*) (logfac(ifac), ifac=1,nfac) 
c      close(10)
c      stop
      

      
      END
