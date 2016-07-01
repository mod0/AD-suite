      program mar2gid
      implicit real*8 (a-h,o-z)
c
c******************************************************************
c This program is based in interface.f. It uses .sinus. files
c and marco result files to construct *.msh and *.res files
c to be read by gid.
c It also produces a fort.30 file to be read by vigie...
c******************************************************************
c
       PARAMETER (imax=3,jmax=16000,numax=9,kmax=4,nelemmax=81000,
     &           nfacmax=8400,nnufrmax=9)
       
       REAL*8    coor(imax,jmax),q(jmax,numax)
       INTEGER isfac(jmax), ifold(nfacmax), ns2, nfac2, isold(jmax)
       INTEGER isnew(jmax)
       INTEGER nu(kmax,nelemmax),logfac(nfacmax),nufac(3,nfacmax)

       character*6  gidname
       character*10 gimsh,gires,eltype
       character*4  wopos
       integer
     .   ndime,nnode,npoin,nelem,nnodb,imate,
     .   idime,inode,ipoin,ielem,inodb,
     .   ivolu
       
c       
c*********************************************************************
c imax = dimension en espace
c nsmax = nombre maximal de noeuds
c numax = nombre maximal d'inconnues a l'interieur du domaine
c kmax = nombre maximal de noeuds par element
c nelemmax = nombre maximal d'elements a l'interieur du domaine
c nfacmax = nombre maximal d'elements sur la surface
c nnufrmax = nombre maximal d'inconnues sur la surface
c*********************************************************************
c      
c*********************************************************************
c Donnees utiles pour le format VIGIE
c*********************************************************************
c
        nblock=1
        itype=0
        nnu=9
        nnufr=9
c        nnu=8
c        nnufr=8
        npatch=0
        nelemtypes=1
        ndim=3
c
c*********************************************************************
c Lecture du maillage
c*********************************************************************
c
        unit1=60
c             
        OPEN(unit1)
c
        READ(unit1,*) nnodes,nelem,nfac
        write (6,*)nnodes,nelem,nfac
c
        READ(unit1,*) ((coor(i,is),i=1,3),is=1,nnodes)
c
        READ(unit1,*) ((nu(i,ie),i=1,4),ie=1,nelem)
c
        READ(unit1,*) (logfac(ifac),ifac=1,nfac)
c
        READ(unit1,*) ((nufac(i,ifac),i=1,3),ifac=1,nfac)
c
        PRINT*, 'Ns = ', nnodes
        PRINT*, 'Nt = ', nelem
        PRINT*, 'Nf = ', nfac
c       
        close(60)
c
c     Computing the coordinates extrema
c
      xmin                         = 1.0e+16
      ymin                         = 1.0e+16
      zmin                         = 1.0e+16
      xmax                         =-1.0e+16
      ymax                         =-1.0e+16
      zmax                         =-1.0e+16
c
      DO 100 is=1,nnodes
         xmin                      = MIN(xmin, coor(1,is)) 
         ymin                      = MIN(ymin, coor(2,is)) 
         zmin                      = MIN(zmin, coor(3,is)) 
         xmax                      = MAX(xmax, coor(1,is)) 
         ymax                      = MAX(ymax, coor(2,is)) 
         zmax                      = MAX(zmax, coor(3,is)) 
100   CONTINUE 
c
      WRITE(6, *) xmin, xmax
      WRITE(6, *) ymin, ymax
      WRITE(6, *) zmin, zmax
c
c**** Cas de l'aile : on cherche a extraire l'aile du maillage 3D
c     Elle est reperee pour logfac=-2
c
      DO 90 ifac = 1,nfac
c
         IF (logfac(ifac) .EQ. 2) THEN
c
            IF ((coor(3,nufac(1,ifac)) .GT. zmin) .OR. 
     &           (coor(3,nufac(2,ifac)) .GT. zmin) .OR. 
     &           (coor(3,nufac(3,ifac)) .GT. zmin)) THEN 
c
               logfac(ifac)        =-2
c
            ENDIF
c
         ENDIF
c
 90   CONTINUE
c
c*********************************************************************
c Recherche du nombre d'elements sur la peau de l'avion : nfac2
c*********************************************************************
c
        nfac2 = 0
c
        do ifac = 1,nfac
           if ((logfac(ifac).eq.200).or.(logfac(ifac).eq.-2).or.
     $          (logfac(ifac).eq.7)) then
c
            nfac2 = nfac2+1
c
            ifold(nfac2)           = ifac
c
            isfac(nufac(1,ifac))   = 1
            isfac(nufac(2,ifac))   = 1     
            isfac(nufac(3,ifac))   = 1
c
         endif
c
         end do
c
c*********************************************************************
c Recherche du nombre de noeuds sur la peau de l'avion : ns2
c*********************************************************************
c
      ns2                          = 0
c
      DO 40 is=1,nnodes
         IF (isfac(is) .EQ. 1) THEN 
            ns2                    = ns2 + 1
            isold(ns2)             = is
            isnew(is)              = ns2
         ENDIF
40    CONTINUE
c
      PRINT*,' '
      PRINT*,'Ns2   = ', ns2
      PRINT*,'Nfac2 = ', nfac2
c  
c*********************************************************************
c Lecture de la solution q(ns,9) dans le fichier fort.20
c     --> masse volumique
c     --> vitesse suivant x
c     --> vitesse suivant y
c     --> vitesse suivant z
c     --> Energie
c     --> Pression
c     --> Nombre de Mach
c     --> Distributions de pression Cp
c     --> Entropie
c*********************************************************************
c
c        nunit=20
c        OPEN(nunit)
cc
c        READ(nunit,113) ((q(i,n),i=1,nnodes),n=1,nnu)
cc
cc      do n=1,nnu
cc         do i=1,nnodes
cc            q(i,n) = 0.
cc         end do
cc      end do
cc
c113   FORMAT(e15.8)
c
c********************************************************************
c ECRITURE DANS FORT.30 DES DONNEES AU FORMAT HDB
c********************************************************************

      
      
         nnodesfr=nfac

         write(6,'(a)') 'Problem name: '
         read (5,*) gidname


c...     Gid msh file:         
c

         gimsh= gidname//'.msh'
         gires= gidname//'.res'
         open (10,file=gimsh,status='unknown')
         open (20,file=gires,status='unknown')

c...     Volume mesh and general coordinates:

         
         ndime=imax
         npoin=nnodes
         imate=1
         nnode=4
         eltype='Tetrahedra'           
         
         write(10,*)
     .     'MESH ',gidname,' dimension ',ndime,
     .     ' Elemtype ',eltype, ' Nnode ',nnode
         
         write(10,*) 'coordinates'
         do ipoin=1,npoin
           write(10,*) ipoin,(coor(idime,ipoin),idime=1,3)
         end do
         write(10,*) 'end coordinates'
         write(10,*) 'elements'

         do ielem=1,nelem
           write(10,*) ielem,(nu(inode,ielem),inode=1,nnode),imate            
         end do
         write(10,*) 'end elements'
         
c...     Surface mesh:

c         imate=2
c         nnodb=3
c         eltype='Triangle'           
c
c         write(10,*)
c     .     'MESH ',' dimension ',ndime,
c     .     ' Elemtype ',eltype, ' Nnode ',nnodb
c         
c
c         write(10,*) 'coordinates'
c         write(10,*) 'end coordinates'
c         write(10,*) 'elements'
c         do ifac=1,nfac
c           write(10,*) ifac,(nufac(inodb,ifac),inodb=1,nnodb),imate            
c         end do
c         write(10,*) 'end elements'
c
         close (10)         
         
         END
