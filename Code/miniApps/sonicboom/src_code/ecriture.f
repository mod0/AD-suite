      SUBROUTINE ECRITURE(ctrl,nnunit)
C
C*** Ecriture au format ZVIS de la nouvelle coque.
C
      INCLUDE 'Paramopt3D.h'
C
      INTEGER IS, JT, I1, I2, I3, nnunit, isp, i
      REAL*8 XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN, COORPBAR(3,NNSP)
      REAL*8 ctrl(nnsp)

C

C [llh]      write(6,*)
C [llh]     .  'Subroutine ECRITURE returned.. writting only in GID files'
      
      return

C [llh]      PRINT*,'Ecriture.f: nnunit=',nnunit
      do isp = 1,nsp
         do i = 1,3
            coorpbar(i,isp) = coorp(i,isp) + ctrl(isp)*vnon(i,isp)
         end do
      end do


      XMAX = - 1000000
      XMIN =   1000000
      YMAX = - 1000000
      YMIN =   1000000
      ZMAX = - 1000000
      ZMIN =   1000000
c
      DO 30 IS= 1,NSP
             XMIN = MIN( XMIN,COORPBAR(1,IS))
             XMAX = MAX( XMAX,COORPBAR(1,IS))
             YMIN = MIN( YMIN,COORPBAR(2,IS))
             YMAX = MAX( YMAX,COORPBAR(2,IS))
             ZMIN = MIN( ZMIN,COORPBAR(3,IS))
             ZMAX = MAX( ZMAX,COORPBAR(3,IS))
30    CONTINUE
C
      PRINT*,'Ecriture.f:',
     $'CALCUL DES MIN ET DES MAX DES NOUVELLES COORDONNEES'
      PRINT*,' DE LA COQUE :'
      WRITE(6, *) ' '
      PRINT*, ' MAXI X = ',XMAX, ' MINI X = ',XMIN
      PRINT*, ' MAXI Y = ',YMAX, ' MINI Y = ',YMIN
      PRINT*, ' MAXI Z = ',ZMAX, ' MINI Z = ',ZMIN
      WRITE(6, *) ' '
C
c

      xman= xmin
      
      write(6,*)
     .  'Subroutine ECRITURE returned.. writting only in GID files'
      
      return
      

      nnunit=37
      print*,'Ecriture.f: sauvegarde sur fort.',nnunit
c
      OPEN(unit=nnunit,ACCESS='APPEND')
            REWIND(nnunit)
            WRITE(nnunit,*) ' couleurs '
            WRITE(nnunit,*) ' 4 '
            WRITE(nnunit,*) '  0   .7   .8   1.  '
            WRITE(nnunit,*) '  1   .8   .9   1.  '
            WRITE(nnunit,*) '  2   0.   1.   1.  '
            WRITE(nnunit,*) '  3   0.   0.   1.  '
            WRITE(nnunit,*) ' 2 '
            WRITE(nnunit,*) ' 0 '
            WRITE(nnunit,*) '   4  0. 0. 0. '
            WRITE(nnunit,*) ' 128  1. 1. 1. '
            WRITE(nnunit,*) ' 1 '
            WRITE(nnunit,*) ' 129 0. 0. 0. '
            WRITE(nnunit,*) ' 255 1. 0. 0. '
            WRITE(nnunit,*) ' domaine '
            WRITE(nnunit,*)  XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
c
            WRITE(nnunit,*) ' triangles '
            WRITE(nnunit,*) NSP
c
            DO 12 IS=1, NSP
                WRITE(nnunit,78) COORPBAR(1,IS),COORPBAR(2,IS),
     $              COORPBAR(3,IS)
12          CONTINUE

            WRITE(nnunit,*) NTP

            DO 13 JT=1, NTP
              I1 =  NUP(1,JT)-1
              I2 =  NUP(2,JT)-1
              I3 =  NUP(3,JT)-1
               WRITE(nnunit,79) I1,I2,I3
13          CONTINUE
c
        CLOSE(nnunit)
78       FORMAT (3x,E12.6,2x,E12.6,2x,E12.6)
79       FORMAT (I6,2x,I6,2x,I6)
c
      return
      end
