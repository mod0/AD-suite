      SUBROUTINE AIRAGNORMAG
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER IS, I, IZ, NZ, IZT, IZG, NVL, IZGT,isp
C
c-------------------------------------------------------------------------
c..........initialisation des aires :                                    |
c                                                                        |
c             gairesp(is)=airesp(is) ---> aires des cellules is, is=1,nsp|
c             gairesp(iz)=0        ---> aires des zones iz=is+1,nzma1    |
c                                   (initialisation a zero)              |
c                                                                        |
c.........initialisation des normales aux celles grossieres              |
c                                                                        |
c             gvncoq(3,is)=vncoq(3,is) ---> normales aux cellules fines  |
c             gvncoq(3,iz)=0          ---> normales aux cellules         |
c                                         grossieres                     |
c-------------------------------------------------------------------------
c
c     Initialisation des normales grossieres GVNCOQ et des aires grossieres 
c     GAIRESP
c
      do 7 isp=1,nsp
         gairesp(isp)=airesp0(isp)
         do i=1,3
            gvncoq(i,isp)=vno(i,isp)
         end do
7     continue
c
      do 8 iz=nsp+1,nzma1
         gairesp(iz)=0.
         do i=1,3
            gvncoq(i,iz)=0.
         end do
8     continue
C
      nz=0
C
      DO 1000 NVL = 1,NIVG
C
c----------------------------------------------------------------------
c    CALCUL DES AIRES ET DE NORMALES DES ZONES DU NIVEAU NVL+1        |
c    =========================================================        |
c                                                                     |
c       Pour iz=1,nbz(nvl+1) :                                        |
c                                                                     |
c           gairesp(nz+nbz(nvl)+iz)=somme des aires des zones iz du   |
c                                niveau nvl incluses dans iz          |
c           gvncoq(3,nz+nbz(nvl)+iz)=somme des normales des zones iz  |
c                                   du niveau nvl incluses dans iz    |
c----------------------------------------------------------------------
c
         do 250 iz=1,nbz(nvl)
            izt=nz+iz
            izg=nuzo(izt)
            izgt=nz+nbz(nvl)+izg
            gairesp(izgt)=gairesp(izgt)+gairesp(izt)
            do i=1,3
               gvncoq(i,izgt)=gvncoq(i,izgt)+gvncoq(i,izt)
            end do
250      continue
c
         nz=nz+nbz(nvl)
C
1000  CONTINUE
c
      RETURN
      END
