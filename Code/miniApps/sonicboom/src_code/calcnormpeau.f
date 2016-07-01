      SUBROUTINE CALCNORMPEAU(VN_CELL,COORD)
C
C*** Cette procedure calcule :
C     * Les normales en chaque triangle de la coque (stockees dans vnp(3,ntp),
c       tableau mis en COMMON)
C     * Les normales en chaque cellule duale de la coque (stockees dans
C       vn_cell(3,nsp))
C
C*** ATTENTION : sur le maillage 3D, les normales sont orientees vers 
C                l'interieur de la coque. En ne considerant que la coque, on
C                aurait tendance a les orienter vers l'exterieur. Pour rester
C                rationnel, on les oriente vers l'interieur (meme sens que
C                les Vnfac).
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c
      real*8 x(3), y(3), z(3), vnx, vny, vnz
      integer jtp, l, isp, i
      real*8 somment(3), sommeni(3)
      real*8 vn_cell(3,nnsp), coord(3,nnsp)
c
      real*8 airx, airy, airz, airfac, diff
      real*8  vx, vy, vz, psca1, psca2
      real*8  sumx, sumy, sumz
c
c*** 1ere etape : calcul des normales en chaque triangle
c
C$AD II-LOOP
      DO 1 jtp = 1,ntp
         Do 2 l = 1,3
            isp = nup(l,jtp)
            x(l)=coord(1,isp)
            y(l)=coord(2,isp)
            z(l)=coord(3,isp)
2        CONTINUE

c     Calcul du produit vectoriel entre 2 vecteurs-aretes d'un triangle
c     La norme du vecteur resultant est l'aire du triangle.
c                         --->    --->    
c     Produit vectoriel : A1A2 /\ A1A3
c
         vnx= - y(1)*(z(3)-z(2)) + y(2)*(z(3)-z(1)) - y(3)*(z(2)-z(1))
C
         vny= + x(1)*(z(3)-z(2)) - x(2)*(z(3)-z(1)) + x(3)*(z(2)-z(1))
C
         vnz= - x(1)*(y(3)-y(2)) + x(2)*(y(3)-y(1)) - x(3)*(y(2)-y(1))

         vnp(1,jtp)= - 0.5*vnx
         vnp(2,jtp)= - 0.5*vny
         vnp(3,jtp)= - 0.5*vnz
c
c*** On ne passe pas par les verifications, on va directement au calcul
c    de la normale en les noeuds de la coque.
c
         goto 1000
c
c         open(90)
c
C    ON VERIFIE QUE C'EST NORMAL AUX COTES DE LA FACE :
C    --------------------------------------------------
C
c         VX= X(2) - X(1)
c         VY= Y(2) - Y(1)
c         VZ= Z(2) - Z(1)
c         PSCA1= VNP(1,jtp)*VX + VNP(2,jtp)*VY + VNP(3,jtp)*VZ
c
c         if (PSCA1.gt.1.e-07) then
c         write(90,*) 'Produit scalaire Normale.A1A2', jtp, ' ', PSCA1
c         endif
C
c         VX= X(3) - X(1)
c         VY= Y(3) - Y(1)
c         VZ= Z(3) - Z(1)
c         PSCA2= VNP(1,jtp)*VX + VNP(2,jtp)*VY + VNP(3,jtp)*VZ
c
c         if (PSCA2.gt.1.e-07) then         
c         write(90,*) 'Produit scalaire Normale.A1A3', jtp, ' ', PSCA2
c         endif
C
C    ON VERIFIE QUE LA NORME EST BIEN LA SURFACE DE LA FACE :
C    --------------------------------------------------------
C
c         AIRTP(JTP)= SQRT(VNP(1,jtp)**2 + VNP(2,jtp)**2 + VNP(3,jtp)**2)
c
c         write(90,*) 'Calcul de la norme :',jt,' ', AIRTP(JT)
C
c         AIRZ= 0.5*( X(1)*Y(2) - X(2)*Y(1) + X(2)*Y(3) -
c     1          X(3)*Y(2) + X(3)*Y(1) - X(1)*Y(3) )
c         AIRX= 0.5*( Z(1)*Y(2) - Z(2)*Y(1) + Z(2)*Y(3) -
c     1          Z(3)*Y(2) + Z(3)*Y(1) - Z(1)*Y(3) )
c         AIRY= 0.5*( X(1)*Z(2) - X(2)*Z(1) + X(2)*Z(3) -
c     1          X(3)*Z(2) + X(3)*Z(1) - X(1)*Z(3) )
c         AIRFAC= SQRT( AIRX*AIRX + AIRY*AIRY + AIRZ*AIRZ )
c
c         write(90,*) 'Calcul de la surface du triangle :',jtp,' ', AIRFAC
C
c         DIFF= AIRFAC - AIRTP(JTP)
c
c         if (DIFF.gt.1.e-07) then
c         write(90,*) 'Calcul de la difference :', jtp, ' ', DIFF
c         endif
c
 1000   continue
c
1       CONTINUE 
C
C  INTEGRALE DES NORMALES EXTERIEURES
C  ----------------------------------
c
c      SUMX= 0.
c      SUMY= 0.
c      SUMZ= 0.
c      DO 300 jtp = 1 , ntp
c         SUMX= SUMX + VNP(1,jtp)
c         SUMY= SUMY + VNP(2,jtp)
c         SUMZ= SUMZ + VNP(3,jtp)
c300   CONTINUE
c
c      if (sumx.gt.1.e-07) write(90,*) 'sumx = ',sumx
c      if (sumy.gt.1.e-07) write(90,*) 'sumy =',sumy
c      if (sumz.gt.1.e-01) write(90,*) 'sumz =',sumz

c      WRITE( 90, * )'       INT. DES NORMALES EXT. : ',
c     &             SUMX, ' ', SUMY, ' ', SUMZ
c
c
c**** 2eme etape : on deduit vn_cell
c
c
c     newnormale(i) = (Somme sur les triangles ayant pour sommet i des
c                     normales a ces triangles)/3. 
c
c     Initialisation a 0 des nouvelles normales vnorp
c     -----------------------------------------------
c

c
C$AD II-LOOP
        DO 4 i = 1,3
           DO 3 isp = 1,nsp
              vn_cell(i,isp)=0.
 3         CONTINUE
 4      CONTINUE
c
C$AD II-LOOP
      DO 5 l = 1,3
         DO 101 jtp = 1,ntp
            vn_cell(1,nup(l,jtp)) = vn_cell(1,nup(l,jtp)) + vnp(1,jtp)/3.
            vn_cell(2,nup(l,jtp)) = vn_cell(2,nup(l,jtp)) + vnp(2,jtp)/3.
            vn_cell(3,nup(l,jtp)) = vn_cell(3,nup(l,jtp)) + vnp(3,jtp)/3.
 101     CONTINUE
 5    CONTINUE

c
c     Test sur le fait que la somme des normales aux triangles sur tous les
c     triangles doit etre nulle, ainsi que la somme des normales en chaque
c     noeud.
c
c      DO 7 i = 1,3
c         somment(i)=0
c         sommeni(i)=0
c7     CONTINUE
c
c      DO 8 jtp = 1,ntp
c         DO 9 i = 1,3
c           somment(i) = somment(i) + vnp(i,jtp)
c9        CONTINUE
c8     CONTINUE
c
c      DO 10 isp = 1,nsp
c         DO 11 i = 1,3
c           sommeni(i) = sommeni(i) + vn_cell(i,isp) 
c11       CONTINUE
c10    CONTINUE
c
c      do i = 1,2
c       if (sommeni(i).gt.1.e-07) write(90,*) 'somme is',sommeni(i)
c      end do
c      if (sommeni(3).gt.1.-01) write(90,*) 'somme jt',sommeni(3)

c      write(90,*) 'Somme des vncoq =',(sommeni(i),i=1,3)
c      write(6,*) 'Somme des normales aux triangles ='
c      write(6,*) (somment(i),i=1,3)
c
c REMARQUE : si nous travaillons sur l'aile M6, la somme en z ne fait
c            pas 0, car la geometrie n'est pas fermee et le plan de symetrie
c            est en z.
c            si nous travaillons sur le falcon, la somme en y ne fait
c            pas 0, car la geometrie n'est pas fermee et le plan de symetrie
c            est en y.
c            
      return
      end
