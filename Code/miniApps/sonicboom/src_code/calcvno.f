

      SUBROUTINE CALCVNO
C
C*** Cette procedure calcule :
C     * Les normales en chaque triangle de la coque (stockees dans vnp(3,ntp),
c       tableau mis en COMMON)
C     * Les normales en chaque cellule duale de la coque initiale (stockees 
c       dans VNO)
C
C*** ATTENTION : sur le maillage 3D, les normales sont orientees vers 
C                l'interieur de la coque. En ne considerant que la coque, on
C                aurait tendance a les orienter vers l'exterieur. Pour rester
C                rationnel, on les oriente vers l'interieur (meme sens que
C                les Vnfac.
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c
      real*8 x(3), y(3), z(3),oldl(3),newl(3),vnnew,vnold,
     .  coorpbar(3,nnsp),vnx, vny, vnz,vnpno,vnpfa,norj,nork
      integer jtp, l, isp, ktp,i,ino
c
c*** 1ere etape : calcul des normales en chaque triangle
c
      DO 1 jtp = 1,ntp
        Do 2 l = 1,3
          isp = nup(l,jtp)
          x(l)=coorp(1,isp)
          y(l)=coorp(2,isp)
          z(l)=coorp(3,isp)
    2   CONTINUE
        
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
    1 CONTINUE 
C
c     Initialisation a 0 des normales vno
c     -----------------------------------
c
      DO  isp = 1,nsp
        DO  i = 1,3
          vno(i,isp)=0.d0
          vnon(i,isp)=0.d0
          coorpbar(i,isp)=0.d0
        end do
      end do
      
      
c
      
      vnpno=1.d0
      vnpfa=1.d0
      DO isp = 1,nsp
        DO jtp = 1,ntp
          DO  l = 1,3
            IF (nup(l,jtp).eq.isp) THEN
              
              vno(1,isp) = vnpno*vno(1,isp) + vnpfa*vnp(1,jtp)/3.d0
              vno(2,isp) = vnpno*vno(2,isp) + vnpfa*vnp(2,jtp)/3.d0
              vno(3,isp) = vnpno*vno(3,isp) + vnpfa*vnp(3,jtp)/3.d0

            ENDIF
          end do
        end do
      end do

      if (intzz(2).eq.0) return                             ! correct normals?

c
c    Correcting conflictive normals:
c
c    1. Acute angles normals
c      
      
      DO isp = 1,nsp
        DO jtp = 1,ntp
          DO l = 1,3
            IF (nup(l,jtp).eq.isp) THEN

              vnpfa= 1.d0
              vnpfa= vnon(1,isp)*vnp(1,jtp)
     .          +vnon(2,isp)*vnp(2,jtp)+vnon(3,isp)*vnp(3,jtp)

              ino=0

              if (vnpfa.lt.0.d0 .and. ino.eq.0) then
                
                vnx= vnpfa
                vnpno=1.d0
                vnpfa=0.d0
                
                if (coorpbar(1,isp).gt.-0.5d0) then
                  
                  ktp=int(coorpbar(1,isp))

                  nork= vnp(1,ktp)*vnp(1,ktp)
     .              +vnp(2,ktp)*vnp(2,ktp)+vnp(3,ktp)*vnp(3,ktp)
                  norj= vnp(1,jtp)*vnp(1,jtp)
     .              +vnp(2,jtp)*vnp(2,jtp)+vnp(3,jtp)*vnp(3,jtp)

                  nork=sqrt(nork)
                  norj=sqrt(norj)

c                  norj=1.d0
c                  nork=1.d0

c                  write (6,*) ktp,jtp,isp,vnx/nork/norj

                  vnon(1,isp) = vnp(1,jtp)/norj+vnp(1,ktp)/nork
                  vnon(2,isp) = vnp(2,jtp)/norj+vnp(2,ktp)/nork
                  vnon(3,isp) = vnp(3,jtp)/norj+vnp(3,ktp)/nork
                  
                  coorpbar(1,isp)= -1.d0
                  
                end if                
                
              else if (coorpbar(1,isp).lt.-0.5d0) then
                vnpno=1.d0
                vnpfa=0.d0
              else
                vnpno=1.d0
                vnpfa=1.d0                                                  
                coorpbar(1,isp)= dble(jtp)
              end if
                
              vnon(1,isp) = vnpno*vnon(1,isp) + vnpfa*vnp(1,jtp)/3.d0
              vnon(2,isp) = vnpno*vnon(2,isp) + vnpfa*vnp(2,jtp)/3.d0
              vnon(3,isp) = vnpno*vnon(3,isp) + vnpfa*vnp(3,jtp)/3.d0

            ENDIF
          end do
        end do
      end do
      
      do isp=1,nsp
        if (coorpbar(1,isp).lt.-0.5d0) then

          oldl(1)=vno(1,isp)
          oldl(2)=vno(2,isp)
          oldl(3)=vno(3,isp)
          vnold= sqrt(oldl(1)*oldl(1)+oldl(2)*oldl(2)+oldl(3)*oldl(3))
          
          newl(1)=vnon(1,isp)
          newl(2)=vnon(2,isp)
          newl(3)=vnon(3,isp)
          vnnew= sqrt(newl(1)*newl(1)+newl(2)*newl(2)+newl(3)*newl(3))

c          vnold=
c     .      oldl(1)*newl(1)+oldl(2)*newl(2)+oldl(3)*newl(3)/vnold
          
          vno(1,isp)= vnold*newl(1)
          vno(2,isp)= vnold*newl(2)
          vno(3,isp)= vnold*newl(3)

          
          
          
        end if
      end do

      
      
      end
