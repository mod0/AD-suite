









      SUBROUTINE LP(kniv,deltaz,delta)
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      REAL*8 D
      integer i, izt, iz, nzt, nzg, nz, nv, kniv
      real*8 nvl, mv, nov, mov, izv, izgt
      real*8 sompoids, poids, poids1
      real*8 scal
C
C     Entree : deltaz(1:nbz(kniv)).
C     ------   On lui fait correspondre un tableau Delta(nzma1), car on
C              a besoin d'aller sur tous les niveaux.
C              Kniv=niveau sur lequel on va travailler.
C     Sortie : tableau delta(nzma1) = LP deltaZ
C     ------
C     But de la procedure: Calculer LP deltaZ
C     --------------------
C
C     L est un operateur de moyenne, pondere par des poids P_ij tels que :
C                -> ->       ->       ->
c     P_ij = max(Ni.Nj,0) ou Ni (resp Nj) est la normale a la cellule i 
c     (resp j). Les normales ont ete stockees dans gvncoq, dans la procedure
c     Zonag.f
c      
c     P est l'operateur de prolongement de la grille grossiere vers la grille
c     fine, sous la forme d'une injection canonique.
c      
      nzt = 0
c
      do i = 1, nivg
         nzt = nzt + nbz(i)
      end do
c
      Do iz = 1,nbz(niveau)
            delta(iz) = 0.
            deltanew(iz) = 0.
      End Do
c
      Do iz = 1+nbz(niveau), nzt
            delta(iz) = 0.
            deltanew(iz) = 0.
      End Do
c
      nz = 0
      nzg = nbz(niveau)
c
c      *****************************************************
c      **** Transfert de DeltaZ vers l'empilement Delta ****
c      *****************************************************
c
      DO 100 NVL = 1,KNIV - 1
c                                                  
          nz = nz + nbz(nvl)
          nzg = nzg + nbz(nvl+1)
c
100    CONTINUE
c
       Do 150 izt = 1+nz,nzg
          delta(izt) = deltaZ(izt-nz)
 150   continue
c
c             **********************************
c             ****   Calcul de LP delta   ****
c             **********************************
c
       DO 200 NVL = KNIV-1,1,-1
c
          nz = nz - nbz(nvl)
          nzg = nzg - nbz(nvl+1)
c
          Do 75 izt = 1+nz,nz+nbz(nvl) ! Calcul de P delta
             izgt = nzg +nuzo(izt)
                delta(izt) = delta(izgt)
75        Continue
c
         Do 95 izt = 1+nz,nz + nbz(nvl) ! Calcul de LP delta
c
c**** Calcul de :
c
c           ___
c           \     p_ij delta_j
c           /
c           --- 
c         j=V(i)U{i}
c        -----------------------
c              ___
c              \     p_ij
c              /
c              ---
c           j=V(i)U{i}
c  
            scal = gvncoq(1,izt)**2 + gvncoq(2,izt)**2
     $           + gvncoq(3,izt)**2
            poids = max(scal,0.)
            sompoids = poids
               d = delta(izt)*poids
            Do 110 nv = 1,nvoi(izt)
               izv = nuvoi(izt,nv) + nz
               scal = gvncoq(1,izt)*gvncoq(1,izv) + 
     &                 gvncoq(2,izt)*gvncoq(2,izv) +
     &                 gvncoq(3,izt)*gvncoq(3,izv)
               poids = max(scal,0.)
               sompoids = sompoids + poids
                  d = d + delta(izv)*poids
110         Continue
c
            if (sompoids.ne.0.) then
               d = d/sompoids
               deltanew(izt) = (1.-theta)*delta(izt) + theta*d
            else
            print*,'Attention a la cellule ',izt,sompoids
            endif
c
95    Continue
c
         Do 160 izt = 1+nz,nz+nbz(nvl)
               delta(izt) = deltanew(izt)
160     Continue
c
200   CONTINUE
c
      RETURN
      END
