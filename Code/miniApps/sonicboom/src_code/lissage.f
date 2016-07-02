

      SUBROUTINE LISSAGE(kniv)
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      REAL*8 D
      integer i, izt, iz, nzt, nzg, nz, nv, kniv
      integer nvl, mv, nov, mov, izv, izgt
      real*8 sompoids, poids, poids1
      real*8 scal
C
C     Entree : Grad(1:nsp).
C     ------   On lui fait correspondre un tableau Delta(nzma1), car on
C              a besoin d'aller sur tous les niveaux.
C              Kniv=niveau sur lequel on va travailler.
C     Sortie : tableau delta(nzma1) = LPP*L*Grad
C     ------
C     But de la procedure: Calculer LPP*L*Grad
C     --------------------
C
C     L est un operateur de moyenne, pondere par des poids P_ij tels que :
C                -> ->       ->       ->
c     P_ij = max(Ni.Nj,0) ou Ni (resp Nj) est la normale a la cellule i 
c     (resp j). Les normales ont ete stockees dans gvncoq, dans la procedure
c     Zonag.f
c      
c     L* est son adjoint calcule avec la norme L2.
c
c     P est l'operateur de prolongement de la grille grossiere vers la grille
c     fine, sous la forme d'une injection canonique.
c      
c     P* est l'adjoint de P dans L2, il est pondere par les aires (gairesp) 
c     des cellules.
c
      nzt = 0
c
      do i = 1, nivg
         nzt = nzt + nbz(i)
      end do
c
      Do iz = 1,nbz(niveau)
            delta(iz) = Grad(iz)
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
c             ********************************
c             ****   Calcul de P*L*Grad   ****
c             ********************************
c
      DO 100 NVL = 1,KNIV - 1
c                                                  
         DO 50 izt = 1+nz,nbz(nvl)+nz ! Calcul de L*Grad
c
            d = 0.
c
c**** On calcule d'abord :
c                         _____
c                         \        P_ij delta_j
c                          \       ------------
c                           \       __
c                           /       \   Pjk
c                          /        /
c                         /         --
c                         ----     k =V(j)U{j}
c                        j = V(i)
c
            do 40 nv = 1,nvoi(izt)
               nov = nuvoi(izt,nv) + nz
               scal = gvncoq(1,nov)**2 + gvncoq(2,nov)**2 +
     $                          gvncoq(3,nov)**2  
               poids = max(scal,0.d0)
               sompoids = poids
               do 30 mv = 1,nvoi(nov)
                  mov=nuvoi(nov,mv)+nz
                  scal = gvncoq(1,nov)*gvncoq(1,mov) + gvncoq(2,nov)*
     $             gvncoq(2,mov) + gvncoq(3,nov)*gvncoq(3,mov) 
                  poids = max(scal,0.d0)
                  sompoids = sompoids + poids
30             continue
c
               scal = gvncoq(1,nov)*gvncoq(1,izt) + gvncoq(2,nov)*
     $             gvncoq(2,izt) + gvncoq(3,nov)*gvncoq(3,izt)
               poids = max(scal,0.d0)
               if (sompoids.ne.0.) then
                  d = d + delta(nov)*poids/sompoids
               else
                  print*,'Attention a la cellule ',nov,sompoids
               endif
c
40          continue
c
c**** On calcule ensuite :
c
c             P_ii delta_i
c             ------------
c               __
c               \   P_ij
c               /   
c               --
c             k=V(i)U{i}
c
            scal = gvncoq(1,izt)**2 + gvncoq(2,izt)**2 +
     $                          gvncoq(3,izt)**2
            poids1 = max(scal,0.d0)
            sompoids = poids1
            do 10 nv = 1,nvoi(izt)
               nov = nuvoi(izt,nv)+nz
               scal = gvncoq(1,nov)*gvncoq(1,izt) + gvncoq(2,nov)*
     $             gvncoq(2,izt) + gvncoq(3,nov)*gvncoq(3,izt)
               poids = max(scal,0.d0)
               sompoids = sompoids + poids
10          continue
c
c**** D'ou finalement :
c                                               _____
c                                               \        P_ij delta_j
c                                                \       ------------
c                                                 \       __
c                                                 /       \   Pjk
c   deltanew_i = (1-theta)* delta_i + theta *    /        /
c                                               /         --
c                                               ----     k =V(j)U{j}
c                                             j=V(i)U{i}
c
            if (sompoids.ne.0.) then
               d = d + delta(izt)*poids1/sompoids
               deltanew(izt) = (1.-theta)*delta(izt) + theta*d
5           continue
            else
            print*,'Attention a la cellule ',izt,sompoids
            endif
c
50       continue
c
         Do 15 izt = 1+nz,nbz(nvl)+nz
               delta(izt) = deltanew(izt)
15       Continue
c
         Do 35 izt = 1+nz,nz+nbz(nvl) ! Calcul de P*L*Grad
            izgt = nzg+nuzo(izt)
               delta(izgt) = delta(izgt) + delta(izt)*gairesp(izt)
35        Continue
c
          Do 55 izgt = 1+nzg,nzg+nbz(nvl+1)
                delta(izgt) = delta(izgt)/gairesp(izgt)
55        Continue
c
          nz = nz + nbz(nvl)
          nzg = nzg + nbz(nvl+1)
c
100    CONTINUE
c
c             **********************************
c             ****   Calcul de LPP*L*Grad   ****
c             **********************************
c
       DO 200 NVL = KNIV-1,1,-1
c
          nz = nz - nbz(nvl)
          nzg = nzg - nbz(nvl+1)
c
          Do 75 izt = 1+nz,nz+nbz(nvl) ! Calcul de PP*L*Grad
             izgt = nzg +nuzo(izt)
                delta(izt) = delta(izgt)
75        Continue
c
         Do 95 izt = 1+nz,nz + nbz(nvl) ! Calcul de LPP*L*Grad
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
            poids = max(scal,0.d0)
            sompoids = poids
               d = delta(izt)*poids
            Do 110 nv = 1,nvoi(izt)
               izv = nuvoi(izt,nv) + nz
               scal = gvncoq(1,izt)*gvncoq(1,izv) + 
     &                 gvncoq(2,izt)*gvncoq(2,izv) +
     &                 gvncoq(3,izt)*gvncoq(3,izv)
               poids = max(scal,0.d0)
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
