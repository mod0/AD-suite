

      SUBROUTINE LISPRES(Presref)
C
c     Dans le cas ou l'on veut calculer la pression desiree par lissage
c     NPRES=1 et sans transpiration.
c
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      REAL*8 D
      integer izt, nv
      integer izv
      real*8 sompoids, poids, poids1
      real*8 scal, Presref(nnsp)
C
C     Entree : Pref(1:nsp).
C     ------   
C              
C     Sortie : tableau PdesP(nsp) = LPref
C     ------
C
C     L est un operateur de moyenne, pondere par des poids P_ij tels que :
C                -> ->       ->       ->
c     P_ij = max(Ni.Nj,0) ou Ni (resp Nj) est la normale a la cellule i 
c     (resp j). Les normales ont ete stockees dans gvncoq, dans la procedure
c     Zonag.f
c      
c             **********************************
c             ****   Calcul de LPref      ****
c             **********************************
c
         Do 95 izt = 1,nsp
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
               d = Presref(izt)*poids
            Do 110 nv = 1,nvoi(izt)
               izv = nuvoi(izt,nv) 
               scal = gvncoq(1,izt)*gvncoq(1,izv) + 
     &                 gvncoq(2,izt)*gvncoq(2,izv) +
     &                 gvncoq(3,izt)*gvncoq(3,izv)
               poids = max(scal,0.d0)
               sompoids = sompoids + poids
                  d = d + Presref(izv)*poids
110         Continue
c
            if (sompoids.ne.0.) then
               d = d/sompoids
               deltanew(izt) = (1.-theta)*Presref(izt) + theta*d
            else
            print*,'Attention a la cellule ',izt,sompoids
            endif
c
95    Continue
c
         Do 160 izt = 1,nsp
               Presref(izt) = deltanew(izt)
160     Continue
c
      RETURN
      END
