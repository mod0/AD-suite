     SUBROUTINE critere(grad,sortie,kniv)
c
c*** Critere de resolution par la methode d additivite
c
c*** Separation des hautes frequences et des basses frequences
c
c*** Exemple avec 2 niveaux :
c    ------------------------
c     n+1   n              1  2            1 2 
c    U   = U - ( rho1 (1- P  P *)g + rho2 P P *g )
c                          2  1            2 1
c
c*** Entree : grad (gradient), kniv=nvlg (niveau grossier)
c*** Sortie : sortie (addition des frequences), g(niveaux) 
c***              (les directions de descente sur chaque niveau)

      include 'paramopt.h'

      dimension grad(nns),solu(nzma2),sortie(nzma2),g(10,nns)

      do iz=1,nbz(niveau)
         if (logfr(iz).eq.0) then
         g(niveau,iz)=grad(iz)
         solu(iz)=grad(iz)
         else
         g(niveau,iz)=0.
         solu(iz)=0.
         endif
         sortie(iz)=0.
      end do

      DO 1 niv=2,kniv
         call evolag(niv,solu)
         do 3 iz=1,nbz(niveau)
            if (logfr(iz).eq.0) then
            g(niv-1,iz)=g(niv-1,iz)-solu(iz)
            g(niv,iz)=solu(iz)
            endif
3        continue
1     CONTINUE

      Do 10 niv=1,kniv
         do 11 iz=1,nbz(niveau)
            if (logfr(iz).eq.0) then
            sortie(iz)=sortie(iz) + rho(niv)*g(niv,iz)
            endif
11       continue
10    Continue

      return
      end
