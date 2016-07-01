      SUBROUTINE PRESCIBLE(CTRL)
C
c     Cette procedure calcule la pression desiree dans le cas de l'aile
c     M6 et du Falcon, mais pas dans le cas de la tuyere 3D.
c
c     NPRES=1 --> Lissage de la pression calculee sans transpiration
c     NPRES=0 --> Calcul de la pression par transpiration simulant une
c                 incidence.
c
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c
      INTEGER ISP, IS, nunit
      REAL*8 presref(nnsp), ctrl(nnsp)
c
      IF (NPRES.EQ.1) THEN !! Il aura fallu mettre ITRANS=0
C
c     On calcule la pression (Presref) et on la lisse (Pdesp = L.Presref)
c     ITRANS = 0
c
         CALL ETAT(Ctrl,ktmax)
c
         do isp = 1,nsp
            is = node2D3D(isp)
            Presref(isp) = gam1*(ua(5,is) - 0.5*(ua(2,is)**2 +
     $           ua(3,is)**2 + ua(4,is)**2)/ua(1,is))
         end do
c
         CALL LISPRES(Presref)
c
         do isp = 1,nsp
            is = node2D3D(isp)
            Pdesp(is) = Presref(isp)
         end do   
c     
c     Ecriture dans fort.22 de la pression desiree afin qu'elle
c     soit visualisee avec VIGIE.
c
c      nunit = 27
c      CALL ResPres(nunit)
C
      ELSE
c
c     Calcul de la pression desiree avec une incidence
c     simulee par transpiration (0.1 DEGRES) : ITRANS=1, contr=0, npres=0
c
         CALL ETAT(Ctrl,ktmax)
c
         do isp = 1,nsp
            is = node2D3D(isp)
            Pdesp(is) = gam1*(ua(5,is) - 0.5*(ua(2,is)**2 +
     $           ua(3,is)**2 + ua(4,is)**2)/ua(1,is))
         end do
c
      ENDIF
C
c      CALL ECRITVIGIE
c
      if (npres.eq.0) npres = 1
      if (itrans.eq.0) itrans = 1
c
c      print*,'npres = ',npres,' itrans = ',itrans
c
c      stop
c
      RETURN
      END
