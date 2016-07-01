



      SUBROUTINE CONDDIR
C
c     Introduction des conditions de Dirichlet sur les noeuds a l'infini.
c
      include 'Param3D.h'
C
      INTEGER IFAC, i, is
C
      Do ifac = 1,nfac
         if (logfac(ifac).eq.4) then
            do i = 1,3
               logfr(nsfac(i,ifac)) = 1
            end do
         endif
      End do
C
      RETURN
      END
