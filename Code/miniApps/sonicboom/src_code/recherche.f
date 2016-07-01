
      SUBROUTINE RECHERCHE
C
c*** Verification de l'orientation des normales aux faces
c
      include 'Param3D.h'
c
      INTEGER JT,IFAC,KS,ISOPP,NBR,NSSAUVE,NLOC,K
      REAL*8 GX,GY,GZ,GPX,GPY,GPZ,PRODSCAL
C
      NBR=0
      NLOC=0
c
c     NLOC = Compte les faces du maillage
c     NBR = Compte les faces dont les normales sont mal orientees
C
      DO JT = 1,NT
c
c     Recherche de tetraedres ayant une face frontiere
c
         IF ((LOGFR(NU(1,JT)).GT.0).OR.
     $       (LOGFR(NU(2,JT)).GT.0).OR.
     $       (LOGFR(NU(3,JT)).GT.0).OR.
     $       (LOGFR(NU(4,JT)).GT.0)) THEN 
            DO IFAC = 1,NFAC
c
c     Recherche des noeuds du tetraedre se trouvant sur la frontiere
c
               IF ((NSFAC(1,IFAC).EQ.NU(1,JT)).OR.
     $             (NSFAC(1,IFAC).EQ.NU(2,JT)).OR.
     $             (NSFAC(1,IFAC).EQ.NU(3,JT)).OR.
     $             (NSFAC(1,IFAC).EQ.NU(4,JT))) THEN
                  IF ((NSFAC(2,IFAC).EQ.NU(1,JT)).OR.
     $             (NSFAC(2,IFAC).EQ.NU(2,JT)).OR.
     $             (NSFAC(2,IFAC).EQ.NU(3,JT)).OR.
     $             (NSFAC(2,IFAC).EQ.NU(4,JT))) THEN
                     IF ((NSFAC(3,IFAC).EQ.NU(1,JT)).OR.
     $             (NSFAC(3,IFAC).EQ.NU(2,JT)).OR.
     $             (NSFAC(3,IFAC).EQ.NU(3,JT)).OR.
     $             (NSFAC(3,IFAC).EQ.NU(4,JT))) THEN
C
                        NLOC = NLOC +1
                        ISOPP = 0
C
                        DO KS = 1,4
c
c     Recherche du noeud ne se trouvant pas sur le bord
c
                           IF ((NU(KS,JT).NE.NSFAC(1,IFAC)).AND.
     $                          (NU(KS,JT).NE.NSFAC(2,IFAC)).AND.
     $                          (NU(KS,JT).NE.NSFAC(3,IFAC))) THEN
                              ISOPP = NU(KS,JT)
                           ENDIF
                        END DO
C
                        IF (ISOPP.EQ.0) PRINT*,'ATTENTION'
C
c     Calcul des coordonnees du barycentre de la face frontiere
c
                        GX = (COOR(1,NSFAC(1,IFAC)) +
     $                        COOR(1,NSFAC(2,IFAC)) +
     $                        COOR(1,NSFAC(3,IFAC)))/3.D0
C
                        GY = (COOR(2,NSFAC(1,IFAC)) +
     $                        COOR(2,NSFAC(2,IFAC)) +
     $                        COOR(2,NSFAC(3,IFAC)))/3.D0
C
                        GZ = (COOR(3,NSFAC(1,IFAC)) +
     $                        COOR(3,NSFAC(2,IFAC)) +
     $                        COOR(3,NSFAC(3,IFAC)))/3.D0
C
c                       ------------------------->
c     Calcul du vecteur Bar,Noeud pas sur la face
c
                        GPX = COOR(1,ISOPP) - GX
                        GPY = COOR(2,ISOPP) - GY
                        GPZ = COOR(3,ISOPP) - GZ
C
c                                ----------------->  ------------------------>
c     Calcul du produit scalaire Normale a la face . Bar,Noeud pas sur la face
c
                        PRODSCAL = VNFAC(1,IFAC)*GPX +
     $                     VNFAC(2,IFAC)*GPY + VNFAC(3,IFAC)*GPZ
C
c     Si le produit scalaire est positif, c'est que la normale est mal dirigee.
c     On inverse alors 2 noeuds de la face.
c
               IF (PRODSCAL.GT.0) THEN
                           NBR = NBR + 1
c                  PRINT*,'NBR = ',NBR,' IFAC = ',IFAC,' JT = ',JT
                  NSSAUVE = NSFAC(2,IFAC)
                  NSFAC(2,IFAC) = NSFAC(3,IFAC)
                  NSFAC(3,IFAC) = NSSAUVE
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      END DO
      ENDIF     
      END DO
C
C [llh]      PRINT*,'NLOC = ',NLOC,' NBR = ',NBR
C
      RETURN
      END
