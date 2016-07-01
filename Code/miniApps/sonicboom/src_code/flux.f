C     -----------------------------------------------
      SUBROUTINE FLUX( NN, NNT, NNSG, ORDRE, NU, NUBO, DX, DY, 
     $                 DZ, VNOCL, COOR, UA, GAMMA, FLUEUL )
C     -----------------------------------------------
C
C     La procedure FLUX dirige le calcul des flux 
C     explicites convectifs.
C     Elle  appelle (suivant  le  schema  choisi)  la
C     procedure chargee du  calcul  de la fonction de
C     flux de Van  Leer. Les  indirections
C     sont  resolues au  moment de  l'appel de 
C     cette procedure  (pour  faciliter la derivation 
C     de la fonction de flux). La boucle sur les seg-
C     -ments du maillage  est donc dans  la procedure
C     FLUX, de meme que  l'assemblage  des flux  con-
C     -vectifs.
C
C     La molecule ordre 1 est  tres simple,  elle est
C     constituee des deux extremites du segment.
C
C     Le segment sur lequel  on calcule le  flux  est 
C     le segment [1,2].
C
C     -----------------------------------------------
C     Description de la liste d'appel
C     -----------------------------------------------
C
C     ENTREE :
C
C     NS : Nombre de sommets du maillage.
C     NSEG : Nombre de segments du maillage.
C     ORDRE : Ordre de la methode.
C     NU : Tableau  d'entiers  contenant les indices 
C          des points constituant chaque triangle.
C     NUBO : Tableau d'entiers contenant les indices
C          des points constituant chaque segment.
C     COOR : Tableau de reels  contenant les coordon-
C          -nees (x,y et z) de chacun des sommets du
C          maillage.
C     UA : Tableau de reels contenant les valeurs des
C          inconnues aux points du maillage.
C     GAMMA : Rapport  des  chaleurs   specifiques  a
C          pression et volume constants.
C
C     SORTIE :
C
C     FLUEUL : Tableau  de reels contenant  les  flux
C          explicites convectifs laminaires.
C
C     -----------------------------------------------
C     Debut des declarations
C     -----------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER NEQUATION, NBLOC
      PARAMETER( NEQUATION = 5 )
      PARAMETER( NBLOC = 2 )
C     
C     variables d'appel
      INTEGER NN, NNT, NNSG, ORDRE
      INTEGER NU(4,NNT), NUBO(2,NNSG)
      REAL*8 COOR(3,NN), GAMMA, VNOCL(3,NNSG)
      REAL*8 UA(NEQUATION,NN), FLUEUL(NEQUATION,NN)
      REAL*8 DX(NEQUATION,NN), DY(NEQUATION,NN), DZ(NEQUATION,NN)
C
C     Variables locales
C
C     Tableaux de travail
      REAL*8 VARLOC(NBLOC,NEQUATION), COORLOC(NBLOC,3)
      REAL*8 GRADLOC(NBLOC,NEQUATION,3)
      REAL*8 FLULOC(NEQUATION)
      INTEGER ILOCAL(NBLOC)
C
C     Indices de boucle
      INTEGER IS, NSG, K
C
C     Divers
      INTEGER NUBO1, NUBO2
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
C     Initialisation des flux convectifs
C
      DO 200 IS = 1, NN, 1
         DO 100 K = 1, NEQUATION, 1
            FLUEUL(K,IS) = 0.
 100     CONTINUE
 200  CONTINUE
C
C     Boucle d'integration sur les segments du maillage.
C
      DO 1000 NSG = 1, NNSG, 1
C
         NUBO1 = NUBO(1,NSG)
         NUBO2 = NUBO(2,NSG)
C
C        Determination des coordonnees et des variables locales.
C        On resout toutes les indirections avant l'appel a 
C        VANLEER.
C
         CALL GTOLOCAL( NEQUATION, NN, NNT, NBLOC, NU, COOR, UA,
     $                  DX,DY,DZ,NUBO1, NUBO2, ILOCAL, COORLOC, 
     $                  VARLOC,GRADLOC)
C
C        Appel du flux de Van Leer pour le segment courant
C
c  debug
c      print*,'nsg,vnloc',nsg,vnocl(1,nsg),VARLOC(1,1),
c     $ GRADLOC(1,1,1),GAMMA, NBLOC, NEQUATION,COORLOC(1,1)
c
         CALL VANLEER( NEQUATION, NBLOC, ORDRE, COORLOC, VNOCL(1,NSG),
     $                 VARLOC(1,1), VARLOC(1,2), VARLOC(1,3), 
     $                 VARLOC(1,4), VARLOC(1,5), GRADLOC, 
     $                 GAMMA, FLULOC )
c
c   debug
c      print*,'flux: fluloc(1)=',fluloc(1)
C
C        Affectation des flux 
C
         DO 800 K = 1, NEQUATION, 1
            FLUEUL(K,NUBO1) = FLUEUL(K,NUBO1) + FLULOC(K)
            FLUEUL(K,NUBO2) = FLUEUL(K,NUBO2) - FLULOC(K)
  800    CONTINUE
C
 1000 CONTINUE
C
      RETURN
      END
