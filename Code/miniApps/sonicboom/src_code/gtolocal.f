
C     -----------------------------------------------
      SUBROUTINE GTOLOCAL( NEQUATION, NN, NNT, NBLOC, NU, COOR, UA,
     $                     DX,DY,DZ,NUBO1, NUBO2, ILOCAL, COORLOC, 
     $                     VARLOC,GRADLOC )
C     -----------------------------------------------
C
C     La procedure de localisation a pour but, etant
C     donnes un segment, ses extremites, ses triangles
C     amont-aval et ses triangles adjacents de fournir
C     les coordonnees des points de la molecule 
C     "papillon" et les variables conservatives
C     correspondantes dans un tableau local.
C
C     -----------------------------------------------
C     Debut des declarations
C     -----------------------------------------------
C
      IMPLICIT NONE
C
C     inclusion du header
C
C     variables d'appel
      INTEGER NEQUATION, NN, NNT, NBLOC
      INTEGER NU(4,NNT)
      REAL*8 COOR(3,NN), UA(NEQUATION,NN)
      REAL*8 DX(NEQUATION,NN), DY(NEQUATION,NN), DZ(NEQUATION,NN)
      REAL*8 GRADLOC(NBLOC,NEQUATION,3)
      INTEGER NUBO1, NUBO2
      INTEGER ILOCAL(NBLOC)
      REAL*8 COORLOC(NBLOC,3), VARLOC(NBLOC,NEQUATION)
C
C     Variables locales
C
C     Tableaux de travail
CJMM$$      INTEGER NEXT(3), OTHER(3)
C
C     Indices de boucle
      INTEGER K, L
C
C     Divers
CJMM$$      INTEGER NTVA1, NTVA2
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
CJMM$$      NEXT(1) = 2
CJMM$$      NEXT(2) = 3
CJMM$$      NEXT(3) = 1
CJMM$$C
CJMM$$      OTHER(1) = 3
CJMM$$      OTHER(2) = 1
CJMM$$      OTHER(3) = 2
C
C     D'abord les deux points du segment considere.
C
C     Geometrie :
C
      ILOCAL(1) = NUBO1
      ILOCAL(2) = NUBO2
C
      COORLOC(1,1) = COOR(1,NUBO1)
      COORLOC(1,2) = COOR(2,NUBO1)
      COORLOC(1,3) = COOR(3,NUBO1)
      COORLOC(2,1) = COOR(1,NUBO2)
      COORLOC(2,2) = COOR(2,NUBO2)
      COORLOC(2,3) = COOR(3,NUBO2)
c
      GRADLOC(1,1,1) = DX(1,NUBO1)
      GRADLOC(1,1,2) = DY(1,NUBO1)
      GRADLOC(1,1,3) = DZ(1,NUBO1)
      GRADLOC(1,2,1) = DX(2,NUBO1)
      GRADLOC(1,2,2) = DY(2,NUBO1)
      GRADLOC(1,2,3) = DZ(2,NUBO1)
      GRADLOC(1,3,1) = DX(3,NUBO1)
      GRADLOC(1,3,2) = DY(3,NUBO1)
      GRADLOC(1,3,3) = DZ(3,NUBO1)
      GRADLOC(1,4,1) = DX(4,NUBO1)
      GRADLOC(1,4,2) = DY(4,NUBO1)
      GRADLOC(1,4,3) = DZ(4,NUBO1)
      GRADLOC(1,5,1) = DX(5,NUBO1)
      GRADLOC(1,5,2) = DY(5,NUBO1)
      GRADLOC(1,5,3) = DZ(5,NUBO1)
c
      GRADLOC(2,1,1) = DX(1,NUBO2)
      GRADLOC(2,1,2) = DY(1,NUBO2)
      GRADLOC(2,1,3) = DZ(1,NUBO2)
      GRADLOC(2,2,1) = DX(2,NUBO2)
      GRADLOC(2,2,2) = DY(2,NUBO2)
      GRADLOC(2,2,3) = DZ(2,NUBO2)
      GRADLOC(2,3,1) = DX(3,NUBO2)
      GRADLOC(2,3,2) = DY(3,NUBO2)
      GRADLOC(2,3,3) = DZ(3,NUBO2)
      GRADLOC(2,4,1) = DX(4,NUBO2)
      GRADLOC(2,4,2) = DY(4,NUBO2)
      GRADLOC(2,4,3) = DZ(4,NUBO2)
      GRADLOC(2,5,1) = DX(5,NUBO2)
      GRADLOC(2,5,2) = DY(5,NUBO2)
      GRADLOC(2,5,3) = DZ(5,NUBO2)
C     
C     Variables hydrodynamiques
C     
      DO 400 K = 1, NEQUATION, 1
         VARLOC(1,K) = UA(K,NUBO1)
         VARLOC(2,K) = UA(K,NUBO2)
  400 CONTINUE
C     
CJMM$$C     Maintenant les extensions pour l'ordre 2, c'est a dire
CJMM$$C     les valeurs correspondant aux points 3,4,5 et 6 de la
CJMM$$C     molecule "papillon".
CJMM$$C     
CJMM$$C     Il n'y a pas d'exceptions, en effet il est de la 
CJMM$$C     responsabilite de TOTU de fournir dans TOUS les cas
CJMM$$C     un triangle amont ou un triangle aval.
CJMM$$C     
CJMM$$C     Classement des points dans les variables locales
CJMM$$C     
CJMM$$C     On remplit completement ces variables locales, c'est a
CJMM$$C     dire que l'on delivre les 6 variables et les 8 points
CJMM$$C     meme lorsque l'on fait de l'ordre 1; dans ce cas, c'est 
CJMM$$C     la procedure d'interpolation qui fait la difference.
CJMM$$C     
CJMM$$      DO 700 L = 1, 3, 1
CJMM$$         IF ( NU(L,TAMONT) .EQ. NUBO1 ) THEN
CJMM$$C           
CJMM$$            NTVA1 = NU(NEXT(L),TAMONT)
CJMM$$            NTVA2 = NU(OTHER(L),TAMONT)
CJMM$$C           
CJMM$$            DO 500 K = 1, NEQUATION, 1
CJMM$$               VARLOC(5,K) = UA(NTVA1,K)
CJMM$$               VARLOC(6,K) = UA(NTVA2,K)
CJMM$$  500       CONTINUE
CJMM$$C           
CJMM$$            COORLOC(5,1) = COOR(1,NTVA1)
CJMM$$            COORLOC(5,2) = COOR(2,NTVA1)
CJMM$$            COORLOC(6,1) = COOR(1,NTVA2)
CJMM$$            COORLOC(6,2) = COOR(2,NTVA2)
CJMM$$C
CJMM$$            ILOCAL(5) = NTVA1
CJMM$$            ILOCAL(6) = NTVA2
CJMM$$C           
CJMM$$         ENDIF
CJMM$$C        
CJMM$$         IF ( NU(L,TAVAL) .EQ. NUBO2 ) THEN
CJMM$$C           
CJMM$$            NTVA1 = NU(NEXT(L),TAVAL)
CJMM$$            NTVA2 = NU(OTHER(L),TAVAL)
CJMM$$C           
CJMM$$            DO 600 K = 1, NEQUATION, 1
CJMM$$               VARLOC(3,K) = UA(NTVA1,K)
CJMM$$               VARLOC(4,K) = UA(NTVA2,K)
CJMM$$  600       CONTINUE
CJMM$$C           
CJMM$$            COORLOC(3,1) = COOR(1,NTVA1)
CJMM$$            COORLOC(3,2) = COOR(2,NTVA1)
CJMM$$            COORLOC(4,1) = COOR(1,NTVA2)
CJMM$$            COORLOC(4,2) = COOR(2,NTVA2)
CJMM$$C           
CJMM$$            ILOCAL(3) = NTVA1
CJMM$$            ILOCAL(4) = NTVA2
CJMM$$C           
CJMM$$         ENDIF
CJMM$$  700 CONTINUE
CJMM$$C     
CJMM$$C     On rajoute les points 7 et 8 necessaires au calcul
CJMM$$C     des normales. Le point 7 est toujours present. Si le
CJMM$$C     point 8 n'existe pas, on donne des coordonnees fictives,
CJMM$$C     qui sont celles du milieu du segment [NUBO1,NUBO2].
CJMM$$C     
CJMM$$      COORLOC(7,1) = COOR(1,NUSO1)
CJMM$$      COORLOC(7,2) = COOR(2,NUSO1)
CJMM$$C     
CJMM$$      IF ( NUSO2 .EQ. 0 ) THEN
CJMM$$         COORLOC(8,1) = ( COOR(1,NUBO1) + COOR(1,NUBO2) )/ 2.
CJMM$$         COORLOC(8,2) = ( COOR(2,NUBO1) + COOR(2,NUBO2) )/ 2.
CJMM$$      ELSE
CJMM$$         COORLOC(8,1) = COOR(1,NUSO2)
CJMM$$         COORLOC(8,2) = COOR(2,NUSO2)
CJMM$$      ENDIF
C     
      RETURN
      END
