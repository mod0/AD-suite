
      SUBROUTINE FULLVCYCLE(KNIV)
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER KNIV
C
C*********************************************************************
C                  METHODE MULTINIVEAU FULL-V-CYCLES                 C
C*********************************************************************
C
C   Exemple ou nbvc = 2
C   (nbvc est le nombre de fois ou un meme V-Cycle est repete)
C
C                                      O           O
C                                     / \         / \
C                    O       O       /   O       /   O
C                   / \     / \     /     \     /     \
C          O   O   /   O   /   O   /       O   /       O
C         / \ / \ /     \ /     \ /         \ /         \
C      O-O   O   O       O       O           O           O ...........
C
C      |-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|->  ITERATIONS
C      1 2 3 4 5 6 7 8 9 ...........
C
C      |-|-------|---------------|-----------------------...........
C       2    2        2 cycles
C
C
C Lors des nbvc premieres iterations, on travaille sur le niveau grossier NIVG
C
C De l'iteration nbvc+1 a l'iteration 3*nbvc, on alterne sur les niveaux
c NIVG et NIVG-1
C
      IF ((ITOPT.GT.nbvc).AND.(ITOPT.LE.3*nbvc)) THEN
         KNIV = KNIV - 1
         IF (MOD(NBVC,2).EQ.0) THEN
            IF (MOD(ITOPT,2).EQ.0) KNIV = NIVG
         ENDIF
         IF (MOD(NBVC,2).EQ.1) THEN
            IF (MOD(ITOPT,2).EQ.1) KNIV = NIVG       
         ENDIF     
      ENDIF
C
      IF (NIVG-2.NE.NIVEAU) THEN
C
C****                 CAS DE L'AILE M6 ET DU FALCON               *****C
C                     =============================
C
C De l'iteration 3*NBVC+1 a l'iteration 6*NBVC, on alterne sur les niveaux 
C NIVG-2, NIVG-1 et NIVG (Dents de Scie).
C A partir de l'iteration 6*NBVC+1, on alterne du niveau fin NIVEAU=1 au niveau
C grossier NIVG (toujours en dents de scie).
C
         IF ((ITOPT.GT.3*NBVC).AND.(ITOPT.LE.6*NBVC)) THEN
            IF (MOD(ITOPT,3).EQ.1) KNIV = NIVG - 2
            IF (MOD(ITOPT,3).EQ.2) KNIV = NIVG - 1
            IF (MOD(ITOPT,3).EQ.0) KNIV = NIVG 
         ENDIF
C
         IF (ITOPT.GT.6*NBVC) THEN
            KNIV = KNIV + 1
            IF (KNIV.EQ.NIVG+1) KNIV = NIVEAU
         ENDIF
C
      ENDIF
C
      IF (NIVG-2.EQ.NIVEAU) THEN
C
C****                       CAS DE LA TUYERE    
C                           ================           
C
C On ne travaille qu'avec 3 niveaux :
C
C                    O       O   
C                   / \     / \  
C          O   O   /   O   /   O 
C         / \ / \ /     \ /     \ 
C      O-O   O   O       O       O ...............
C
         IF (ITOPT.GT.3*NBVC) THEN
            KNIV = KNIV + 1
            IF (KNIV.EQ.NIVG+1) KNIV = NIVEAU
         ENDIF
      ENDIF
C
      RETURN
      END
