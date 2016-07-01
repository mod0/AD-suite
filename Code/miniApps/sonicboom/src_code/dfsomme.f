
C     -----------------------------------------------
      SUBROUTINE DFSOMME( NEQUATION, ALPHA, BETA, ROT, DFD, DFX,
     $                    NORME, DMAT )
C     -----------------------------------------------
C
C     CETTE PROCEDURE CALCULE LA SOMME DES DERIVEES DU FLUX CENTRE 
C                            _
C     ET DU FLUX DECENTRE : DF
C
C     ON EN DEDUIT ALORS LA MATRICE IMPLICITE DES FLUX DE VAN LEER :
C                                 -1   _ _
C               DF ( W ) = {N} * R  * DF(W) * R
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
      INTEGER NEQUATION
      REAL*8 ALPHA, BETA, ROT(NEQUATION,NEQUATION)
      REAL*8 DFD(NEQUATION,NEQUATION), DFX(NEQUATION,NEQUATION)
      REAL*8 DMAT(NEQUATION,NEQUATION),NORME
C
C     Variables locales
C
C     Tableaux de travail
      REAL*8 DF(5,5), MTMP(5,5)
C
C     Indices de boucle
      INTEGER IA, IB, IC
C
C     Divers
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
      DO 200 IA = 1, NEQUATION, 1
         DO 100 IB = 1, NEQUATION, 1
C
C           Decentrage
C
           DF(IA,IB) = ALPHA*DFD(IA,IB) + BETA*DFX(IA,IB)
C
C           mise a zero de DMAT et de MTMP
C
            DMAT(IA,IB) = 0.
            MTMP(IA,IB) = 0.
C
  100    CONTINUE
c$$$         MTMP(IA,IA) = 1.
  200 CONTINUE
C
C     Retour aux coordonnees (X,Y,Z) par rotation inverse
C
C                   -1  _ _
C     DF(W) = |N | R   DF(W) R
C
C     ou N est le vecteur normal a la facette de la cellule
C     contenant le point courant, R la rotation de l'espace,
C     _           _
C     W = R(W) et F est (dans ce sous programme) represente 
C     par la variable DF.
C
C                      _ _
C     On commence par DF(W) R
C
      DO 500 IA = 1, NEQUATION, 1
         DO 400 IB = 1, NEQUATION, 1
            DO 300 IC = 1, NEQUATION, 1
               MTMP(IA,IB) = MTMP(IA,IB) + DF(IA,IC)* ROT(IC,IB)
  300       CONTINUE
  400    CONTINUE
  500 CONTINUE
C
C                              -1  _ _
C     Maintenant, on effectue R  [DF(W) R], en remarquant que
C      -1    t
C     R   = R       (poil au nez)
C
      DO 800 IA = 1, NEQUATION, 1
         DO 700 IB = 1, NEQUATION, 1
            DO 600 IC = 1, NEQUATION, 1
               DMAT(IA,IB) = DMAT(IA,IB) + NORME*ROT(IC,IA)*MTMP(IC,IB)
c$$$               DMAT(IA,IB) = DMAT(IA,IB) + ROT(IC,IA)*ROT(IC,IB)
  600       CONTINUE
  700    CONTINUE
  800 CONTINUE
C
c$$$      DO 1000 IA = 1, NEQUATION, 1
c$$$         DO 900 IB = 1, NEQUATION, 1
c$$$            MTMP(IA,IB) = MTMP(IA,IB) - DMAT(IA,IB)
c$$$            IF (ABS(MTMP(IA,IB)) .GT. 1.E-6) THEN
c$$$               PRINT *, MTMP(IA,IB)
c$$$            ENDIF
c$$$ 900     CONTINUE
c$$$ 1000 CONTINUE
C
      RETURN
      END
