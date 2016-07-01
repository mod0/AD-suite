C     -----------------------------------------------
      SUBROUTINE MATVL( NEQUATION, NNdim, NSGMAXdim,NUBO,
     &          VNOCL, GAMMA, UA, DIAG,ZM,nn,nseg)
C     -----------------------------------------------
C
C     Construction de la matrice implicite : Flux Van-Leer
C
C     -----------------------------------------------
C     Debut des declarations
C     -----------------------------------------------
C
      IMPLICIT NONE
C
C     variables d'appel
C
      INTEGER NEQUATION, NN, NSGMAX, NSEG
      integer nndim,nsgmaxdim
      REAL*8 VNOCL(3,NSGMAXdim), UA(NEQUATION,NNdim)
      REAL*8 ZM(NEQUATION,NEQUATION,2,NSGMAXdim),
     $DIAG(NNdim,NEQUATION,NEQUATION)
      INTEGER NUBO(2,NSGMAXdim)
      REAL*8 GAMMA
C
C     Variables locales
C
C     Tableaux de travail
      REAL*8 DFD(5,5), DFX(5,5), DMAT(5,5,2), ROT(5,5)
C
C     Indices de boucle
      INTEGER NSG, IA, IB, is
C
C     Divers
      REAL*8 GAM1, GAM3, GAM4, GAM5, GAM6
      REAL*8 C, UN, VN, WN, NORME
      REAL*8 AA, BB, CC, XMM, XMP,AA1,BB1
      INTEGER NUBO1, NUBO2
c
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
C     CONSTANTES DU CALCUL
C
C     GAM = 1.4 , GAM1 = GAM - 1 , GAM3 = GAM - 3
C     GAM4 = GAM*GAM1, GAM5 = 1./GAM , GAM6 = 1. / (GAM*GAM-1)
C
C     -----------------------------------------------
C
      GAM1= GAMMA - 1. 
      GAM3= GAMMA - 3.
      GAM4= GAMMA*GAM1
      GAM5= 1./GAMMA
      GAM6= 1./( GAMMA**2 - 1. )
c
C     Boucle sur les segments
C     =======================
C
      DO 30 NSG = 1 , NSEG , 1
C        
C        1/ Quantites en NUBO1
C        ========================
C
         NUBO1= NUBO(1,NSG)
         NUBO2= NUBO(2,NSG)
C
c*** Debut du test de validation sur le calcul des derivees.
c    On enleve la rotation
c
         CALL ROTATION( NEQUATION,UA(1,NUBO1), UA(2,NUBO1), UA(3,NUBO1),
     $                  UA(4,NUBO1), UA(5,NUBO1), GAMMA, VNOCL(1,NSG),
     $                  UN, VN, WN, C, NORME, ROT  )

C
C        Calcul des valeurs servant a determiner le decentrage.
C                                         
C        Calcul de DFX jacobien exact de F en NUBO1
C
C        DF(I,J) =  DFJ / DWI
C
         XMM = 0.5*( UN/C - 1. )
         XMP = 0.5*( UN/C + 1. )
c
         CALL DFEXACT( NEQUATION, UA(1,NUBO1), UN, VN, WN,
     $                 UA(NEQUATION,NUBO1), GAMMA, GAM1, GAM3, DFX )
C
C        Calcul des derivees de F-  (DECENTRAGE)
C
C        DF(I,J) = DFI/DWJ
C
         CALL DFDECENTRE( NEQUATION, UA(1,NUBO1), UN, VN, WN,
     $                    UA(NEQUATION,NUBO1), GAMMA, GAM1,
     $                    GAM4, GAM5, GAM6, XMM, XMP, C, DFD )
c     
c**** Fin du test en nubo1
C
C        Calcul de DF(NUBO1) :      +
C        DF = (AA + BB)*DFX - AA*DFD
c
C        Avec DFX jacobien exact de F et DFD jacobien decentre de F+
C        et :
C
C     AA = 0 et BB = 0 si UN/C >= 1
C     AA = 1 et BB = 0 si |UN/C| < 1
C     AA = 0 et BB = 1 si UN/C <= -1
C
C        DF(I,J) =  DFJ / DWI
C
         IF (UN/C.ge.1.) THEN
            AA = 0.
            BB = 0.
         ELSE
            IF (ABS(UN/C).lt.1.) THEN
               AA = 1.
               BB = 0.
            ELSE
               IF (UN/C.le.-1.) THEN
                  AA = 0.
                  BB = 1.
               ENDIF
            ENDIF
         ENDIF
c

         AA1= -AA
         BB1= BB+AA
         CALL DFSOMME(NEQUATION,AA1,BB1,ROT,DFD,DFX,NORME,
     $                DMAT(1,1,1))
c
C        2/ Quantites en NUBO2
C        ========================
C
C        Rotation (U,V,W) --> (UN,VN,WN) et variables auxilliaires
C
         CALL ROTATION(NEQUATION, UA(1,NUBO2), UA(2,NUBO2), UA(3,NUBO2),
     $                  UA(4,NUBO2), UA(5,NUBO2), GAMMA, VNOCL(1,NSG),
     $                  UN, VN, WN, C, NORME, ROT)
C                                        
C        CALCUL DE DFX JACOBIEN EXACT DE F EN NUBO2
C
C        DF(I,J) =  DFJ / DWI
C
         XMM = 0.5*( UN/C - 1. )
         XMP = 0.5*( UN/C + 1. )
c
         CALL DFEXACT( NEQUATION, UA(1,NUBO2), UN, VN, WN,
     $                 UA(NEQUATION,NUBO2), GAMMA, GAM1, GAM3, DFX )
C
C        Calcul des derivees de F+  (decentrage)
C
C        DFD(I,J) = DFDI/DWJ
C
         CALL DFDECENTRE( NEQUATION, UA(1,NUBO2), UN, VN, WN, 
     $                    UA(NEQUATION,NUBO2), GAMMA, GAM1,
     $                    GAM4, GAM5, GAM6, XMM, XMP, C, DFD )
c
C        Calcul de DF(NUBO2) :          
C                   + 
C        DF = AA*DFD + BB*DFX
C
C        Avec DFX Jacobien exact de F et DFD Jacobien decentre de F+
C
C     AA = 0 et BB = 1 si UN/C >= 1
C     AA = 1 et BB = 0 si |UN/C| < 1
C     AA = 0 et BB = 0 si UN/C <= -1
C
C        DF(I,J) =  DFJ / DWI
C
         IF (UN/C.ge.1.) THEN
            AA = 0.
            BB = 1.
         ELSE
            IF (ABS(UN/C).lt.1.) THEN
               AA = 1.
               BB = 0.
            ELSE
               IF (UN/C.le.-1.) THEN
                  AA = 0.
                  BB = 0.
               ENDIF
            ENDIF
         ENDIF
c
         CALL DFSOMME(NEQUATION,AA,BB,ROT,DFD,DFX,NORME,DMAT(1,1,2))
C
C        On remplit la matrice:
C        DIAG = termes diagonaux
C        ZM   = termes extra-diagonaux.
C
         DO 20 IA = 1, NEQUATION, 1
            DO 10 IB = 1, NEQUATION, 1
c
               DIAG(NUBO1,IA,IB) = DIAG(NUBO1,IA,IB) - DMAT(IA,IB,1)
               DIAG(NUBO2,IA,IB) = DIAG(NUBO2,IA,IB) + DMAT(IA,IB,2)
c
               ZM(IA,IB,1,NSG)   = ZM(IA,IB,1,NSG) + DMAT(IA,IB,1)
               ZM(IA,IB,2,NSG)   = ZM(IA,IB,2,NSG) - DMAT(IA,IB,2)
c
10         CONTINUE
20      CONTINUE
C
30    CONTINUE
c
c   debug
c        is=2203
c        print*,'fin matvl: is= 2203 diag=',
c     $DIAG(IS,1,1),DIAG(IS,2,2),DIAG(IS,3,3),DIAG(IS,4,4)
c
      RETURN
      END
