C     PARAMETRES ET COMMON POUR LA PHASE D'OPTIMISATION :
C     ---------------------------------------------------
C
C     VARIABLES CONCERNANT LE MAILLAGE DE PEAU :
C
C     NSP : NOMBRE EXACT DE NOEUDS 
C     NTP : NOMBRE EXACT DE TRIANGLES
C     NSEGP : NOMBRE EXACT DE SEGMENTS
C     COORP : COORDONNEES DES NOEUDS
C     NUP : NUMERO DES SOMMETS DES TRIANGLES
C     NUBOP : NUMEROS DES EXTREMITES D'UN SEGMENT DONNE
C     VNO : NORMALE INITIALE A UNE CELLULE ISP DONNEE EN SON CENTRE
C     VNP : NORMALE A UN TRIANGLE TP DONNE
C     AIRTP : AIRE D'UN TRIANGLE
C     AIRESP : AIRE D'UNE CELLULE DE LA COQUE
C     AIRESP0 : AIRE D'UNE CELLULE DE LA COQUE INITIALE
C     NODE2D3D : TAB. DE CORRESPONDANCE ENTRE COQUE ET MAILLAGE SUR LES NOEUDS
C     FACE2D3D :   "            "         "   COQUE ET MAILLAGE SUR LES FACES
C     ITRANS   : CONDITIONS DE TRANSPIRATION
C
C     curva: (normal part of the) curvature tensor, defined on the coque
C     curno: curva's norm
C
      INTEGER NNSP, NNTP, NNSEGP
      PARAMETER (NNSP=3219, NNTP=6368, NNSEGP= 9586)
C                     ******    AILE M6 a 15460 noeuds
      INTEGER NSP, NTP, NSEGP
      REAL*8 COORP(3,NNSP), VNO(3,NNSP), VNP(3,NNTP), AIRTP(NNTP)
      REAL*8 AIRESP(NNSP),AIRESP0(NNSP) , VNON(3,NNSP)
      INTEGER NUBOP(2,NNSEGP), NUP(3,NNTP), LOGFRP(NNSP)
      INTEGER TMPT(10000), NODE2D3D(NNSP), FACE2D3D(NNTP)
      INTEGER ITRANS
C
      COMMON/DONMAIL1/COORP, VNO, VNP, AIRTP, AIRESP, AIRESP0, VNON
      COMMON/DONMAIL2/NSP, NTP, NSEGP, NUBOP, NUP, LOGFRP
      COMMON/CORRES2D3D/TMPT, NODE2D3D, FACE2D3D
      COMMON/TRANSPI/ITRANS

      real*8 curva(3,3,nnsp),curno(nnsp)
      common/zzsurcurv/curva,curno
C
C     VARIABLES CONCERNANT LA PARAMETRISATION HIERARCHIQUE :
Cg
C          nzma1 = somme des nombres de zones de tous les niveaux,
c          -----
c          nvma  = nb. maxi de niveaux de grilles,
c          ----
c          nsgma = somme de tous les segments de tous les niveaux
C          -----
C          nvoimax = nombre maximum de voisins pour une cellule
c          -------
c          nuvoi(iz,k) = numero de la kieme zone voisine iz   
c          -----------
c          nvoi(iz)    = nb. total de zones voisines a iz     
c          --------
c          nbz(nvl) = nombre de zones sur le niveau nvl
c          --------
c          nuzo(iz) = numero courant de la cellule ou se trouve iz
c          --------
c          gairesp(iz) = aire de la cellule iz, pour tous les niveaux
c          ----------
c          gvncoq(3,iz) = normale a la cellule iz, pour tous les niveaux
c          -----------
      INTEGER NZMA1, NSGMA, NVMA, NVOIMAX
      PARAMETER (NZMA1=2*NNSP, NSGMA=5*NNSP, NVMA=10, NVOIMAX=50)
C
      INTEGER NUVOI(NZMA1,NVOIMAX), NVOI(NZMA1), NUZO(NZMA1), NBZ(NVMA)
      REAL*8 GAIRESP(NZMA1), GVNCOQ(3,NZMA1)
C
      COMMON/ZONES1/NUVOI, NVOI, NUZO, NBZ
      COMMON/ZONES2/GAIRESP, GVNCOQ
C
      REAL*8 THETA
      COMMON/LIS/THETA

C
C     VARIABLES POUR L''ACCROISSEMENT DE LA FORME :
C
      REAL*8 DELTA(NZMA1), DELTANEW(NZMA1), DELTAZ(NNSP)
C
      COMMON/ACCROIS/DELTA, DELTANEW, DELTAZ
C 
      REAL*8 GRAD(NZMA1),COUT, COUTP, DJDGAM(NNSP), GRADZ(NNSP)
      REAL*8 DJDGAMZ(NNSP)
C
      COMMON/OPTIM/GRAD,COUT,COUTP,DJDGAM,GRADZ,DJDGAMZ

c...  variables zz

      real*8 reazz(30),xman,dixol,faczz(30)
      common/zzvari/reazz,xman,dixol,faczz       

      integer iliss,isuct,imvmsh,kmov,kmvcou,imvng
      common/ilissage/iliss,isuct,imvmsh,kmov,kmvcou,imvng





