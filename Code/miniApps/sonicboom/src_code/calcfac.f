      SUBROUTINE CALCFAC
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  VERIFICATION SUR LES VOLUMES ET LES NORMALES AUX FACETTES
C
C  PARAMETRES D'ENTREE :
C  ---------------------
C  NS      : NOMBRE DE SOMMETS
C  NT      : NOMBRE DE TETRAEDRES
C  NFAC    : NOMBRE DE FACES FRONTIERES
C  
      INCLUDE 'Param3D.h'
C      
      REAL*8 VOL,SUMX,SUMY,SUMZ
      INTEGER IS,JT,IFAC
C
C           ........... SOMME DES VOLUMES ..........
C
      VOL= 0.
      DO 210 IS = 1 , NS , 1
         VOL= VOL + VOLS(IS)
210   CONTINUE
c
C [llh]      WRITE(6,*) 'Verification sur la somme des normales exterieurs :'
C [llh]      WRITE(6,*) '*************************************************** '
C [llh]      WRITE(6, *) ' '
C [llh]      WRITE(6,*) 'SOMME DES VOLUMES EN LES CELLULES: ',VOL
C
      VOL= 0.
      DO 220 JT = 1 , NT , 1
         VOL= VOL + VOLT(JT)
220   CONTINUE
C [llh]      WRITE(6,*) 'SOMME DES VOLUMES EN LES TETRAEDRES: ',VOL
c
C [llh]      WRITE(6, *) ' '
C
C  INTEGRALE DES NORMALES EXTERIEURES
C  ==================================
C
      SUMX= 0.
      SUMY= 0.
      SUMZ= 0.
      DO 300 IFAC = 1 , NFAC , 1
         SUMX= SUMX + VNFAC(1,IFAC)
         SUMY= SUMY + VNFAC(2,IFAC)
         SUMZ= SUMZ + VNFAC(3,IFAC)
300   CONTINUE
C [llh]      WRITE(6,*) 'Integrale des normales exterieures :'
C [llh]      WRITE(6, *) 'Suivant X :',SUMX
C [llh]      WRITE(6, *) 'Suivant Y :',SUMY
C [llh]      WRITE(6, *) 'Suivant Z :',SUMZ
C [llh]      WRITE(6, *) ' '
C [llh]      WRITE(6, *) ' '
C
      RETURN
      END
