

      SUBROUTINE MAIL3D2D(numfic)
c
      include 'Param3D.h'
      include 'Paramopt3D.h'
c
C     Creation des tableaux de correspondance entre les maillages 3D
C     et 2D/3D (coque) :
C     NODE2D3D(1:NSP), NODE3D2D(1:NS), FACE2D3D(1:NTP), FACE3D2D(1:NFAC)
C
C     CREATION DU NOMBRE DE   NOEUDS  SUR LA COQUE : NSP
C        "          "       TRIANGLES     "        : NTP
C     CREATION D'UN NOUVEAU TABLEAU DE COORDONNEES SUR LA COQUE : COORP(3,NSP)
C        "                             SOMMETS            "     : NUP(3,NTP)
C        "                             LOGFR              "     : LOGFRP(NSP)
C     ECRITURE DANS FORT.88 DES DONNEES DU MAILLAGE
C
      integer i, is, jt, ifac, k, idrap(nsmax),numfic,nunit
      real*8 XMAXI, XMINI, YMAXI, YMINI, ZMAXI, ZMINI
      integer I1, I2, I3
c
      DO 10 IS=1,NS
         NODE3D2D(IS)=0
         IDRAP(IS) = 0
10    CONTINUE
C
      DO 15 IFAC = 1,NFAC
         FACE3D2D(IFAC) = 0
15    CONTINUE
C
      NTP = 0
      NSP = 0
c
      DO 20 IFAC=1,NFAC
C
           IF ((LOGFAC(IFAC).EQ.200).OR.(LOGFAC(IFAC).EQ.-2)
     $        .OR.(LOGFAC(IFAC).EQ.7)) THEN
             NTP = NTP + 1
             FACE2D3D(NTP) = IFAC
             FACE3D2D(IFAC) = NTP
             IDRAP(NSFAC(1,IFAC) ) = 1
             IDRAP(NSFAC(2,IFAC) ) = 1
             IDRAP(NSFAC(3,IFAC) ) = 1
             DO K = 1,3
                IF (NODE3D2D(NSFAC(K,IFAC)).EQ.0) THEN
                   NSP = NSP+1
                   LOGFRP(NSP)  = 2
                   NODE2D3D(NSP) = NSFAC(K,IFAC)
                   NODE3D2D(NSFAC(K,IFAC)) = NSP
                   
                ENDIF
             END DO
         ENDIF
C
20    CONTINUE

      XMAXI = - 1000000
      XMINI =   1000000
      YMAXI = - 1000000
      YMINI =   1000000
      ZMAXI = - 1000000
      ZMINI =   1000000
c
      DO 30 IS= 1,NS
         IF(IDRAP(IS) .EQ. 1) THEN
             COORP(1,NODE3D2D(IS))=COOR(1,IS)
             COORP(2,NODE3D2D(IS))=COOR(2,IS)
             COORP(3,NODE3D2D(IS))=COOR(3,IS)
             XMINI = MIN( XMINI,COORP(1,NODE3D2D(IS)))
             XMAXI = MAX( XMAXI,COORP(1,NODE3D2D(IS)))
             YMINI = MIN( YMINI,COORP(2,NODE3D2D(IS)))
             YMAXI = MAX( YMAXI,COORP(2,NODE3D2D(IS)))
             ZMINI = MIN( ZMINI,COORP(3,NODE3D2D(IS)))
             ZMAXI = MAX( ZMAXI,COORP(3,NODE3D2D(IS)))
         ENDIF
30    CONTINUE
C

C [llh]      WRITE(6, *) '              4- Creation du maillage de peau : '
C [llh]      WRITE(6, *) '              ********************************* '
C [llh]      WRITE(6, *) ' '
c
C [llh]      PRINT*,'LE NOMBRE DE NOEUDS SUR LA COQUE EST :',NSP
C [llh]      PRINT*,'LE NOMBRE DE TRIANGLES SUR LA COQUE EST :',NTP
C
C [llh]      PRINT*,'CALCUL DES MIN ET DES MAX DES COORDONNEES DE LA COQUE :'
C [llh]      WRITE(6, *) ' '
C [llh]      PRINT*, ' MAXI X = ',XMAXI, ' MINI X = ',XMINI
C [llh]      PRINT*, ' MAXI Y = ',YMAXI, ' MINI Y = ',YMINI
C [llh]      PRINT*, ' MAXI Z = ',ZMAXI, ' MINI Z = ',ZMINI
C [llh]      WRITE(6, *) ' '
C [llh]      WRITE(6, *) ' '
C
      DO 40 IFAC=1,NFAC
         IF (FACE3D2D(IFAC).NE.0) THEN
             NUP(1,FACE3D2D(IFAC))= NODE3D2D(NSFAC(1,IFAC))
             NUP(2,FACE3D2D(IFAC))= NODE3D2D(NSFAC(2,IFAC))
             NUP(3,FACE3D2D(IFAC))= NODE3D2D(NSFAC(3,IFAC))
         ENDIF
40    CONTINUE
c
            REWIND(numfic)
            WRITE(numfic,*) ' couleurs '
            WRITE(numfic,*) ' 4 '
            WRITE(numfic,*) '  0   .7   .8   1.  '
            WRITE(numfic,*) '  1   .8   .9   1.  '
            WRITE(numfic,*) '  2   0.   1.   1.  '
            WRITE(numfic,*) '  3   0.   0.   1.  '
            WRITE(numfic,*) ' 2 '
            WRITE(numfic,*) ' 0 '
            WRITE(numfic,*) '   4  0. 0. 0. '
            WRITE(numfic,*) ' 128  1. 1. 1. '
            WRITE(numfic,*) ' 1 '
            WRITE(numfic,*) ' 129 0. 0. 0. '
            WRITE(numfic,*) ' 255 1. 0. 0. '
            WRITE(numfic,*) ' domaine '
            WRITE(numfic,*)  XMINI,XMAXI,YMINI,YMAXI,ZMINI,ZMAXI
c
c 111        format(6(E12.6,1x))
c
            WRITE(numfic,*) ' triangles '
            WRITE(numfic,*) NSP
c
            DO 12 IS=1, NSP
                WRITE(numfic,78) COORP(1,IS),COORP(2,IS),COORP(3,IS)
12          CONTINUE

            WRITE(numfic,*) NTP

            DO 13 JT=1, NTP
              I1 =  NUP(1,JT)-1
              I2 =  NUP(2,JT)-1
              I3 =  NUP(3,JT)-1
               WRITE(numfic,79) I1,I2,I3
13          CONTINUE
c
C=== ECRITURE DU MAILLAGE DANS FORT.88
c
            NUNIT=88
            OPEN(NUNIT,FORM='FORMATTED')
            write(NUNIT,80) NSP,NTP
            write(NUNIT,81) ((COORP(I,IS),I=1,3),IS=1,NSP)
            write(NUNIT,82) ((NUP(K,JT),K=1,3), JT=1,NTP)
            write(NUNIT,83) (LOGFRP(is),is=1,NSP)


78       FORMAT (3x,E12.6,2x,E12.6,2x,E12.6)
79       FORMAT (I6,2x,I6,2x,I6)
80       FORMAT (2i5)      
81       FORMAT (8E15.8)
82       FORMAT (8i5)
83       FORMAT (i5)
c
      return
      end
