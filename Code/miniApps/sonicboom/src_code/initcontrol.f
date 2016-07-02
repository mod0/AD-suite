      SUBROUTINE INITCONTROL(CTRL)
C
c*** Initialisation du controle
c
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER IS, ISP, KVAR, I, nunit
      REAL*8 PI, CTRL(NNSP),coorpbar(3,nnsp),x,ddx,y,ray
      REAL*8 NORMALIS
C
      IF (ITRANS.EQ.1) THEN
c
         IF (COEFM1.EQ.10) THEN ! Cas de la tuyere
c
            IF (contr.eq.1) THEN
c
C     Definition du controle CTRL pour transpirer la tuyere a
c     partir de la conduite a section constante.
c
               PI = acos(-1.0)
               ray=2.99719
c
               DO ISP = 1,NSP
c
                  IS = NODE2D3D(ISP) 
                  IF (((COOR(1,IS).GE.0.).AND.(COOR(1,IS).LE.2.)).AND.
     $                 (COOR(2,IS).EQ.5.)) THEN
                     x = 0.5*coor(1,is) + 0.5
                     if (abs(x-1.0).le.0.5) then
                     ddx=x-1.0
                     y=-ray+0.042+sqrt(ray*ray-ddx*ddx)
                     CTRL(ISP) = - 2.*y
                     endif
                     
c                     CTRL(ISP) = 4.75 + 0.25 *SIN(PI*(COOR(1,IS)
c     $                                    +0.5)) - 5.
                  ELSE
                     PRINT*,'ATTENTION DANS INITCONTROL'
                     CTRL(ISP) = 0.
                  ENDIF
               END DO
c
            ELSE
c
               IF (contr.eq.2) THEN
C
c Initialisation du controle de forme CTRL a 0.
c
                  DO ISP = 1,NSP
                     CTRL(ISP) = 0.
                  END DO
C      
               ENDIF
c
            ENDIF
c
         ELSE !! On a toujours ITRANS = 1 mais COEFM1<>10
c
c Cas de l'aile d'avion ou du falcon (contr n'est egal ni a 1 ni a 2)
c
            DO ISP = 1,NSP
               CTRL(ISP) = 0.
            END DO
c
         ENDIF !! Fin du cas ou ITRANS=1
c
      ELSE  !! Cas de non transpiration ITRANS=0 quel que soit COEFM1
C
         DO ISP = 1,NSP
            CTRL(ISP) = 0.
         END DO
C
      ENDIF
c
      nunit=99
C [llh]      print*,'InitControl.f: nsp=',nsp,'nunit=',nunit
      CALL ECRITURE(CTRL,nunit)
c
      RETURN
      END
