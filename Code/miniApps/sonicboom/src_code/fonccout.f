
      SUBROUTINE FONCCOUT
C
C**** Calcul de la fonctionnelle cout J(W). Elle est definie uniquement
C     sur la coque et elle depend de la pression :
C            nsp
C            --  1                                                   2
C     J(W) = \   -  Aire(isp) * (Pression(is) - Pression_desiree(is)) 
C            /   3
C            --
C          isp=1
C     is=node2d3d(isp)
C
C     La fonctionnelle est stockee dans le scalaire Cout
C
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
C
      INTEGER ISP, IS
      REAL*8 PRESSION
C
      COUT = 0.
C
      DO 10 ISP = 1,NSP
         IS = NODE2D3D(ISP)
         PRESSION = GAM1 * (UA(5,IS) - 0.5*(UA(2,IS)**2 + UA(3,IS)**2
     $              + UA(4,IS)**2)/UA(1,IS))
      COUT = COUT + AIRESP0(ISP)/3. *(PRESSION - PDESP(IS))**2
 10   CONTINUE
C
      ncf = ncf + 1
c
      WRITE(6, *) ' '
      WRITE(6, *) '*** LA FONCTIONNELLE COUT VAUT : ',Cout
      WRITE(6, *) ' '
c
      OPEN(7,ACCESS='APPEND')
      WRITE(7,*)ITOPT,COUT
      CLOSE(7)
C
      RETURN
      END
