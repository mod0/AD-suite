
      SUBROUTINE RESPRES(nunit)
c---------------------------------------------------------------------
c   Ecriture de la pression lissee desiree au format Vigie
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------   
      INTEGER nps   , is, nunit
      REAL*8    pmax  , pmin
      REAL*8    pm
c
      pmax                         =-1.0e+30
      pmin                         = 1.0e+30
c
      nps                          = 0
c
      DO 10 is=1,ns
c
         pm                        = Pdesp(is)
         pmax                      = MAX(pmax, pm)
         pmin                      = MIN(pmin, pm)
c
10    CONTINUE
c
c      Ecriture dans nunit des solutions au format hdb
c      pour la visualisation avec VIGIE.
c
         Do is = 1,ns
           WRITE(nunit,113) Pdesp(is)
         End Do
c
113   FORMAT(e15.8)
c
      RETURN
      END
