      SUBROUTINE CONDDIRFLUX
C
c     Introduction des conditions de Dirichlet:
c     si logfac(ifac)=5 (et logfr(is)=5)   ou   si logfr(is)=1
c
c     On met les flux explicites a 0 en ces noeuds.
c
      include 'Param3D.h'
C
      INTEGER IS, I
C
C$AD II-LOOP
      DO IS = 1,NS
        IF (LOGFR(IS).EQ.1 .or. LOGFR(IS).EQ.5) THEN
          DO I = 1,5
            CE(I,IS) = 0.d0
          END DO
        ENDIF
      END DO

      END
