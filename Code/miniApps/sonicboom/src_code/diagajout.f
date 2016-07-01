      SUBROUTINE DIAGAJOUT
C
      INCLUDE 'Param3D.h'
C
      INTEGER IA,IS
C
C     Ajout a la diagonale de la matrice implicite du terme 
c     (Id/ Delta t* Volume_Seg).
C
      DO IA = 1, 5
czz         DO IS=1,NSMAX
         DO IS=1,NS
            DIAG(IS,IA,IA)  = DIAG(IS,IA,IA) + VOLS(IS)/DTL(IS)
         END DO
      END DO
C
      RETURN
      END
