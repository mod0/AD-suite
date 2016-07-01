
C     -----------------------------------------------
      SUBROUTINE FDECFEX ( CMULT, GAMMA, RHO, UN, VN, WN, P, E, FM )
C     -----------------------------------------------
C
C     Calcule tous les decentrages dans tous les sens ;-)
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
C
      REAL*8 GAMMA, CMULT
      REAL*8 RHO, E, P
      REAL*8 FM(5), C, UN, VN, WN
C
C     Variables locales
C
C     Tableaux de travail
C
C     Indices de boucle
C
C     Divers
      REAL*8 ACOEFF, BCOEFF
      REAL*8 GAM1, GAM2, GAM3, GAM4
C
C     Procedures et Fonctions
C
C     -----------------------------------------------
C     Fin des declarations
C     -----------------------------------------------
C
      C = SQRT( GAMMA * P / RHO )
C
      GAM1 = GAMMA - 1.
      GAM2= 1./( GAMMA**2 - 1. )
      GAM3= 1./GAMMA
      GAM4= 1./GAM1
C
      ACOEFF = 0.5*( 1. - CMULT* SIGN( 1., (UN/C - CMULT) ) )
      BCOEFF = 0.5*( 1. + CMULT* SIGN( 1., (UN/C + CMULT) ) )
c debug
c      print*,'fdecfex1:', C , ACOEFF,BCOEFF ,cmult,rho,un
C
      FM(1)= BCOEFF * ( (1. - ACOEFF) * RHO*UN +
     $                  CMULT / 4.* ACOEFF * RHO *
     $                  C * ( UN/C + CMULT )**2 )
c debug
c       print*,'fdecfex2:',FM(1)
C
      FM(2)= BCOEFF * ( (1. - ACOEFF) * 
     $                  ( RHO* UN**2 + P ) +
     $                  CMULT/4.* ACOEFF * RHO * C*
     $                  (UN/C + CMULT )**2 *
     $                  ( GAM1 * UN + CMULT * 2.* C )* GAM3 )
C
      FM(3)= BCOEFF * ( (1. - ACOEFF) *
     $                  RHO * UN * VN +
     $                  CMULT/ 4.* ACOEFF* RHO * C*
     $                  ( UN/C + CMULT )**2 * VN )
C
      FM(4)= BCOEFF * ( (1. - ACOEFF) *
     $                  RHO * UN * WN +
     $                  CMULT/ 4.* ACOEFF* RHO * C*
     $                  ( UN/C + CMULT )**2 * WN )
C
      FM(5)= BCOEFF * ( ( 1. - ACOEFF ) *
     $                  ( E + P ) * UN +
     $                  CMULT /8.* ACOEFF* RHO * C*
     $                  ( UN/C + CMULT )**2 *
     $                  ( VN**2 + WN**2 + ( GAM1*UN + CMULT * 2.* C)**2*
     $                    GAM2 ) )
C
c debug
c       print*,'fdecfex5:',FM(1)

      RETURN
      END
