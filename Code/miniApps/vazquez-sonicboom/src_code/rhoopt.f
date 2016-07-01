      SUBROUTINE RHOOPT(CTRL,fic)
c-----------------------------------------------------------
C     Calcul du pas optimal dans une
c     recherche monodimensionnelle
c      Algorithme de type RoAlg7 du livre de J. Cea
c-----------------------------------------------------------
c Entrees :
c  Rinit: une estimation du pas Ro
c  Train: valeur initiale de la fonctionnelle (R=0)
c  Itopt: numero d'iteration de gradient
c
c Ce programme appelle une routine fcnelj(R,FJ) qui
c a un pas R fait correspondre un scalaire FJ
c FJ (R) = fonctionnelle( v - R Grad )
c les autres variables de fcnel passent en common.
c-----------------------------------------------------------
      include 'Param3D.h'
      include 'Paramopt3D.h'
c
      INTEGER nunit, kpas, fic, kdiv,kpamax
      REAL*8 rk, rk0, q0, q1, q2, q4
      REAL*8 rk1, rk2, rk4, FJ, anum, bnum, ctrl(nnsp)
c
      kdiv = 0
      Rk = Roinit
      Rk0 = 0.
c
c
C *** On prend ici comme cout de reference le minimum entre
C     le cout calcule avec res. syst. lin. et le cout calcule 
C     avec moins d'iterations 
C  
      PRINT*,'======================================= '
      PRINT*,' '
      PRINT*,'         ENTREE DANS RHOOPT'
      PRINT*,'         =================='
      PRINT*,' '
      PRINT*,'======================================= '
      PRINT*,'Iteration ',itopt
      WRITE(6, *) ' '
c

      kpamax=20
      
      if (itopt.eq.1 .or. imvng.eq.1) then
c...    Initialize steepest descent: either because the cycle starts (itopt=1)
c...     or because bruno has just moved the mesh (imvng=1)
        CALL fcnelj(CTRL,Rk0,FJ,fic)
        COUTP = FJ
        write (6,*) 'holaaaa' ,roinit, coutp, cout,imvng
      endif
      kpas=1
c
      

      Q0 = min(COUT,COUTP)
c
      FJ = Q0
      Rk1 = Rk
c
10    CONTINUE
c

      CALL fcnelj(CTRL,Rk1,FJ,fic)
      Q1 = FJ
c

      IF (Q1 .lt. Q0)  THEN
c                      ****
c       1. Descente au premier pas
c       **************************
c
      PRINT*,' '
      print*,'  RHOOPT',itopt,':      DESCENTE,  pas',kpas
      PRINT*,'======================================= '
      PRINT*,' '
         print*,'Q0 =',Q0,' Q1 = ',Q1
c
           Rk2 = 2.0*Rk1
c           Rk2 = Roinit+Rk1
           CALL fcnelj(CTRL,Rk2,FJ,fic)
           Q2 = FJ
c
20         CONTINUE
c
           IF (Q2 .lt. Q1 .and. kpas.le.kpamax) THEN
c                          ====
c            1.1 La descente continue
c            ========================
c
      kpas=kpas+1
c
      PRINT*,' '
      print*,'  RHOOPT',itopt,':     LA DESCENTE CONTINUE,  pas',kpas
      PRINT*,'======================================= '
      PRINT*,' '
              print*,'Q1 = ',Q1,' Q2 = ',Q2
              print*,'Rk0 =',Rk0,' Rk1 = ',Rk1,' Rk2 = ',Rk2
c
               Rk0 = Rk1
               Rk1 = Rk2
               Rk2 = 2.*Rk1
c               Rk2 = Roinit+Rk1
               Q0 = Q1
               Q1 = Q2
c
               print*,'Q0 = ',Q0,' Q1 = ',Q1
c
               CALL fcnelj(CTRL,Rk2,FJ,fic)
               Q2 = FJ
               GOTO 20
c
            ELSE
c           ====
c            1.2 Fin de descente
c            ===================
c
c                  1.2.1 Interpolation 
c
      PRINT*,' '
      print*,'        RHOOPT',itopt,': FIN DE DESCENTE,  pas',kpas
      PRINT*,'======================================= '
               print*,'Q2 = ',Q2,' Q1 = ',Q1
c
               if (kpas.gt.kpamax) then

                 write(6,*) 'Maximal iteration step reached= ',kpas

               end if
               
             bnum = Rk0*Rk0*(Q1-Q2) -Q0*(Rk1*Rk1-Rk2*Rk2)
     &              + Rk1*Rk1*Q2 - Rk2*Rk2*Q1
             anum = Q0*(Rk1-Rk2) -Rk0*(Q1-Q2) + Q1*Rk2 -Rk1*Q2
c
             print*,'bnum = ',bnum,' anum = ',anum
c
             Rk4 = -0.5*bnum/anum
c
             print*,'Rk4 = ',Rk4
c
             CALL fcnelj(CTRL,Rk4,FJ,fic)
             Q4 = FJ
c
c                  1.2.2 Controle Interpolation
c On choisit le parametre qui rend la fonctionnelle 
c la plus petit entre le minimum dichotomique
c et le minimum interpole
c
      PRINT*,' '
      print*,'        RHOOPT',itopt,': CONTROLE INTERPOLATION'
      PRINT*,'======================================= '
             print*,'Q4 = ',Q4,' Q1 = ',Q1
c
             if (Q4 .lt. Q1) then
                 Roopt = Rk4
                 Roinit = Roopt
                 COUTP = Q4
                             else
                 Roopt = Rk1
                 Roinit = Roopt
                 COUTP = Q1
             endif
c
c 
             open(34,ACCESS='APPEND')
             write(34,100) itopt,ncf,Roinit,COUTP
 100         FORMAT(2i5,2e15.8)
             close(34)

      ENDIF
c     ====
c     
      ELSE
c     ****
c            2. premier pas trop grand
c            *************************
c
      PRINT*,' '
      print*,'        RHOOPT',itopt,':  MONTEE'
      PRINT*,'======================================= '
         print*,'Q0 = ',Q0,' Q1 = ',Q1
c



           Rk1 = 0.5*Rk1
c           Rk1 = 0.1*Rk1
           kdiv = kdiv +1
           if(kdiv.ge.6)then
              print*,'RHOOPT: Mauvaise direction de descente:stop'
C [llh Nan]              call zzoutsur(ctrl,1)         
              stop
              endif
           GOTO 10
      ENDIF
c     *****
c
      RETURN
      END
