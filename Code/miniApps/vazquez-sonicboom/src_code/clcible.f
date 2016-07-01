      SUBROUTINE CLCIBLE(CTRL)
C
c     Cette procedure calcule la portance CL desiree dans le cas de l'aile
c     M6 et du Falcon, mais pas dans le cas de la tuyere 3D.
c     On demarre avec itrans = 1 et CTRL=0. tetacdcl = 3.06 degres.
c
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'

      INTEGER ISP, IS, nunit
      REAL*8 ctrl(nnsp), Pres(nsmax), pression,cdaux,cdpct
      REAL*8 cx, cy, cz, cyz, vncoq(3,nnsp), denominateur,wingsur

      INTEGER ii,iwhico
      REAL*8 sol(nsmax,9),coefr(5),vcarre,denom,coqtri,coqsur
c      REAL*8 cd,cl

      CX=0.
      CY=0.
      CZ=0.



      CALL ETAT(Ctrl,ktmax)

c
c*****************
      if (ncont.eq.0) then
c******************
c
      coefr(1)                     = rhoref
      coefr(2)                     = rhoref*vref
      coefr(3)                     = rhoref*vref
      coefr(4)                     = rhoref*vref
      coefr(5)                     = pref

      DO is = 1,ns
         sol(is,1) = ua(1,is)
         sol(is,2) = ua(2,is)/ua(1,is)
         sol(is,3) = ua(3,is)/ua(1,is)
         sol(is,4) = ua(4,is)/ua(1,is)
         sol(is,5) = ua(5,is)
         sol(is,6) = gam1*(sol(is,5)*coefr(5) - 0.5*sol(is,1)*(sol(is,2)
     &        *sol(is,2) + sol(is,3)*sol(is,3) + sol(is,4)*sol(is,4)))
         sol(is,7) = SQRT( (sol(is,2)*sol(is,2)+sol(is,3)*sol(is,3)+
     &               sol(is,4)*sol(is,4))/(gam*sol(is,6)/sol(is,1)))
         vcarre = sol(is,2)**2+sol(is,3)**2+sol(is,4)
         denom = roin*0.5*vcarre
         sol(is,8) = (pin - sol(is,6))/denom
         sol(is,9) = (sol(is,6)/pin)*(roin/sol(is,1))**gam - 1.
      END DO
c
c  Sauvegarde de la solution dans fort.21
c
      print *,'Sauvegarde solution-cible dans fort.21'

      Do ii = 1,9
         Do is = 1,ns
           WRITE(21,113) sol(is,ii)
         End Do
      End Do

113   FORMAT(e15.8)

      write(59,112) kt, t

112    FORMAT(i6,2x,e12.5)

      endif
c
C********************
c*********************
c
      
      CALL NORMCOQ(CTRL,VNCOQ)

      do is = 1,ns
         pdesp(is) = 0.
      end do


      coqsur=0.d0
      Do isp = 1,nsp
        is = node2d3d(isp)
        pression = gam1*(ua(5,is) - 0.5*(ua(2,is)**2 + ua(3,is)**2 +
     .    ua(4,is)**2)/ua(1,is))
        pdesp(is) = pression
        
        cx = cx + (vncoq(1,isp)*pression)
        cy = cy + (vncoq(2,isp)*pression)
        cz = cz + (vncoq(3,isp)*pression)
        
        coqtri= vncoq(1,isp)*vncoq(1,isp)
     .    + vncoq(2,isp)*vncoq(2,isp)
     .    + vncoq(3,isp)*vncoq(3,isp)        
        coqsur= coqsur+sqrt(coqtri)

      end do
      
      Denominateur = roin*(uxin**2+uyin**2+uzin**2)/2.      !dynamic pressure (at inf) 

      wingsur=faczz(30)
      if (wingsur.lt.0.d0) wingsur=faczz(25)                ! <<---- a constant number...

      wingsur=coqsur                                        ! <<---- ... or the shell surface

      CX = CX/Denominateur/wingsur
      CY = CY/Denominateur/wingsur
      CZ = CZ/Denominateur/wingsur
c
c    Calcul des coefficients de trainee CDTARGET et de portance CLTARGET
c
      iwhico= int(faczz(6))                                 
      cyz= cy
      if (iwhico.eq.3) cyz= cz

      
      cltarget= -sin(tetacdcl)*CX + cos(tetacdcl)*cyz
      cdaux   =  cos(tetacdcl)*CX + sin(tetacdcl)*cyz

      
      
      write(6,*)
      write(6,*) 'Target values of the Aerodynamic coefficients:'      
      write(6,*)

      cdpct=1.d0
      if (igicc.eq.2) then
        cdpct   = cdtarget
        cdtarget= cdtarget*cdaux
        write(6,*) '     CDTar = ',cdtarget,
     .    ' (  ',cdpct,'  of the original one) '
        write(6,*) '     CLTar = ',cltarget
        write(6,*)
      else        
        write(6,*) '     CDTar = ',cdtarget
        write(6,*) '     CLTar = ',cltarget
        write(6,*)        
      end if
        
      npres = 2
      itrans = 1

      RETURN
      END
