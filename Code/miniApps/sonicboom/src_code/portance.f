      SUBROUTINE PORTANCE(CTRL,CCD,CCL)
c---------------------------------------------------------------------   
c  Calcul des coefficients de portance et de trainee
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c---------------------------------------------------------------------   
c
      REAL*8 cx, cy, cz, cyz,ccd, ccl, vncoq(3,nnsp), ctrl(nnsp)
      REAL*8 denominateur, pression,wingsur,coqsur,coqtri
      INTEGER isp, is,iwhico
c
      CX=0.d0
      CY=0.d0
      CZ=0.d0
      iwhico= int(faczz(6))                                 
c
c    Calcul des normales a la peau
c
      CALL NORMCOQ(CTRL,VNCOQ)
c
      coqsur=0.d0
      Do 100 isp = 1,nsp
         is = node2d3d(isp)
         pression = gam1*(ua(5,is) - 0.5*(ua(2,is)**2 + ua(3,is)**2 +
     $              ua(4,is)**2)/ua(1,is))
c
         cx = cx + (vncoq(1,isp)*pression)
         cy = cy + (vncoq(2,isp)*pression)
         cz = cz + (vncoq(3,isp)*pression)

         coqtri= vncoq(1,isp)*vncoq(1,isp)
     .     + vncoq(2,isp)*vncoq(2,isp)
     .     + vncoq(3,isp)*vncoq(3,isp)

         coqsur= coqsur+sqrt(coqtri)
c
 100  continue
c
      Denominateur = roin*(uxin**2+uyin**2+uzin**2)/2.      !dynamic pressure (at inf) 

      wingsur=faczz(30)
      if (wingsur.lt.0.d0) wingsur=faczz(25)                ! <<---- a constant number...

      wingsur=coqsur                                        ! <<---- ... or the shell surface
      
      CX = CX/Denominateur/wingsur
      CY = CY/Denominateur/wingsur
      CZ = CZ/Denominateur/wingsur

c
c
c    Calcul des coefficients de trainee CD et de portance CL
c
      cyz= cy
      if (iwhico.eq.3) cyz= cz
      
      ccl= -sin(tetacdcl)*cx + cos(tetacdcl)*cyz
      ccd=  cos(tetacdcl)*cx + sin(tetacdcl)*cyz
c
c      write(6,*) 'CD = ',CCD,' CL = ',CCL
c

      
      
      END
