

      SUBROUTINE MAIL2D
c
      include 'Paramopt3D.h'
c
      integer jaretp(nnsp,nvoimax,2),nor(3),nex(3),nop(3)
      integer i, is, jt, k, j, kv, is1, is2, is3, inew
      integer ig0, nub1, nub2, nunit, iseg
      real*8 eps

c     LECTURE DE LA GEOMETRIE sur fort.88
c     -----------------------------------
c      nunit=88
c      rewind(nunit)
c      read(nunit,12) nsp,ntp
c      read(nunit,13) ((coorp(i,is),i=1,3),is=1,nsp)
c      read(nunit,14) ((nup(k,jt),k=1,3),jt=1,ntp)
c      read(nunit,15) (logfrp(is),is=1,nsp)
12    format(2i5)   
13    format(8e15.8)
14    format(8i5)
15    format(i5)
c
c     Calcul du nombre de segments (nsegp) ainsi que des extremites
c     des segments (nubop(1,.) et nubop(2,.)) a l'aide des tableaux
c     jaretp(nnsp,kv,2)

c     Initialisations :
c     -----------------

      nor(1)=1
      nex(1)=2
      nop(1)=3
c
      nor(2)=2
      nex(2)=3
      nop(2)=1
c
      nor(3)=3
      nex(3)=1
      nop(3)=2
c


      do j=1,nnsegp
        do k=1,2
           nubop(k,j)=0.
        end do
      end do

      do 11 is=1,nsp
         do 6 kv=1,nvoimax
            jaretp(is,kv,1)=0
            jaretp(is,kv,2)=0
6        continue
11    continue
c
      eps=1./3.
      nsegp=0
c                               -- boucle / elements
      do 1000 jt=1,ntp
c                               -- boucle / aretes 
         do 500 k=1,3
            is1=nup(nor(k),jt)
            is2=nup(nex(k),jt)
            is3=nup(nop(k),jt)
c                               -- inew=1 si nouvelle arete,2 sinon
      inew=2
c                               -- is2 est-il voisin de is1 ?
      do 100 kv=1,nvoimax
         if (jaretp(is1,kv,1).eq.0) goto 110
         if (jaretp(is1,kv,1).eq.is2) then
            ig0 = jaretp(is1,kv,2)
            nub1= iabs(nubop(1,ig0))
            nubop(1,ig0)=-nubop(1,ig0)
            nub2 = nubop(2,ig0)
            if (is1.eq.nub1.and.is2.eq.nub2) then
               eps=abs(eps)
            else
               if (is1.eq.nub2.and.is2.eq.nub1) then
                  eps=-abs(eps)
               else
                  print *,' ERREUR SEG2 ',is1,is2,nub1,nub2
               endif
            endif
            goto 200
         endif
100   continue
c                               -- is1 est-il voisin de is2 ?
110   do 120 kv=1,nvoimax
         if (jaretp(is2,kv,1).eq.0) goto 130
         if (jaretp(is2,kv,1).eq.is1) then
            ig0=jaretp(is2,kv,2)
            nub1=iabs(nubop(1,ig0))
            nubop(1,ig0)=-nubop(1,ig0)
            nub2=nubop(2,ig0)
            if (is1.eq.nub1.and.is2.eq.nub2) then
               eps=abs(eps)
            else
               if (is1.eq.nub2.and.is2.eq.nub1) then
                  eps=-abs(eps)
               else
                  print *,' ERREUR SEG2 ',is1,is2,nub1,nub2
               endif
            endif
            goto 200
         endif
120   continue
c                          -- L'arete [is1,is2] n'exite pas : on l'ajoute
      kv=kv+1
      if (kv.gt.nvoimax) then
         print *,' augmenter nb de voisins max ',kv
      endif
c
130   nsegp=nsegp+1

      inew=1
c
      if (nsegp.gt.nnsegp) then
         print *,' augmenter nnsegp = ',nnsegp
      endif
      jaretp(is2,kv,1)=is1
      jaretp(is2,kv,2)=nsegp
      nubop(1,nsegp)=-is1
      nubop(2,nsegp)=is2
      ig0=nsegp
      eps=abs(eps)
c
200   continue

500   continue
1000  continue
c
      do 2000 iseg=1,nsegp
         nub1=nubop(1,iseg)
         nub2=nubop(2,iseg)
         if (nub1.lt.0) then
            nubop(1,iseg)=-nubop(1,iseg)
         endif
         if (nub2.lt.0) nubop(2,iseg)=-nubop(2,iseg)
2000  continue
c
C [llh]      print *,'Le nombre de segments sur la coque est egal a:',nsegp
C [llh]      WRITE(6, *) ' '
C [llh]      WRITE(6, *) ' '
c
      return
      end
