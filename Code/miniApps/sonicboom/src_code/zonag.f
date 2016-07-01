
         
      SUBROUTINE ZONAG
c
c*** Cette procedure calcule, a partir d'un maillage, le nombre
c    de zones grossieres, les segments grossiers, le nombre de
c    voisins (et leur numero) d'une cellule en question.
c
      include 'Param3D.h'
      include 'Paramopt3D.h'
c
c--------------------------------------------------------------c
c             parametres d'entree :                            c
c             ---------------------                            c
c                  nsp, ntp, nsegp, nubop, airesp, vncoq       c
c                                                              c
c             parametres de sortie :                           c
c             ----------------------                           c
c                  nubot, nbz, nsegt, nuvoi, nvoi, nuzo        c
c--------------------------------------------------------------c
c
c     Variables locales :
c     -------------------
      integer jaretmg(nzma1,nvoimax,2),ip(2),nuzo1(nzma1)
      integer nubot(nsgma,2),nsgvoi(nnsp,nvoimax)
      integer nsegt(nvma)
      integer k, iz, iseg, nvl, nvg, nsg, nz, ik, ivz
      integer isegt, iz1, iz2, iz1t, iz2t, iv, iza, izb, isegzt
      integer izt, n, ivoiz, ivzt, nubo1, nubo2
      integer nubot1, nubot2, kv, isegz, nuz1, nuz2, nsgt
      integer ntest, nvlplu
c

      intzz(30)=0
      
      data ip/2,1/

c...........................initialisation.............................
c
c          nzma1 = somme des nombres de zones de tous les niveaux,
c          -----
c          nvma  = nb. maxi de niveaux de grilles,
c          ----
c          nsgma = somme de tous les segments de tous les niveaux,
c          -----
c.....................................................................
c
      do 1 k=1,nvoimax
         do 1 iz=1,nzma1
            nuvoi(iz,k)=0
1     continue
c
      do 2 iz=1,nzma1
        nvoi(iz)=0
        nuzo(iz)=0
2     continue
c
      do 3 k=1,nvma
         nbz(k)=0
         nsegt(k)=0
3     continue
c
c
c                  ********************************
c                  *                              *
c                  *  INITIALISATION GRILLE FINE  *
c                  *                              *
c                  ********************************
c
c----------------------------------------------------------------------
c...........nbz(1)=nb zones du 1er niveau                             |
c...........nsegt(1)=nb de segments du 1er niveau                     |
c----------------------------------------------------------------------
c
      nbz(1)=nsp
      nsegt(1)=nsegp
c
c----------------------------------------------------------------------
c...........nubot(iseg,k)=nubop(k,iseg) k=1,2 et iseg=1,nsegt(1)       |
c----------------------------------------------------------------------
c
      do 5 iseg=1,nsegp
         nubot(iseg,1)=nubop(1,iseg)
         nubot(iseg,2)=nubop(2,iseg)
5     continue
      do 6 iseg=nsegp+1,nsgma
         nubot(iseg,1)=0
         nubot(iseg,2)=0
6     continue
c
C [llh]      WRITE(6, *) '5-Methode d"agglomeration des volumes sur la coque: '
C [llh]      WRITE(6, *) '*************************************************** '
C [llh]      WRITE(6, *) ' '
C [llh]      print *,'Niveau le plus fin = 1, nombre de segments =',nsegp
C [llh]      WRITE(6, *) ' '
c
c                  ***********************************
c                  *                                 *
c                  *  FIN INITIALISATION GRILLE FINE *
c                  *                                 *
c                  ***********************************
c
c----------------------------------------------------------------------
c............nvg correspond au niveau le plus grossier                |
c----------------------------------------------------------------------
c
      nvg=niveau+nivg-1
c
C [llh]      print *,' Niveau le plus grossier = ',nivg
C [llh]      WRITE(6, *) ' '
c
c----------------------------------------------------------------------
c............initialisation a zero des parametres de position :       |
c                                                                     |
c       *  pour atteindre les elements du 1er niveau (niveau fin) :   |
c                       nsg,nz = zero                                 |
c       *  pour atteindre les elements du niveau k, k > 1 :           |
c                 nsg = somme des nsegt(j), j=1,k-1                   |
c                 nz  = somme des nbz(j) , j=1,k-1                    |
c----------------------------------------------------------------------
      nsg=0
      nz=0
      ik=-1 
c
c                  ***************************************************
c                  *                                                 *
c                  *  BOUCLE SUR LE NB. DE NIVEAUX GROSSIERS = NVG-1 *
c                  *                                                 *
c                  ***************************************************
c
c
c
      DO 1000 nvl=1,nvg-1
c     ===================
c
c----------------------------------------------------------------------
c                CONSTRUCTION DES TABLEAUX NUVOI ET NVOI              |
c                =======================================              |
c                                                                     |
c        Pour la construction du niveau nvl+1 :                       |
c        ------------------------------------                         |
c                                                                     |
c           _ nz   = somme des nbz(k)  , k=1,nvl-1                    |
c           _ nsg  = somme des nsegt(k) , k=1,nvl-1                   |
c                                                                     |
c             Pour k=1,kmax(=nvoimax), iz=1,nbz(nvl) :                |
c                                                                     |
c           _ nuvoi(nz+iz,k) = numero de la kieme zone voisine iz     |
c                                 (niveau nvl)                        |
c           _ nvoi(nz+iz)    = nb. total de zones voisines a iz       |
c                                 (niveau nvl)                        |
c                (boucle sur les segments niveau nvl)                 |
c----------------------------------------------------------------------
c
         do 100 iseg=1,nsegt(nvl)
            isegt=nsg+iseg
            do 100 k=1,2
               iz1=nubot(isegt,k)
               iz2=nubot(isegt,ip(k))
               iz1t=nz+iz1
               iz2t=nz+iz2
               if (nvoi(iz1t).eq.0) then
                  nuvoi(iz1t,1)=iz2
                  jaretmg(iz1t,1,1)=iz2 
                  jaretmg(iz1t,1,2)=iseg 
                  nvoi(iz1t)=1
               else
                  do 99 iv=1,nvoi(iz1t)
                     if (nuvoi(iz1t,iv).eq.iz2) goto 100
99                continue
                  nvoi(iz1t)=nvoi(iz1t)+1
                  nuvoi(iz1t,nvoi(iz1t))=iz2
                  jaretmg(iz1t,nvoi(iz1t),1)=iz2 
                  jaretmg(iz1t,nvoi(iz1t),2)=iseg 
               endif
100      continue
c
c----------------------------------------------------------------------
c              CONSTRUCTION DES ZONES DU NIVEAU NVL+1                 |
c              ======================================                 |
c                                                                     |
c       Pour iz=1,nbz(nvl) :                                          |
c                                                                     |
c       _ nuzo(nz+iz) = numero de la zone niveau nvl+1 dans laquelle  |
c                       est incluse la zone iz du niveau nvl          |
c       _ izg         = compteur de zones du niveau nvl+1             |
c                                                                     |
c      On construit les differentes zones de la facon suivante :      |
c         nuzo(nz+iz) = izg +1,                                       |
c         Pour tous les voisins ivz de iz (niveau nvl) :              |
c             nuzo(nz+ivz) = nuzo(nz+iz) si nuzo(nz+ivz) valait zero  |
c                   (ie. on inclut ivz dans la nouvelle zone)         |
c                            passage au voisin suivant sinon          |
c                            (ie. ivz deja precedemment incluse dans  |
c                                   une nouvelle zone                 |
c----------------------------------------------------------------------
c
c--- Balayage symetrique

          ik=ik*(-1) 
          if (ik.eq.1) then 
            iza=1 
            izb=nbz(nvl) 
          else 
            iza=nbz(nvl) 
            izb=1 
          endif 

         do 200 iz=iza,izb,ik 
            izt=nz+iz
            if (nuzo(izt).ne.0) goto 200
            nbz(nvl+1)=nbz(nvl+1)+1
            nuzo(izt)=nbz(nvl+1)
c
c----------------------------------------------------------------------
c    n = compteur local de zones voisines a iz qui vont contribuer    |
c        a construire la nouvelle zone no izg+1                       |
c                                                                     |
c    Le traitement du cas n=0 suivant n'est plus pris en compte  :    |
c        Si n = 0, alors la nvelle zone se reduit a iz, ce qu'on veut |
c                  eviter.                                            |
c        Dans ce cas, on supprime cette nvelle zone en diminuant de 1 |
c        le compteur de zones et en affectant a iz le numero de zone  |
c        de son premier voisin                                        |
c----------------------------------------------------------------------
            n=0
            do 199 ivz=1,nvoi(izt)
              ivoiz=nuvoi(izt,ivz)
              ivzt=nz+ivoiz
              if (nuzo(ivzt).ne.0) goto 199
              nuzo(ivzt)=nbz(nvl+1)
              n=n+1
199         continue
c
c    Le traitement du cas n=0 suivant n'est plus pris en compte :
c              if (n.eq.0) then
c                 ivz1=nuvoi(izt,1)
c                 ivzt1=nz+ivz1
c                 nuzo(izt)=nuzo(ivzt1)
c                 nbz(nvl+1)=nbz(nvl+1)-1
c              endif
c
200      continue

        if (ik.eq.-1) then 
           do 4000 iz=1,nbz(nvl) 
             izt=nz+iz 
             nuzo1(izt)=nbz(nvl+1)+1-nuzo(izt) 
4000        continue 
           do 5000 iz=1,nbz(nvl) 
             izt=nz+iz 
             nuzo(izt)=nuzo1(izt) 
5000        continue 
         endif 
c
C [llh]      WRITE(6, *) 'Sur le niveau',nvl+1,'   le nombre de zones vaut :'
C [llh]      WRITE(6, *) ' nbz(',nvl+1,') = ',nbz(nvl+1) 
C [llh]      WRITE(6, *) ' '
c
c----------------------------------------------------------------------
c             CONSTRUCTION DES SEGMENTS                               |
c             =========================                               |
c                                                                     |
c       Pour chaque segment iseg=1,nsegt(nvl) ,on a                   |
c                                                                     |
c         _ nubo1,nubo2 correspondant aux numeros de zones niveau nvl |
c                       des extremites du segment iseg                |
c         - ou bien nubo1 et nubo2 sont dans 2 zones niveau nvl+1     |
c                   differentes et alors :                            |
c             * on conserve ce segment pour le niveau nvl+1 pour      |
c               representer la separation de ces 2 zones s'il n'y a   |
c               pas deja eu de creation de ce representant et dans    |
c               ce cas on incremente de 1 le compteur de segments du  |
c               niveau nvl+1, le nouveau segment cree a pour          |
c               extremites les numeros de ces 2 zones.                |
c                                                                     |
c         _ ou bien nubo1 et nubo2 sont dans une meme zone du niveau  |
c                   nvl+1 et alors on passe au segment iseg du niveau |
c                   nvl suivant.                                      |
c                                                                     |
c   VERSION RAPIDE DE CALCUL DES SEGMENTS GROSSIERS  26/10/88         |
c           avec les segments voisins aux zones ...                   |
c                                                                     |
c----------------------------------------------------------------------
c         
         do 41 kv=1,nvoimax
            do 42 iz=1,nbz(nvl+1)-1
               nsgvoi(iz,kv)=0
42    continue
41    continue
C
c                   -- Boucle sur les segments du niveau fin --
c
         DO 300 iseg=1,nsegt(nvl)
            isegt=nsg+iseg
            nubo1=nubot(isegt,1)
            nubo2=nubot(isegt,2)
            nubot1=nz+nubo1
            nubot2=nz+nubo2
            iz1=nuzo(nubot1)
            iz2=nuzo(nubot2)
            if (iz1.eq.iz2) goto 300
c
c                    -- Recherche du segment suivant --
c
            iz = iz1
298         kv = 0
  299       kv = kv+1
            if (kv.gt.intzz(30)) intzz(30)=kv
            if (kv.gt.nvoimax) print*,'Attention augmenter nvoimax!'
            isegz = nsgvoi(iz,kv)
            if (isegz.eq.0) then
               nsgvoi(iz,kv)=nsegt(nvl+1)+1
               if (iz.eq.iz1) then
                  iz=iz2
                  goto 298
               else
                  goto 310
               endif
            endif
            isegzt=nsg+nsegt(nvl)+isegz
            nuz1=nubot(isegzt,1)
            nuz2=nubot(isegzt,2)
            if (iz1.eq.nuz1.and.iz2.eq.nuz2) then
                  ntest=1
            else
                  if (iz1.eq.nuz2.and.iz2.eq.nuz1) then
                      ntest=1
                  else
                     ntest=0
                  endif
            endif
c
c                      -- Le segment existe deja --
c
            if (ntest.eq.1) goto 300
            goto 299
c
c                      -- Creation du segment suivant --
c
310       continue
            nsegt(nvl+1)=nsegt(nvl+1)+1
            nsgt=nsg+nsegt(nvl)+nsegt(nvl+1)
            nubot(nsgt,1)=iz1
            nubot(nsgt,2)=iz2
300      continue
c
C [llh]      WRITE(6, *) 'Sur le niveau',nvl+1,'  le nombre de segments vaut :'
C [llh]      WRITE(6, *) 'nsegt(',nvl+1,') = ',nsegt(nvl+1) 
C [llh]      WRITE(6, *) ' '
C [llh]      WRITE(6, *) ' '
c
c----------------------------------------------------------------------
c NOUVELLES VALEURS DES PARAMETRES DE POSITION POUR LE NIVEAU SUIVANT |
c----------------------------------------------------------------------
c
         nz=nz+nbz(nvl)
         nsg=nsg+nsegt(nvl)
c
      if(nbz(nvl+1).eq.1)then
         nvlplu=nvl+1-1
         print *,'Niveau trop grossier'
         print *,'Changer le nombre de niveaux a',nvlplu
         stop
         endif
c
1000  CONTINUE
c     ========

         do 114 iseg=1,nsegt(nvg)
            isegt=nsg+iseg
            do 114 k=1,2
               iz1=nubot(isegt,k)
               iz2=nubot(isegt,ip(k))
               iz1t=nz+iz1
               iz2t=nz+iz2
               if (nvoi(iz1t).eq.0) then
                  nuvoi(iz1t,1)=iz2
                  jaretmg(iz1t,1,1)=iz2 
                  jaretmg(iz1t,1,2)=iseg 
                  nvoi(iz1t)=1
               else
                  do 109 iv=1,nvoi(iz1t)
                     if (nuvoi(iz1t,iv).eq.iz2) goto 114
109                continue
                  nvoi(iz1t)=nvoi(iz1t)+1
                  nuvoi(iz1t,nvoi(iz1t))=iz2
                  jaretmg(iz1t,nvoi(iz1t),1)=iz2 
                  jaretmg(iz1t,nvoi(iz1t),2)=iseg 
               endif
114       continue
c


          return
      end
