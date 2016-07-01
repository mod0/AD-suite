      SUBROUTINE INIM3D
c---------------------------------------------------------------------
c Initialises the physical solution depending on the problem under 
c consideration   
c---------------------------------------------------------------------
      INCLUDE 'Param3D.h'
c---------------------------------------------------------------------
c     Local variables definition
      INTEGER is, ivar,iinex,idumy,itake
      INTEGER is1, is2, is3, if11, if2, ifac,ihola,id,idi,idf
      character*4 wopos
      REAL*8    xcoor, q(nsmax,8),uval(5),riter
c
      DO 5 is=1,ns
         adh1(is)                  = 1.0
5     CONTINUE  
c
c*** Initialisation de l'etat-adjoint
c
      DO is = 1,ns
         Do ivar = 1,5
            Piadj(ivar,is) = 0.
         End Do
      End Do
c
      IF (ncont .EQ. 0) THEN
c     
         IF (coefm1 .EQ. 1) THEN
c
      WRITE(6, *) 'On part d"une solution uniforme '
      WRITE(6, *) ' '
c
c           Shock tube problem
c
            DO 10 is=1,ns
c
               xcoor               = coor(1,is)
c
               IF (xcoor .LT. 0.5) THEN
c
                  ua(1,is)         = roin/rhoref
                  ua(2,is)         = roin*uxin/(rhoref*vref)
                  ua(3,is)         = roin*uyin/(rhoref*vref)
                  ua(4,is)         = roin*uzin/(rhoref*vref)
                  ua(5,is)         = (0.5*(uxin*uxin + uyin*uyin +
     &                                     uzin*uzin)*roin + 
     &                                pin/gam1)/pref
c
               ELSE
c
                  ua(1,is)         = roout/rhoref
                  ua(2,is)         = roout*uxout/(rhoref*vref)
                  ua(3,is)         = roout*uyout/(rhoref*vref)
                  ua(4,is)         = roout*uzout/(rhoref*vref)
                  ua(5,is)         = (0.5*(uxout*uxout + uyout*uyout +
     &                                     uzout*uzout)*roout + 
     &                                pout/gam1)/pref
c
               ENDIF
c
10          CONTINUE  
c
         ENDIF
c
         IF ((coefm1 .EQ. 2) .OR. (coefm1 .EQ. 5) .OR.
     $        (coefm1 .EQ. 10)) THEN
c
c           External flow 
c
            DO 40 is=1,ns
c
               ua(1,is)            = roin/rhoref
               ua(2,is)            = roin*uxin/(rhoref*vref)
               ua(3,is)            = roin*uyin/(rhoref*vref)
               ua(4,is)            = roin*uzin/(rhoref*vref)
               ua(5,is)            = (0.5*(uxin*uxin + uyin*uyin + 
     &                                     uzin*uzin)*roin + 
     &                                pin/gam1)/pref
c
40          CONTINUE

            iinex=intzz(1)
            
            if (iinex.eq.2) then
c
c...  Modify some values from an external gid file:              
c              
              write(6,*) 'Reading INIT.RES (GiD results file)'
                                                
              open (11,file='init.res',status='old')

              ihola=0

              
              do while (ihola.eq.0) 

                read(11,1000,end=500)
     .            wopos,idumy,riter,idumy,idumy,idumy

                itake = 0
                idi   = 1
                idf   = 1
                if (wopos.eq.'VELO') then                  
                  itake = 1
                  idi   = 2
                  idf   = 4
                else if (wopos.eq.'DENS') then
                  itake = 1
                  idi   = 1 
                  idf   = idi
                else if (wopos.eq.'ENER') then
                  itake = 1
                  idi   = 5
                  idf   = idi
                end if

                ivar=0
                
                do is=1,ns
                  read(11,*) idumy,(uval(id),id=idi,idf) 
                  
                  if (logfr(is).eq.5 .and.
     .              coor(3,is).gt.-600.0d0) then
                    ivar=ivar+1
                    do id=idi,idf                    
                      if (itake.eq.1)
     .                  ua(id,is)=uval(id)                    
                    end do
                  end if

                end do                                        

              end do

  500         close(11)

              
              
            end if



            
         ENDIF                
c
         IF (ivis .EQ. 1) THEN 
c
            DO 25 if11=1,nf11
c
               ifac                = noe1(if11)
c
               is1                 = nsfac(1,ifac)
               is2                 = nsfac(2,ifac)
               is3                 = nsfac(3,ifac)
c
               ua(2,is1)           = 0.0
               ua(3,is1)           = 0.0
               ua(4,is1)           = 0.0
               ua(5,is1)           = ua(1,is1)*tbrd
c
               ua(2,is2)           = 0.0
               ua(3,is2)           = 0.0
               ua(4,is2)           = 0.0
               ua(5,is2)           = ua(1,is2)*tbrd
c
               ua(2,is3)           = 0.0
               ua(3,is3)           = 0.0
               ua(4,is3)           = 0.0
               ua(5,is3)           = ua(1,is3)*tbrd
c
               adh1(is1)           = 0.0
               adh1(is2)           = 0.0
               adh1(is3)           = 0.0
c
25          CONTINUE
c
         ENDIF
c
                           
      ELSE
c
c        Reinitialization strategy
c        Reading initial physical solution in global representation
c
c****** Quand le programme se termine normalement, il s'affiche dans fort.58
c       La solution finale de l'ecoulement.
c
c         REWIND(58)
c
c         DO 30 is=1,ns
c            READ(58, 115) is1, ua(1,is), ua(2,is), ua(3,is),
c     &                         ua(4,is), ua(5,is)
c30       CONTINUE
c         READ(58, 112) kt0, t0
c
c         CLOSE(58)
c
c****** Dans le fichier fort.21, il y a la solution de l'ecoulement pour
c       calculer la pression cible et la trainee et portance cibles. Dans
c       le fichier fort.59, il s'affiche le numero de l'iteration en temps
c       et le temps physique.
c  
         read(21,113) ((q(is,ivar),is = 1,ns),ivar = 1,8)
c
113      FORMAT(e15.8)
c
         do is = 1,ns
               ua(1,is) = q(is,1)
               ua(2,is) = q(is,2)*q(is,1)
               ua(3,is) = q(is,3)*q(is,1)
               ua(4,is) = q(is,4)*q(is,1)
               ua(5,is) = q(is,5)
         end do
c
         READ(59, 112) kt0, t0
c
c         REWIND(57)
c         READ(57,110) som1, dro1
c         CLOSE(57)
c
c110      FORMAT(2e14.7)
c
112      FORMAT(i6,2x,e12.5)
115      FORMAT(i6,5(e15.8,1x))

         
         WRITE(6, *) 'On part d"une solution non uniforme'
         WRITE(6, *) 'L"iteration en temps de depart vaut :',kt0
         WRITE(6, *) ' '
c
      ENDIF
c
 1000 format(a15,i5,e12.5,3i5)
 1200 format(1x,i8,3(2x,e12.5))
 2000 format(17h MESH   dimension,i3,26h Elemtype  Triangle  Nnode,i3)
 3000 format(17h MESH   dimension,i3,28h Elemtype  Tetrahedra  Nnode,i3)

      
      do is=1,ns
        if (logfr(is).eq.1 .or. logfr(is).eq.4) then
          
          ua(1,is)            = roin/rhoref
          ua(2,is)            = roin*uxin/(rhoref*vref)
          ua(3,is)            = roin*uyin/(rhoref*vref)
          ua(4,is)            = roin*uzin/(rhoref*vref)
          ua(5,is)            = (0.5*(uxin*uxin + uyin*uyin + 
     &      uzin*uzin)*roin + 
     &      pin/gam1)/pref
          
        else if (iinex.eq.2) then
          
          ua(2,is)            = ua(2,is)*ua(1,is)
          ua(3,is)            = ua(3,is)*ua(1,is)
          ua(4,is)            = ua(4,is)*ua(1,is)
          
        end if
        
      end do

      
      END
