

      subroutine stmat2zm(vecto,matri,nma1,nma2,nma3,nma4,nmat,nve)

c---------------------------------------------------------------------   
c 
c     stmat(nsegs)        ---> vecto(nve)
c     zm(5,5,2,nsgmax)    ---> matri(nma1,nma2,nma3,nma4)
c 
c 
c     where: 
c 
c             nve   =   nsegs = 50*nsgmax
c 
c             nma1  =   5
c             nma2  =   5
c             nma3  =   2
c             nma4  =   nsgmax
c 
c 
c---------------------------------------------------------------------   

      implicit none

      integer*4
     .  nve,nma1,nma2,nma3,nma4,iseg,j,id1,nmat
      real*8
     .  vecto(nve), matri(nma1,nma2,nma3,nma4)

      
      DO iseg=1,nmat                                        ! NOTE: remember nma4 = nsgmax...!
        
        
        DO j=1,5
          
          id1           = 50*(iseg - 1) + j
          
c                     stmat(id1)    = stmat(id1)    + dmatf(1,j)
c                     stmat(id1+5)  = stmat(id1+5)  + dmatf(2,j)
c                     stmat(id1+10) = stmat(id1+10) + dmatf(3,j)
c                     stmat(id1+15) = stmat(id1+15) + dmatf(4,j)
c                     stmat(id1+20) = stmat(id1+20) + dmatf(5,j)

          matri( j ,   1 ,   1 ,iseg)    = vecto(id1)    
          matri( j ,   2 ,   1 ,iseg)    = vecto(id1+5)  
          matri( j ,   3 ,   1 ,iseg)    = vecto(id1+10) 
          matri( j ,   4 ,   1 ,iseg)    = vecto(id1+15) 
          matri( j ,   5 ,   1 ,iseg)    = vecto(id1+20) 
          
        end do

        
        DO j=1,5

          id1           = 50*(iseg - 1) + 25 + j

c                     stmat(id1)    = stmat(id1)    + dmatf(1,j)
c                     stmat(id1+5)  = stmat(id1+5)  + dmatf(2,j)
c                     stmat(id1+10) = stmat(id1+10) + dmatf(3,j)
c                     stmat(id1+15) = stmat(id1+15) + dmatf(4,j)
c                     stmat(id1+20) = stmat(id1+20) + dmatf(5,j)


          matri( j ,   1 ,   2 ,iseg)    = vecto(id1)    
          matri( j ,   2 ,   2 ,iseg)    = vecto(id1+5)  
          matri( j ,   3 ,   2 ,iseg)    = vecto(id1+10) 
          matri( j ,   4 ,   2 ,iseg)    = vecto(id1+15) 
          matri( j ,   5 ,   2 ,iseg)    = vecto(id1+20) 

        end do
        
        
      end do
      
      END
