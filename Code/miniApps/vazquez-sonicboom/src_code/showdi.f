      subroutine showdi
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c---------------------------------------------------------------------

      integer*4 nvmaxi,ndmaxi,is

      nvmaxi=0
      do is=1,ns
        if (nbvoi(is).gt.nvmaxi) nvmaxi=nbvoi(is)
      end do
      
      
      write (6,*) ' --->>>> RECOMMENDED DIMENSIONS : '
      write (6,*) '================================ '
      write (6,1000) '   nsmax   =  ', ns + 3,
     .  '   ---> actual value =',nsmax
      write (6,1000) '   ntmax   =  ', nt + 3,
     .  '   ---> actual value =',ntmax
      write (6,1000) '   nfcmax  =  ', nfac+ 3,
     .  '   ---> actual value =',nfcmax
      write (6,1000) '   nsgmax  =  ', nseg+ 3 ,
     .  '   ---> actual value =',nsgmax
      write (6,1000) '   nvmax   =  ', nvmaxi+ 3,
     .  '   ---> actual value =',nvmax
      write (6,1000) '   nnsp    =  ', nsp + 3 ,
     .  '   ---> actual value =',nnsp
      write (6,1000) '   nntp    =  ', ntp + 3,
     .  '   ---> actual value =',nntp
      write (6,1000) '   nnsegp  =  ', nsegp+ 3,
     .  '   ---> actual value =',nnsegp
      write (6,1000) '   nvoimax =  ', intzz(30)+ 3,
     .  '   ---> actual value =',nvoimax
      write (6,1000) '================================ '

 1000 format(a,i10,a,i10,a)
      

 
      end
      
