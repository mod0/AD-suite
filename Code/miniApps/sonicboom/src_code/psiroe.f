      subroutine psiroe(ctrl,ctrlno)
c---------------------------------------------------------------------   
c 
c 
c 
c     ce     -->  psi (ctrl,ua)
c     ctrl   -->  gamma
c     ua     -->  w
c 
c 
c 
c---------------------------------------------------------------------   
      INCLUDE 'Param3D.h'
      INCLUDE 'Paramopt3D.h'
c---------------------------------------------------------------------

      integer*4 is,ia
      real*8    ctrlno,ctrl(nnsp)


      DO  is=1,ns
        DO  ia=1,5
          dx(ia,is)            = 0.0
          dy(ia,is)            = 0.0
          dz(ia,is)            = 0.0
        END DO
      END DO

      DO is=1,nsmax
         ce(1,is)                  = 0.d0
         ce(2,is)                  = 0.d0
         ce(3,is)                  = 0.d0
         ce(4,is)                  = 0.d0
         ce(5,is)                  = 0.d0
      END DO

c...  Compute the hermitian nodal gradients
      CALL GRADNOD

c...  Compute the convective fluxes
      CALL FLUROE

      
c...  Compute the boundary conditions

      CALL VCURVM(ctrlno)      
      IF (itrans.eq.1 .and. ctrlno.gt.1.0d-10)
     .  CALL TRANSPIRATION(CE,CTRL)        
      CALL CONDDIRFLUX        
      
      END
