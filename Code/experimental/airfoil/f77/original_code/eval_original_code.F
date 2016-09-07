ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     c
c     Nonlinear airfoil lift calculation                           c
c     c
c     Copyright Devendra Ghate and Mike Giles, 2005                c
c     but can be freely used with due acknowledgement              c
c     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      program airfoil
      implicit none

c     directory to read all input data from
      character(len = 100) :: data_directory
c     directory to write results into
      character(len = 100) :: results_directory

c     max number of iterations, minimum iterations
c     between successive lift computation
      integer :: maxiter, minitergap
      parameter(maxiter = 50000, minitergap = 100)

c     dimension for independent, dependent, and the parameter  
      integer :: n_dim, m_dim, p_dim
c     input variables
      double precision, dimension(1) :: x_ind
c     output variables
      double precision, dimension((maxiter/minitergap)) :: y_dep
c     parameter variables
      double precision, dimension(11) :: param

      call get_filepaths(data_directory, results_directory)

c     Not applicable for the airfoil test case as F77 does not
c     support dynamic memory allocation
      call get_independent_size( n_dim, data_directory )
      call get_dependent_size(   m_dim, data_directory )
      call get_parameter_size(   p_dim, data_directory )

c     Not applicable for the airfoil test case as F77 does not
c     support dynamic memory allocation
c     call allocate_independent_variables( n_dim, x       )
c     call allocate_dependent_variables(   m_dim, y       )
c     call allocate_parameter_variables(   p_dim, param   )

      call initialize_parameter_variables(   p_dim, param,
     $     data_directory )
      call initialize_independent_variables( n_dim, x_ind,
     $     data_directory )

      call evaluate_original_code(n_dim, m_dim, p_dim, x_ind, y_dep,
     $     param, data_directory, results_directory)

      call save_dependent_variables( m_dim, y_dep , data_directory )

c     Not applicable for the airfoil test case as F77 does not
c     support dynamic memory allocation
c     call deallocate_independent_variables( n_dim, x )
c     call deallocate_dependent_variables(   m_dim, y )
c     call deallocate_parameter_variables(   p_dim, param )

      contains

      subroutine get_filepaths(data_directory, results_directory)
      integer :: iargs  
      character(len = 100) :: data_directory
      character(len = 100) :: results_directory

c     =================================
c     READ DATA AND RESULTS DIRECTORY
c     =================================
      iargs = iargc()
      if(iargs /= 2) then
         write(*,*) "Incorrect number of arguments passed."
         write(*,*) "Expected data and results directory."
      endif
      call getarg(1, data_directory)
      call getarg(2, results_directory)
      write(*,*) "Data directory: ", data_directory
      write(*,*) "Results directory: ", results_directory  
      end subroutine get_filepaths


      subroutine get_independent_size( n_dim, data_directory )
      integer :: ncid
      integer, intent(out):: n_dim
      character(len = 100) :: data_directory

      n_dim = 1
      end subroutine get_independent_size


      subroutine get_dependent_size( m_dim, data_directory )
      integer :: ncid
      integer, intent(out):: m_dim
      character(len = 100) :: data_directory

      m_dim = maxiter/minitergap
      end subroutine get_dependent_size


      subroutine get_parameter_size( p_dim, data_directory )
      integer :: ncid
      integer, intent(out):: p_dim
      character(len = 100) :: data_directory

      p_dim = 11
      end subroutine get_parameter_size


c     subroutine allocate_independent_variables( n_dim, x)
c     integer, intent(in):: n_dim
c     double precision, dimension(:), allocatable, intent(out):: x
c     User-Application specific
c     ===========================
c     N/A
c     Standard AD-Suite Interface
c     ===========================
c     allocate( x(n_dim) )
c     end subroutine allocate_independent_variables
c     
c     
c     subroutine allocate_dependent_variables( m_dim, y)
c     integer, intent(in):: m_dim
c     double precision, dimension(:), allocatable, intent(out):: y
c     User-Application specific
c     ===========================
c     N/A
c     Standard AD-Suite Interface
c     ===========================
c     allocate( y(m_dim) )
c     end subroutine allocate_dependent_variables
c     
c     
c     subroutine allocate_parameter_variables( p_dim, param)
c     integer, intent(out):: p_dim
c     double precision, dimension(:), allocatable, intent(out):: param
c     User-Application specific
c     ===========================
c     p_dim = 9                
c     Standard AD-Suite Interface
c     ===========================
c     allocate( param(p_dim) )
c     end subroutine allocate_parameter_variables


      subroutine initialize_independent_variables(n_dim, x,
     $     data_directory)
      integer :: i
      integer, intent(in):: n_dim
      character(len = 100) :: data_directory
      double precision, dimension(n_dim), intent(inout):: x
c     User-Application specific
c     ===========================
c     x(1) = 3.0 * datan(1.0d0)/45.0d0       ! alpha
      open(1, file=trim(adjustl(data_directory))//'x.dat',
     $     status='old')

      do i = 1, n_dim
         read(1,*,err=999,end=999) x(i)
      enddo
 999  close(1)      
c     Standard AD-Suite Interface
c     =========================== 
c     N/A
      end subroutine initialize_independent_variables


      subroutine initialize_parameter_variables(p_dim, param,
     $     data_directory)
      integer :: i
      integer, intent(in):: p_dim
      character(len = 100) :: data_directory
      double precision, dimension(p_dim), intent(inout):: param
c     User-Application specific
c     ===========================
c     param(1) = 1.4d0               ! gam
c     param(2) = 0.4d0               ! gm1
c     param(3) = 0.9d0               ! cfl
c     param(4) = 0.05d0              ! eps
c     param(5) = 0.4d0               ! mach
c     param(6) = 1.0d0               ! p
c     param(7) = 1.0d0               ! r
c     param(8) = 0.4732863826479693  ! u
c     param(9) = 2.612               ! e
c     param(10) = 20000              ! iter
c     param(11) = 100                ! lift iter
      open(1, file=trim(adjustl(data_directory))//'parameters1.dat',
     $     status='old')

      do i = 1, p_dim
         read(1,*,err=999,end=999) param(i)
      enddo
 999  close(1)      
c     Standard AD-Suite Interface
c     ===========================
c     N/A
      end subroutine initialize_parameter_variables


      subroutine evaluate_original_code(n_dim, m_dim, p_dim, x_ind,
     $     y_dep, param, data_directory, results_directory)

c     ===========================================================
c     FORMAL ARGUMENTS
c     ===========================================================

      character(len = 100) :: data_directory
      character(len = 100) :: results_directory

c     dimension for independent, dependent, and the parameter  
      integer ::  n_dim, m_dim, p_dim
c     input variables
      double precision, dimension(n_dim) :: x_ind
c     output variables
      double precision, dimension(m_dim) :: y_dep
c     parameter variables
      double precision, dimension(p_dim) :: param

c     ===========================================================
c     SIMULATION VARIABLES
c     ===========================================================

      integer  ::  maxnode,     maxcell,       maxedge
      parameter (maxnode=721801,maxcell=720000+1,maxedge=1441400)
c     
      integer :: ecell(2,maxedge), boun(maxedge),
     &     edge(2,maxedge), cell(4,maxcell)

      double precision ::  x(2,maxnode), q(4,maxcell), qold(4,maxcell),
     &     adt(maxcell), res(4,maxcell)

c     
      integer :: in, in1, in2, in3, in4, ic, ic1, ic2, ie, ipde,
     &     iter, niter, nitergap, k, nnode, ncell, nedge

      double precision ::  lift, rms, xt, yt, p, r, u, e

c     ===============================================
c     ALLOCATE SPACE FOR INDIVIDUAL INDEPENDENTS
c     ===============================================

c     ===============================================
c     INCLUDE VARIABLES IN COMMON BLOCK
c     ===============================================
#include "./src/const.inc"
c     

c     ===============================================
c     READ INDEPENDENT VARIABLES
c     ===============================================

      alpha = x_ind(1)

c     =====================================================
c     INTIALIZE SCENARIO AND APPLICATION SPECIFIC VARIABLES
c     =====================================================

c     ===============================================
c     set universal constants
c     ===============================================
      gam     =  param(1)       ! gam
      gm1     =  param(2)       ! gm1
      cfl     =  param(3)       ! cfl
      eps     =  param(4)       ! eps
      mach    =  param(5)       ! mach
      p       =  param(6)       ! p
      r       =  param(7)       ! r
      u       =  param(8)       ! u
      e       =  param(9)       ! e
      niter    = int(param(10)) ! iter
      nitergap = int(param(11)) ! lift iter gap

c     =========================================================
c     read in grid and flow data
c     =========================================================
      call read_grid(maxnode, maxcell, maxedge, nnode, ncell, nedge,
     $     x, cell, edge, ecell, boun, data_directory)


      call read_flow(maxcell, ncell, p, r, u, e, q, data_directory)

c     =========================================================
c     OVER-RIDE PARAMETER VARIABLES, UPDATE VALUES ETC.
c     =========================================================

c     =====================================================
c     rotate grid to specified angle of attack
c     =====================================================
      do in = 1, nnode
         xt = x(1,in)
         yt = x(2,in)
         x(1,in) = cos(alpha)*xt + sin(alpha)*yt
         x(2,in) = cos(alpha)*yt - sin(alpha)*xt
      enddo
c     


c     ===========================
c     EXECUTE CODE
c     ===========================


c     =========================================================
c     main time-marching loop
c     =========================================================
      
      do iter = 1, niter
c     =========================================================
c     save old flow solution
c     =========================================================
         do ic = 1, ncell-1
            do ipde = 1,4
               qold(ipde,ic) = q(ipde,ic)
            enddo
         enddo
c     =========================================================
c     predictor/corrector update loop
c     =========================================================
         do k = 1, 2
            do ic = 1, ncell
               do ipde = 1, 4
                  res(ipde,ic) = 0.0d0
               enddo
            enddo
c     ========================================================
c     calculate area/timstep
c     ========================================================
            do ic = 1, ncell-1
               in1 = cell(1,ic)
               in2 = cell(2,ic)
               in3 = cell(3,ic)
               in4 = cell(4,ic)
               call time_cell(x(1,in1),x(1,in2),x(1,in3),x(1,in4),
     &              q(1,ic),adt(ic))
            enddo

            adt(ncell) = 0.0d0
c     =======================================================
c     flux evaluation loop 
c     =======================================================
            do ie = 1, nedge
               in1 = edge(1,ie)
               in2 = edge(2,ie)
               ic1 = ecell(1,ie)
               ic2 = ecell(2,ie)     
               if(boun(ie).eq.0) then
                  call flux_face(x(1,in1),x(1,in2),q(1,ic1),q(1,ic2), 
     &                 adt(ic1),adt(ic2),res(1,ic1),res(1,ic2))
               elseif (boun(ie).eq.1) then
                  call flux_wall(x(1,in1),x(1,in2),q(1,ic2),res(1,ic2))
               elseif (boun(ie).eq.2) then
                  stop 'Problem'
               endif
            enddo
c     =====================================================
c     flow field update
c     =====================================================
            rms = 0.0d0
            
            do ic = 1, ncell-1
               do ipde = 1, 4
                  q(ipde,ic) = qold(ipde,ic) - res(ipde,ic)/adt(ic)
                  rms = rms + (q(ipde,ic)-qold(ipde,ic))**2
               enddo
            enddo
         enddo

         rms = dsqrt(rms/ncell)
c     ===================================================
c     print iteration history, including lift calculation
c     ===================================================
         if(mod(iter, nitergap).eq.0) then
            lift = 0.0d0    
            do ie = 1,nedge
               if(boun(ie).eq.1) then
                  in1 = edge(1,ie)
                  in2 = edge(2,ie)
                  ic2 = ecell(2,ie)
                  call lift_wall(x(1,in1),x(1,in2),q(1,ic2),lift)
               endif
            enddo
            write(*,'(i5,e13.6,2e23.16)') iter, rms, lift

            y_dep(iter/nitergap) = lift
         endif
      enddo

c     ===========================
c     WRITE RESULTS
c     ===========================
      call write_flow(ncell, q, results_directory)

c     ===========================
c     WRITE DEPENDENT VARIABLES
c     ===========================
c     N/A
c     ==========================
c     DEALLOCATE VARIABLES
c     ===========================
c     N/A
      end subroutine evaluate_original_code


      subroutine save_dependent_variables(m_dim, y, results_directory)
      integer :: i
      integer, intent(in):: m_dim
      character(len = 100) :: results_directory
      double precision, dimension(m_dim), intent(in):: y
c     Standard AD-Suite Interface
c     ===========================  
c     N/A
c     Standard AD-Suite Interface
c     ===========================  
      open(1, file=trim(adjustl(results_directory))//'y.dat',
     $     status='old')
      do i = 1, m_dim
         write(1,*,err=999) y(i)
      enddo
 999  close(1)
      end subroutine save_dependent_variables


c     subroutine deallocate_independent_variables( n_dim, x)
c     integer, intent(in):: n_dim
c     double precision, dimension(:), allocatable, intent(inout):: x
c     User-Application specific
c     ===========================
c     N/A
c     Standard AD-Suite Interface
c     ===========================  
c     deallocate( x )
c     end subroutine deallocate_independent_variables
c     
c     
c     subroutine deallocate_dependent_variables( m_dim, y)
c     integer, intent(in):: m_dim
c     double precision, dimension(:), allocatable, intent(inout):: y
c     User-Application specific
c     ===========================
c     N/A
c     Standard AD-Suite Interface
c     ===========================  
c     deallocate( y )
c     end subroutine deallocate_dependent_variables
c     
c     
c     subroutine deallocate_parameter_variables( p_dim, param)
c     integer, intent(in):: p_dim
c     double precision, dimension(:), allocatable, intent(inout):: param
c     User-Application specific
c     ===========================
c     N/A
c     Standard AD-Suite Interface
c     ===========================  
c     deallocate( param )
c     end subroutine deallocate_parameter_variables

      end program
