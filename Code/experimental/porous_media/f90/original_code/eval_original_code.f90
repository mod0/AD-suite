program runspe10
  use parameters
  use utils
  use matrix
  use finitevolume
  use simulation
  use netcdf
  use netcdfwrapper
  implicit none  
  ! dimension for independent, dependent, and the parameter  
  integer :: n_dim, m_dim, p_dim
  ! input variables
  double precision, dimension(:), allocatable :: x
  ! output variables
  double precision, dimension(:), allocatable :: y
  ! parameter variables
  double precision, dimension(:), allocatable :: param

  call allocate_independent_variables( n_dim, x     )
  call allocate_dependent_variables(   m_dim, y     )
!!!!!!!!!!!!!!! Do we now remove the following statement?
  call allocate_parameter_variables(   p_dim, param )
!!!!!!!!!!!!!!! Do we now remove the following statement?
  call initialize_parameter_variables(   p_dim, param )
  call initialize_independent_variables( n_dim, x     )

  call evaluate_original_code(n_dim, m_dim, p_dim, x, y, param)

  call save_dependent_variables( m_dim, y )

  call deallocate_independent_variables( n_dim, x )
  call deallocate_dependent_variables(   m_dim, y )
!!!!!!!!!!!!!!! Do we now remove the following statement?
  call deallocate_parameter_variables(    p_dim, param )

  return
contains

  subroutine allocate_independent_variables( n_dim, x)
    integer, intent(out):: n_dim
    double precision, dimension(:), allocatable, intent(out):: x
    ! User-Application specific
    ! ===========================
    n_dim = 6
    ! Standard AD-Suite Interface
    ! ===========================
    allocate( x(n_dim) )
  end subroutine allocate_independent_variables

  subroutine allocate_dependent_variables( m_dim, y)
    integer, intent(out):: m_dim
    double precision, dimension(:), allocatable, intent(out):: y
    ! User-Application specific
    ! ===========================
    m_dim = 1
    ! Standard AD-Suite Interface
    ! ===========================
    allocate( y(m_dim) )
  end subroutine allocate_dependent_variables

!!!!!!!!!!!!!!! Do we now remove this?
  subroutine allocate_parameter_variables( p_dim, param)
    integer, intent(out):: p_dim
    double precision, dimension(:), allocatable, intent(out):: param
    ! User-Application specific
    ! ===========================
    ! N/A
    ! Standard AD-Suite Interface
    ! ===========================
    allocate( param(p_dim) )
  end subroutine allocate_parameter_variables

  subroutine initialize_independent_variables(n_dim, x)
    integer, intent(in):: n_dim
    double precision, dimension(:), allocatable, intent(inout):: x
    ! User-Application specific
    ! ===========================
    x(1) = -1.0d0
    x(2) =  0.0d0
    x(3) = +1.0d0
    x(4) = 0.0001d0
    x(5) = 0.0001d0
    x(6) = 0.0001d0
    ! Standard AD-Suite Interface
    ! =========================== 
    ! N/A
  end subroutine initialize_independent_variables

!!!!!!!!!!!!!!! Do we now remove this?
  subroutine initialize_parameter_variables(p_dim, param)
    integer, intent(in):: p_dim
    double precision, dimension(:), allocatable, intent(inout):: param
    ! User-Application specific
    ! ===========================
    ! N/A
    ! Standard AD-Suite Interface
    ! ===========================
    ! N/A
  end subroutine initialize_parameter_variables

  subroutine evaluate_original_code(n_dim, m_dim, p_dim, x, y, param)
    integer, intent(in) :: n_dim, m_dim, p_dim
    double precision, dimension(n_dim), intent(in) :: x
    double precision, dimension(m_dim), intent(inout) :: y
    double precision, dimension(p_dim), intent(in) :: param
    character(len = 100) :: data_directory
    character(len = 100) :: results_directory
    integer :: iargs  
    integer :: i, j, k, n_dof
    integer :: nx, ny, nz
    integer :: nd, st, pt
    double precision, dimension(:), allocatable :: Q
    double precision, dimension(:), allocatable :: S
    double precision, dimension(:, :, :), allocatable :: P
    double precision, dimension(:, :, :, :), allocatable :: V
    double precision :: totaloil
    double precision, dimension(:), allocatable :: Tt
    double precision, dimension(:, :), allocatable :: Pc
    double precision, dimension(:), allocatable :: mu
    double precision, dimension(:), allocatable :: sigma
    n_dof = n_dim/2
    allocate( mu(n_dof)    )
    allocate( sigma(n_dof) )    
    ! READ INDEPENDENT VARIABLES
    ! ===========================
    mu(1)    = x(1)
    mu(2)    = x(2)
    mu(3)    = x(3)
    sigma(1) = x(4)
    sigma(2) = x(5)
    sigma(3) = x(6)
    ! =================================
    ! READ DATA AND RESULTS DIRECTORY
    ! =================================
    iargs = iargc()
    if(iargs /= 2) then
       write(*,*) "Incorrect number of arguments passed. Expected data and results directory"
    endif
    call getarg(1, data_directory)
    call getarg(2, results_directory)
    write(*,*) "Data directory: ", data_directory
    write(*,*) "Results directory: ", results_directory  
    ! =======================================================
    ! INITIALIZE SCENARIO AND APPLICATION SPECIFIC VARIABLES
    ! =======================================================
    call initialize_scenario(data_directory, nx, ny, nz, st, pt, nd, solver_inner, solver_outer)
    call allocate_arrays(nx, ny, nz, st, pt, nd, Q, S, P, V, Tt, Pc)
    call allocate_shared_arrays(nx, ny, nz, st, pt, nd)
    call initialize_arrays(Q, S, P, V, Tt, Pc)
    call initialize_shared_arrays(data_directory, nx, ny, nz, st, pt, nd)
    verbose = .false.   
    ! EXECUTE CODE
    ! ===========================
    call wrapper(nx, ny, nz, nd, n_dof, pt, st, mu, sigma, Q, S, P, V, Tt, Pc, totaloil)
    ! ===========================
    call write_results(results_directory, nx, ny, nz, st, pt, nd, n_dof, mu, sigma, Tt, Pc, totaloil)
    ! WRITE DEPENDENT VARIABLES
    ! ===========================
    y(1) = totaloil
    ! ===========================
    call deallocate_arrays(Q, S, P, V, Tt, Pc)
    call deallocate_shared_arrays()
    deallocate(mu)
    deallocate(sigma)
  end subroutine evaluate_original_code

  subroutine save_dependent_variables(m_dim, y)
    integer, intent(in):: m_dim
    double precision, dimension(:), allocatable, intent(in):: y
    ! Standard AD-Suite Interface
    ! ===========================  
    ! N/A
    ! Standard AD-Suite Interface
    ! ===========================  
  end subroutine save_dependent_variables

  subroutine deallocate_independent_variables( n_dim, x)
    integer, intent(in):: n_dim
    double precision, dimension(:), allocatable, intent(inout):: x
    ! User-Application specific
    ! ===========================
    ! N/A
    ! Standard AD-Suite Interface
    ! ===========================  
    deallocate( x )
  end subroutine deallocate_independent_variables

  subroutine deallocate_dependent_variables( m_dim, y)
    integer, intent(in):: m_dim
    double precision, dimension(:), allocatable, intent(inout):: y
    ! User-Application specific
    ! ===========================
    ! N/A
    ! Standard AD-Suite Interface
    ! ===========================  
    deallocate( y )
  end subroutine deallocate_dependent_variables

!!!!!!!!!!!!!!! Do we now remove this?
  subroutine deallocate_parameter_variables( p_dim, param)
    integer, intent(in):: p_dim
    double precision, dimension(:), allocatable, intent(inout):: param
    ! User-Application specific
    ! ===========================
    ! N/A
    ! Standard AD-Suite Interface
    ! ===========================  
    ! N/A
  end subroutine deallocate_parameter_variables

  subroutine initialize_parameters(data_directory, nx, ny, nz, st, pt, nd, solver_inner, solver_outer)
    implicit none

    integer :: ncid                                                  ! netcdf file handle

    character(len = *) :: data_directory                             ! directory location of parameters 
    integer :: varid                                                 ! generic variable id
    integer :: nc_chunksize                                          ! chunk size for reading

    integer :: nx, ny, nz
    integer :: st, pt, nd
    integer :: solver_inner, solver_outer  

    ! Open parameters file
    call ncopen(trim(adjustl(data_directory))//"/param1.nc", &
         NF90_NOWRITE, ncid)

    ! read NX, NY, and NZ
    call ncread(ncid, "NX", nx)
    call ncread(ncid, "NY", ny)
    call ncread(ncid, "NZ", nz)

    ! read St, Pt, ND
    call ncread(ncid, "St",  st)
    call ncread(ncid, "Pt",  pt)
    call ncread(ncid, "ND",  nd)

    ! initialize all constants insize parameters file.
    ! read ir_const
    call ncread(ncid, "ir_const", ir)

    ! compute ir
    ir = ir * nx * ny * nz

    ! read hX, hY, hZ
    call ncread(ncid, "hX", hx)
    call ncread(ncid, "hY", hy)
    call ncread(ncid, "hZ", hz)

    ! compute V
    vol = hx * hy * hz

    ! read vw, vo, swc, sor
    call ncread(ncid, "vw", vw)
    call ncread(ncid, "vo", vo)
    call ncread(ncid, "swc", swc)
    call ncread(ncid, "sor", sor)

    ! read solver_parameters
    call ncread(ncid, "solver_inner", solver_inner)    
    call ncread(ncid, "solver_outer", solver_outer)    

    ! Close parameters file
    call ncclose(ncid)
  end subroutine initialize_scenario
end program runspe10
