program runspe10
  use parameters
  use utils
  use matrix
  use finitevolume
  use simulation
  use netcdf
  use netcdfwrapper
  implicit none  

  ! directory to read all input data from
  character(len = 100) :: data_directory
  ! directory to write results into
  character(len = 100) :: results_directory


  ! dimension for independent, dependent, and the parameter  
  integer :: n_dim, m_dim, p_dim
  ! input variables
  double precision, dimension(:), allocatable :: x
  ! output variables
  double precision, dimension(:), allocatable :: y
  ! parameter variables
  double precision, dimension(:), allocatable :: param

  call get_filepaths(data_directory, results_directory)

  call get_independent_size( n_dim, data_directory )
  call get_dependent_size(   m_dim, data_directory )

  call allocate_independent_variables( n_dim, x       )
  call allocate_dependent_variables(   m_dim, y       )
  call allocate_parameter_variables(   p_dim, param   )


  call initialize_parameter_variables(   p_dim, param )
  call initialize_independent_variables( n_dim, x,   data_directory )

  call evaluate_original_code(n_dim, m_dim, p_dim, x, y, param, &
       data_directory, results_directory)

  call save_dependent_variables( m_dim, y       , data_directory )

  call deallocate_independent_variables( n_dim, x )
  call deallocate_dependent_variables(   m_dim, y )
  call deallocate_parameter_variables(   p_dim, param )

  return
contains

  subroutine get_filepaths(data_directory, results_directory)
    integer :: iargs  
    character(len = 100) :: data_directory
    character(len = 100) :: results_directory

    ! =================================
    ! READ DATA AND RESULTS DIRECTORY
    ! =================================
    iargs = iargc()
    if(iargs /= 2) then
       write(*,*) "Incorrect number of arguments passed. Expected data and results directory."
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

    call ncopen(trim(adjustl(data_directory))//"x.nc", NC_NOWRITE, ncid)
    call ncread(ncid, "n_dim", n_dim)
    call ncclose(ncid)
  end subroutine get_independent_size

  subroutine get_dependent_size( m_dim, data_directory )
    integer :: ncid
    integer, intent(out):: m_dim
    character(len = 100) :: data_directory

    call ncopen(trim(adjustl(data_directory))//"y.nc", NC_NOWRITE, ncid)
    call ncread(ncid, "m_dim", m_dim)
    call ncclose(ncid)
  end subroutine get_dependent_size

  subroutine allocate_independent_variables( n_dim, x)
    integer, intent(in):: n_dim
    double precision, dimension(:), allocatable, intent(out):: x
    ! User-Application specific
    ! ===========================
    ! N/A
    ! Standard AD-Suite Interface
    ! ===========================
    allocate( x(n_dim) )
  end subroutine allocate_independent_variables

  subroutine allocate_dependent_variables( m_dim, y)
    integer, intent(in):: m_dim
    double precision, dimension(:), allocatable, intent(out):: y
    ! User-Application specific
    ! ===========================
    ! N/A
    ! Standard AD-Suite Interface
    ! ===========================
    allocate( y(m_dim) )
  end subroutine allocate_dependent_variables

  subroutine allocate_parameter_variables( p_dim, param)
    integer, intent(out):: p_dim
    double precision, dimension(:), allocatable, intent(out):: param
    ! User-Application specific
    ! ===========================
    p_dim = 1 ! parameter is not used through this interface in this app.
    ! Standard AD-Suite Interface
    ! ===========================
    allocate( param(p_dim) )
  end subroutine allocate_parameter_variables

  subroutine initialize_independent_variables(n_dim, x, data_directory)
    integer :: ncid
    integer, intent(in):: n_dim
    character(len = 100) :: data_directory
    double precision, dimension(:), allocatable, intent(inout):: x
    ! User-Application specific
    ! ===========================
    x = 0.0d0
    ! Standard AD-Suite Interface
    ! =========================== 
    call ncopen(trim(adjustl(data_directory))//"x.nc", NC_NOWRITE, ncid)
    call ncread(ncid, "x", x)
    call ncclose(ncid)
  end subroutine initialize_independent_variables

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

  subroutine evaluate_original_code(n_dim, m_dim, p_dim, x, y, param, &
	 data_directory, results_directory)
    integer, intent(in) :: n_dim, m_dim, p_dim
    double precision, dimension(n_dim), intent(in) :: x
    double precision, dimension(m_dim), intent(inout) :: y
    double precision, dimension(p_dim), intent(in) :: param

    character(len = 100) :: data_directory
    character(len = 100) :: results_directory

    integer :: i, j, k, n_dof
    integer :: nx, ny, nz
    integer :: nd, st, pt

    double precision :: totaloil
    double precision, dimension(:), allocatable :: Q
    double precision, dimension(:), allocatable :: S
    double precision, dimension(:, :, :), allocatable :: P
    double precision, dimension(:, :, :, :), allocatable :: V
    double precision, dimension(:), allocatable :: Tt
    double precision, dimension(:, :), allocatable :: Pc
    double precision, dimension(:), allocatable :: mu
    double precision, dimension(:), allocatable :: sigma

    ! ==========================================
    ! ALLOCATE SPACE FOR INDIVIDUAL INDEPENDENTS
    ! ==========================================
    n_dof = n_dim/2
    allocate( mu(n_dof)    )
    allocate( sigma(n_dof) )    
   
    ! ===========================
    ! READ INDEPENDENT VARIABLES
    ! ===========================
    mu    = x(1  :  n_dof)
    sigma = x(n_dof + 1 :)

    ! =======================================================
    ! INITIALIZE SCENARIO AND APPLICATION SPECIFIC VARIABLES
    ! =======================================================
    call initialize_size_param_vars(nx, ny, nz, st, pt, nd, data_directory)

    call allocate_vector_param_vars(nx, ny, nz)
    call allocate_vector_sim_vars(nx, ny, nz, st, pt, nd, Q, S, P, V, Tt, Pc)

    call initialize_scalar_param_vars(data_directory)
    call initialize_vector_param_vars(nx, ny, nz, data_directory)
    call initialize_vector_sim_vars(nx, ny, nz, Q, S, P, V, Tt, Pc, data_directory)

    ! ======================================================
    ! OVER-RIDE PARAMETER VARIABLES, UPDATE VALUES ETC.
    ! ======================================================
    verbose = .false.   

    ! ===========================
    ! EXECUTE CODE
    ! ===========================
    call wrapper(nx, ny, nz, nd, n_dof, pt, st, mu, sigma, Q, S, P, V, Tt, Pc, totaloil)

    ! ===========================
    ! WRITE RESULTS
    ! ===========================
    call write_results(results_directory, nx, ny, nz, st, pt, nd, n_dof, mu, sigma, Tt, Pc, totaloil)

    ! ===========================
    ! WRITE DEPENDENT VARIABLES
    ! ===========================
    y(1) = totaloil

    ! ==========================
    ! DEALLOCATE VARIABLES
    ! ===========================
    call deallocate_vector_sim_vars(Q, S, P, V, Tt, Pc)
    call deallocate_vector_param_vars()
    deallocate(mu)
    deallocate(sigma)
  end subroutine evaluate_original_code

  subroutine save_dependent_variables(m_dim, y, data_directory)
    integer :: ncid
    integer, intent(in):: m_dim
    character(len = 100) :: data_directory
    double precision, dimension(:), allocatable, intent(in):: y
    ! Standard AD-Suite Interface
    ! ===========================  
    ! N/A
    ! Standard AD-Suite Interface
    ! ===========================  
    call ncopen(trim(adjustl(data_directory))//"y.nc", NC_WRITE, ncid)
    call ncwrite(ncid, "m_dim", m_dim)
    call ncwrite(ncid, "y", y, (/"m_dim"/))
    call ncclose(ncid)
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

  subroutine deallocate_parameter_variables( p_dim, param)
    integer, intent(in):: p_dim
    double precision, dimension(:), allocatable, intent(inout):: param
    ! User-Application specific
    ! ===========================
    ! N/A
    ! Standard AD-Suite Interface
    ! ===========================  
    deallocate( param )
  end subroutine deallocate_parameter_variables

  subroutine initialize_size_param_vars(nx, ny, nz, st, pt, nd, data_directory)
    integer :: ncid
    integer :: nx, ny, nz
    integer :: st, pt, nd
    character(len = *) :: data_directory     

    ! Open parameters file
    call ncopen(trim(adjustl(data_directory))//"/parameters1.nc", &
         NC_NOWRITE, ncid)

    ! read NX, NY, and NZ
    call ncread(ncid, "NX", nx)
    call ncread(ncid, "NY", ny)
    call ncread(ncid, "NZ", nz)

    ! read St, Pt, ND
    call ncread(ncid, "St",  st)
    call ncread(ncid, "Pt",  pt)
    call ncread(ncid, "ND",  nd)

    call ncclose(ncid)
  end subroutine initialize_size_param_vars

  subroutine initialize_scalar_param_vars(data_directory)
    implicit none
    integer :: ncid  
    character(len = *) :: data_directory 

    ! Open parameters file
    call ncopen(trim(adjustl(data_directory))//"/parameters1.nc", &
         NC_NOWRITE, ncid)

    ! initialize all constants insize parameters file.
    ! read ir_const
    call ncread(ncid, "ir_const", ir)

    ! read hX, hY, hZ
    call ncread(ncid, "hX", hx)
    call ncread(ncid, "hY", hy)
    call ncread(ncid, "hZ", hz)

    ! real volume
    call ncread(ncid, "vol", vol)

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
  end subroutine initialize_scalar_param_vars

  subroutine allocate_vector_param_vars(nx, ny, nz)
    integer :: nx, ny, nz

    ! Allocate space for constant permeability/porosity
    allocate(POR(nx * ny * nz), PERM(3, nx, ny, nz))    
  end subroutine allocate_vector_param_vars

  subroutine deallocate_vector_param_vars()
    ! also deallocate Porosity and Permeability
    deallocate(POR, PERM)
  end subroutine deallocate_vector_param_vars

  subroutine initialize_vector_param_vars(nx, ny, nz, data_directory)
    character(len = *) :: data_directory                            
    integer :: nx, ny, nz

    ! Now read the permeabilities and porosities 
    call read_permeability_and_porosity(nx, ny, nz, PERM, POR, data_directory)
  end subroutine initialize_vector_param_vars

  subroutine read_permeability_and_porosity(nx, ny, nz, PERM, POR, data_directory)   
    integer :: ncid                                                  
    integer :: nx, ny, nz
    character(len = *) :: data_directory                             

    double precision, dimension((nx * ny * nz)) :: POR       ! Porosities
    double precision, dimension(nx, ny, nz) :: POR_temp    ! Porosities_temp
    double precision, dimension(3, nx, ny, nz)  :: PERM      ! Permeabilities

    ! initialize porosity and permeability to zero
    PERM = 0.0d0
    POR = 0.0d0
    POR_temp = 0.0d0

    ! read permeability
    call ncopen(trim(adjustl(data_directory))//"/parameters2.nc", &
         NC_NOWRITE, ncid)
    call ncread(ncid, "permeability", PERM)
    call ncclose(ncid)

    ! read porosity
    call ncopen(trim(adjustl(data_directory))//"/parameters3.nc", &
         NC_NOWRITE, ncid)
    call ncread(ncid, "porosity", POR_temp)
    call ncclose(ncid)    

    POR = reshape(POR_temp, (/nx * ny * nz/))

    call mymax_1_0_double(POR, 1.0d-3, POR)
  end subroutine read_permeability_and_porosity

  subroutine allocate_vector_sim_vars(nx, ny, nz, st, pt, nd, Q, S, P, V, Tt, Pc)
    implicit none
    integer :: nx, ny, nz
    integer :: st, pt, nd

    double precision, dimension(:), allocatable          :: Q
    double precision, dimension(:), allocatable          :: S
    double precision, dimension(:, :, :), allocatable    :: P
    double precision, dimension(:, :, :, :), allocatable :: V
    double precision, dimension(:), allocatable          :: Tt
    double precision, dimension(:, :), allocatable       :: Pc

    ! Allocate space for Q, S, P, V
    allocate(Q(nx * ny * nz), &
         S(nx * ny * nz), &
         P(nx, ny, nz), &
         V(3, nx+1, ny+1, nz+1))

    ! Allocate space for Time and Production data 
    allocate(Tt((nd/st) + 1), & 
         Pc(4, (nd/st) + 1))
  end subroutine allocate_vector_sim_vars

  subroutine deallocate_vector_sim_vars(Q, S, P, V, Tt, Pc)
    implicit none

    double precision, dimension(:), allocatable          :: Q
    double precision, dimension(:), allocatable          :: S
    double precision, dimension(:, :, :), allocatable    :: P
    double precision, dimension(:, :, :, :), allocatable :: V
    double precision, dimension(:), allocatable          :: Tt
    double precision, dimension(:, :), allocatable       :: Pc

    ! deallocate variables
    deallocate(Q, S, P, V)
    deallocate(Tt, Pc)
  end subroutine deallocate_vector_sim_vars

  subroutine initialize_vector_sim_vars(nx, ny, nz, Q, S, P, V, Tt, Pc, data_directory)
    integer :: ncid
    integer :: nx, ny, nz
    double precision, dimension(:)          :: Q
    double precision, dimension(:)          :: S
    double precision, dimension(:, :, :)    :: P
    double precision, dimension(:, :, :, :) :: V
    double precision, dimension(:)          :: Tt
    double precision, dimension(:, :)       :: Pc
    double precision, dimension(:,:,:), allocatable :: S_temp

    character(len = *) :: data_directory                             

    ! initialize memory for inflow and saturation
    Q = 0.0d0

    ! read initial saturation from file
    ! Open parameters file
    allocate(S_temp(nx, ny, nz))
    call ncopen(trim(adjustl(data_directory))//"/parameters4.nc", &
         NC_NOWRITE, ncid)
    call ncread(ncid, "saturation", S_temp) 
    call ncclose(ncid)
    S = reshape(S_temp, (/ nx * ny * nz /))
    deallocate(S_temp)

    ! initialize memory for P and V
    P = 0.0d0

    ! note that V vector has an additional
    ! length in each dimension x,y,z
    V = 0.0d0

    Tt = 0.0d0               ! simulation time
    Pc = 0.0d0               ! production data
  end subroutine initialize_vector_sim_vars

  subroutine write_results(results_directory, nx, ny, nz, st, pt, nd, n_dof, mu, sigma, Tt, Pc, totaloil)
    implicit none

    integer :: ncid                                                  
    integer :: nx, ny, nz, st, pt, nd, n_dof
    double precision :: totaloil
    double precision, dimension(:)             :: mu
    double precision, dimension(:)             :: sigma
    double precision, dimension(:)             :: Tt
    double precision, dimension(:, :)          :: Pc
    character(len = *) :: results_directory

    call ncopen(trim(adjustl(results_directory))//"/results_eval_original_code.nc", &
         NC_WRITE, ncid)

    call ncwrite(ncid, "NX", nx)
    call ncwrite(ncid, "NY", ny)
    call ncwrite(ncid, "NZ", nz)
    call ncwrite(ncid, "MU", mu)
    call ncwrite(ncid, "SIGMA", sigma)
    call ncwrite(ncid, "TIME", Tt)
    call ncwrite(ncid, "PROD_CURVES", Pc)
    call ncwrite(ncid, "OIL", totaloil)

    ! Close the file.
    call ncclose(ncid)
  end subroutine write_results
end program runspe10
