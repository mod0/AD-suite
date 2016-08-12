program runspe10
  use parameters
  use utils
  use matrix
  use finitevolume
  use simulation
  use netcdf
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
  call allocate_parameter_variables(   p_dim, param )

  call initialize_parameter_variables(   p_dim, param )
  call initialize_independent_variables( n_dim, x     )

  call evaluate_original_code(n_dim, m_dim, p_dim, x, y, param)

  call save_dependent_variables( m_dim, y )

  call deallocate_independent_variables( n_dim, x )
  call deallocate_dependent_variables(   m_dim, y )
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
  ! ===========================
  iargs = iargc()
  if(iargs /= 2) then
     write(*,*) "Incorrect number of arguments passed. Expected data and results directory"
  endif
  call getarg(1, data_directory)
  call getarg(2, results_directory)
  write(*,*) "Data directory: ", data_directory
  write(*,*) "Results directory: ", results_directory  
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


  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ! The following subroutines should all be in the src folder!!!!!!!!!
  subroutine allocate_shared_arrays(nx, ny, nz, st, pt, nd)
    integer :: nx, ny, nz
    integer :: st, pt, nd

    ! Allocate space for constant permeability/porosity
    allocate(POR(nx * ny * nz), PERM(3, nx, ny, nz))    
  end subroutine allocate_shared_arrays

  subroutine deallocate_shared_arrays()
    ! also deallocate Porosity and Permeability
    deallocate(POR, PERM)
  end subroutine deallocate_shared_arrays

  subroutine initialize_shared_arrays(data_directory, nx, ny, nz, st, pt, nd)
    character(len = *) :: data_directory                             ! directory location of parameters 

    integer :: nx, ny, nz
    integer :: st, pt, nd

    ! Now read the permeabilities and porosities which are global
    call read_permeability_and_porosity(data_directory, nx, ny, nz, PERM, POR)
  end subroutine initialize_shared_arrays

  subroutine allocate_arrays(nx, ny, nz, st, pt, nd, Q, S, P, V, Tt, Pc)
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
         Pc(2, (nd/st) + 1))
  end subroutine allocate_arrays

  subroutine deallocate_arrays(Q, S, P, V, Tt, Pc)
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
  end subroutine deallocate_arrays

  subroutine initialize_arrays(Q, S, P, V, Tt, Pc)
    double precision, dimension(:)          :: Q
    double precision, dimension(:)          :: S
    double precision, dimension(:, :, :)    :: P
    double precision, dimension(:, :, :, :) :: V
    double precision, dimension(:)          :: Tt
    double precision, dimension(:, :)       :: Pc

    ! initialize memory for inflow and saturation
    Q = 0.0d0
    S = 0.0d0

    ! initialize memory for P and V
    P = 0.0d0

    ! note that V vector has an additional
    ! length in each dimension x,y,z
    V = 0.0d0

    Tt = 0.0d0               ! simulation time
    Pc = 0.0d0               ! production data
  end subroutine initialize_arrays

  subroutine initialize_scenario(data_directory, nx, ny, nz, st, pt, nd, solver_inner, solver_outer)
    implicit none
    
    integer :: ncid                                                  ! netcdf file handle

    character(len = *) :: data_directory                             ! directory location of parameters 
    integer :: varid                                                 ! generic variable id
    integer :: nc_chunksize                                          ! chunk size for reading

    integer :: nx, ny, nz
    integer :: st, pt, nd
    integer :: solver_inner, solver_outer  

    ! Open parameters file
    nc_chunksize = 4096
    call iserror(nf90_open(trim(adjustl(data_directory))//"/parameters.nc", &
         NF90_NOWRITE, ncid, nc_chunksize))
    write(*,*) "Chosen chunk size is ", nc_chunksize

    ! read scenario_id
    call iserror(nf90_inq_varid(ncid, "scenario_id", varid))
    call iserror(nf90_get_var(ncid, varid, scenario_id))
    write (*,*) "Scenario: ", scenario_id

    ! read NX, NY, and NZ
    call iserror(nf90_inq_varid(ncid, "NX", varid))
    call iserror(nf90_get_var(ncid, varid, nx))
    call iserror(nf90_inq_varid(ncid, "NY", varid))
    call iserror(nf90_get_var(ncid, varid, ny))
    call iserror(nf90_inq_varid(ncid, "NZ", varid))
    call iserror(nf90_get_var(ncid, varid, nz))
    write (*,*) "(NX, NY, NZ) = ", NX, NY, NZ

    ! read St, Pt, ND
    call iserror(nf90_inq_varid(ncid, "St", varid))
    call iserror(nf90_get_var(ncid, varid, st))
    call iserror(nf90_inq_varid(ncid, "Pt", varid))
    call iserror(nf90_get_var(ncid, varid, pt))
    call iserror(nf90_inq_varid(ncid, "ND", varid))
    call iserror(nf90_get_var(ncid, varid, nd))
    write (*,*) "(St, Pt, ND) = ", st, pt, nd

    ! initialize all constants insize parameters file.
    ! read ir_const
    call iserror(nf90_inq_varid(ncid, "ir_const", varid))
    call iserror(nf90_get_var(ncid, varid, ir))

    ! compute ir
    ir = ir * nx * ny * nz
    write (*,*) "(ir) = ", ir

    ! read hX, hY, hZ
    call iserror(nf90_inq_varid(ncid, "hX", varid))
    call iserror(nf90_get_var(ncid, varid, hx_))
    call iserror(nf90_inq_varid(ncid, "hY", varid))
    call iserror(nf90_get_var(ncid, varid, hy_))
    call iserror(nf90_inq_varid(ncid, "hZ", varid))
    call iserror(nf90_get_var(ncid, varid, hz_))
    write (*,*) "(hX, hY, hZ) = ", hx_, hy_, hz_

    ! compute V
    V_ = hx_ * hy_ * hz_
    write (*,*) "(V) = ", V_

    ! read vw, vo, swc, sor
    call iserror(nf90_inq_varid(ncid, "vw", varid))
    call iserror(nf90_get_var(ncid, varid, vw_))
    call iserror(nf90_inq_varid(ncid, "vo", varid))
    call iserror(nf90_get_var(ncid, varid, vo_))
    call iserror(nf90_inq_varid(ncid, "swc", varid))
    call iserror(nf90_get_var(ncid, varid, swc_))
    call iserror(nf90_inq_varid(ncid, "sor", varid))
    call iserror(nf90_get_var(ncid, varid, sor_))
    write (*,*) "(vw, vo, swc, sor) = ", vw_, vo_, swc_, sor_ 

    ! read solver_parameters
    call iserror(nf90_inq_varid(ncid, "solver_inner", varid))
    call iserror(nf90_get_var(ncid, varid, solver_inner))    
    call iserror(nf90_inq_varid(ncid, "solver_outer", varid))
    call iserror(nf90_get_var(ncid, varid, solver_outer))    

    ! Close parameters file
    call iserror(nf90_close(ncid))
  end subroutine initialize_scenario

  subroutine write_results(results_directory, nx, ny, nz, st, pt, nd, n_dof, mu, sigma, Tt, Pc, totaloil)
    implicit none
    character(len = *) :: results_directory

    integer :: nx, ny, nz, st, pt, nd, n_dof
    double precision :: totaloil, totaloil_mud, totaloil_sigmad
    double precision, dimension(:)             :: mu
    double precision, dimension(:)             :: sigma
    double precision, dimension(:)             :: Tt
    double precision, dimension(:, :)          :: Pc

    ! netCDF variables
    integer :: ncid                                                     ! file handle
    integer :: var_scene_id                                             ! scene variable id
    character(len = *), parameter :: var_scene_name = "Scenario"        ! scene variable name
    integer :: var_nx_id                                                ! nx variable id
    character(len = *), parameter :: var_nx_name = "NX"                 ! nx variable name
    integer :: var_ny_id                                                ! ny variable id
    character(len = *), parameter :: var_ny_name = "NY"                 ! ny variable name
    integer :: var_nz_id                                                ! nz variable id
    character(len = *), parameter :: var_nz_name = "NZ"                 ! nz variable name
!     integer :: var_mu_id                                                ! mu variable id
!     character(len = *), parameter :: var_mu_name = "Mu"                 ! mu variable name
!     integer :: var_sigma_id                                             ! sigma variable id
!     character(len = *), parameter :: var_sigma_name = "Sigma"           ! sigma variable name
                                                                        
    integer :: dim_time_id                                              ! time dimension id
    integer :: dim_time_len                                             ! time dimension length
    character(len = *), parameter :: dim_time_name = "Time"             ! time dimension name
                                                                        
    integer :: var_time_id                                              ! time variable id
    character(len = *), parameter :: var_time_name = "Time"             ! time variable name
                                                                        
    integer :: dim_mobility_id                                          ! mobility dimension id
    integer :: dim_mobility_len                                         ! mobility dimension length
    character(len = *), parameter :: dim_mobility_name = "Mobilities"   ! mobility dimension name
                                                                        
    integer :: var_mobility_id                                          ! mobility variable id
    character(len = *), parameter :: var_mobility_name = "WaterOil"     ! mobility variable name
    integer, dimension(2) :: var_mobility_dimids                        ! id of dimensions for mobilities

    integer :: var_oil_id                                               ! oil variable id
    character(len = *), parameter :: var_oil_name = "Oil"               ! oil variable name

    ! Start writing netCDF file having all the computed values
    ! Open file
    call iserror(nf90_create(trim(adjustl(results_directory))//"results_eval_original_code.nc", &
         nf90_clobber, ncid))

    ! setup the sizes of the dimensions
    dim_time_len = (nd/st) + 1
    dim_mobility_len = 2

    ! Define all dimension
    ! Define time dimensions
    call iserror(nf90_def_dim(ncid, dim_time_name, dim_time_len, dim_time_id))
    ! Define mobility dimensions
    call iserror(nf90_def_dim(ncid, dim_mobility_name, dim_mobility_len, dim_mobility_id))

    ! Define scalar scenario input
    call iserror(nf90_def_var(ncid, var_scene_name, NF90_INT, var_scene_id))
    ! Define scalar nx input
    call iserror(nf90_def_var(ncid, var_nx_name, NF90_INT, var_nx_id))
    ! Define scalar ny input
    call iserror(nf90_def_var(ncid, var_ny_name, NF90_INT, var_ny_id))
    ! Define scalar nz input
    call iserror(nf90_def_var(ncid, var_nz_name, NF90_INT, var_nz_id))
    ! Define scalar mu input
!     call iserror(nf90_def_var(ncid, var_mu_name, NF90_DOUBLE, var_mu_id))
    ! Define scalar sigma input
!     call iserror(nf90_def_var(ncid, var_sigma_name, NF90_DOUBLE, var_sigma_id))

    ! Define the time variables. Varid is returned.
    call iserror(nf90_def_var(ncid, var_time_name, NF90_DOUBLE, dim_time_id, var_time_id))

    ! Define the netCDF variables. The dimids array is used to pass the
    ! dimids of the dimensions of the netCDF variables.
    var_mobility_dimids = (/ dim_mobility_id, dim_time_id /)
    call iserror(nf90_def_var(ncid, var_mobility_name, NF90_DOUBLE, var_mobility_dimids &
         , var_mobility_id))

    ! Define scalar oil output
    call iserror(nf90_def_var(ncid, var_oil_name, NF90_DOUBLE, var_oil_id))

    ! End define mode.
    call iserror(nf90_enddef(ncid))

    ! Write scalar inputs
    call iserror(nf90_put_var(ncid, var_scene_id, scenario_id))
    call iserror(nf90_put_var(ncid, var_nx_id, nx))
    call iserror(nf90_put_var(ncid, var_ny_id, ny))
    call iserror(nf90_put_var(ncid, var_nz_id, nz))
!     call iserror(nf90_put_var(ncid, var_mu_id, mu))
!     call iserror(nf90_put_var(ncid, var_sigma_id, sigma))

    ! Write the time data. 
    call iserror(nf90_put_var(ncid, var_time_id, Tt))

    ! Write the data. This will write our mobility and oil data.
    ! The arrays of data are the same size as
    ! the netCDF variables we have defined.
    call iserror(nf90_put_var(ncid, var_mobility_id, Pc))
    call iserror(nf90_put_var(ncid, var_oil_id, totaloil))

    ! Close the file.
    call iserror(nf90_close(ncid))
  end subroutine write_results

  !
  ! This routine opens the permeability and porosity used by
  ! the MATLAB program and uses it for the simulation.
  !
  subroutine read_permeability_and_porosity(data_directory, nx, ny, nz, PERM, POR)   
    integer :: ncid                                                  ! netcdf file handle

    character(len = *) :: data_directory                             ! directory location of parameters 
    integer :: varid                                                 ! generic variable id
    integer :: nc_chunksize                                          ! chunk size for reading

    integer :: i, j, k, l, m
    integer :: nx, ny, nz
    
    double precision, dimension((nx * ny * nz)) :: POR       ! Porosities
    double precision, dimension((nx * ny * nz)) :: POR_temp  ! Porosities_temp
    double precision, dimension(3, nx, ny, nz)  :: PERM      ! Permeabilities

    ! initialize porosity and permeability to zero
    PERM = 0.0d0
    POR = 0.0d0
    POR_temp = 0.0d0

    ! read permeability
        nc_chunksize = 4096
    call iserror(nf90_open(trim(adjustl(data_directory))//"/permeability.nc", &
         NF90_NOWRITE, ncid, nc_chunksize))
    write(*,*) "Chosen chunk size is ", nc_chunksize
    
    call iserror(nf90_inq_varid(ncid, "permeability", varid))
    call iserror(nf90_get_var(ncid, varid, PERM))
    call iserror(nf90_close(ncid))

    ! read porosity
    nc_chunksize = 4096
    call iserror(nf90_open(trim(adjustl(data_directory))//"/porosity.nc", &
         NF90_NOWRITE, ncid, nc_chunksize))
    write(*,*) "Chosen chunk size is ", nc_chunksize
    
    call iserror(nf90_inq_varid(ncid, "porosity", varid))
    call iserror(nf90_get_var(ncid, varid, POR_temp))
    call iserror(nf90_close(ncid))

    call mymax_1_0_double(POR_temp, 1.0d-3, POR)
  end subroutine read_permeability_and_porosity

  ! netCDF Error Check Routine
  subroutine iserror(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Error encountered. Stopped."
    end if
  end subroutine iserror

end program runspe10
