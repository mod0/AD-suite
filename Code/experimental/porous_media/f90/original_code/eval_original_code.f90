program runspe10
  use parameters
  use gnufor2
  use utils
  use matrix
  use finitevolume
  use simulation
  use netcdf

  implicit none

  integer :: i, j, k
  integer :: nx, ny, nz
  integer :: nd, st, pt
  integer :: maxnx, maxny, maxnz

  ! input/intermediate variables
  double precision :: mu, sigma
  double precision, dimension(:), allocatable :: Q
  double precision, dimension(:), allocatable :: S
  double precision, dimension(:, :, :), allocatable :: P
  double precision, dimension(:, :, :, :), allocatable :: V

  ! output variables
  double precision :: totaloil
  double precision, dimension(:), allocatable :: Tt
  double precision, dimension(:, :), allocatable :: Pc

  ! netCDF variables
  integer :: ncid                                                  ! file handle
  integer :: var_scene_id                                          ! scene variable id
  character(len = *), parameter :: var_scene_name = "Scenario"     ! scene variable name
  integer :: var_nx_id                                             ! nx variable id
  character(len = *), parameter :: var_nx_name = "NX"              ! nx variable name
  integer :: var_ny_id                                             ! ny variable id
  character(len = *), parameter :: var_ny_name = "NY"              ! ny variable name
  integer :: var_nz_id                                             ! nz variable id
  character(len = *), parameter :: var_nz_name = "NZ"              ! nz variable name
  integer :: var_mu_id                                             ! mu variable id
  character(len = *), parameter :: var_mu_name = "Mu"              ! mu variable name
  integer :: var_sigma_id                                          ! sigma variable id
  character(len = *), parameter :: var_sigma_name = "Sigma"        ! sigma variable name
  integer :: dim_time_id                                           ! time dimension id
  integer :: var_time_id                                           ! time variable id
  integer :: dim_time_len                                          ! time dimension length
  character(len = *), parameter :: dim_time_name = "Time"          ! time dimension name
  integer :: dim_mobility_id                                       ! mobility dimension id
  integer :: var_mobility_id                                       ! mobility variable id
  integer :: dim_mobility_len                                      ! mobility dimension length
  character(len = *), parameter :: dim_mobility_name = "Mobility"  ! mobility dimension name
  integer, dimension(2) :: dimids                                  ! id of dimensions
  integer :: var_oil_id                                            ! oil variable id
  character(len = *), parameter :: var_oil_name = "Oil"            ! oil variable name

  ! initialize nx, ny, nz
  nx = 10
  ny = 10
  nz = 2

  ! initialize nd, pt, st
  nd = 2000
  pt = 100
  st = 5

  ! initialize maxnx, maxny, maxnz
  maxnx = 60
  maxny = 220
  maxnz = 85

  ! Allocate space for constant permeability/porosity
  allocate(POR(nx * ny * nz), PERM(3, nx, ny, nz))

  ! Allocate space for Q, S, P, V
  allocate(Q(nx * ny * nz), &
       S(nx * ny * nz), &
       P(nx, ny, nz), &
       V(3, nx+1, ny+1, nz+1))

  ! Allocate space for Time and Production data 
  allocate(Tt((nd/st) + 1), & 
       Pc(2, (nd/st) + 1))

  ! Compute the magic number 
  ir = (795.0 * nx * ny * nz) / (maxnx * maxny * maxnz) ! Magic number

  ! initialize simulation output
  Tt = 0.0d0               ! simulation time
  Pc = 0.0d0               ! production data
  totaloil = 0.0d0         ! total oil

  call initialize_scenario(scenario_id, nx, ny, nz, maxnx, maxny, maxnz, mu, sigma, Q, S, P, V)
  call wrapper(nx, ny, nz, nd, pt, st, mu, sigma, Q, S, P, V, Tt, Pc, totaloil)

  ! Start writing netCDF file having all the computed values
  ! Open file
  call iserror(nf90_create(results_eval_original_code, nf90_clobber, ncid))

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
  call iserror(nf90_def_var(ncid, var_mu_name, NF90_DOUBLE, var_mu_id))
  ! Define scalar sigma input
  call iserror(nf90_def_var(ncid, var_sigma_name, NF90_DOUBLE, var_sigma_id))

  ! Define the time variables. Varid is returned.
  call iserror(nf90_def_var(ncid, dim_time_name, NF90_DOUBLE, dim_time_id, var_time_id))

  ! Define the netCDF variables. The dimids array is used to pass the
  ! dimids of the dimensions of the netCDF variables.
  dimids = (/ dim_mobility_id, dim_time_id /)
  call iserror(nf90_def_var(ncid, dim_mobility_name, NF90_DOUBLE, dimids, var_mobility_id))

  ! Define scalar oil output
  call iserror(nf90_def_var(ncid, var_oil_name, NF90_DOUBLE, var_oil_id))

  ! End define mode.
  call iserror(nf90_enddef(ncid))

  ! Write scalar inputs
  call iserror(nf90_put_var(ncid, var_scene_id, scenario_id))
  call iserror(nf90_put_var(ncid, var_nx_id, nx))
  call iserror(nf90_put_var(ncid, var_ny_id, ny))
  call iserror(nf90_put_var(ncid, var_nz_id, nz))
  call iserror(nf90_put_var(ncid, var_mu_id, mu))
  call iserror(nf90_put_var(ncid, var_sigma_id, sigma))

  ! Write the time data. 
  call iserror(nf90_put_var(ncid, var_time_id, Tt))

  ! Write the data. This will write our mobility and oil data.
  ! The arrays of data are the same size as
  ! the netCDF variables we have defined.
  call iserror(nf90_put_var(ncid, var_mobility_id, Pc))
  call iserror(nf90_put_var(ncid, var_oil_id, totaloil))

  ! Close the file.
  call iserror(nf90_close(ncid))

  ! deallocate variables
  deallocate(POR, PERM)
  deallocate(Q, S, P, V)
  deallocate(Tt, Pc)

  return
contains

  subroutine initialize_scenario(scenario_no, nx, ny, nz, maxnx, maxny, maxnz, mu, sigma, Q, S, P, V)
    integer :: scenario_no
    integer :: nx, ny, nz
    integer :: maxnx, maxny, maxnz
    double precision :: mu, sigma
    double precision, dimension((nx * ny * nz)) :: Q
    double precision, dimension((nx * ny * nz)) :: S
    double precision, dimension(nx, ny, nz) :: P
    double precision, dimension(3, nx + 1, ny + 1, nz + 1) :: V

    ! initialize mu
    mu = 0.0d0
    sigma = 1.0d0

    ! initialize memory for inflow and saturation
    Q = 0.0d0
    S = 0.0d0

    ! initialize memory for P and V
    P = 0.0d0

    ! note that V vector has an additional
    ! length in each dimension x,y,z
    V = 0.0d0

    ! Now read the permeabilities and porosities which are global
    call read_permeability_and_porosity(nx, ny, nz, maxnx, maxny, maxnz, PERM, POR)
  end subroutine initialize_scenario

  ! netCDF Error Check Routine
  subroutine iserror(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Error encountered. Stopped."
    end if
  end subroutine iserror

end program runspe10
