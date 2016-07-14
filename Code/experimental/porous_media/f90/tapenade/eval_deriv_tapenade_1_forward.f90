program runspe10
use parameters_d
use gnufor2
use utils
use matrix_d
use finitevolume_d
use simulation_d
use netcdf

implicit none

integer :: i, j, k

! input/intermediate variables
double precision :: mu, sigma
double precision, dimension(N_) :: Q
double precision, dimension(N_) :: S
double precision, dimension(Nx_, Ny_, Nz_) :: P
double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

! output variables
double precision :: totaloil
double precision, dimension((ND/St) + 1) :: Tt
double precision, dimension(2, (ND/St) + 1) :: Pc

! netCDF variables
integer :: ncid                                                  ! file handle
integer :: dim_time_id                                           ! time dimension id
integer :: var_time_id                                           ! time variable id
integer, parameter :: dim_time_len = (ND/St) + 1                 ! time dimension length
character(len = *), parameter :: dim_time_name = "Time"          ! time dimension name
integer :: dim_mobility_id                                       ! mobility dimension id
integer :: var_mobility_id                                       ! mobility variable id
integer, parameter :: dim_mobility_len = 2                       ! mobility dimension length
character(len = *), parameter :: dim_mobility_name = "Mobility"  ! mobility dimension name
integer, dimension(2) :: dimids                                  ! id of dimensions
integer :: var_oil_id                                            ! oil variable id
character(len = *), parameter :: var_oil_name = "Oil"            ! oil variable name

! initialize simulation output
Tt = 0.0d0               ! simulation time
Pc = 0.0d0               ! production data
totaloil = 0.0d0         ! total oil

call initialize_scenario(1, mu, sigma, Q, S, P, V)
call wrapper(mu, sigma, Q, S, P, V, Tt, Pc, totaloil)
 
! Start writing netCDF file having all the computed values
! Open file
call iserror(nf90_create(results_eval_deriv_tapenade_1_fwd, nf90_clobber, ncid))

! Define time dimensions
call iserror(nf90_def_dim(ncid, dim_time_name, dim_time_len, dim_time_id))
! Define mobility dimensions
call iserror(nf90_def_dim(ncid, dim_mobility_name, dim_mobility_len, dim_mobility_id))

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

! Write the time data. 
call iserror(nf90_put_var(ncid, var_time_id, Tt))

! Write the data. This will write our mobility and oil data.
! The arrays of data are the same size as
! the netCDF variables we have defined.
call iserror(nf90_put_var(ncid, var_mobility_id, Pc))
call iserror(nf90_put_var(ncid, var_oil_id, totaloil))

! Close the file.
call iserror(nf90_close(ncid))

return
contains

subroutine initialize_scenario(scenario_no, mu, sigma, Q, S, P, V)
  integer :: scenario_no
  double precision :: mu, sigma
  double precision, dimension(N_) :: Q
  double precision, dimension(N_) :: S
  double precision, dimension(Nx_, Ny_, Nz_) :: P
  double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

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
  call read_permeability_and_porosity(PERM, POR)
end subroutine

! netCDF Error Check Routine
subroutine iserror(status)
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop "Error encountered. Stopped."
  end if
end subroutine iserror

end program runspe10
