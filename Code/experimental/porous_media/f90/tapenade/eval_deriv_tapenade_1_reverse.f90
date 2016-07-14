program runspe10
use parameters_b
use gnufor2
use utils
use matrix_b
use finitevolume_b
use simulation_b
use netcdf

implicit none

integer :: i, j, k

! input/intermediate variables
double precision :: mu, sigma
double precision :: mub, sigmab
double precision, dimension(N_) :: Q
double precision, dimension(N_) :: S
double precision, dimension(Nx_, Ny_, Nz_) :: P
double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

! output variables
double precision :: totaloil
double precision :: totaloilb
double precision, dimension((ND/St) + 1) :: Tt
double precision, dimension(2, (ND/St) + 1) :: Pc

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
integer :: var_mu_id                                                ! mu variable id
character(len = *), parameter :: var_mu_name = "Mu"                 ! mu variable name
integer :: var_sigma_id                                             ! sigma variable id
character(len = *), parameter :: var_sigma_name = "Sigma"           ! sigma variable name
integer :: dim_time_id                                              ! time dimension id
integer :: var_time_id                                              ! time variable id
integer, parameter :: dim_time_len = (ND/St) + 1                    ! time dimension length
character(len = *), parameter :: dim_time_name = "Time"             ! time dimension name
integer :: dim_mobility_id                                          ! mobility dimension id
integer :: var_mobility_id                                          ! mobility variable id
integer, parameter :: dim_mobility_len = 2                          ! mobility dimension length
character(len = *), parameter :: dim_mobility_name = "Mobility"     ! mobility dimension name
integer, dimension(2) :: dimids                                     ! id of dimensions
integer :: var_oil_id                                               ! oil variable id
character(len = *), parameter :: var_oil_name = "Oil"               ! oil variable name
integer :: var_oil_mub_id                                           ! oil_mu variable id
character(len = *), parameter :: var_oil_mub_name = "Oil_mu"        ! oil_mu variable name
integer :: var_oil_sigmab_id                                        ! oil_sigma variable id
character(len = *), parameter :: var_oil_sigmab_name = "Oil_sigma"  ! oil_sigma variable name

! initialize simulation output
Tt = 0.0d0               ! simulation time
Pc = 0.0d0               ! production data
totaloil = 0.0d0         ! total oil

! initialize the adjoint direction
totaloilb = 1.0d0
mub = 0.0d0
sigmab = 0.0d0

! Calling the adjoint function to compute the sensitivities
call initialize_scenario(1, mu, sigma, Q, S, P, V)
call wrapper_b(mu, mub, sigma, sigmab, Q, S, P, V, Tt, Pc, totaloil, totaloilb)

! Calling the function to get meaningful values for Pc, Tt
call initialize_scenario(1, mu, sigma, Q, S, P, V)
call wrapper(mu, sigma, Q, S, P, V, Tt, Pc, totaloil)
 
! Start writing netCDF file having all the computed values
! Open file
call iserror(nf90_create(results_eval_deriv_tapenade_1_rev, nf90_clobber, ncid))

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
! Define scalar oil sensitivity in mu direction
call iserror(nf90_def_var(ncid, var_oil_mub_name, NF90_DOUBLE, var_oil_mub_id))
! Define scalar oil sensitivity in sigma direction
call iserror(nf90_def_var(ncid, var_oil_sigmab_name, NF90_DOUBLE, var_oil_sigmab_id))

! End define mode.
call iserror(nf90_enddef(ncid))

! Write scalar inputs
call iserror(nf90_put_var(ncid, var_scene_id, scenario_id))
call iserror(nf90_put_var(ncid, var_nx_id, Nx_))
call iserror(nf90_put_var(ncid, var_ny_id, Ny_))
call iserror(nf90_put_var(ncid, var_nz_id, Nz_))
call iserror(nf90_put_var(ncid, var_mu_id, mu))
call iserror(nf90_put_var(ncid, var_sigma_id, sigma))

! Write the time data. 
call iserror(nf90_put_var(ncid, var_time_id, Tt))

! Write the data. This will write our mobility and oil data.
! The arrays of data are the same size as
! the netCDF variables we have defined.
call iserror(nf90_put_var(ncid, var_mobility_id, Pc))
call iserror(nf90_put_var(ncid, var_oil_id, totaloil))
call iserror(nf90_put_var(ncid, var_oil_mub_id, mub))
call iserror(nf90_put_var(ncid, var_oil_sigmab_id, sigmab))

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
