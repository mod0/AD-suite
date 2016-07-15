module parameters
  ! FIXED PARAMETERS
  ! grid parameters 
  integer :: scenario_id
  double precision :: hx_, hy_, hz_, V_, ir

  parameter(scenario_id = 1, &
       hx_ = 20.0d0 * 0.3048d0, &    ! step size in x-direction
       hy_ = 10.0d0 * 0.3048d0, &    ! step size in y-direction
       hz_ = 2.0d0 * 0.3048d0, &     ! step size in z-direction
       V_ = hx_ * hy_ * hz_)        ! Volume of each grid cell

  ! fluid parameters
  double precision :: vw_, vo_, swc_, sor_
  parameter(vw_ = 3d-4, &      ! Viscosity of Water
            vo_ = 3d-3, &      ! Viscosity of Oil
            swc_ = 0.2d0, &    ! Saturation of water cut
            sor_ = 0.2d0)      ! Saturation of oil cut

  ! timestepping parameters
  integer :: St_, Pt_, ND_
  parameter(St_ = 5,            &                  ! Max saturation time step
            Pt_ = 100,          &                  ! Pressure time step
            ND_ = 2000)                            ! Number of days in simulation

  ! filenames
  character(*), parameter :: data_directory = "../data/data_1/"
  character(*), parameter :: results_directory = "results/"
  character(*), parameter :: porosity_file = data_directory//"pUr.txt"
  character(*), parameter :: permeability_file = data_directory//"KUr.txt"
  character(*), parameter :: results_eval_original_code = &
       results_directory//"results_eval_original_code.nc"
  character(*), parameter :: results_eval_deriv_tapenade_1_fwd = &
       results_directory//"results_eval_deriv_tapenade_1_forward.nc"
  character(*), parameter :: results_eval_deriv_tapenade_1_rev = &
       results_directory//"results_eval_deriv_tapenade_1_reverse.nc"

  ! PARAMETERS READ FROM FILE
  ! porosity and permeability parameters
  double precision, dimension(:), allocatable :: POR      ! Porosities
  double precision, dimension(:, :, :, :), allocatable :: PERM  ! Permeabilities  

  ! PARAMETERS SET IN DRIVER
  ! linear solver parameters
  logical :: verbose
  integer :: solver_inner, solver_outer

  parameter(verbose = .false.,   &              ! Verbose solver output
       solver_inner = 64,        &              ! Number of inner iterations
       solver_outer = 100000)                   ! Number of outer iterations
end module parameters
