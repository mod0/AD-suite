module parameters
  ! FIXED PARAMETERS
  ! grid parameters 
  integer :: Nx_, Ny_, Nz_, N_
  double precision :: hx_, hy_, hz_, V_, ir
  integer :: maxNx, maxNy, maxNz

  parameter(maxNx=60, maxNy=220, maxNz=85)

  parameter(Nx_ = 10, &              ! Dimension in x-direction
       Ny_ = 10, &                   ! Dimension in y-direction
       Nz_ = 2, &                    ! Dimension in z-direction
       hx_ = 20.0d0 * 0.3048d0, &    ! step size in x-direction
       hy_ = 10.0d0 * 0.3048d0, &    ! step size in y-direction
       hz_ = 2.0d0 * 0.3048d0, &     ! step size in z-direction
       N_ = Nx_ * Ny_ * Nz_, &       ! Total number of grid cells
       V_ = hx_ * hy_ * hz_, &        ! Volume of each grid cell
       ir = (795.0 * Nx_ * Ny_ * Nz_) / (maxNx * maxNy * maxNz)) ! Magic number

  ! fluid parameters
  double precision :: vw_, vo_, swc_, sor_
  parameter(vw_ = 3d-4, &      ! Viscosity of Water
            vo_ = 3d-3, &      ! Viscosity of Oil
            swc_ = 0.2d0, &    ! Saturation of water cut
            sor_ = 0.2d0)      ! Saturation of oil cut

  ! timestepping parameters
  integer :: St, Pt, ND

  parameter(St = 5,            &                  ! Max saturation time step
            Pt = 100,          &                  ! Pressure time step
            ND = 2000)                            ! Number of days in simulation

  ! filenames
  character(*), parameter :: data_directory = "../../data/"
  character(*), parameter :: porosity_file = data_directory//"/shared/pUr.txt"
  character(*), parameter :: permeability_file = data_directory//"/shared/KUr.txt"
 

  ! PARAMETERS READ FROM FILE
  ! porosity and permeability parameters
  double precision, dimension(N_) :: POR      ! Porosities
  double precision, dimension(3, Nx_, Ny_, Nz_) :: PERM  ! Permeabilities  

  ! PARAMETERS SET IN DRIVER
  ! linear solver parameters
  logical :: verbose
  integer :: solver_inner, solver_outer
end module parameters
