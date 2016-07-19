module parameters
  ! FIXED PARAMETERS
  ! grid parameters 
  integer :: scenario_id
  double precision :: hx_, hy_, hz_, V_, ir

  ! fluid parameters
  double precision :: vw_, vo_, swc_, sor_

  ! PARAMETERS READ FROM FILE
  ! porosity and permeability parameters
  double precision, dimension(:), allocatable :: POR      ! Porosities
  double precision, dimension(:, :, :, :), allocatable :: PERM  ! Permeabilities  

  ! PARAMETERS SET IN DRIVER
  ! linear solver parameters
  logical :: verbose
  integer :: solver_inner, solver_outer
end module parameters
