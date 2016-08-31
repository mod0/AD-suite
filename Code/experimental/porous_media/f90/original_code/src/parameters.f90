module parameters
  use netcdfwrapper

  ! grid parameters 
  double precision :: hx, hy, hz, vol, ir

  ! fluid parameters
  double precision :: vw, vo, swc, sor

  ! porosity and permeability parameters
  double precision, dimension(:), allocatable :: POR            ! Porosities
  double precision, dimension(:, :, :, :), allocatable :: PERM  ! Permeabilities  

  ! linear solver parameters
  logical :: verbose
  integer :: solver_inner, solver_outer
end module parameters
