program runspe10
use parameters
use gnufor2
use utils
use matrix
use finitevolume
use simulation

implicit none

integer :: i, j, k
double precision :: mu, sigma
double precision, dimension(N_) :: Q
double precision, dimension(N_) :: S
double precision, dimension(Nx_, Ny_, Nz_) :: P
double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

double precision :: totaloil
double precision, dimension((ND/St) + 1) :: Tt
double precision, dimension(2, (ND/St) + 1) :: Pc


character(len = 30) :: filename


! initialize simulation output
Tt = 0.0d0               ! simulation time
Pc = 0.0d0               ! production data
totaloil = 0.0d0         ! total oil

call initialize_scenario(1, mu, sigma, Q, S, P, V)

call wrapper(mu, sigma, Q, S, P, V, Tt, Pc, totaloil)

write (*,*) mu, sigma, totaloil
 
! ---------------------------------------------------------------------------
! Call GNUPLOT through the interface module.
! Uncomment these plot calls after verifying you have GNUPlot installed.
! ---------------------------------------------------------------------------
! Plot the final solution for the forward trajectory
call plot(Tt, Pc(1,:), Tt, Pc(2,:), terminal='png', filename='WaterOilProductionCurves.png')
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
end program runspe10
