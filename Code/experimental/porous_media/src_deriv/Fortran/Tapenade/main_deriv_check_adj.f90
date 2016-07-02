program runspe10
use fvm_b
use grid_b
use fluid_b
use matrix_b
use gnufor2
use head_b
use utils

!use print_active

implicit none

integer :: i, j, k
integer :: solver_inner, solver_outer
integer, dimension(3) :: solver_inner_values
integer, dimension(1) :: solver_outer_values
double precision :: start_time, end_time
double precision :: ir, Mw, Mo, Mt, mu, mub, sigma, sigmab, totaloil,&
                    totaloilb, dummyTotalOil
double precision, dimension((ND/St) + 1) :: Tt
double precision, dimension(N_) :: S
double precision, dimension(N_) :: Q
double precision, dimension(2, (ND/St) + 1) :: Pc
double precision, dimension(2, (ND/St) + 1) :: Pcb
double precision, dimension(Nx_, Ny_, Nz_) :: P
double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

integer :: eps_values_count
parameter(eps_values_count = 20)
double precision, dimension(eps_values_count) :: eps
double precision, dimension(eps_values_count) :: totaloil_mu_perturbed
double precision, dimension(eps_values_count) :: totaloil_sigma_perturbed
double precision, dimension(eps_values_count) :: finite_difference

character(len = 30) :: filename
! setup step sizes for the experiment
eps = &
(/ 1.000000000000000e-09, 3.359818286283788e-09, &
   1.128837891684688e-08, 3.792690190732253e-08, &
   1.274274985703132e-07, 4.281332398719396e-07, &
   1.438449888287663e-06, 4.832930238571752e-06, &
   1.623776739188721e-05, 5.455594781168514e-05, &
   1.832980710832438e-04, 6.158482110660266e-04, &
   2.069138081114790e-03, 6.951927961775605e-03, &
   2.335721469090121e-02, 7.847599703514607e-02, &
   2.636650898730356e-01, 8.858667904100832e-01, &
   2.976351441631313e+00, 1.000000000000000e+01/)

solver_inner_values = (/64, 128, 256/)
solver_outer_values = (/ 1000000 /)

                        
! Now read the permeabilities and porosities
call read_permeability_and_porosity(PERM, POR)

! compute ir
ir = (795.0 * Nx_ * Ny_ * Nz_) / (maxNx * maxNy * maxNz)

!Q(1:N_:Nx_*Ny_) = ir
!Q(Nx_*Ny_:N_:Nx_*Ny_) = -ir
!Q(1) = ir
!Q(N_) = -ir


! initialize memory for P and V
P = 0.0d0

! note that V vector has an additional
! length in each dimension x,y,z
V = 0.0d0

Tt = 0.0d0               ! simulation time
Pc = 0.0d0               ! production data

! initialize memory for inflow and saturation
Q = 0.0d0
S = 0.0d0

! initialize mu
mu = 0.0d0
sigma = 1.0d0
!

solver_inner = solver_inner_values(2)
solver_outer = solver_outer_values(1)

!
call wrapper(ir, mu, sigma, Q, S, P, V, St, Pt, Tt, ND, Mw, Mo, &
             Mt, Pc, totaloil, solver_inner, solver_outer, verbose = .false.)
! ---------------------------------------------------------------------------
! Call GNUPLOT through the interface module.
! Uncomment these plot calls after verifying you have GNUPlot installed.
! ---------------------------------------------------------------------------
! Plot the final solution for the forward trajectory
call plot(Tt, Pc(1,:), Tt, Pc(2,:), terminal='png', filename='WaterOilProductionCurves.png')

write (*,*) "Total oil is ", totaloil

! Set filename to collect results.
open(unit = 1, file = 'Pc1', &
      form = 'unformatted', access = 'stream')
write(unit=1) Pc(1, :)
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pc2', &
      form = 'unformatted', access = 'stream')
write(unit=1) Pc(2, :)
close(unit=1)

write (*,*) "End of function evaluation"
write (*,*) "Computing the derivative with respect to mu and sigma using adjoint model"
!
! initialize memory for P and V
P = 0.0d0

! note that V vector has an additional
! length in each dimension x,y,z
V = 0.0d0

Tt = 0.0d0               ! simulation time
Pc = 0.0d0               ! production data

! initialize memory for inflow and saturation
Q = 0.0d0
S = 0.0d0

! Set mu and sigma
mu = 0.0
sigma = 1.0

Pcb = 0.0d0
mub = 0.0d0
sigmab = 0.0d0

! Set the adjoint direction
totaloilb = 1.0d0

call wrapper_b(ir, mu, mub, sigma, sigmab, Q, S, P, V, St, Pt, &
&         Tt, ND, Mw, Mo, Mt, Pc, Pcb, totaloil, totaloilb, &
&         solver_inner, solver_outer, verbose = .false.) 
 
write (*,*) "Total oil from adjoint model is ", totaloil

! Set filename to collect results.
open(unit = 1, file = 'Pc1_adj', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pc(1, :)
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pc2_adj', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pc(2, :)
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pcb1_adj', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pcb(1, :)
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pcb2_adj', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pcb(2, :)
close(unit=1)

write (*,*) "mub = ", mub
write (*,*) "sigmab = ", sigmab
write (*,*) "End of adjoint model computation."

! Set filename to collect results.
open(unit = 1, file = 'Pc_porours_point_values', &
      form = 'unformatted', access = 'stream')
write(unit=1) totaloil, mub, sigmab
close(unit=1)

write (*,*) "Computing the effect of perturbation of mu on output of oil"

do i = 1,eps_values_count
  ! initialize memory for P and V
  P = 0.0d0

  ! note that V vector has an additional
  ! length in each dimension x,y,z
  V = 0.0d0

  Tt = 0.0d0               ! simulation time
  Pc = 0.0d0               ! production data

  ! initialize memory for inflow and saturation
  Q = 0.0d0
  S = 0.0d0

  write (*,*) "Perturbing mu by ", eps(i)

  call wrapper(ir, mu + eps(i), sigma, Q, S, P, V, St, Pt, &
               Tt, ND, Mw, Mo, Mt, Pc, totaloil_mu_perturbed(i),&
               solver_inner, solver_outer, verbose = .false.)

  write(filename, '(A20, I2)') "Pc1_mu_perturbed_eps", i
  open(unit = 1, file = filename, &
         form = 'unformatted', access = 'stream')
  write(unit=1) Pc(1,:)
  close(unit=1)

  write(filename, '(A20, I2)') "Pc2_mu_perturbed_eps", i
  open(unit = 1, file = filename, &
         form = 'unformatted', access = 'stream')
  write(unit=1) Pc(2,:)
  close(unit=1)
enddo

open(unit = 1, file = 'Pc_totaloil_mu_perturbed', &
       form = 'unformatted', access = 'stream')
write(unit=1) totaloil_mu_perturbed(:)
close(unit=1)

finite_difference = (totaloil_mu_perturbed - totaloil)/eps
write (*,*) "Finite difference (mu): ", finite_difference

open(unit = 1, file = 'Pc_fd_minus_directional_derivative_mu', &
       form = 'unformatted', access = 'stream')
write(unit=1) abs(finite_difference - mub)
close(unit=1)

write (*,*) "Computing the effect of perturbation of sigma on output of oil"

do i = 1,eps_values_count
  ! initialize memory for P and V
  P = 0.0d0

  ! note that V vector has an additional
  ! length in each dimension x,y,z
  V = 0.0d0

  Tt = 0.0d0               ! simulation time
  Pc = 0.0d0               ! production data

  ! initialize memory for inflow and saturation
  Q = 0.0d0
  S = 0.0d0

  write (*,*) "Perturbing sigma by ", eps(i)

  call wrapper(ir, mu, sigma + eps(i), Q, S, P, V, St, Pt, Tt, &
               ND, Mw, Mo, Mt, Pc, totaloil_sigma_perturbed(i), solver_inner, &
               solver_outer, verbose = .false.)

   write(filename, '(A23, I2)') "Pc1_sigma_perturbed_eps", i
   open(unit = 1, file = filename, &
          form = 'unformatted', access = 'stream')
   write(unit=1) Pc(1,:)
   close(unit=1)

   write(filename, '(A23, I2)') "Pc2_sigma_perturbed_eps", i
   open(unit = 1, file = filename, &
          form = 'unformatted', access = 'stream')
   write(unit=1) Pc(2,:)
   close(unit=1)
enddo

open(unit = 1, file = 'Pc_totaloil_sigma_perturbed', &
       form = 'unformatted', access = 'stream')
write(unit=1) totaloil_sigma_perturbed(:)
close(unit=1)

finite_difference = (totaloil_sigma_perturbed - totaloil)/eps
write (*,*) "Finite difference (sigma): ", finite_difference


open(unit = 1, file = 'Pc_fd_minus_directional_derivative_sigma', &
       form = 'unformatted', access = 'stream')
write(unit=1) abs(finite_difference - sigmab)
close(unit=1)
return
end program runspe10
