program runspe10
use fvm_d
use grid_d
use fluid_d
use matrix_d
use gnufor2
use head_d

!use print_active

implicit none

integer :: i, j, k
integer :: solver_inner, solver_outer
integer, dimension(3) :: solver_inner_values
integer, dimension(1) :: solver_outer_values
double precision :: start_time, end_time
double precision :: ir, Mw, Mo, Mt, mu, mud, sigma, sigmad, totaloil, totaloild, dummyTotalOil
double precision, dimension((ND/St) + 1) :: Tt
double precision, dimension(N_) :: S
double precision, dimension(N_) :: Q
double precision, dimension(2, (ND/St) + 1) :: Pc
double precision, dimension(2, (ND/St) + 1) :: Pcd
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
write (*,*) "Computing the directional derivative with respect to mu using tangent linear model"
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

! Perturb mu and get the total oil
mud = 1.0
sigmad = 0.0
call wrapper_d(ir, mu, mud, sigma, sigmad, Q, S, P, V, St, Pt, &
 &   Tt, ND, Mw, Mo, Mt, Pc, Pcd, dummyTotalOil, totaloild, solver_inner, &
 &   solver_outer, verbose = .false.)

write (*,*) "Total oil from tangent linear model (mu) is ", dummyTotalOil
write (*,*) "Directional derivative of total oil in mu direction is ", totaloild

! Set filename to collect results.
open(unit = 1, file = 'Pc1_mu_tanglin', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pc(1, :)
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pc2_mu_tanglin', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pc(2, :)
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pcd1_mu_tanglin', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pcd(1, :)
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pcd2_mu_tanglin', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pcd(2, :)
close(unit=1)

write (*,*) "End of tangent linear evaluation in mu direction"
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

open(unit = 1, file = 'totaloil_mu_perturbed', &
       form = 'unformatted', access = 'stream')
write(unit=1) totaloil_mu_perturbed(:)
close(unit=1)

finite_difference = (totaloil_mu_perturbed - totaloil)/eps
write (*,*) "Finite difference (mu): ", finite_difference

open(unit = 1, file = 'fd_minus_directional_derivative_mu', &
       form = 'unformatted', access = 'stream')
write(unit=1) abs(finite_difference - totaloild)
close(unit=1)

call plot(eps, abs(finite_difference - totaloild), terminal='png', filename='mu_perturbed_vs_fd.png')

write (*,*) "Computing the directional derivative with respect to sigma using tangent linear model"

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

! Perturb mu and get the total oil
mud = 0.0
sigmad = 1.0
call wrapper_d(ir, mu, mud, sigma, sigmad, Q, S, P, V, St, Pt, &
&   Tt, ND, Mw, Mo, Mt, Pc, Pcd, dummyTotalOil, totaloild, solver_inner, &
&   solver_outer, verbose = .false.)

write (*,*) "Total oil from tangent linear model (sigma) is ", dummyTotalOil
write (*,*) "Directional derivative of total oil in sigma direction is ", totaloild

! Set filename to collect results.
open(unit = 1, file = 'Pc1_sigma_tanglin', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pc(1, :)
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pc2_sigma_tanglin', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pc(2, :)
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pcd1_sigma_tanglin', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pcd(1, :)
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pcd2_sigma_tanglin', &
       form = 'unformatted', access = 'stream')
write(unit=1) Pcd(2, :)
close(unit=1)

write (*,*) "End of tangent linear evaluation in sigma direction"
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

open(unit = 1, file = 'totaloil_sigma_perturbed', &
       form = 'unformatted', access = 'stream')
write(unit=1) totaloil_sigma_perturbed(:)
close(unit=1)

finite_difference = (totaloil_sigma_perturbed - totaloil)/eps
write (*,*) "Finite difference (sigma): ", finite_difference


open(unit = 1, file = 'fd_minus_directional_derivative_sigma', &
       form = 'unformatted', access = 'stream')
write(unit=1) abs(finite_difference - totaloild)
close(unit=1)

call plot(eps, abs(finite_difference - totaloild), terminal='png', filename='sigma_perturbed_vs_fd.png')
return

contains

! FROM: http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90
! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      INTEGER                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)               ! assume the first is the min
      Location = Start                  ! record its position
      DO i = Start+1, End               ! start with next elements
         IF (x(i) < Minimum) THEN       !   if x(i) less than the min?
            Minimum  = x(i)             !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      INTEGER, INTENT(INOUT) :: a, b
      INTEGER                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1                  ! except for the last
         Location = FindMinimum(x, i, Size)     ! find min from this to last
         CALL  Swap(x(i), x(Location))  ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort


end program runspe10
