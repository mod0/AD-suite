program runspe10
use utils
use fvm_d
use grid_d
use fluid_d
use matrix_d
use gnufor2
use head_d


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

integer :: mu_count
integer :: sigma_count
parameter(mu_count = 10, sigma_count = 10)
double precision, dimension(mu_count) :: mus
double precision, dimension(sigma_count) :: sigmas

character(len = 30) :: filename

solver_inner_values = (/ 64, 128, 256/)
solver_outer_values = (/ 1000000 /)
                       
! Now read the permeabilities and porosities
call read_permeability_and_porosity(PERM, POR)

! compute ir
ir = (795.0 * Nx_ * Ny_ * Nz_) / (maxNx * maxNy * maxNz)

solver_inner = solver_inner_values(2)  ! 128
solver_outer = solver_outer_values(1)  ! 1000000

!Q(1:N_:Nx_*Ny_) = ir
!Q(Nx_*Ny_:N_:Nx_*Ny_) = -ir
!Q(1) = ir
!Q(N_) = -ir

mus = (/-1.500000000000000e+00, -1.166666666666667e+00, -8.333333333333334e-01,&
        -5.000000000000000e-01, -1.666666666666667e-01,  1.666666666666667e-01,&
        5.000000000000000e-01,  8.333333333333335e-01,   1.166666666666667e+00,&
        1.500000000000000e+00 /)
        
sigmas = (/ 1,2,3,4,5,6,7,8,9,10 /)


do i = 1, mu_count
 do j = 1, sigma_count
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
    mu = mus(i)
    sigma = sigmas(j)

    totaloil = 0.0d0

    !
    call cpu_time(start_time)
    call wrapper(ir, mu, sigma, Q, S, P, V, St, Pt, Tt, ND, Mw, Mo, Mt, &
    &            Pc, totaloil, solver_inner, solver_outer, verbose = .false.)
    call cpu_time(end_time) 
    write (*,*) mus(i), sigmas(j), totaloil, (end_time - start_time)
 enddo
enddo

! ---------------------------------------------------------------------------
! Call GNUPLOT through the interface module.
! Uncomment these plot calls after verifying you have GNUPlot installed.
! ---------------------------------------------------------------------------
! Plot the final solution for the forward trajectory
! call plot(Tt, Pc(1,:), Tt, Pc(2,:), terminal='png', filename='WaterOilProductionCurves.png')
return
end program runspe10
