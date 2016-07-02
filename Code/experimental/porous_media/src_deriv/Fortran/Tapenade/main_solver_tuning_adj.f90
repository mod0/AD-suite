program runspe10
use fvm_b
use grid_b
use fluid_b
use matrix_b
use gnufor2
use head_b
use utils

implicit none

integer :: i, j, k
integer :: solver_inner, solver_outer
integer, dimension(3) :: solver_inner_values
integer, dimension(1) :: solver_outer_values
double precision :: start_time, end_time
double precision :: ir, Mw, Mo, Mt, mu, mub, sigma, sigmab, totaloil, &
                    totaloilb, dummyTotalOil
double precision, dimension((ND/St) + 1) :: Tt
double precision, dimension(N_) :: S
double precision, dimension(N_) :: Q
double precision, dimension(2, (ND/St) + 1) :: Pc
double precision, dimension(2, (ND/St) + 1) :: Pcb
double precision, dimension(Nx_, Ny_, Nz_) :: P
double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

character(len = 30) :: filename

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
k = 1
do i = 1,3
  do j = 1,1
    if (solver_inner_values(i) <= N_) then
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

      solver_inner = solver_inner_values(i)
      solver_outer = solver_outer_values(j)

      write (*,*) "Function Evaluation: Inner: ", solver_inner, "Outer: ", solver_outer

      !
      call cpu_time(start_time)
      call wrapper(ir, mu, sigma, Q, S, P, V, St, Pt, Tt, ND, Mw, Mo, Mt,&
                   Pc, totaloil, solver_inner, solver_outer, verbose = .true.)
      call cpu_time(end_time)
      ! ---------------------------------------------------------------------------
      ! Call GNUPLOT through the interface module.
      ! Uncomment these plot calls after verifying you have GNUPlot installed.
      ! ---------------------------------------------------------------------------
      ! Plot the final solution for the forward trajectory
      call plot(Tt, Pc(1,:), Tt, Pc(2,:), terminal='png', filename='WaterOilProductionCurves.png')

      write (*,*) "Total oil is ", totaloil
      
      write(filename, '(A17, I1)') "Pc1_solver_tuning", k
      ! Set filename to collect results.
      open(unit = 1, file = filename, &
            form = 'unformatted', access = 'stream')
      write(unit=1) Pc(1, :)
      close(unit=1)

      write(filename, '(A17, I1)') "Pc2_solver_tuning", k
      ! Set filename to collect results.
      open(unit = 1, file = filename, &
            form = 'unformatted', access = 'stream')
      write(unit=1) Pc(2, :)
      close(unit=1)

      write (*,*)  "End of function evaluation. Elapsed time is", (end_time - start_time)
     
      k = k + 1
    endif
  enddo
enddo


k = 1
do i = 1,3
  do j = 1,1
    if (solver_inner_values(i) <= N_) then
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
      
      solver_inner = solver_inner_values(i)
      solver_outer = solver_outer_values(j)
      
      write (*,*) "Derivative Evaluation: Inner: ", solver_inner, "Outer: ", solver_outer

      call cpu_time(start_time)
      
      call wrapper_b(ir, mu, mub, sigma, sigmab, Q, S, P, V, St, Pt, &
&         Tt, ND, Mw, Mo, Mt, Pc, Pcb, dummyTotalOil, totaloilb, &
&         solver_inner, solver_outer, verbose = .true.)
      
      call cpu_time(end_time)
      
      write (*,*) "Total oil (computed in adjoint mode) is ", dummyTotalOil
      write (*,*) "Directional derivative of total oil in mu direction is ", mub
      write (*,*) "Directional derivative of total oil in sigma direction is ", sigmab

      ! Set filename to collect results.
      write(filename, '(A7, I1)') "Pc1_adj", k
      open(unit = 1, file = filename, &
            form = 'unformatted', access = 'stream')
      write(unit=1) Pc(1, :)
      close(unit=1)

      ! Set filename to collect results.
      write(filename, '(A7, I1)') "Pc2_adj", k
      open(unit = 1, file = filename, &
            form = 'unformatted', access = 'stream')
      write(unit=1) Pc(2, :)
      close(unit=1)

      ! Set filename to collect results.
      write(filename, '(A8, I1)') "Pcb1_adj", k
      open(unit = 1, file = filename, &
            form = 'unformatted', access = 'stream')
      write(unit=1) Pcb(1, :)
      close(unit=1)

      ! Set filename to collect results.
      write(filename, '(A8, I1)') "Pcb2_adj", k
      open(unit = 1, file = filename, &
            form = 'unformatted', access = 'stream')
      write(unit=1) Pcb(2, :)
      close(unit=1)
      
      write (*,*)  "End of derivative evaluation. Elapsed time is", (end_time - start_time)
     
      k = k + 1
    endif
  enddo
enddo

return
end program runspe10
