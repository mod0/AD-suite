module head
  use grid
  use matrix
  use fvm

  implicit none
  integer :: St, Pt, ND

  parameter(St = 5,            &                  ! Max saturation time step
            Pt = 100,          &                  ! Pressure time step
            ND = 2000)                            ! Number of days in simulation
contains

!
! This routine opens the permeability and porosity used by
! the MATLAB program and uses it for the simulation.
!
subroutine read_permeability_and_porosity(PERM, POR)
    integer :: i, j, k, l, m

    double precision, dimension(N_) :: POR      ! Porosities
    double precision, dimension(3, Nx_, Ny_, Nz_) :: PERM  ! Permeabilities

    double precision, dimension(Nx_, Ny_, Nz_) :: P
    double precision, dimension(maxNx * maxNy * maxNz) :: pUr
    double precision, dimension(3 * maxNx, maxNy * maxNz) :: KUr
    double precision, dimension(3 * maxNx * maxNy * maxNz) :: KUrl

    integer, dimension(Nx_ * Ny_ * Nz_) :: Pindices
    integer, dimension(3 * Nx_ * Ny_ * Nz_) :: Kindices

    ! initialize porosity and permeability to zero
    PERM = 0.0d0
    POR = 0.0d0

    ! read KUr
    open(1,file='KUr.txt',status='old')
    read(1,*) ((KUr(i,j), j=1,maxNy * maxNz), i=1,3 * maxNx)
    close(1)

    ! reshape 2 dimension to 1 dimension
    call myreshape_2_1(KUr, KUrl)

    ! select according to specified dimension
    m = 0
    do l = 1, Nz_
        do k = 1,Ny_
            do j = 1,Nx_
                do i = 1,3
                    m = m + 1
                    Kindices(m) = ((l - 1) * (maxNx * maxNy * 3) &
                                  + (k - 1) * (maxNx * 3) &
                                  + 3 * (j-1) + i)
                end do
            end do
        end do
    end do

    ! then reshape 1 dimension to 4 dimension (hack for time being)
    call myreshape_1_4(KUrl(Kindices), PERM)

    ! read KUr
    open(1,file='pUr.txt',status='old')
    read(1,*) (pUr(i), i=1,maxNx * maxNy * maxNz)
    close(1)

    m = 0
    do k = 1,Nz_
        do j = 1,Ny_
            do i = 1,Nx_
                m = m + 1
                Pindices(m) = ((k - 1) * (maxNx * maxNy) &
                              + (j - 1) * (maxNx) + i)
            end do
        end do
    end do

    !POR = max(pUr(Pindices), 1.0d-3)
    call mymax_1_0_double(pUr(Pindices), 1.0d-3, POR)
end subroutine read_permeability_and_porosity


!
! Initialize inflow and outflow.
!
subroutine init_flw_trnc_norm_xin_pt_out(ir, mu, sigma, Q)
    !!use gnufor2
    integer :: i, j
    double precision, dimension(Nx_) :: idx
    double precision :: x, pi, pdf, mass
    double precision :: mu, sigma, ir
    double precision, dimension(N_) :: Q
    double precision, dimension(Nx_) :: Q_x

    ! value of pi
    pi = 3.14159265358979323d0

    !initialize the total mass to 0
    mass = 0.0d0
    Q_x = 0.0d0

    ! Note that the portion of the  Standard Normal distribution between
    ! -3sigma/2 to 3sigma/2 is assumed to fit the 1..Nx where sigma is 1
    do i = 1, Nx_
        ! get the real x coordinate
        x = -1.5d0 + ((i - 1) * 3.0d0)/(Nx_ - 1)    ! Mapping x = [-1.5, 1.5] to Nx_ dimension

        ! Now use mu and sigma to find the pdf value at x
        pdf = 1.0d0/(sigma * sqrt(2.0d0 * pi)) * exp(-(((x - mu)/sigma)**2.0d0)/2.0d0)

        ! set the value at the index equal to the pdf value at that point
        Q_x(i) = pdf

        ! increment the mass by the value of the pdf
        mass = mass + pdf

        ! index to test initialization by plot
        idx(i) = i * 1.0
    end do

    ! now rescale all the entities
    Q_x = Q_x/mass * ir

    ! Assign Q_x to Q
    j = 1
    do i = 1, Nx_* Ny_, Ny_
      Q(i) = Q_x(j)
      j = j + 1
    end do

    ! now set the output
    Q(N_) = -ir

   !---------------------------------------------------------------------------
   ! Call GNUPLOT through the interface module.
   ! Uncomment these plot calls after verifying you have GNUPlot installed.
   !---------------------------------------------------------------------------
   ! Plot the Q and check if it is correct.
   !call plot(idx, Q_x, terminal='png', filename='inflow.png')
end subroutine init_flw_trnc_norm_xin_pt_out


!
! This subroutine simulates the reservoir
! model.
!
subroutine simulate_reservoir(Q, S, P, V, St, Pt, Tt, ND, Mw, Mo, Mt, Pc, oil, &
                              solver_inner, solver_outer, verbose)
    use grid
    use fluid
    logical :: verbose
    integer :: solver_inner, solver_outer
    integer :: i, j, k, ND, St, Pt
    double precision :: Mw, Mo, Mt, tempoil1, tempoil2
    double precision ::  oil
    double precision, dimension((ND/St) + 1) :: Tt
    double precision, dimension(N_) :: S
    double precision, dimension(N_) :: Q
    double precision, dimension(2, (ND/St) + 1) :: Pc
    double precision, dimension(Nx_, Ny_, Nz_) :: P
    double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

    S = swc_                            ! initial saturation

    Pc(1, 1) = 0.0d0                    ! initial production
    Pc(2, 1) = 1.0d0
    Tt(1) = 0.0d0                       ! initial time.

    tempoil1 = 0.0d0
    tempoil2 = 0.0d0
    k = 1
    do i = 1, ND/Pt
        call Pres(S, Q, P, V, solver_inner, solver_outer, verbose) ! Pressure solver
!        write (*,*) 'PRESSURE========================================================'
         do j = 1, Pt/St
!            write(*,*) "Outer:", i, "Inner:", j
            k = k + 1
            call NewtRaph(S, V, Q, St, solver_inner, solver_outer, verbose) ! Solve for saturation
            call RelPerm(S(N_), Mw, Mo)             ! Mobilities in well-block

            Mt = Mw + Mo

            Tt(k) = 1.0d0 * k * St
            Pc(1,k) = Mw/Mt
            Pc(2,k) = Mo/Mt

            call update_oil(Pc, k, St, tempoil1, tempoil2)
            tempoil1 = tempoil2
        end do
    end do

    oil = tempoil2
end subroutine simulate_reservoir


subroutine update_oil(Pc, k, St, oilin, oilout)
  integer :: St, k
  double precision ::  oilin
  double precision ::  oilout
  double precision, dimension(2, (ND/St) + 1) :: Pc

  oilout = oilin +  Pc(2, k) * St                     ! Reimann sum
end subroutine update_oil

subroutine wrapper(ir, mu, sigma, Q, S, P, V, St, Pt, Tt, ND, Mw, Mo, &
&                  Mt, Pc, oil, solver_inner, solver_outer, verbose)
  use grid
  use fluid
  integer :: ND, St, Pt
  logical :: verbose
  integer :: solver_inner, solver_outer
  double precision :: mu, sigma, ir
  double precision :: Mw, Mo, Mt
  double precision :: oil
  double precision, dimension((ND/St) + 1) :: Tt
  double precision, dimension(N_) :: S
  double precision, dimension(N_) :: Q
  double precision, dimension(2, (ND/St) + 1) :: Pc
  double precision, dimension(Nx_, Ny_, Nz_) :: P
  double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

  call init_flw_trnc_norm_xin_pt_out(ir, mu, sigma, Q)
  call simulate_reservoir(Q, S, P, V, St, Pt, Tt, ND, Mw, Mo, Mt, Pc, &
&                         oil, solver_inner, solver_outer, verbose)

end subroutine wrapper
end module head
