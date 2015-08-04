program runspe10
use fvm
use grid
use fluid
use matrix
!use gnufor2
!use print_active

implicit none

integer :: St, Pt, ND
double precision :: ir, Mw, Mo, Mt, mu, sigma, oil
double precision, dimension(:), pointer :: Tt, S, Q
double precision, dimension(:,:), pointer :: Pc
double precision, dimension(:,:,:), pointer :: P
double precision, dimension(:,:,:,:), pointer :: V

St = 5                              ! Max saturation time step
Pt = 100                            ! Pressure time step
ND = 2000                           ! Number of days in simulation

! setup memory for P and V
call zeros(Nx_, Ny_, Nz_, P)

! note that V vector has an additional
! length in each dimension x,y,z
call zeros(3, Nx_ + 1, Ny_ + 1, Nz_ + 1, V)

call zeros((ND/St) + 1, Tt)               ! simulation time
call zeros(2, (ND/St) + 1, Pc)            ! production data

! setup memory for inflow and saturation
call zeros(N_, Q)
call zeros(N_, S)

! Initialize permeability and porosity for testing
call zeros(3, Nx_, Ny_, Nz_, K_)
call zeros(N_, Por_)
call inputKP()                      ! Now read the permeabilities and porosities

! compute ir
ir = (795.0 * Nx_ * Ny_ * Nz_) / (maxNx * maxNy * maxNz)

!Q(1:N_:Nx_*Ny_) = ir
!Q(Nx_*Ny_:N_:Nx_*Ny_) = -ir
!Q(1) = ir
!Q(N_) = -ir

!$openad independent(mu)
!$openad independent(sigma)
!$openad independent(sigma)

! initialize mu
mu = 0.0d0
sigma = 1.0d0

call inflow_truncated_normal_x_outflow_point(Q, ir, mu, sigma)
call simulate_reservoir()

!! Set filename to collect results.
!open(unit = 1, file = 'Pc1', &
!        form = 'unformatted', access = 'stream')
!
!! write the production curve 1
!write(unit=1) Pc(1, :)
!
!! close files
!close(unit=1)
!
!! Set filename to collect results.
!open(unit = 1, file = 'Pc2', &
!        form = 'unformatted', access = 'stream')
!
!! write the production curve 2
!write(unit=1) Pc(2, :)
!
!! close files
!close(unit=1)

!---------------------------------------------------------------------------
! Call GNUPLOT through the interface module.
! Uncomment these plot calls after verifying you have GNUPlot installed.
!---------------------------------------------------------------------------
! Plot the final solution for the forward trajectory
call plot(Tt, Pc(1,:), Tt, Pc(2,:), terminal='png', filename='WaterOilCut.png')


!call print_array(Pc, 1,0,1,0,output)

call free_mat(Tt)
call free_mat(S)
call free_mat(Q)
call free_mat(Pc)
call free_mat(P)
call free_mat(V)

contains


!
! This routine opens the permeability and porosity used by
! the MATLAB program and uses it for the simulation.
!
subroutine inputKP()
    integer :: i, j, k, l, m

    double precision, dimension(maxNx * maxNy * maxNz) :: pUr
    double precision, dimension(3 * maxNx, maxNy * maxNz) :: KUr
    double precision, dimension(3 * maxNx * maxNy * maxNz) :: KUrl

    integer, dimension(Nx_ * Ny_ * Nz_) :: Pindices
    integer, dimension(3 * Nx_ * Ny_ * Nz_) :: Kindices

    ! read KUr
    open(1,file='KUr.txt',status='old')
    read(1,*) ((KUr(i,j), j=1,maxNy * maxNz), i=1,3 * maxNx)
    close(1)

    ! reshape 2 dimension to 1 dimension
    call myreshape(KUr, KUrl)
    ! then reshape 1 dimension to 4 dimension (hack for time being)
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

    call myreshape(KUrl(Kindices), K_)

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

    Por_ = max(pUr(Pindices), 1.0d-3)
end subroutine inputKP


!
! Initialize inflow and outflow.
!
subroutine inflow_truncated_normal_x_outflow_point(Q, ir, mu, sigma)
!    use gnufor2
    integer :: i
!    double precision, dimension(Nx_) :: idx
    double precision :: x, pi, pdf, mu, sigma, ir, mass
    double precision, dimension(:), target :: Q
    double precision, dimension(:), pointer :: Q_x

    ! set Q_x to be only the values along the first plane and x direction
    Q_x => Q(1:Nx_ * Ny_:Ny_)

    ! value of pi
    pi = 3.14159265358979323d0

    !initialize the total mass to 0
    mass = 0.0d0

    ! Note that the portion of the  Standard Normal distribution between
    ! -3sigma/2 to 3sigma/2 is assumed to fit the 1..Nx
    do i = 1, Nx_
        ! get the real x coordinate
        x = -1.5d0 + ((i - 1) * 3.0d0)/(Nx_ - 1)

        ! Now use mu and sigma to find the pdf value at x
        pdf = 1.0d0/(sigma * sqrt(2.0d0 * pi)) * exp(-(((x - mu)/sigma)**2.0d0)/2.0d0)

        ! set the value at the index equal to the pdf value at that point
        Q_x(i) = pdf

        ! increment the mass by the value of the pdf
        mass = mass + pdf

!        ! index to test initialization by plot
!        idx(i) = i * 1.0
    end do

    ! now rescale all the entities
    Q_x = Q_x/mass * ir

    ! now set the output
    Q(N_) = -ir

!    !---------------------------------------------------------------------------
!    ! Call GNUPLOT through the interface module.
!    ! Uncomment these plot calls after verifying you have GNUPlot installed.
!    !---------------------------------------------------------------------------
!    ! Plot the Q and check if it is correct.
!    call plot(idx, Q_x, terminal='png', filename='inflow.png')
end subroutine inflow_truncated_normal_x_outflow_point


!
! This subroutine simulates the reservoir
! model.
!
subroutine simulate_reservoir()
    integer :: i, j, k

    S = swc_                            ! initial saturation

    Pc(1, 1) = 0.0d0                    ! initial production
    Pc(2, 1) = 1.0d0
    Tt(1) = 0.0d0                       ! initial time.

    oil = 0.0d0

    k = 1
    do i = 1, ND/Pt
        call Pres(S, Q, P, V)                       ! Pressure solver

        do j = 1, Pt/St
            k = k + 1
            call NewtRaph(S, V, Q, St * 1.0d0)      ! Solve for saturation
            call RelPerm(S(N_), Mw, Mo)             ! Mobilities in well-block

            Mt = Mw + Mo

            Tt(k) = 1.0d0 * k * St
            Pc(1,k) = Mw/Mt
            Pc(2,k) = Mo/Mt

            oil = oil + Pc(2, k) * St                     ! Reimann sum
        end do
    end do
end subroutine simulate_reservoir

end program runspe10
