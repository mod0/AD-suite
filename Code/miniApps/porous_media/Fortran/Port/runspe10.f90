program runspe10
use fvm
use grid
use fluid
use matrix
use gnufor2
use print_active

implicit none

integer :: i, j, k, St, Pt, ND
double precision :: ir, Mw, Mo, Mt
double precision, dimension(:), pointer :: Tt, S, Q
double precision, dimension(:,:), pointer :: Pc
double precision, dimension(:,:,:), pointer :: P
double precision, dimension(:,:,:,:), pointer :: V

ir = (795.0 * Nx_ * Ny_ * Nz_) / (maxNx * maxNy * maxNz)

call zeros(N_, Q)
!Q(1:N_:Nx_*Ny_) = ir
!Q(Nx_*Ny_:N_:Nx_*Ny_) = -ir
Q(1) = ir
Q(N_) = -ir

! setup P and V
call zeros(Nx_, Ny_, Nz_, P)

! note that V vector has an additional
! length in each dimension x,y,z
call zeros(3, Nx_ + 1, Ny_ + 1, Nz_ + 1, V)

! Initialize permeability and porosity for testing
call zeros(3, Nx_, Ny_, Nz_, K_)
call zeros(N_, Por_)
call inputKP()                      ! Now read the permeabilities and porosities

St = 5                              ! Max saturation time step
Pt = 100                            ! Pressure time step
ND = 2000                           ! Number of days in simulation

call ones(N_, S)
S = S * swc_

call zeros((ND/St) + 1, Tt)               ! simulation time
call zeros(2, (ND/St) + 1, Pc)            ! production data

Pc(1, 1) = 0.0d0                    ! initial production
Pc(2, 1) = 1.0d0
Tt(1) = 0.0d0                       ! initial time.

k = 1
do i = 1, ND/Pt
    !call print_array(S,1,0,output)
    !call print_array(Q,1,0,output)
    call Pres(S, Q, P, V)                       ! Pressure solver
    !call print_array(P,1,0,1,0,1,0,output)
    !call print_array(V,1,0,1,0,1,0,1,0,output)
    do j = 1, Pt/St
        k = k + 1

        print *, "Timestep:", k


        call NewtRaph(S, V, Q, St * 1.0d0)      ! Solve for saturation

        !call print_array(S,1,0,output)

        call RelPerm(S(N_), Mw, Mo)             ! Mobilities in well-block

        Mt = Mw + Mo

        Tt(k) = 1.0d0 * k * St
        Pc(1,k) = Mw/Mt
        Pc(2,k) = Mo/Mt
    end do
end do

! Set filename to collect results.
open(unit = 1, file = 'Pc1', &
        form = 'unformatted', access = 'stream')

! write the production curve 1
write(unit=1) Pc(1, :)

! close files
close(unit=1)

! Set filename to collect results.
open(unit = 1, file = 'Pc2', &
        form = 'unformatted', access = 'stream')

! write the production curve 2
write(unit=1) Pc(2, :)

! close files
close(unit=1)


! Set filename to collect results.
open(unit = 1, file = 'Pc12', &
        form = 'unformatted', access = 'stream')

! write the both production curves
write(unit=1) Pc

! close files
close(unit=1)


! open file to write saturation data
open(unit=1, file='Sat', &
          form = 'unformatted', access='stream')
! Write saturation data
write(unit=1) S

! close file
close(unit=1)


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

    !print *, K_
    !print *, Por_
end subroutine inputKP

end program runspe10
