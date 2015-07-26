program runspe10
use fvm
use grid
use fluid
use matrix

implicit none

integer :: i, j, k, St, Pt, ND
double precision :: ir, Mw, Mo, Mt
double precision, dimension(:), pointer :: Tt, S, Q
double precision, dimension(:,:), pointer :: Pc
double precision, dimension(:,:,:), pointer :: P
double precision, dimension(:,:,:,:), pointer :: V

ir = 795 * (Nx_ * Ny_ / (Nx_ * Ny_ * Nz_))

call zeros(N_, Q)
Q(1:N_:Nx_*Ny_) = ir
Q(Nx_*Ny_:N_:Nx_*Ny_) = -ir

! setup P and V
call zeros(Nx_, Ny_, Nz_, P)

! note that V vector has an additional
! length in each dimension x,y,z
call zeros(3, Nx_ + 1, Ny_ + 1, Nz_ + 1, V)

! Initialize permeability and porosity for testing
call zeros(3, Nx_, Ny_, Nz_, K_)
call zeros(N_, Por_)
call inputKP()


St = 5                              ! Max saturation time step
Pt = 100                            ! Pressure time step
ND = 2000                           ! Number of days in simulation


call ones(N_, S)
S = S * swc_

call zeros(ND/St, Tt)               ! simulation time
call zeros(2, ND/St, Pc)            ! production data

Pc(1, 1) = 0.0d0                    ! initial production
Pc(2, 1) = 1.0d0
Tt(1) = 0.0d0                       ! initial time.

k = 1
do i = 1, ND/Pt
    call Pres(S, Q, P, V)                       ! Pressure solver
    do j = 1, Pt/St
        k = k + 1

        print *, "Time step: ", k

        call NewtRaph(S, V, Q, St * 1.0d0)      ! Solve for saturation
        call RelPerm(S(N_), Mw, Mo)             ! Mobilities in well-block

        Mt = Mw + Mo

        Tt(k) = k * St
        Pc(1,k) = Mw/Mt
        Pc(2,k) = Mo/Mt
    end do
end do

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
    integer, dimension(3 * Nx_ * Ny_ * Nz_) :: Kindices
    double precision, dimension(3366000) :: KUrl
    double precision, dimension(180,18700) :: KUr
    double precision, dimension(1122000) :: pUr

    ! read KUr
    open(1,file='KUr.txt',status='old')
    read(1,*) ((KUr(i,j), j=1,18700), i=1,180)
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
                    Kindices(m) = ((l-1) * (Nx_ * Ny_ * 3) + (k - 1) * (Nx_ * 3) &
                                    + 3 * (j-1) + i)
                end do
            end do
        end do
    end do

    call myreshape(KUrl(Kindices), K_)

    ! read KUr
    open(1,file='pUr.txt',status='old')
    read(1,*) (pUr(i), i=1,1122000)
    close(1)

    Por_ = max(pUr((/(((((l-1) *(Nx_ * Ny_) + (k - 1) * Nx_ &
                + j), j=1,Nx_), k=1,Ny_), l=1,Nz_)/)), 1.0d-3)
end subroutine inputKP

end program runspe10
