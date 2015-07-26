module fvm
use grid
use fluid
use matrix
use linsolve
use mathutil

implicit none

interface RelPerm
    module procedure RelPerm_scalar
    module procedure RelPerm_vector
end interface RelPerm
contains

!
! Performs Newton Raphson to solve for saturations
!
subroutine NewtRaph(S, V, Q, St)
    implicit none
    double precision :: St                             ! Maximum saturation Time Step
    double precision, dimension(:) :: S, Q
    double precision, dimension(:,:,:,:) :: V

    integer :: i, j, it
    type(spmat) :: A, B, dG
    logical :: converged
    double precision :: dt, dsn
    double precision, dimension(N_) :: S_copy, S_iter_copy, dtx, fi, fw, Mw, &
                                       Mo, dMw, dMo, dF, G, dS, bfw

    ! nullify the fields of A, B, dG
    nullify(A%row_index, A%col_index, A%values, &
            B%row_index, B%col_index, B%values, &
            dG%row_index, dG%col_index, dG%values)

    ! not yet converged
    converged = .false.

    ! Assemble system matrix
    call GenA(V, Q, A)

    ! copy S over
    S_copy = S

    ! set scaling factor
    it = 0.0d0

    do while(.not. converged)
        dt = St/(2**it)
        dtx = dt/(V_ * Por_)
        fi = max(Q, 0.0d0) * dtx

        call spmat_multiply(A, dtx, B, "POS")

        i = 0

        ! This loop is very badly implemented in the original MATLAB code
        do while (i < 2**it)
            j = 0
            i = i + 1
            dsn = 1.0d0
            S_iter_copy = S

            do while (dsn > 1.0d-3 .and. j < 10)
                call RelPerm(S, Mw, Mo, dMw, dMo)
                dF = dMw/(Mw + Mo) - (Mw/((Mw + Mo)**2) * (dMw + dMo))

                call spmat_multiply(B, dF, dG, "PRE")
                call addx_diagonal(dG, -1.0d0, 0)

                fw = Mw / (Mw + Mo)
                call spmat_multiply(B, fw, bfw, "PRE")

                G = S - S_iter_copy - bfw - fi

                call solve(dG, G, dS)
                S = S + dS
                call dnrm2(dS, N_, dsn)

                j = j + 1
            end do

            if (dsn > 1.0d-3) then
                i = 2**it               ! Breaks out of while loop.
                S = S_copy
            end if
        end do

        if (dsn < 1.0d-3) then
            converged = .true.
        else
            it = it + 1
        end if
    end do

    call free_mat(A)
    call free_mat(B)
    call free_mat(dG)
end subroutine NewtRaph

!
! Pressure Solver
!
subroutine Pres(S, Q, P, V)
    double precision, dimension(:) :: S, Q
    double precision, dimension(:,:,:) :: P
    double precision, dimension(:,:,:,:) :: V

    double precision, dimension(3 * N_) :: M
    double precision, dimension(3, Nx_, Ny_, Nz_) :: KM

    call RelPerm(S, M(1 : 3 * N_ : 3), M(2 : 3 * N_ : 3))

    M(1 : 3 * N_ : 3) = M(1 : 3 * N_ : 3) + M(2 : 3 * N_ : 3)
    M(2 : 3 * N_ : 3) = M(1 : 3 * N_ : 3)
    M(3 : 3 * N_ : 3) = M(1 : 3 * N_ : 3)

    call myreshape(M, KM)

    ! point-wise multiply
    KM = KM * K_

    call tpfa(KM, Q, P, V)
end subroutine


!
! Relative Permeabilities
!
subroutine RelPerm_vector(S, Mw, Mo, dMw, dMo)
    implicit none
    double precision, dimension(:) :: S, Mw, Mo
    double precision, dimension(:), optional :: dMw, dMo

    double precision, dimension(size(S, 1)) :: S_

    S_ = (S - swc_)/(1.0d0 - swc_ - sor_)   ! rescale saturation
    Mw = S_**2/vw_
    Mo = (1 - S_)**2/vo_

    if (present(dMo) .and. present(dMw)) then
        dMw = 2 * S_/vw_/(1 - swc_ - sor_)
        dMo = -2 * (1 - S_)/vo_/(1 - swc_ - sor_)
    endif
end subroutine RelPerm_vector

!
! Relative Permeabilities
!
subroutine RelPerm_scalar(S, Mw, Mo, dMw, dMo)
    implicit none
    double precision :: S, Mw, Mo
    double precision, optional :: dMw, dMo

    double precision :: S_

    S_ = (S - swc_)/(1.0d0 - swc_ - sor_)   ! rescale saturation
    Mw = S_**2/vw_
    Mo = (1 - S_)**2/vo_

    if (present(dMo) .and. present(dMw)) then
        dMw = 2 * S_/vw_/(1 - swc_ - sor_)
        dMo = -2 * (1 - S_)/vo_/(1 - swc_ - sor_)
    endif
end subroutine RelPerm_scalar

!
! Generate A matrix
!
subroutine GenA(V, Q, A)
    implicit none

    type(spmat) :: A
    double precision, dimension(:) :: Q
    double precision, dimension(:,:,:,:) :: V  ! V has an extra length
                                               ! across each x, y, z

    ! the matrix containing the diagonal entries
    double precision, dimension(N_, 7) :: diags

    ! reshape arrays first
    call myreshape(V(3,1:Nx_, 1:Ny_, 2:Nz_ + 1), diags(:, 1)) ! z2
    call myreshape(V(2,1:Nx_, 2:Ny_ + 1, 1:Nz_), diags(:, 2)) ! y2
    call myreshape(V(1,2:Nx_ + 1, 1:Ny_, 1:Nz_), diags(:, 3)) ! x2
    call myreshape(V(1,1:Nx_, 1:Ny_, 1:Nz_), diags(:, 5)) ! x1
    call myreshape(V(2,1:Nx_, 1:Ny_, 1:Nz_), diags(:, 6)) ! y1
    call myreshape(V(3,1:Nx_, 1:Ny_, 1:Nz_), diags(:, 7)) ! z1

    diags(:, 1) = max(diags(:,1), 0.0d0)
    diags(:, 2) = max(diags(:,2), 0.0d0)
    diags(:, 3) = max(diags(:,3), 0.0d0)
    diags(:, 5) = -min(diags(:,5), 0.0d0)
    diags(:, 6) = -min(diags(:,6), 0.0d0)
    diags(:, 7) = -min(diags(:,7), 0.0d0)
    diags(:, 4) = min(Q, 0.0d0) - diags(:, 5) - diags(:, 3) &
                                - diags(:, 6) - diags(:, 2) &
                                - diags(:, 7) - diags(:, 1)

    call spdiags(diags, (/ -Nx_ * Ny_, -Nx_, -1, 0, 1, Nx_, Nx_ * Ny_ /), &
                 N_, N_, A)
end subroutine GenA


!
! Two point flux approximation.
!
subroutine tpfa(K, Q, P, V)
    implicit none
    double precision, dimension(:) :: Q
    double precision, dimension(:,:,:) :: P
    double precision, dimension(:,:,:,:) :: V, K


    ! local variables
    double precision :: tx_, ty_, tz_
    double precision, dimension(:,:,:), pointer :: TX, TY, TZ

    ! the matrix containing the diagonal entries
    double precision, dimension(N_, 7) :: diags

    ! solution to the linear system
    double precision, dimension(N_) :: u

    ! point-wise inverse of permeability
    double precision, dimension(3,Nx_,Ny_,Nz_) :: L

    ! sparse matrix
    type(spmat) :: A
    nullify(A%row_index)
    nullify(A%col_index)
    nullify(A%values)

    ! get the point-wise inverse of the permeability matrix
    L = 1.0d0/K

    tx_ = 2.0d0 * hy_ * hz_ / hx_
    ty_ = 2.0d0 * hx_ * hz_ / hy_
    tz_ = 2.0d0 * hy_ * hx_ / hz_

    call zeros(Nx_ + 1, Ny_, Nz_, TX);
    call zeros(Nx_, Ny_ + 1, Nz_, TY);
    call zeros(Nx_, Ny_, Nz_ + 1, TZ);

    ! Compute transmissibilities by averaging harmonically
    TX(2:Nx_,:,:) = tx_/(L(1, 1:Nx_ - 1, :, :) + L(1, 2:Nx_, :, :))
    TY(:,2:Ny_,:) = ty_/(L(2, :, 1:Ny_ - 1, :) + L(2, :, 2:Ny_, :))
    TZ(:,:,2:Nz_) = tz_/(L(3, :, :, 1:Nz_ - 1) + L(3, :, :, 2:Nz_))

    call myreshape(-TX(1:Nx_,:,:), diags(:, 5))          ! -x1
    call myreshape(-TY(:,1:Ny_,:), diags(:, 6))          ! -y1
    call myreshape(-TZ(:,:,1:Nz_), diags(:, 7))          ! -z1
    call myreshape(-TX(2:Nx_ + 1,:,:), diags(:, 3))      ! -x2
    call myreshape(-TY(:,2:Ny_ + 1,:), diags(:, 2))      ! -y2
    call myreshape(-TZ(:,:,2:Nz_ + 1), diags(:, 1))      ! -z2

    ! Assemble discretization matrix
    diags(:, 4) = -(diags(:,1) + diags(:,2) + diags(:,3) &
                    + diags(:,5) + diags(:,6) + diags(:,7))

    call spdiags(diags, (/ -Nx_ * Ny_, -Nx_, -1, 0, 1, Nx_, Nx_ * Ny_ /), &
                 N_, N_, A)

    ! Increment the 1,1 element of A
    call addx_elem(A, K_(1,1,1,1) + K_(2,1,1,1) + K_(3,1,1,1), 1, 1)

    ! solve the linear system
    call solve(A, Q, u)

    ! reshape the solution
    call myreshape(u, P)

    ! V.x
    V(1, 2:Nx_, 1:Ny_, 1:Nz_) = (P(1:Nx_ - 1, :, :) - P(2:Nx_, :, :)) * TX(2:Nx_,:,:)
    ! V.y
    V(2, 1:Nx_, 2:Ny_, 1:Nz_) = (P(:, 1:Ny_ - 1, :) - P(:, 2:Ny_, :)) * TY(:,2:Ny_,:)
    ! V.z
    V(3, 1:Nx_, 1:Ny_, 2:Nz_) = (P(:, :, 1:Nz_ - 1) - P(:, :, 2:Nz_)) * TZ(:,:,2:Nz_)

    ! free matrices
    call free_mat(A)
    call free_mat(TX)
    call free_mat(TY)
    call free_mat(TZ)
end subroutine tpfa

end module fvm
