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

integer :: output
parameter(output = 0)

contains

!
! Performs Newton Raphson to solve for saturations
!
subroutine NewtRaph1(S, V, Q, St)
    use print_active
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
    it = 0

    ! print A
    !call disp_spmat(A, output)

    do while(.not. converged)
        dt = St/(2**it)
        dtx = dt/(V_ * Por_)
        fi = max(Q, 0.0d0) * dtx

        !call print_array(dtx,1,0,output)

        call spmat_multiply(A, dtx, B, "POS")

        !call disp_spmat(B, output)

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

                !call print_array(dF,1,0,output)

                call spmat_multiply(B, dF, dG, "PRE")

                !call disp_spmat(dG, output)

                call add_x(dG, -1.0d0, 0)

                !call disp_spmat(dG, output)

                fw = Mw / (Mw + Mo)
                call spmat_multiply(B, fw, bfw, "PRE")

                !call print_array(bfw,1,0,output)

                G = S - S_iter_copy - bfw - fi

                !call print_array(G,1,0,output)

                call solve(dG, G, dS)
                S = S + dS

                !call print_array(dS,1,0,output)

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
end subroutine NewtRaph1


!
! Performs Newton Raphson to solve for saturations
!
subroutine NewtRaph2(S, V, Q, St)
    use print_active
    implicit none
    double precision :: St                             ! Maximum saturation Time Step
    double precision, dimension(:) :: S, Q
    double precision, dimension(:,:,:,:) :: V
    integer :: i, j, it

    integer :: arows, acols, annz
    integer, dimension(:), pointer :: arow_index
    integer, dimension(:), pointer :: acol_index
    double precision, dimension(:), pointer :: avalues

    integer :: brows, bcols, bnnz
    integer, dimension(:), pointer :: brow_index
    integer, dimension(:), pointer :: bcol_index
    double precision, dimension(:), pointer :: bvalues

    integer :: dgrows, dgcols, dgnnz
    integer, dimension(:), pointer :: dgrow_index
    integer, dimension(:), pointer :: dgcol_index
    double precision, dimension(:), pointer :: dgvalues

    logical :: converged
    double precision :: dt, dsn
    double precision, dimension(N_) :: S_copy, S_iter_copy, dtx, fi, fw, Mw, &
                                       Mo, dMw, dMo, dF, G, dS, bfw

    ! nullify the fields of A, B, dG
    nullify(arow_index, acol_index, avalues, &
            brow_index, bcol_index, bvalues, &
            dgrow_index, dgcol_index, dgvalues)

    ! not yet converged
    converged = .false.

    ! Assemble system matrix
    call GenA2(V, Q, arows, acols, annz, arow_index, acol_index, avalues)

    ! copy S over
    S_copy = S

    ! set scaling factor
    it = 0

    do while(.not. converged)
        dt = St/(2**it)
        dtx = dt/(V_ * Por_)
        fi = max(Q, 0.0d0) * dtx

        ! Matrix-diagonal matrix product
        call spmat_multiply(arows, acols, annz, arow_index, acol_index, avalues, dtx, &
                            brows, bcols, bnnz, brow_index, bcol_index, bvalues, "POS")

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

                call spmat_multiply(brows, bcols, bnnz, brow_index, bcol_index, bvalues, dF, &
                                    dgrows, dgcols, dgnnz, dgrow_index, dgcol_index, dgvalues, "PRE")

                call add_x(dgrows, dgcols, dgnnz, dgrow_index, dgcol_index, dgvalues, -1.0d0, 0)

                fw = Mw / (Mw + Mo)
                call spmat_multiply(brows, bcols, bnnz, brow_index, bcol_index, bvalues, fw, bfw, "PRE")

                G = S - S_iter_copy - bfw - fi

                call solve(dgrows, dgcols, dgnnz, dgrow_index, dgcol_index, dgvalues, G, dS)
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

    call free_mat(arow_index, acol_index, avalues)
    call free_mat(brow_index, bcol_index, bvalues)
    call free_mat(dgrow_index, dgcol_index, dgvalues)
end subroutine NewtRaph2


!
! Pressure Solver
!
subroutine Pres(S, Q, P, V)
    use print_active
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

    ! Use tpfa if sparse matrix user defined type is desired
    ! Use tpfa2 if the type has to be flattented out.
    call tpfa2(KM, Q, P, V)
end subroutine


!
! Relative Permeabilities
!
subroutine RelPerm_vector(S, Mw, Mo, dMw, dMo)
    use print_active
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
    use print_active
    implicit none

    type(spmat) :: A
    double precision, dimension(:) :: Q
    double precision, dimension(:,:,:,:) :: V  ! V has an extra length
                                               ! across each x, y, z

    ! the matrix containing the diagonal entries
    double precision, dimension(N_, 7) :: diags

    !call print_array(V, 1,0,1,0,1,0,1,0,output)

    ! reshape arrays first
    call myreshape(V(3,1:Nx_, 1:Ny_, 2:Nz_ + 1), diags(:, 1)) ! z2
    call myreshape(V(2,1:Nx_, 2:Ny_ + 1, 1:Nz_), diags(:, 2)) ! y2
    call myreshape(V(1,2:Nx_ + 1, 1:Ny_, 1:Nz_), diags(:, 3)) ! x2
    call myreshape(V(1,1:Nx_, 1:Ny_, 1:Nz_), diags(:, 5)) ! x1
    call myreshape(V(2,1:Nx_, 1:Ny_, 1:Nz_), diags(:, 6)) ! y1
    call myreshape(V(3,1:Nx_, 1:Ny_, 1:Nz_), diags(:, 7)) ! z1

    !call print_array(diags,1,0,1,0,output)

    diags(:, 1) = max(diags(:,1), 0.0d0)
    diags(:, 2) = max(diags(:,2), 0.0d0)
    diags(:, 3) = max(diags(:,3), 0.0d0)
    diags(:, 5) = -min(diags(:,5), 0.0d0)
    diags(:, 6) = -min(diags(:,6), 0.0d0)
    diags(:, 7) = -min(diags(:,7), 0.0d0)
    diags(:, 4) = min(Q, 0.0d0) - diags(:, 5) - diags(:, 3) &
                                - diags(:, 6) - diags(:, 2) &
                                - diags(:, 7) - diags(:, 1)

    !call print_array(diags,1,0,1,0,output)

    ! This can be sped up by passing 3 arrays having rows, cols and diagind
    ! this can be done because diag positions are fixed.
    !TODO: Have a variant of spdiags which will instead of writing it in
    !the spmat type, it will write it in separate arrays sent to it.
    call spdiags(diags, (/ -Nx_ * Ny_, -Nx_, -1, 0, 1, Nx_, Nx_ * Ny_ /), &
                 N_, N_, A)
end subroutine GenA


!
! Generate A matrix
!
subroutine GenA2(V, Q, arows, acols, annz, arow_index, acol_index, avalues)
    use print_active
    implicit none

    integer :: arows, acols, annz
    integer, dimension(:),pointer :: arow_index
    integer, dimension(:), pointer :: acol_index
    double precision, dimension(:), pointer :: avalues
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

    ! This can be sped up by passing 3 arrays having rows, cols and diagind
    ! this can be done because diag positions are fixed.
    !TODO: Have a variant of spdiags which will instead of writing it in
    !the spmat type, it will write it in separate arrays sent to it.

    arows = N_
    acols = N_
    call spdiags(diags, (/ -Nx_ * Ny_, -Nx_, -1, 0, 1, Nx_, Nx_ * Ny_ /), &
                 arows, acols, annz, arow_index, acol_index, avalues)
end subroutine GenA2

!
! Two point flux approximation.
!
subroutine tpfa(K, Q, P, V)
    use print_active
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


    !TODO: Have a variant of spdiags which will instead of writing it in
    !the spmat type, it will write it in separate arrays sent to it.
    call spdiags(diags, (/ -Nx_ * Ny_, -Nx_, -1, 0, 1, Nx_, Nx_ * Ny_ /), &
                 N_, N_, A)


    ! Increment the 1,1 element of A
    call add_x(A, K_(1,1,1,1) + K_(2,1,1,1) + K_(3,1,1,1), 1, 1)

    ! solve the linear system
    ! Pass the rows_index, cols_index, values separately.
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


!
! Two point flux approximation.
!
subroutine tpfa2(K, Q, P, V)
    use print_active
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
    integer :: arows, acols, annz
    integer, dimension(:), pointer :: arow_index
    integer, dimension(:), pointer :: acol_index
    double precision, dimension(:), pointer :: avalues
    nullify(arow_index)
    nullify(acol_index)
    nullify(avalues)

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


    !TODO: Have a variant of spdiags which will instead of writing it in
    !the spmat type, it will write it in separate arrays sent to it.
    arows = N_
    acols = N_
    call spdiags(diags, (/ -Nx_ * Ny_, -Nx_, -1, 0, 1, Nx_, Nx_ * Ny_ /), &
                 arows, acols, annz, arow_index, acol_index, avalues)


    ! Increment the 1,1 element of A
    call add_x(arows, acols, annz, arow_index, acol_index, avalues, &
                    K_(1,1,1,1) + K_(2,1,1,1) + K_(3,1,1,1), 1, 1)

    ! solve the linear system
    ! Pass the rows_index, cols_index, values separately.
    call solve(arows, acols, annz, arow_index, acol_index, avalues, Q, u)

    ! reshape the solution
    call myreshape(u, P)

    ! V.x
    V(1, 2:Nx_, 1:Ny_, 1:Nz_) = (P(1:Nx_ - 1, :, :) - P(2:Nx_, :, :)) * TX(2:Nx_,:,:)
    ! V.y
    V(2, 1:Nx_, 2:Ny_, 1:Nz_) = (P(:, 1:Ny_ - 1, :) - P(:, 2:Ny_, :)) * TY(:,2:Ny_,:)
    ! V.z
    V(3, 1:Nx_, 1:Ny_, 2:Nz_) = (P(:, :, 1:Nz_ - 1) - P(:, :, 2:Nz_)) * TZ(:,:,2:Nz_)

    ! free matrices
    call free_mat(arow_index, acol_index, avalues)
    call free_mat(TX)
    call free_mat(TY)
    call free_mat(TZ)
end subroutine tpfa2

end module fvm
