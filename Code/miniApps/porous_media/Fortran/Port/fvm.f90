module fvm
    use grid
    use fluid
    use matrix
    use linsolve

    implicit none
contains
    subroutine tpfa(Q, P, V)
        implicit none
        double precision, dimension(:), pointer :: Q
        double precision, dimension(:,:,:), pointer :: P
        double precision, dimension(:,:,:,:), pointer :: V

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
        L = 1.0d0/K_

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
        V(1, 2:Nx_, :, :) = (P(1:Nx_ - 1, :, :) - P(2:Nx_, :, :)) * TX(2:Nx_,:,:)
        ! V.y
        V(2, :, 2:Ny_, :) = (P(:, 1:Ny_ - 1, :) - P(:, 2:Ny_, :)) * TY(:,2:Ny_,:)
        ! V.z
        V(3, :, :, 2:Nz_) = (P(:, :, 1:Nz_ - 1) - P(:, :, 2:Nz_)) * TZ(:,:,2:Nz_)

        ! free matrices
        call free_mat(A)
        call free_mat(TX)
        call free_mat(TY)
        call free_mat(TZ)
    end subroutine tpfa
end module fvm
