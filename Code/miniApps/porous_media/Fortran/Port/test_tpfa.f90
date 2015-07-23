program test_tpfa
    use grid
    use fluid
    use matrix
    use fvm

    integer :: dealloc_err
    double precision, dimension(:), pointer :: Q
    double precision, dimension(:,:,:), pointer :: P
    double precision, dimension(:,:,:,:), pointer :: V

    ! setup permeabilities
    call ones(3, Nx_, Ny_, Nz_, K_)

    ! setup Flow?
    call zeros(N_, Q)

    Q(1) = 1
    Q(N_) = -1

    call tpfa(Q, P, V)

    ! free all allocated memory

end program
