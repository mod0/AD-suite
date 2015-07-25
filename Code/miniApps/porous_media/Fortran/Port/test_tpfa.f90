program test_tpfa
    use grid
    use fluid
    use matrix
    use fvm

    implicit none

    character(len = 10) :: filename
    integer :: dealloc_err, p_filenum
    double precision, dimension(:), pointer :: Q
    double precision, dimension(:,:,:), pointer :: P
    double precision, dimension(:,:,:,:), pointer :: V

    p_filenum = 11

    ! Set filename to collect results.
    write(filename, '(A6, I2.2)') "output", p_filenum
    open(unit = p_filenum, file = filename, &
            form = 'unformatted', access = 'stream')

    ! setup permeabilities
    call ones(3, Nx_, Ny_, Nz_, K_)

    ! setup Flow?
    call zeros(N_, Q)

    ! setup P and V
    call zeros(Nx_, Ny_, Nz_, P)

    ! note that V vector has an additional
    ! length in each dimension x,y,z
    call zeros(3, Nx_ + 1, Ny_ + 1, Nz_ + 1, V)

    Q(1) = 1
    Q(N_) = -1

    call tpfa(K_, Q, P, V)

    ! write the solution
    write(unit=p_filenum) P

    ! close files
    close(unit=p_filenum)

    ! free all allocated memory
    call free_mat(Q)
    call free_mat(P)
    call free_mat(V)
    call free_mat(K_)
end program
