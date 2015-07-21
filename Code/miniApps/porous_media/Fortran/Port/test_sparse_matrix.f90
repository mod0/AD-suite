program test_sparse_matrix
    use sparse_matrix
    implicit none

    integer, dimension(3) :: idiags
    double precision, dimension(6,3) :: B
    type(spmat) :: matA, matB, matC, matD

    ! instantiate B
    B(1,1) = 1
    B(2,1) = 2
    B(3,1) = 3
    B(4,1) = 4
    B(5,1) = 5
    B(6,1) = 6

    B(1,2) = 7
    B(2,2) = 8
    B(3,2) = 9
    B(4,2) = 10
    B(5,2) = 11
    B(6,2) = 12

    B(1,3) = 13
    B(2,3) = 14
    B(3,3) = 15
    B(4,3) = 16
    B(5,3) = 17
    B(6,3) = 18

    ! instantiate idiags
    idiags(1) = -2
    idiags(2) = 0
    idiags(3) = 2

    ! construct sparse matrices
    call spdiags(B, idiags, 6, 6, matA)

    ! display matrices
    call disp_spmat(matA)

    ! free sparse matrices
    call free_spmat(matA)
end program

