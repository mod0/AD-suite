program test_sparse_matrix
    use matrix
    implicit none

    integer, dimension(3) :: idiags
    double precision, dimension(6,3) :: B
    type(spmat) :: matA, matB, matC, matD, matE, matF, matG, matH

! Get the results from running spdiags_matlabtest.m to compare implementation
! with MATLAB. Matches for all the below cases.

    ! instantiate B
    B(1,1) = 1
    B(2,1) = 2
    B(3,1) = 3
    B(4,1) = 4
    B(5,1) = 5
    B(6,1) = 6

    B(1,2) = 37
    B(2,2) = 38
    B(3,2) = 39
    B(4,2) = 30
    B(5,2) = 31
    B(6,2) = 32

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

    nullify(matA%row_index, matA%col_index, matA%values)
    nullify(matB%row_index, matB%col_index, matB%values)
    nullify(matC%row_index, matC%col_index, matC%values)
    nullify(matD%row_index, matD%col_index, matD%values)
    nullify(matE%row_index, matE%col_index, matE%values)
    nullify(matF%row_index, matF%col_index, matF%values)
    nullify(matG%row_index, matG%col_index, matG%values)
    nullify(matH%row_index, matH%col_index, matH%values)

    ! construct sparse matrices
    call spdiags(B, idiags, 6, 6, matA)
    call spdiags(B, idiags, 6, 5, matB)
    call spdiags(B, idiags, 5, 6, matC)
    call spdiags(B, idiags, 5, 5, matD)
    call spdiags(B, idiags, 5, 4, matE)
    call spdiags(B, idiags, 4, 5, matF)
    call spdiags(B, idiags, 2, 2, matG) ! Expect it will throw error or crash
                                        ! Interestingly code handles this case.
                                        ! So no.
    call spdiags(B, idiags, 1, 1, matH) ! Expect it will throw error or crash
                                        ! Interestingly code handles this case.
                                        ! So no.

    ! display matrices
    write (*,*) "Matrix A:"
    call disp_spmat(matA)
    write (*,*) "Matrix B:"
    call disp_spmat(matB)
    write (*,*) "Matrix C:"
    call disp_spmat(matC)
    write (*,*) "Matrix D:"
    call disp_spmat(matD)
    write (*,*) "Matrix E:"
    call disp_spmat(matE)
    write (*,*) "Matrix F:"
    call disp_spmat(matF)
    write (*,*) "Matrix G:"
    call disp_spmat(matG)
    write (*,*) "Matrix H:"
    call disp_spmat(matH)


    ! free sparse matrices
    call free_mat(matA)
    call free_mat(matB)
    call free_mat(matC)
    call free_mat(matD)
    call free_mat(matE)
    call free_mat(matF)
    call free_mat(matG)
    call free_mat(matH)
end program

