program test_linsolve
    use matrix
    use linsolve
    implicit none

    type(spmat) :: matA
    integer, dimension(3) :: idiags
    double precision, dimension(6) :: x, u
    double precision, dimension(6,3) :: B

! Get the results from running spdiags_matlabtest.m to compare implementation
! with MATLAB. Matches for all the below cases.

    ! instantiate B
    B(1,1) = 1.0d0
    B(2,1) = 2.0d0
    B(3,1) = 3.0d0
    B(4,1) = 4.0d0
    B(5,1) = 5.0d0
    B(6,1) = 6.0d0

    B(1,2) = 37.0d0
    B(2,2) = 38.0d0
    B(3,2) = 39.0d0
    B(4,2) = 30.0d0
    B(5,2) = 31.0d0
    B(6,2) = 32.0d0

    B(1,3) = 13.0d0
    B(2,3) = 14.0d0
    B(3,3) = 15.0d0
    B(4,3) = 16.0d0
    B(5,3) = 17.0d0
    B(6,3) = 18.0d0

    ! instantiate idiags
    idiags(1) = -2
    idiags(2) = 0
    idiags(3) = 2

    nullify(matA%row_index, matA%col_index, matA%values)

    call spdiags(B, idiags, 6, 6, matA)

    call disp_spmat(matA)

    u(1) = 1.7d0
    u(2) = 3.5d0
    u(3) = 7.7d0
    u(4) = -2.11d0
    u(5) = 4.3d0
    u(6) = 15.8d0

    call solve(matA, u, x)

    print *, x
end program test_linsolve
