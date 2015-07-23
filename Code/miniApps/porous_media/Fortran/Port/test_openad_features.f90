program test_openad_features
    double precision :: x
    double precision, dimension(2,2) :: A
    double precision, dimension(2,2) :: B
    double precision, dimension(2,2) :: C
    double precision, dimension(4) :: D

!$openad independent(x)
!$openad dependent(D)

    x = 4.3456d0

    A(1,1) = 1.0d0 * x
    A(2,1) = 2.0d0 * x

    A(1,2) = 3.0d0 * x
    A(2,2) = 4.0d0 * x

    B(1,1) = 1.0d0
    B(2,1) = 2.0d0

    B(1,2) = 3.0d0
    B(2,2) = 4.0d0

    C = A + B

    ! Open AD does not handle reshape.
    !D = reshape(C, (/4/))

    ! Open AD does not handle the array literal with implied do loop either.
    ! D = (/((C(i,j), i = 1,2,1), j = 1,2,1)/)

    call myreshape_2_1(C, D)
!    D(1) = C(1,1)
!    D(2) = C(2,1)
!    D(3) = C(2,2)
!    D(4) = C(1,2)


    !write (*,*), D
contains

    subroutine myreshape_2_1(amatrix, bmatrix)
        integer :: i, j, k
        double precision, dimension(:,:) :: amatrix
        double precision, dimension(:) :: bmatrix

        k = 0

        do i = 1, size(amatrix, 1)
            do j = 1, size(amatrix, 2)
                k = k + 1
                bmatrix(k) = amatrix(i, j)
            end do
        end do
    end subroutine


end program test_openad_features
