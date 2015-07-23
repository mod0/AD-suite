program test_fortran_wholearray_ops

    double precision, dimension(3) :: A
    double precision, dimension(3) :: B
    double precision, dimension(3, 3) :: C
    double precision, dimension(3, 3) :: D
    double precision, dimension(3, 3) :: E
    double precision, dimension(3, 3) :: F

    A(1) = 1.0d0
    A(2) = 2.0d0
    A(3) = 3.0d0

    ! point-wise inverse: literal constants, constants and scalar variables are
    ! conformable with all arrays
    B = 1.0d0/A

    C(1,1) = 1.0d0
    C(2,1) = 2.0d0
    C(3,1) = 3.0d0

    C(1,2) = 4.0d0
    C(2,2) = 5.0d0
    C(3,2) = 6.0d0

    C(1,3) = 7.0d0
    C(2,3) = 8.0d0
    C(3,3) = 9.0d0

    D = 1.0d0/C


    E(1,1) = 1.5d0
    E(2,1) = 2.5d0
    E(3,1) = 3.5d0

    E(1,2) = 4.5d0
    E(2,2) = 5.5d0
    E(3,2) = 6.5d0

    E(1,3) = 7.5d0
    E(2,3) = 8.5d0
    E(3,3) = 9.5d0

    F = 1.0d0/(C + E + D)

    write (*,*) A, B, C, D, E, F

end program
