module linsolvewadjoint

    public :: solve

    interface solve
        module procedure linsolvewadjoint_stub(A, b, x)
    end interface

contains

    subroutine linsolvewadjoint_stub(A, b, x)
        use matrix

        !$openad xxx oad_template_linsolve.f90

        type(spmat) :: A
        double precision, dimension(:) :: b
        double precision, dimension(:) :: x
        double precision :: sum
        integer :: i

        sum = 0.0d0
        do i = 1, A%nnz
            sum = sum + A%values(i)
        end do

        x = b/sum
    end subroutine linsolvewadjoint_stub

end module
