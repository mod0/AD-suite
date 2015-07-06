! stub for module LU

module conj_gradStub_mod

  public :: solve

  interface solve
     module procedure conj_gradStub
  end interface

contains 
      subroutine conj_gradStub (x, b, A)
    !$openad xxx template oad_template_conj_grad.f90
!-----------------------------

!      KRISHNA: THIS IS THE LINEAR SYSTEM SOLVER
        use stream_vel_variables
        real(8), intent(inout), dimension(n) :: x
        real(8), intent(in), dimension(n) :: b
        real(8), intent(in), dimension(n,3) :: A
        integer :: i
        integer :: j
        do i=1,n
         x(i) = A(i,1)+A(i,2)+A(i,3)+b(i)
        enddo

      end subroutine conj_gradStub

end module


