module conj_grad_mod

  public :: solve

  interface solve
     module procedure conj_grad
!     module procedure adconj_grad
  end interface

contains 


!-----------------------------
        subroutine conj_grad (x, b, A)
!-----------------------------

!      KRISHNA: THIS IS THE LINEAR SYSTEM SOLVER

        use stream_vel_variables

        real(8), intent(inout), dimension(n) :: x
        real(8), intent(in), dimension(n) :: b
        real(8), intent(in), dimension(n,3) :: A
        
        real(8) :: alpha, beta, dp1, dp2, res_init, resid
        integer :: i, ii, k_iter

        r(:) = 0.
        r_old(:) = 0.
        p(:) = 0.
        p_old(:) = 0.
        k_iter = 0
        x(:) = 0.

        ! r_0 = b - A(x_0)
        do i=1,n
         r(i) = b(i) - A(i,2) * x(i)
         if (i.gt.1) then
          r(i) = r(i) - A(i,1) * x(i-1)
         endif
         if (i.lt.n) then
          r(i) = r(i) - A(i,3) * x(i+1)
         endif
        enddo

        dp1 = 0.
        do i=1,n
         dp1 = dp1 + r(i)*r(i)
        enddo
        res_init = sqrt (dp1)
        resid = res_init

        p(:) = r(:)

        do ii=1,200*n

        if (resid.gt. 1e-10*res_init) then

         k_iter = k_iter + 1

         ! A(p_ii)
         do i=1,n
          ax(i) = A(i,2) * p(i)
          if (i.gt.1) then
           ax(i) = ax(i) + A(i,1) * p(i-1)
          endif
          if (i.lt.n) then
           ax(i) = ax(i) + A(i,3) * p(i+1)
          endif
         enddo


         ! alpha_ii = r'r / p'Ap
         dp1 = 0.
         dp2 = 0.
         do i=1,n
          dp1 = dp1 + r(i)*r(i)
          dp2 = dp2 + p(i)*ax(i)
         enddo
         alpha = dp1/dp2

         r_old(:) = r(:)

         x(:) = x(:) + alpha*p(:)
         r(:) = r(:) - alpha*ax(:)

         ! size of resid
         dp1 = 0.
         do i=1,n
          dp1 = dp1 + r(i)*r(i)
         enddo
         resid = sqrt(dp1)

         ! beta_ii = r'r / p'Ap
         dp1 = 0.
         dp2 = 0.
         do i=1,n
          dp1 = dp1 + r(i)*r(i)
          dp2 = dp2 + r_old(i)*r_old(i)
         enddo
         beta = dp1/dp2

         p(:) = r(:) + beta * p(:)          

        endif
        enddo


       end subroutine conj_grad
end module
