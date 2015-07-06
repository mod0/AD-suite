module conj_grad_ad_mod

  public :: adsolve

  interface adsolve
     module procedure adconj_grad
  end interface

contains 


       subroutine adconj_grad( x_state, b, b_state, x, a, ada )
!      will be called as
!      subroutine adconj_grad_2( x, adx, b, adb, a, ada )
!      KRISHNA: THIS IS THE MANUAL ADJOINT OF THE LINEAR SYSTEM SOLVER

        use stream_vel_variables_passive
        use conj_grad_mod

        real(8), intent(inout), dimension(n) :: b, b_state
        real(8), intent(inout), dimension(n) :: x, x_state
        real(8), intent(inout), dimension(n,3) :: A, adA
        
        real(8) :: alpha, beta, dp1, dp2, res_init, resid
        integer :: i, ii, k_iter
        real(8), dimension(n) :: x_in


        r(:) = 0.
        r_old(:) = 0.
        p(:) = 0.
        p_old(:) = 0.
        k_iter = 0
        x_state = 0.

        x_in = x

        call solve (x_state, b_state, A)

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

        do ii=1,20*n

        if (resid.gt. 1e-6*res_init) then

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



        do i=1,n
         if (i.gt.1) then
          ada (i,1) = ada(i,1) - x(i) * x_state(i-1)
         endif
         ada (i,2) = ada(i,2) - x(i) * x_state(i)
         if (i.lt.n) then
          ada (i,3) = ada(i,3) - x(i) * x_state(i+1)
         endif

        enddo

        

        x = x_in + x


        b(:) = 0._8

       end subroutine adconj_grad
end module
