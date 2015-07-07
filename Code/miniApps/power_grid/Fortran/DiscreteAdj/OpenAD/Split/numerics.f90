module numerics
    implicit none

    ! Constants
    double precision :: pi, pmax, omegaB, omegaS, H, D, phiS, beta, c, t0, &
                        tend, tf, tcl, dt, perturb

    ! Shared
    double precision, save :: pm, quad_f, quad_b

    ! Set the parameter values
    parameter(                          &
        pi = 3.141592653589793d0,       &
        pmax = 2.0877d0,                &
        omegaB = 120.0d0 * pi,          &
        omegaS = 1d0,                   &
        H = 5.0d0,                      &
        D = 5.0d0,                      &
        phiS = 1.0d0,                   &
        beta = 2.0d0,                   &
        c = 10000.0d0,                  &
        t0 = 0.0d0,                     &
        tend = 10.0d0,                  &
        tf = 0.095d0,                   &  ! trying to evade floating pt issues
        tcl = 0.205d0,                  &  ! trying to evade floating pt issues
        dt = 0.01d0,                    &
        perturb = 1.0d0)

contains

    !
    ! Wrapper method which obtains the cost function and gradient
    ! used by the upper level optimization routine.
    !
    subroutine get_cost_function_and_gradient(f, g, tlen)
        integer :: tlen, i
        double precision, intent(out) :: f
        double precision, dimension(:), intent(out) :: g

        ! -------------------------------------------------------------------------
        ! declare variables for the integration scheme - use fixed time stepping
        ! may be use the SAVE attribute to avoid repeated allocation/deallocation
        ! -------------------------------------------------------------------------
        double precision, dimension(2) :: x
        double precision, dimension(2) :: lm
        double precision, dimension(tlen) :: tout_f
        double precision, dimension(tlen) :: tout_b
        double precision, dimension(tlen,2) :: yout_f
        double precision, dimension(tlen,2) :: yout_b

#ifdef USE_OPENAD
!$openad independent(pm)
!$openad dependent(f)
#endif
        ! Initialize the time array.
        tout_f(1) = t0
        do i = 2, tlen - 1
            tout_f(i) = tout_f(i - 1) + dt
        enddo
        tout_f(tlen) = tend

        ! Reverse the forward time array for adjoint use.
        tout_b = tout_f(tlen:1:-1)

        ! print the current parameter value
        ! print *, "pm is ", pm

        ! Initialize the x value
        x(1) = asin(perturb*pm/pmax)
        x(2) = perturb*1.0d0

        ! Compute the function value.
        call power_grid_cost_function(x, tout_f, yout_f, f)

#if !defined(USE_OPENAD)
        ! call the function that computes the sensitivites using
        ! continuous adjoints.

        ! Initialize the lambda value for the adjoint integration
        lm = (/0.0d0, 0.0d0/)

        ! Compute the gradient value.
        call power_grid_cost_gradient(lm, tout_b, yout_b, tout_f, yout_f, g)
#endif
    end subroutine get_cost_function_and_gradient


    !
    ! The subroutine computes power grid cost function
    ! This subroutine should call the time integrator
    ! in the forward mode along with the mode being
    ! set to "FWD".
    !
    subroutine power_grid_cost_function(x, tout_f, yout_f, f)
        implicit none

        double precision, dimension(2), intent(inout) :: x
        double precision, dimension(:), intent(inout) :: tout_f
        double precision, dimension(:,:), intent(out) :: yout_f
        double precision, intent(out) :: f

        call crank_nicolson(x, tout_f, yout_f, "FWD")

        f = -pm + quad_f
        !print *, "The cost function is ", f
    end subroutine power_grid_cost_function

#if !defined(USE_OPENAD)
    !
    ! The subroutine computes the sensitivity of cost
    ! function with respect to the parameter p
    ! This subroutine should call the time integrator
    ! in the adjoint mode along with the mode being
    ! set to "ADJ".
    !
    subroutine power_grid_cost_gradient(x, tout_b, yout_b, tout_f, yout_f, g)
        implicit none
        ! The gradient can be computed using continuous adjoints
        ! or by discrete adjoints

        double precision, dimension(2), intent(inout) :: x
        double precision, dimension(:), intent(inout) :: tout_b
        double precision, dimension(:,:), intent(out) :: yout_b
        double precision, dimension(:), intent(inout) :: tout_f
        double precision, dimension(:,:), intent(inout) :: yout_f
        double precision, dimension(:), intent(out) :: g

        call crank_nicolson(x, tout_b, yout_b, "ADJ", tout_f, yout_f)

        g(1) = -1.0d0 + quad_b &
                + 1.0d0/sqrt(pmax**2 - pm**2) * yout_b(size(yout_b, 1), 1)

        !print *, "The gradient is ", g(1)
    end subroutine power_grid_cost_gradient
#endif

    !
    ! Perform time integration by Crank Nicolson scheme.
    !
    subroutine crank_nicolson(x0, tout, yout, mode, tout_traj, yout_traj)
        implicit none

        logical :: converged
        integer :: max_conv_iter, iter, idx
        double precision :: dh, t, rcond
        double precision, dimension(:), intent(inout) :: x0
        integer, dimension(2) :: ipvt
        double precision, dimension(2) :: x_n_1, x_n, dx_n, &
                                          f_n_1, f_n, wz

        double precision :: conv_err
        double precision, dimension(2, 2) :: J_n, I_n
        double precision, dimension(:), intent(inout) :: tout
        double precision, dimension(:, :), intent(out) :: yout
        double precision, dimension(:), intent(inout), optional :: tout_traj
        double precision, dimension(:, :), intent(inout), optional :: yout_traj
        character(len = 3) :: mode

        ! set the maximum number of iterations (newton) before convergence.
        max_conv_iter = 10

        ! get the first time step
        t = tout(1)

        ! set the solution at the first time step as starting value
        yout(1, :) = x0

#if !defined(USE_OPENAD)
        if (mode == "FWD") then
            ! get the step size depending on whether we are going forward or
            ! backward. (FWD)
            dh = dt

            ! Initialize the forward quadrature value
            call forward_ode_quad(t, x0, "INIT")
        elseif (mode == "ADJ") then
            ! Check whether optional arguments are present
            if(.not.(present(tout_traj) .and. present(yout_traj))) then
                stop "The solution trajectories are not present for &
                &the adjoint mode"
            endif

            ! get the step size depending on whether we are going forward or
            ! backward. (FWD)
            dh = -dt

            ! Initialize the backward quadrature value
            call adjoint_ode_quad(t, x0, "INIT")
        endif
#else
        ! get the step size depending on whether we are going forward or
        ! backward. (FWD)
        dh = dt

        ! Initialize the forward quadrature value
        call forward_ode_quad(t, x0, "INIT")
#endif
        ! read the first time instant x value
        x_n_1 = x0

        ! start newton iteration with previous time instant value
        ! this assignment was initially inside the newton iteration
        ! loop. It has however been pulled out as it is duplicating
        ! the work done inside the norm check if-block
        x_n = x_n_1

        ! have an identity matrix handy
        call identity_matrix(I_n)

        ! compute all the other x values based on
        ! solving the nonlinear system (semi-implicit)
        ! using newton iteration.
        do idx = 2, size(tout)
#if !defined(USE_OPENAD)
            ! compute the f value at the previous time instant
            ! and initial previous x value.
            if (mode == "FWD") then
                call forward_ode_rhs(t, x_n_1, f_n_1)
            elseif (mode == "ADJ") then
                call adjoint_ode_rhs(t, x_n_1, tout_traj, yout_traj, f_n_1)
            endif
#else
            ! compute the f value at the previous time instant
            ! and initial previous x value.
            call forward_ode_rhs(t, x_n_1, f_n_1)
#endif
            ! get the new t value
            t = tout(idx)

            ! set converged to false
            converged = .false.

            ! try to converge within max_conv_iter
            do iter = 1, max_conv_iter
                if(.not. converged) then
#if !defined(USE_OPENAD)
                    ! get the function value and the jacobian
                    ! at the unconverged x and new t
                    if (mode == "FWD") then
                        call forward_ode_rhs(t, x_n, f_n)
                        call forward_ode_jac(t, x_n, J_n)
                    elseif (mode == "ADJ") then
                        call adjoint_ode_rhs(t, x_n, tout_traj, yout_traj, f_n)
                        ! Note that the jacobian in the adjoint mode is
                        ! independent of the current iterate. However,
                        ! to keep things symmetric between fwd and adj,
                        ! we are recomputing the jacobian each time.
                        ! This is wasteful, yes. But its a 2x2 system
                        ! in this particular case.
                        call adjoint_ode_jac(t, tout_traj, yout_traj, J_n)
                    endif
#else
                    ! get the function value and the jacobian
                    ! at the unconverged x and new t
                    call forward_ode_rhs(t, x_n, f_n)
                    call forward_ode_jac(t, x_n, J_n)
#endif

                    ! construct the left hand side matrix
                    J_n =  (I_n - dh/2 * J_n)

                    ! construct the right hand side vector
                    dx_n = x_n - x_n_1 - dh/2 * (f_n + f_n_1)

                    ! solve the system for dx_n
                    ! first factorize the matrix
                    ! call dgeco(J_n, 2, 2, ipvt, rcond, wz)

                    ! check the conditioning of the matrix
                    if (1.0d0 + rcond .eq. 1.0d0) then
                        !stop "Matrix is badly conditioned."
                        !print *,"Matrix is badly conditioned."
                    endif

                    ! solve the system
                    ! call dgesl(J_n, 2, 2, ipvt, dx_n, 0)
                    call linsolve2(J_n, dx_n, dx_n)

                    ! update the iterate
                    x_n = x_n - dx_n

                    ! check whether the right hand side vector with
                    ! new iterate has converged
                    ! Note that moving the if block down may change
                    ! the optimal slightly.
                    call dnrm2(dx_n, 2, conv_err)

                    if (conv_err < 1.0d-8) then
                        ! print *, "Converged at t = ", t
                        x_n_1 = x_n
                        f_n_1 = f_n
                        converged = .true.
                    endif
                else
                    ! do nothing.
                endif
            enddo

            ! update the new value in the solution trajectory
            yout(idx, :) = x_n


            ! Perform quadrature along
#if !defined(USE_OPENAD)
            if (mode == "FWD") then
                ! Iterate the forward quadrature value
                call forward_ode_quad(t, x_n, "ITER")
            elseif (mode == "ADJ") then
                ! Iterate the backward quadrature value
                call adjoint_ode_quad(t, x_n, "ITER")
            endif
#else
            call forward_ode_quad(t, x_n, "ITER")
#endif
        enddo


#if !defined(USE_OPENAD)
        ! Finalize the quadrature term.
        if (mode == "FWD") then
            ! Finalize the forward quadrature value
            call forward_ode_quad(t, x_n, "DONE")
        elseif (mode == "ADJ") then
            ! Finalize the backward quadrature value
            call adjoint_ode_quad(t, x_n, "DONE")
        endif
#else
        call forward_ode_quad(t, x_n, "DONE")
#endif
    end subroutine crank_nicolson

    !
    ! Subroutine initializes the values of I to
    ! be a square/non-square (meaningless) identity
    ! matrix.
    !
    subroutine identity_matrix(I_n)
        integer :: idx, jdx
        double precision, dimension(:,:), intent(out) :: I_n

        do idx = 1, size(I_n,1)
            do jdx = 1, size(I_n, 2)
                if (idx == jdx) then
                    I_n(idx,jdx) = 1
                else
                    I_n(idx,jdx) = 0
                endif
            enddo
        enddo
    end subroutine identity_matrix

    !
    ! The RHS function of the forward ODE
    !
    subroutine forward_ode_rhs(t, x, y)
        implicit none

        double precision, intent(inout) :: t
        double precision, dimension(2), intent(inout) :: x
        double precision, dimension(2), intent(out) :: y
        double precision :: phi
        double precision :: omega
        double precision :: pmax_

        phi = x(1)
        omega = x(2)

        if (t > tf .AND. t <= tcl) then
            pmax_ = 0.0d0
        else
            pmax_ = pmax
        endif

        y(1) = omegaB * (omega - omegaS)
        y(2) = omegaS * (pm - pmax_ * sin(phi) - D * (omega - omegaS))/(2*H)
    end subroutine forward_ode_rhs

    !
    ! The JAC function of the forward ODE
    !
    subroutine forward_ode_jac(t, x, J)
        implicit none

        double precision, intent(inout) :: t
        double precision, dimension(2), intent(inout) :: x
        double precision, dimension(2,2), intent(out) :: J
        double precision :: phi
        double precision :: omega
        double precision :: pmax_

        phi = x(1)
        omega = x(2)

        if (t > tf .AND. t <= tcl) then
            pmax_ = 0.0d0
        else
            pmax_ = pmax
        endif

        J(1, 1) = 0.0d0
        J(1, 2) = omegaB
        J(2, 1) = -omegaS * pmax_ * cos(phi) / (2.0d0 * H)
        J(2, 2) = -D * omegaS / (2.0d0 * H)
    end subroutine forward_ode_jac

    !
    ! The subroutine to compute the quadrature in FWD mode
    !
    subroutine forward_ode_quad(t, x, state)
        implicit none

        character(len = 4) :: state
        double precision, intent(inout) :: t
        double precision, dimension(2), intent(inout) :: x
        double precision :: phi, omega, new_h
        double precision,save :: prev_h

        phi = x(1)
        omega = x(2)

        if (state == "ITER") then
            ! perform trapezoidal integration using previous height
            ! and new height.
            new_h =  c * max(0.0d0, (phi - phiS))**beta
            quad_f = quad_f + 0.5d0 * dt * (prev_h + &
                       new_h)
            ! update the prev height.
            prev_h = new_h
        elseif (state == "INIT") then
            quad_f = 0.0d0
            prev_h = c * max(0.0d0, (phi - phiS))**beta
        elseif (state == "DONE") then
            !print *, "The forward quadrature value is ", quad_f
        endif
    end subroutine forward_ode_quad


#if !defined(USE_OPENAD)
    !
    ! The RHS function of the adjoint ODE
    ! Uses the trajectory of the forward solution to
    ! compute adjoint RHS function
    !
    subroutine adjoint_ode_rhs(t, x, tout_f, yout_f, y)
        implicit none

        integer :: idx
        double precision :: phi
        double precision, intent(inout) :: t
        double precision, dimension(2,2) :: J
        double precision, dimension(2) :: forcing
        double precision, dimension(2), intent(inout) :: x
        double precision, dimension(:), intent(inout) :: tout_f
        double precision, dimension(:, :), intent(inout) :: yout_f
        double precision, dimension(2), intent(out) :: y

        ! find the index of the current t
        idx = findindex(t, tout_f)

        ! get the jacobian of the forward ode
        call forward_ode_jac(t, yout_f(idx, :), J)

        ! get phi at the given time instant
        phi = yout_f(idx, 1)

        ! compute the forcing term
        forcing = (/ c * beta * max(0.0d0, &
                    (phi - phiS))**(beta - 1.0d0), 0.0d0 /)

        ! now compute y
        y = matmul(transpose(-J), x) - forcing
    end subroutine

    !
    ! The jacobian of the adjoint ODE
    ! Uses the trajectory of the forward solution to
    ! compute the jacobian at a given `t`
    !
    subroutine adjoint_ode_jac(t, tout_f, yout_f, J)
        implicit none

        integer :: idx
        double precision, intent(inout) :: t
        double precision, dimension(2,2), intent(out) :: J
        double precision, dimension(:), intent(inout) :: tout_f
        double precision, dimension(:, :), intent(inout) :: yout_f

        ! find the index of the current t
        idx = findindex(t, tout_f)

        ! get the jacobian of the forward ode
        call forward_ode_jac(t, yout_f(idx, :), J)

        J = transpose(-J)
    end subroutine

    !
    ! The subroutine to compute the quadrature in ADJ mode
    ! Uses the trajectory of the forward solution to compute
    ! the quadrature term
    !
    subroutine adjoint_ode_quad(t, x, state)
        implicit none

        character(len = 4) :: state
        double precision, intent(inout) :: t
        double precision, dimension(2), intent(inout) :: x
        double precision :: phi, omega, new_h
        double precision,save :: prev_h

        phi = x(1)
        omega = x(2)

        if (state == "ITER") then
            ! perform trapezoidal integration using previous height
            ! and new height.
            new_h = dot_product((/0.0d0, omegaS/(2.0d0*H)/), x)
            quad_b = quad_b + 0.5d0 * dt * (prev_h + &
                       new_h)
            ! update the previous height.
            prev_h = new_h
        elseif (state == "INIT") then
            quad_b = 0.0d0
            prev_h = dot_product((/0.0d0, omegaS/(2.0d0*H)/), x)
        elseif (state == "DONE") then
            !print *, "The adjoint quadrature value is ", quad_b
        endif
    end subroutine adjoint_ode_quad

    !
    ! Function returns the first index when traversing
    ! from lowest index (1) to highest index (length)
    ! whose value matches the desired value
    !
    integer function findindex(t, array)
        implicit none

        double precision, intent(inout) :: t
        double precision, dimension(:), intent(inout) :: array

        integer :: idx, default_index

        default_index = 0
        do idx = 1, size(array, 1)
            if(array(idx) == t) then
                default_index = idx
                exit
            endif
        enddo

        findindex = default_index
    end function findindex
#endif

    !
    ! Subroutine computes the 2 norm of the vector
    ! adapted from the original function version
    ! of the corresponding blas routine from
    ! NETLIB
    !
    subroutine dnrm2(v, len_v, n)
        integer :: i, len_v
        double precision :: n
        double precision, dimension(len_v) :: v
        double precision :: scale

        n = 0.0d0
        scale = 0.0d0

        do  i = 1, len_v
            scale = max(scale, abs(v(i)))
        enddo

        if (scale .eq. 0.0d0) then
            n = 0.0d0
        else
            do i = 1, len_v
                n = n + (v(i)/scale)**2
            enddo

            n = scale*sqrt(n)
        endif
    end subroutine dnrm2

!====================== The end of dnrm2 ===============================

    !
    ! Solve a 2x2 linear system
    !
    subroutine linsolve2(A, b, x)
        double precision, dimension(2,2), intent(inout) :: A
        double precision, dimension(2), intent(inout) :: b
        double precision, dimension(2), intent(out) :: x

        double precision, dimension(2) :: b_copy
        double precision :: b_copy_row
        double precision, dimension(2,2) :: A_copy
        double precision, dimension(2):: A_copy_row

        double precision :: factor

        ! Maintain a copy of A
        A_copy = A
        b_copy = b

        ! Swap the rows if the abs(A(1,1)) < abs(A(2,1))
        ! Not checking for zero coefficients
        if(dabs(A_copy(1,1)) < dabs(A_copy(2,1))) then
            A_copy_row = A_copy(2,:)
            A_copy(2,:) = A_copy(1,:)
            A_copy(1,:) = A_copy_row

            b_copy_row = b(2)
            b_copy(2) = b_copy(1)
            b_copy(1) = b_copy_row
        endif

        factor = (A_copy(2,1)/A_copy(1,1))

        A_copy(2, 1) = A_copy(2, 1) -  factor * A_copy(1, 1)
        A_copy(2, 2) = A_copy(2, 2) -  factor * A_copy(1, 2)
        b_copy(2) = b_copy(2) - factor * b_copy(1)

        x(2) = b_copy(2) / A_copy(2,2)
        x(1) = (b_copy(1) - (A_copy(1,2) * x(2)))/A_copy(1,1)
    end subroutine
!====================== The end of linsolve2 ===========================

end module numerics
