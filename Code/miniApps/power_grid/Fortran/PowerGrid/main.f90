program power_grid
    use power_grid_constants_and_shared
    use print_active
    implicit none

    ! -------------------------------------------------------------------------
    ! declare variables for the main optimization routine
    ! -------------------------------------------------------------------------
    ! we are using L-BFGS-B(v.3) from Nocedal et al.

    ! n is the number of parameters - pm
    ! m is the maximum number of limited memory corrections.
    integer :: np, nx, node, m, max_opt_iter, tlen, alloc_state
    parameter (np = 1, nx = 2, node = 2, m = 10, max_opt_iter = 20)

    ! Declare the variables needed by the code.
    ! Some of these variables are from LBFGS-B,
    ! rest are either passed between methods or
    ! shared by methods.
    ! Shared (globally accessible) : pm
    ! Passed: x, tout_f, tout_b, yout_f, yout_b
    character*60     :: task, csave
    character*10      :: filename
    logical          :: lsave(4)
    integer          :: i, iprint, nbd(np), iwa(3*np), isave(44)
    double precision :: f, factr, pgtol, &
                        l(np), u(np), g(np), dsave(29), &
                        wa(2*m*np + 5*np + 11*m*m + 8*m)

    ! -------------------------------------------------------------------------
    ! declare variables for the integration scheme - use fixed time stepping
    ! -------------------------------------------------------------------------
    double precision, dimension(nx) :: x
    double precision, dimension(node) :: lm
    double precision, dimension(:), allocatable :: tout_f
    double precision, dimension(:), allocatable :: tout_b
    double precision, dimension(:,:), allocatable :: yout_f
    double precision, dimension(:,:), allocatable :: yout_b

    ! -------------------------------------------------------------------------
    ! declare variables to hold file name numbers
    ! -------------------------------------------------------------------------
    integer :: tout_f_filenum, yout_f_filenum, tout_b_filenum, yout_b_filenum

    ! Check that the number of parameters does not exceed one.
    ! The formulation of the problem has to be changed otherwise.
    if (np > 1) then
        stop "The number of parameters is greater than one."
    endif

    ! -------------------------------------------------------------------------
    ! initialize variables/constants for the main optimization routine
    ! -------------------------------------------------------------------------
    ! We wish to have output at every iteration.
    iprint = 1

    ! We specify the tolerances in the stopping criteria.
    factr = 1.0d+1
    pgtol = 0.0d0

    ! Set the bound on the parameter p
    nbd(1) = 3
    u(1) = 1.1

    ! -------------------------------------------------------------------------
    ! initialize variables/constants for time integration scheme - fixed t.step
    ! -------------------------------------------------------------------------
    tlen = int((tend - t0)/dt) + 1            ! Verify that the tlen is correct
                                              ! Number of endpoints of intervals



    ! Set the filename numbers
    tout_f_filenum = 11
    yout_f_filenum = 12
    tout_b_filenum = 13
    yout_b_filenum = 14

    ! Now allocate those many end points.
    allocate(tout_f(tlen), yout_f(tlen, nx), &
             tout_b(tlen), yout_b(tlen, nx), stat = alloc_state)

    ! Check if allocation was successful.
    if (alloc_state /= 0) then
        stop "Unable to allocate memory for tout and yout"
    endif

    ! Initialize the time array.
    tout_f(1) = t0
    do i = 2, tlen - 1
        tout_f(i) = tout_f(i - 1) + dt
    enddo
    tout_f(tlen) = tend

    ! Reverse the forward time array for adjoint use.
    tout_b = tout_f(tlen:1:-1)

    ! Set the initial value of the parameter pm
    pm = 0.4

    ! Set the task variable
    task = "START"

    ! Iterate maximum number of iterations calling the optimization routine.
    do i = 1, max_opt_iter
        ! Set filename to collect results.
!        write(filename, '(A6, I2, I2)') "output", tout_f_filenum, i
!        open(unit = tout_f_filenum, file = filename)
        write(filename, '(A6, I2, I02)') "output", yout_f_filenum, i
        open(unit = yout_f_filenum, file = filename)
!        write(filename, '(A6, I2, I2)') "output", tout_b_filenum, i
!        open(unit = tout_b_filenum, file = filename)
        write(filename, '(A6, I2, I02)') "output", yout_b_filenum, i
        open(unit = yout_b_filenum, file = filename)

        call setulb(np, m, pm, l, u, nbd, f, g, factr, pgtol, wa, iwa, task, &
                    iprint, csave,lsave,isave,dsave)

        if (task(1:2) == "FG") then
            ! Initialize the x value
            x = (/asin(perturb*pm/pmax), perturb*1.0/)

            ! Compute the function value.
            call power_grid_cost_function(x, tout_f, yout_f, f)

#ifdef USE_OPENAD
            ! call the OPENAD generated discrete adjoint to get the
            ! gradient used for optimization.
#else
            ! call the function that computes the sensitivites using
            ! continuous adjoints.

            ! Initialize the lambda value for the adjoint integration
            lm = (/0.0d0, 0.0d0/)

            ! Compute the gradient value.
            call power_grid_cost_gradient(lm, tout_b, yout_b, tout_f, yout_f, g)
#endif

            ! write the solutions to the file.
            call write_array2(yout_f, 1, 0, 1, &
                    0, yout_f_filenum)
            call write_array2(yout_b, 1, 0, 1, &
                    0, yout_b_filenum)
        elseif (task(1:5) == "NEW_X") then
            if (i < max_opt_iter) then
                ! continue the iteration at the new point
                print *, "Continuing with new starting value. pm is ", pm
                cycle
            else
                ! stop the iteration. print current values
                print *, "The value of pm is ", pm
                exit
            endif
        elseif (task(1:4) == "CONV") then
            ! stop the iteration. print current values
            print *, "Convergence to within tolerance achieved pm is ", pm
            exit
        elseif (task(1:4) == "ABNO") then
            ! stop the iteration. print the current values
            ! also print the reason for abnormal termination
            print *, "Optimization terminated abnormally. pm is ", pm
            exit
        elseif (task(1:5) == "ERROR") then
            ! stop the iteration. print the current values
            ! also print the reason for abnormal termination
            print *, "Optimization encountered an error. pm is ", pm
            exit
        else
            print *, "The task is ", task
        endif

        close(yout_f_filenum)
        close(yout_b_filenum)
    enddo
contains

    !
    ! The subroutine computes power grid cost function
    ! This subroutine should call the time integrator
    ! in the forward mode along with the mode being
    ! set to "FWD".
    !
    subroutine power_grid_cost_function(x, tout_f, yout_f, f)
        use power_grid_constants_and_shared
        implicit none

        double precision, dimension(2), intent(in) :: x
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:,:), intent(out) :: yout_f
        double precision, intent(out) :: f

        call crank_nicolson(x, tout_f, yout_f, "FWD")

        f = -pm + quad_f
        print *, "The cost function is ", f
    end subroutine power_grid_cost_function

    !
    ! The subroutine computes the sensitivity of cost
    ! function with respect to the parameter p
    ! This subroutine should call the time integrator
    ! in the adjoint mode along with the mode being
    ! set to "ADJ".
    !
    subroutine power_grid_cost_gradient(x, tout_b, yout_b, tout_f, yout_f, g)
        use power_grid_constants_and_shared
        implicit none
        ! The gradient can be computed using continuous adjoints
        ! or by discrete adjoints

        double precision, dimension(2), intent(in) :: x
        double precision, dimension(:), intent(in) :: tout_b
        double precision, dimension(:,:), intent(out) :: yout_b
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:,:), intent(out) :: yout_f
        double precision, dimension(:), intent(out) :: g

        call crank_nicolson(x, tout_b, yout_b, "ADJ", tout_f, yout_f)

        g(1) = -1 + quad_b + dot_product( (/1/sqrt(pmax**2 - pm**2) , 0.0d0/), &
                            yout_b(size(yout_b, 1), :))
        print *, "The gradient is ", g(1)
    end subroutine power_grid_cost_gradient

    !
    ! Perform time integration by Crank Nicolson scheme.
    !
    subroutine crank_nicolson(x0, tout, yout, mode, tout_traj, yout_traj)
        use power_grid_constants_and_shared
        implicit none

        integer :: max_conv_iter, iter, idx
        double precision :: dh, t, rcond
        double precision, dimension(:), intent(in) :: x0
        integer, dimension(2) :: ipvt
        double precision, dimension(2) :: x_n_1, x_n, dx_n, &
                                          f_n_1, f_n, wz
        double precision, dimension(2, 2) :: J_n, I_n
        double precision, dimension(:), intent(in) :: tout
        double precision, dimension(:, :), intent(out) :: yout
        double precision, dimension(:), intent(in), optional :: tout_traj
        double precision, dimension(:, :), intent(in), optional :: yout_traj
        character(len = 3) :: mode

        ! add reference to external dnrm2 function
        interface
            function dnrm2(n,x,incx)
                integer :: dnrm2
                integer :: n,incx
                double precision :: x(n)
            end function dnrm2
        end interface

        ! set the maximum number of iterations (newton) before convergence.
        max_conv_iter = 10

        ! get the first time step.
        t = tout(1)

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

        ! read the first time instant x value
        x_n_1 = x0

        ! have an identity matrix handy
        call identity_matrix(I_n)

        ! compute all the other x values based on
        ! solving the nonlinear system (semi-implicit)
        ! using newton iteration.
        do idx = 2, size(tout)
            ! compute the f value at the previous time instant
            ! and initial previous x value.
            if (mode == "FWD") then
                call forward_ode_rhs(t, x_n_1, f_n_1)
            elseif (mode == "ADJ") then
                call adjoint_ode_rhs(t, x_n_1, tout_traj, yout_traj, f_n_1)
            endif

            ! get the new t value
            t = tout(idx)

            ! start newton iteration with previous time instant value
            ! this assignment is not needed - it can be pulled out
            ! of the loop.
            x_n = x_n_1

            ! try to converge within max_conv_iter
            do iter = 1, max_conv_iter
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

                ! construct the left hand side matrix
                J_n =  (I_n - dh/2 * J_n)

                ! construct the right hand side vector
                dx_n = x_n - x_n_1 - dh/2 * (f_n + f_n_1)

                ! check whether the right hand side vector with
                ! new iterate has converged
                if (dnrm2(2,dx_n, 1) < 1.0e-8) then
                    ! print *, "Converged at t = ", t
                    x_n_1 = x_n
                    exit
                endif

                ! solve the system for dx_n
                ! first factorize the matrix
                call dgeco(J_n, 2, 2, ipvt, rcond, wz)

                ! check the conditioning of the matrix
                if (1.0d0 + rcond .eq. 1.0d0) then
                    stop "Matrix is badly conditioned."
                endif

                ! solve the system
                call dgesl(J_n, 2, 2, ipvt, dx_n, 0)

                ! update the iterate
                x_n = x_n - dx_n
            enddo

            ! update the new value in the solution trajectory
            yout(idx, :) = x_n

            ! Perform quadrature along
            if (mode == "FWD") then
                ! Iterate the forward quadrature value
                call forward_ode_quad(t, x_n, "ITER")
            elseif (mode == "ADJ") then
                ! Iterate the backward quadrature value
                call adjoint_ode_quad(t, x_n, "ITER")
            endif
        enddo

        ! Finalize the quadrature term.
        if (mode == "FWD") then
            ! Finalize the forward quadrature value
            call forward_ode_quad(t, x_n, "DONE")
        elseif (mode == "ADJ") then
            ! Finalize the backward quadrature value
            call adjoint_ode_quad(t, x_n, "DONE")
        endif

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
        use power_grid_constants_and_shared
        implicit none

        double precision, intent(in) :: t
        double precision, dimension(2), intent(in) :: x
        double precision, dimension(2), intent(out) :: y
        double precision :: phi
        double precision :: omega
        double precision :: pmax_

        phi = x(1)
        omega = x(2)

        if (t > tf .AND. t <= tcl) then
            pmax_ = 0
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
        use power_grid_constants_and_shared
        implicit none

        double precision, intent(in) :: t
        double precision, dimension(2), intent(in) :: x
        double precision, dimension(2,2), intent(out) :: J
        double precision :: phi
        double precision :: omega
        double precision :: pmax_

        phi = x(1)
        omega = x(2)

        if (t > tf .AND. t <= tcl) then
            pmax_ = 0
        else
            pmax_ = pmax
        endif

        J(1, 1) = 0.0d0
        J(1, 2) = omegaB
        J(2, 1) = -omegaS * pmax_ * cos(phi) / (2 * H)
        J(2, 2) = -D * omegaS / (2 * H)
    end subroutine forward_ode_jac

    !
    ! The subroutine to compute the quadrature in FWD mode
    !
    subroutine forward_ode_quad(t, x, state)
        use power_grid_constants_and_shared
        implicit none

        character(len = 4), intent(in) :: state
        double precision, intent(in) :: t
        double precision, dimension(2), intent(in) :: x
        double precision :: phi, omega, new_h
        double precision,save :: prev_h

        phi = x(1)
        omega = x(2)

        if (state == "ITER") then
            ! perform trapezoidal integration using previous height
            ! and new height.
            new_h =  c * max(0.0d0, (phi - phiS))**beta
            quad_f = quad_f + 0.5 * dt * (prev_h + &
                       new_h)
            ! update the prev height.
            prev_h = new_h
        elseif (state == "INIT") then
            quad_f = 0
            prev_h = c * max(0.0d0, (phi - phiS))**beta
        elseif (state == "DONE") then
            print *, "The forward quadrature value is ", quad_f
        endif
    end subroutine forward_ode_quad

    !
    ! The RHS function of the adjoint ODE
    ! Uses the trajectory of the forward solution to
    ! compute adjoint RHS function
    !
    subroutine adjoint_ode_rhs(t, x, tout_f, yout_f, y)
        use power_grid_constants_and_shared
        implicit none

        integer :: idx
        double precision :: phi
        double precision, intent(in) :: t
        double precision, dimension(2,2) :: J
        double precision, dimension(2) :: forcing
        double precision, dimension(2), intent(in) :: x
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:, :), intent(in) :: yout_f
        double precision, dimension(2), intent(out) :: y

        ! find the index of the current t
        idx = findindex(t, tout_f)

        ! get the jacobian of the forward ode
        call forward_ode_jac(t, yout_f(idx, :), J)

        ! get phi at the given time instant
        phi = yout_f(idx, 1)

        ! compute the forcing term
        forcing = (/ c * beta * max(0.0d0, (phi - phiS))**(beta - 1), 0.0d0 /)

        ! now compute y
        y = matmul(transpose(-J), x) - forcing
    end subroutine

    !
    ! The jacobian of the adjoint ODE
    ! Uses the trajectory of the forward solution to
    ! compute the jacobian at a given `t`
    !
    subroutine adjoint_ode_jac(t, tout_f, yout_f, J)
        use power_grid_constants_and_shared
        implicit none

        integer :: idx
        double precision, intent(in) :: t
        double precision, dimension(2,2), intent(out) :: J
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:, :), intent(in) :: yout_f

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
        use power_grid_constants_and_shared
        implicit none

        character(len = 4), intent(in) :: state
        double precision, intent(in) :: t
        double precision, dimension(2), intent(in) :: x
        double precision :: phi, omega, new_h
        double precision,save :: prev_h

        phi = x(1)
        omega = x(2)

        if (state == "ITER") then
            ! perform trapezoidal integration using previous height
            ! and new height.
            new_h = dot_product((/0.0d0, omegaS/(2*H)/), x)
            quad_b = quad_b + 0.5 * dt * (prev_h + &
                       new_h)
            ! update the previous height.
            prev_h = new_h
        elseif (state == "INIT") then
            quad_b = 0
            prev_h = dot_product((/0.0d0, omegaS/(2*H)/), x)
        elseif (state == "DONE") then
            print *, "The adjoint quadrature value is ", quad_b
        endif
    end subroutine adjoint_ode_quad

    !
    ! Function returns the first index when traversing
    ! from lowest index (1) to highest index (length)
    ! whose value matches the desired value
    !
    integer function findindex(t, array)
        implicit none

        double precision, intent(in) :: t
        double precision, dimension(:), intent(in) :: array

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
end program power_grid
