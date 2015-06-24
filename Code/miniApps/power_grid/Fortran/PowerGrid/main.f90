#define fwd_ode_rhs_fn()
#define adj_ode_rhs_fn()
#define fwd_ode_rhs_jac()
#define adj_ode_rhs_jac()
#define fwd_ode_op_fn()
#define adj_ode_op_fn()
#define time_integrator(x, to, yo) crank_nicolson(x, to, yo)
program power_grid
    use power_grid_constants
    implicit none

    ! -------------------------------------------------------------------------
    ! declare variables for the main optimization routine
    ! -------------------------------------------------------------------------
    ! we are using L-BFGS-B(v.3) from Nocedal et al.

    ! n is the number of parameters - pm
    ! m is the maximum number of limited memory corrections.
    integer :: np, nx, m, maxiter, tlen, alloc_state
    parameter (np = 1, nx = 2, m = 5, maxiter = 5)

    ! Declare the variables needed by the code.
    ! Some of these variables are from LBFGS-B,
    ! rest are either passed between methods or
    ! shared by methods.
    ! Shared (globally accessible) : pm
    ! Passed: x, tout_f, tout_b, yout_f, yout_b
    character*60     :: task, csave
    logical          :: lsave(4)
    integer          :: i, iprint, nbd(np), iwa(3*np), isave(44)
    double precision :: pm, f, factr, pgtol, quad_f, quad_b, &
                        l(np), u(np), g(np), dsave(29), &
                        wa(2*m*np + 5*np + 11*m*m + 8*m)

    ! -------------------------------------------------------------------------
    ! declare variables for the integration scheme - use fixed time stepping
    ! -------------------------------------------------------------------------
    double precision, dimension(nx) :: x
    double precision, dimension(:), allocatable :: tout_f
    double precision, dimension(:), allocatable :: tout_b
    double precision, dimension(:,:), allocatable :: yout_f
    double precision, dimension(:,:), allocatable :: yout_b

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
    factr = 1.0d+7
    pgtol = 1.0d-5

    ! Set the bound on the parameter p
    nbd(1) = 3
    u(1) = 1.1

    ! -------------------------------------------------------------------------
    ! initialize variables/constants for time integration scheme - fixed t.step
    ! -------------------------------------------------------------------------
    tlen = int((tend - t0)/dt) + 1            ! Verify that the tlen is correct
                                              ! Number of endpoints of intervals

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
    task = 'START'

    ! Iterate maximum number of iterations calling the optimization routine.
    do i = 1, maxiter

        ! Initialize the x value
        x = (/asin(perturb*pm/pmax), perturb*1.0/)

        call setulb(np, m, pm, l, u, nbd, f, g, factr, pgtol, wa, iwa, task, &
                    iprint, csave,lsave,isave,dsave)

        if (task(1:2) == 'FG') then
            ! Compute the function value.
            call power_grid_cost_function(x, tout_f, yout_f, f)

            ! Compute the gradient value.
#ifdef USE_OPENAD
            ! call the OPENAD generated discrete adjoint to get the
            ! gradient used for optimization.
#else
            ! call the function that computes the sensitivites using
            ! continuous adjoints.
#endif

        elseif (task(1:5) == 'NEW_X') then
            if (i < maxiter) then
                ! continue the iteration at the new point

            else
                ! stop the iteration. print current values

            endif
        elseif (task(1:4) == 'CONV') then
            ! stop the iteration. print current values

        elseif (task(1:4) == 'ABNO') then
            ! stop the iteration. print the current values
            ! also print the reason for abnormal termination

        elseif (task(1:5) == 'ERROR') then
            ! stop the iteration. print the current values
            ! also print the reason for abnormal termination

        else

        endif
    enddo
contains

    !
    ! The subroutine computes power grid cost function
    ! This subroutine should call the time integrator
    ! in the forward mode along with the mode being
    ! set to 'FWD'.
    !
    subroutine power_grid_cost_function(x, tout_f, yout_f, f)
        use power_grid_constants
        implicit none

        double precision, dimension(2), intent(in) :: x
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:,:), intent(out) :: yout_f
        double precision, intent(out) :: f

        call crank_nicolson(x, tout_f, yout_f, 'FWD')

        f = -pm + quad_f
    end subroutine power_grid_cost_function

    !
    ! The subroutine computes the sensitivity of cost
    ! function with respect to the parameter p
    !
    subroutine power_grid_cost_gradient()
        use power_grid_constants
        implicit none
        ! The gradient can be computed using continuous adjoints
        ! or by discrete adjoints

    end subroutine power_grid_cost_gradient

    !
    ! Perform time integration by Crank Nicolson scheme.
    !
    subroutine crank_nicolson(x0, tout, yout, mode)
        use power_grid_constants
        implicit none

        double precision, dimension(:), intent(in) :: x0
        double precision, dimension(:), intent(in) :: tout
        double precision, dimension(:, :), intent(out) :: yout
        character(len=3) :: mode

        if (mode == 'FWD') then
            ! Call ForwardODERHS, ForwardODEJac, ForwardODEQuad
        elseif (mode == 'ADJ') then
            ! Call AdjointODERHS, AdjointODEJac, AdjointODEQuad
        endif
    end subroutine crank_nicolson

    !
    ! The RHS function of the forward ODE
    !
    subroutine forward_ode_rhs(t, x, y)
        use power_grid_constants
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
        use power_grid_constants
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
        use power_grid_constants
        implicit none

        double precision, intent(in) :: t
        double precision, dimension(2), intent(in) :: x
        character(len = 4), intent(in) :: state
        double precision :: phi
        double precision :: omega

        double precision,save :: d_quad_f

        phi = x(1)
        omega = x(2)

        if (state == "ITER") then
            ! perform trapezoidal integration using previous height
            ! and current height.
            quad_f = quad_f + 0.5 * dt * (d_quad_f + &
                       c * max(0.0d0, (phi - phiS))**beta)
            ! update the current height.
            d_quad_f = c * max(0.0d0, (phi - phiS))**beta
        elseif (state == "INIT") then
            quad_f = 0
            d_quad_f = c * max(0.0d0, (phi - phiS))**beta
        elseif (state == "DONE") then
            print *, "The forward quadrature value is ", quad_f
        endif
    end subroutine forward_ode_quad

    !
    ! The RHS function of the adjoint ODE
    !
    subroutine adjoint_ode_rhs(t, x, tout_f, yout_f, y)
        use power_grid_constants
        implicit none

        integer :: i
        double precision :: phi
        double precision, intent(in) :: t
        double precision, dimension(2,2) :: J
        double precision, dimension(2) :: forcing
        double precision, dimension(2), intent(in) :: x
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:, :), intent(out) :: yout_f
        double precision, dimension(2), intent(out) :: y

        ! find the index of the current t
        i = findindex(t, tout_f)

        ! get the jacobian of the forward ode
        call forward_ode_jac(t, yout_f(i, :), J)

        ! get phi at the given time instant
        phi = yout_f(i, 1)


        ! compute the forcing term
        forcing = (/ c * beta * max(0.0d0, (phi - phiS))**(beta - 1), 0.0d0 /)

        ! now compute y
        y = matmul(transpose(-J), x) - forcing
    end subroutine

    !
    ! The jacobian of the adjoint ODE
    !
    subroutine adjoint_ode_jac(t, x, tout_f, yout_f, J)
        use power_grid_constants
        implicit none

        integer :: i
        double precision, intent(in) :: t
        double precision, dimension(2), intent(in) :: x
        double precision, dimension(2,2), intent(out) :: J
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:, :), intent(out) :: yout_f

        ! find the index of the current t
        i = findindex(t, tout_f)

        ! get the jacobian of the forward ode
        call forward_ode_jac(t, yout_f(i, :), J)

        J = transpose(-J)
    end subroutine

    !
    ! Function returns the first index when traversing
    ! from lowest index (1) to highest index (length)
    ! whose value matches the desired value
    !
    integer function findindex(t, array)
        implicit none

        double precision, intent(in) :: t
        double precision, dimension(:), intent(in) :: array

        integer :: i, default_index

        default_index = 0
        do i = 1, size(array, 1)
            if(array(i) == t) then
                default_index = i
                exit
            endif
        enddo

        findindex = default_index
    end function findindex
end program power_grid
