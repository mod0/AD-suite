program power_grid
    use print_active
    use gnufor2
    use numerics
#ifdef USE_OPENAD
    use OAD_active
    use w2f__types
    use OAD_rev
    use OAD_tape
#endif
    implicit none

    ! -------------------------------------------------------------------------
    ! declare variables for the main optimization routine
    ! -------------------------------------------------------------------------
    ! we are using L-BFGS-B(v.3) from Nocedal et al.

    ! n is the number of parameters - pm
    ! m is the maximum number of limited memory corrections.
    integer :: np, m, max_opt_iter, tlen
    parameter (np = 1, m = 10, max_opt_iter = 50)

    ! Declare the variables needed by the code.
    ! Some of these variables are from LBFGS-B,
    ! rest are either passed between methods or
    ! shared by methods.
    ! Shared (globally accessible) : pm
    ! Passed: x, tout_f, tout_b, yout_f, yout_b
    character*60     :: task, csave
    character*10     :: filename
    logical          :: lsave(4)
    integer          :: i, iprint, nbd(np), iwa(3*np), isave(44)
    double precision :: f, factr, pgtol, &
                        l(np), u(np), g(np), dsave(29), &
                        wa(2*m*np + 5*np + 11*m*m + 8*m)

#ifdef USE_OPENAD
    type(active) :: f_OAD
#endif

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
    ! Check the documentation of LBFGS for more information
    ! about these parameters.
    factr = 1.0d+3              ! factr higher -> looser, lower -> stricter
    pgtol = 0.0d0               ! Projected gradient tolerance

    ! Set the bound on the parameter p
    nbd(1) = 3                  ! only upper-bound on the variable
    u(1) = 1.1d0                ! the value of the upper bound.


    ! -------------------------------------------------------------------------
    ! initialize variables/constants for time integration scheme - fixed t.step
    ! -------------------------------------------------------------------------
    tlen = int((tend - t0)/dt) + 1            ! Verify that the tlen is correct
                                              ! Number of endpoints of intervals

    ! Set the initial value of the parameter pm
#if !defined(USE_OPENAD)
    pm = 0.4d0
#else
    pm%v = 0.4d0
#endif

    ! Set the task variable
    task = "START"

    ! Iterate maximum number of iterations calling the optimization routine.
    do i = 1, max_opt_iter
#if !defined(USE_OPENAD)
        call setulb(np, m, pm, l, u, nbd, f, g, factr, pgtol, wa, iwa, task, &
                    iprint, csave,lsave,isave,dsave)
#else
        call setulb(np, m, pm, l, u, nbd, f_OAD%v, pm%d, factr, pgtol, wa, iwa, task, &
                    iprint, csave,lsave,isave,dsave)
#endif
        if (task(1:2) == "FG") then
#if !defined(USE_OPENAD)
            call get_cost_function_and_gradient(f, g, tlen)
#else
            pm%d = 0.0d0
            f_OAD%d = 1.0d0

            !Call function in Tape mode
            our_rev_mode%arg_store=.FALSE.
            our_rev_mode%arg_restore=.FALSE.
            our_rev_mode%res_store=.FALSE.
            our_rev_mode%res_restore=.FALSE.
            our_rev_mode%plain=.FALSE.
            our_rev_mode%tape=.TRUE.
            our_rev_mode%adjoint=.FALSE.
            call get_cost_function_and_gradient(f_OAD, g, tlen)
            !Call function in Adjoint mode
            our_rev_mode%arg_store=.FALSE.
            our_rev_mode%arg_restore=.FALSE.
            our_rev_mode%res_store=.FALSE.
            our_rev_mode%res_restore=.FALSE.
            our_rev_mode%plain=.FALSE.
            our_rev_mode%tape=.FALSE.
            our_rev_mode%adjoint=.TRUE.
            call get_cost_function_and_gradient(f_OAD, g, tlen)
#endif
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
    enddo

    ! stop the iteration. print current values
    print *, "Maximum number of iterations reached. The value of pm is ", pm
end program power_grid
