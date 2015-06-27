module power_grid_constants_and_shared
    implicit none
    save

    ! Constants
    double precision :: pi, pmax, omegaB, omegaS, H, D, phiS, beta, c, t0, &
                        tend, tf, tcl, dt, perturb

    ! Shared
    double precision :: pm, quad_f, quad_b

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
        tf = 0.1d0,                     &
        tcl = 0.2d0,                    &
        dt = 0.01d0,                    &
        perturb = 1.0d0)

end module power_grid_constants_and_shared
