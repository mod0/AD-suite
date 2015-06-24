module power_grid_constants
    implicit none
    save

    double precision :: pi, pmax, omegaB, omegaS, H, D, phiS, beta, c, t0, &
                        tend, tf, tcl, dt, perturb

    ! Set the parameter values
    parameter(                      &
        pi = 3.14159265358979323,   &
        pmax = 2.0877,              &
        omegaB = 120 * pi,          &
        omegaS = 1,                 &
        H = 5.0,                    &
        D = 5.0,                    &
        phiS = 1.0,                 &
        beta = 2,                   &
        c = 10000,                  &
        t0 = 0,                     &
        tend = 10.0,                &
        tf = 0.1,                   &
        tcl = 0.2,                  &
        dt = 0.01,                  &
        perturb = 1.0)

end module power_grid_constants
