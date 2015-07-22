module fluid
    double precision :: vw, vo, swc, sor

    parameter(vw = 3d-4, &      ! Viscosity
              vo = 3d-3, &      ! Viscosity
              swc = 0.2d0, &    ! Saturation
              sor = 0.2d0)      ! Saturation

end module fluid
