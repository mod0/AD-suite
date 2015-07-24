module fluid
    double precision :: vw_, vo_, swc_, sor_

    parameter(vw_ = 3d-4, &      ! Viscosity
              vo_ = 3d-3, &      ! Viscosity
              swc_ = 0.2d0, &    ! Saturation
              sor_ = 0.2d0)      ! Saturation

end module fluid
