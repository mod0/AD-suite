module fluid
    double precision :: vw_, vo_, swc_, sor_

    parameter(vw_ = 3d-4, &      ! Viscosity of Water
              vo_ = 3d-3, &      ! Viscosity of Oil
              swc_ = 0.2d0, &    ! Saturation of water cut
              sor_ = 0.2d0)      ! Saturation of oil cut

end module fluid
