module grid
    integer :: Nx_, Ny_, Nz_, N_
    double precision :: hx_, hy_, hz_, V_
    double precision, dimension(:,:,:,:), pointer :: K_  ! Permeabilities

    parameter(Nx_ = 32, &                  ! Dimension in x-direction
              Ny_ = 32, &                  ! Dimension in y-direction
              Nz_ = 1, &                  ! Dimension in z-direction
              hx_ = 1.0d0/Nx_, &    ! step size in x-direction
              hy_ = 1.0d0/Ny_, &    ! step size in y-direction
              hz_ = 1.0d0/Nz_, &    ! step size in z-direction
              N_ = Nx_ * Ny_ * Nz_, &       ! Total number of grid cells
              V_ = hx_ * hy_ * hz_)         ! Volume of each grid cell

end module grid
