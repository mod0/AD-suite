module grid
    integer :: Nx_, Ny_, Nz_
    double precision :: hx_, hy_, hz_, V_
    double precision, dimension(:,:,:,:), pointer :: K_  ! Permeabilities

    parameter(Nx_ = 60, &                   ! Dimension in x-direction
              Ny_ = 220, &                  ! Dimension in y-direction
              Nz_ = 1, &                    ! Dimension in z-direction
              hx_ = 20.0d0 * 0.3048d0, &    ! step size in x-direction
              hy_ = 10.0d0 * 0.3048d0, &    ! step size in y-direction
              hz_ = 2.0d0 * 0.3048d0, &     ! step size in z-direction
              N_ = Nx_ * Ny_ * Nz_, &       ! Total number of grid cells
              V_ = hx_ * hy_ * hz_)         ! Volume of each grid cell

end module grid