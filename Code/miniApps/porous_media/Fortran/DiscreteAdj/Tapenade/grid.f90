module grid
    integer :: Nx_, Ny_, Nz_, N_
    double precision :: hx_, hy_, hz_, V_
    integer :: maxNx, maxNy, maxNz
    parameter(maxNx=60, maxNy=220, maxNz=85)

    parameter(Nx_ = 4, &                   ! Dimension in x-direction
              Ny_ = 4, &                  ! Dimension in y-direction
              Nz_ = 2, &                    ! Dimension in z-direction
              hx_ = 20.0d0 * 0.3048d0, &    ! step size in x-direction
              hy_ = 10.0d0 * 0.3048d0, &    ! step size in y-direction
              hz_ = 2.0d0 * 0.3048d0, &     ! step size in z-direction
              N_ = Nx_ * Ny_ * Nz_, &       ! Total number of grid cells
              V_ = hx_ * hy_ * hz_)         ! Volume of each grid cell

    double precision, dimension(N_) :: Por_      ! Porosities
    double precision, dimension(3, Nx_, Ny_, Nz_) :: K_  ! Permeabilities 
    
    integer, dimension(7) :: diags_
    parameter(diags_ = (/ -Nx_ * Ny_, -Nx_, -1, 0, 1, Nx_, Nx_ * Ny_ /))
end module grid
