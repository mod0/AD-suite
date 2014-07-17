!-----------------------------
        module stream_vel_variables
!-----------------------------

        integer, parameter :: n = 79                ! grid size
        integer, parameter :: n_nl = 60            ! # nonlin iterations (mat solves)
        integer, parameter :: n_timesteps = 200         
        real(8), parameter :: Lx = 79.e3            ! domain length
        real(8), parameter :: dt = 0.04            ! length of timestep (years)
        real(8), parameter :: tol = 1.e-7    !tolerance of fwd loop
        real(8), parameter :: adjtol = 1.e-7 !tolerance of adjoint loop
        real(8), parameter :: ep_glen = 1.e-7       ! reg. parameter for glen's law
        real(8), parameter :: eps = 1.e-5           
        real(8), parameter :: Aglen = 5.0002e-17           ! glen's law constant
        real(8), parameter :: nglen =3.            ! glen's law exponent
        real(8), parameter :: g = 9.81              ! 
        real(8), parameter :: rhoi = 910.
        real(8), parameter :: rhow = 1035.
        real(8), parameter :: R_bed = -900.         
        real(8), parameter :: beta_const = 5.   ! \beta for sliding law
        real(8), parameter :: h_left = 1050.        
        real(8), parameter :: h_right = 1050.
       
        
        real(8), parameter :: PI = 3.14159265358979323844d0
        real(8), dimension (n) :: r, r_old, p, p_old, ax, x_old
        real(8), dimension (n,3) :: tridiag_0

        real(8) :: dx

        end module stream_vel_variables
