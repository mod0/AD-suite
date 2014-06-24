!-----------------------------
        module stream_vel_variables
!-----------------------------

        integer, parameter :: n = 79                ! grid size
        integer, parameter :: n_nl = 20            ! # nonlin iterations (mat solves)
        real(8), parameter :: Lx = 79.e3            ! domain length
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
! stub for module LU

module conj_gradStub_mod

  public :: solve

  interface solve
     module procedure conj_gradStub
  end interface

contains 
      subroutine conj_gradStub (x, b, A)
    !$openad xxx template oad_template_conj_grad.f90
!-----------------------------

!      KRISHNA: THIS IS THE LINEAR SYSTEM SOLVER
        use stream_vel_variables
        real(8), intent(inout), dimension(n) :: x
        real(8), intent(in), dimension(n) :: b
        real(8), intent(in), dimension(n,3) :: A
        integer :: i
        integer :: j
        do i=1,n
         x(i) = A(i,1)+A(i,2)+A(i,3)+b(i)
        enddo

      end subroutine conj_gradStub

end module


!$TAF SUBROUTINE CONJ_GRAD ADNAME  = adCONJ_GRAD
!$TAF SUBROUTINE CONJ_GRAD INPUT   = 1,2,3
!$TAF SUBROUTINE CONJ_GRAD OUTPUT  = 1
!$TAF SUBROUTINE CONJ_GRAD ACTIVE  = 1,2,3
!$TAF SUBROUTINE CONJ_GRAD DEPEND  = 1,2,3




!-----------------------------
        subroutine stream_vel_init (h, beta)
!-----------------------------

        use stream_vel_variables

        real(8), intent(inout), dimension(n) :: h, beta

        integer :: i

        dx = Lx / real(n)

        do i=1,n
          beta (i) = beta_const
          h(i) = h_left + (h_right-h_left)/Lx * (i-0.5)*dx
        enddo

        end subroutine stream_vel_init


!-----------------------------
        subroutine stream_vel (u, bb, fc)
!-----------------------------



        use stream_vel_variables
        use conj_gradStub_mod

        real(8), intent(inout), dimension(n+1) :: u
        real(8), intent(inout), dimension(n) :: bb
        real(8), intent(inout) :: fc

        real(8), dimension(n,3) :: A
        real(8), dimension(n) :: f, beta_fric, h, beta_0, h0, utmp, b, nu
        real(8) :: fend
        integer :: i,j

!$openad INDEPENDENT(u)
!$openad INDEPENDENT(bb)
!$openad DEPENDENT(fc)

!$TAF INIT tape_inner = static, n_nl
        
        call stream_vel_init (h0, beta_0)
        beta_fric = beta_0
        h=h0+bb
   
!         print *, bb(n)

        call stream_vel_taud (h, f, fend)

        u = 0.

!------ driving stress -------------

        do i=1,n
         b(i) = -dx * f(i)
         if (i.lt.n) then
          b(i) = b(i) - f(i+1) * dx
         endif
        enddo
        
        b(n) = b(n) + fend

!-----------------------------------
        do i=1, n_nl

!$TAF STORE u = tape_inner

         call stream_vel_visc (h, u, nu)                 ! update viscosities
              
         call stream_assemble (nu, beta_fric, A)         ! assemble tridiag matrix
                                                         ! this represents discretization of
                                                         !  (nu^(i-1) u^(i)_x)_x - \beta^2 u^(i) = f

         utmp = 0.
         call solve (utmp, b, A)                     ! solve linear system for new u

         do j=1,n
          u(j+1) = utmp(j)                               ! effectively apply boundary condition u(1)==0
         enddo


 
        enddo

        fc=0.
        do i=2,n+1
         fc = fc + u(i) * u(i)
        enddo

        end subroutine stream_vel

!-----------------------------
        subroutine stream_vel_taud (h, f, fend)
!-----------------------------

        use stream_vel_variables

        real(8), intent(in), dimension(n) :: h
        real(8), intent(inout), dimension(n) :: f
        real(8), intent(inout) :: fend

        integer :: i

        do i=1,n
         if ((i.gt.1).and.(i.lt.n)) then
          f (i) = rhoi * g * h(i) * (h(i+1)-h(i-1))/2./dx
         elseif (i.eq.1) then
          f (i) = rhoi * g * h(i) * (h(i+1)-h(i))/dx
         elseif (i.eq.n) then
          f (i) = rhoi * g * h(i) * (h(i)-h(i-1))/dx
         endif
        enddo


        fend = .5 * (rhoi * g * (h(n))**2 - rhow * g * R_bed**2)

        


        end subroutine stream_vel_taud

!-----------------------------
        subroutine stream_vel_visc (h, u, nu)
!-----------------------------

        use stream_vel_variables

        real(8), intent(in), dimension(n) :: h
        real(8), intent(in), dimension(n+1) :: u
        real(8), intent(inout), dimension(n) :: nu

        real(8) :: ux, tmp
        integer :: i

        do i=1,n
         ux = (u(i+1)-u(i)) / dx
         tmp = ux**2 + ep_glen**2
         nu(i) = .5 * h(i) * Aglen**(-1./nglen) * tmp ** ((1-nglen)/2./nglen)
        enddo

        end subroutine stream_vel_visc
 
!-----------------------------
        subroutine stream_assemble (nu, beta_fric, A)
!-----------------------------

        use stream_vel_variables

        real(8), intent(in), dimension(n) :: nu, beta_fric
        real(8), intent(inout), dimension(n,3) :: A
        integer :: i

        do i = 1,n
         A(i,2) = 4*nu(i)/dx + dx/3. * beta_fric(i)**2
         if (i.gt.1) then
          A(i,1) = -4*nu(i)/dx + dx/6. * beta_fric(i)**2
!           print *, A(i,1)
         endif
         if (i.lt.n) then
          A(i,2) = A(i,2) + 4*nu(i+1)/dx + dx/3. * beta_fric(i+1)**2
          A(i,3) = -4*nu(i+1)/dx + dx/6. * beta_fric(i+1)**2
         endif
        enddo
                
        end subroutine stream_assemble


