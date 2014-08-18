!$TAF SUBROUTINE CONJ_GRAD ADNAME  = adCONJ_GRAD
!$TAF SUBROUTINE CONJ_GRAD INPUT   = 1,2,3
!$TAF SUBROUTINE CONJ_GRAD OUTPUT  = 1
!$TAF SUBROUTINE CONJ_GRAD ACTIVE  = 1,2,3
!$TAF SUBROUTINE CONJ_GRAD DEPEND  = 1,2,3


!-----------------------------
        subroutine stream_vel_timedep (h, u, bb, fc)
!-----------------------------
        use stream_vel_variables

        real(8), intent(inout), dimension(n) :: h, bb
        real(8), intent(out), dimension(n+1) :: u
        real(8), intent(inout) :: fc

        real(8), dimension(n) :: h0, beta_0, beta_fric
        integer :: nt,i

!$openad INDEPENDENT(u)
!$openad INDEPENDENT(bb)
!$openad DEPENDENT(fc)
        call stream_vel_init (h0, beta_0)
        beta_fric = beta_0
        h=h0+bb
        dx = Lx / real(n)

        call stream_vel (u,h,beta_fric)

        call stream_vel_timedep_stage(h, u, beta_fric)

        fc=0.
        do i=2,n+1
         fc = fc + u(i) * u(i)
        enddo

        end subroutine stream_vel_timedep

!-----------------------------
        subroutine stream_vel_timedep_stage (h, u, beta_fric)
!-----------------------------
!$openad xxx template ad_revolve_template.F
        use stream_vel_variables
        real(8), intent(inout), dimension(n) :: h
        real(8), intent(inout), dimension(n+1) :: u
        real(8), intent(inout), dimension(n) :: beta_fric
        integer :: nt

        do nt=1,n_timesteps
         call OpenAD_forward_step (h, u, beta_fric)
        enddo


        end subroutine stream_vel_timedep_stage

!-----------------------------
subroutine OpenAD_forward_step (h, u, beta_fric)
!-----------------------------
        use stream_vel_variables
        real(8), intent(inout), dimension(n) :: h
        real(8), intent(inout), dimension(n+1) :: u
        real(8), intent(inout), dimension(n) :: beta_fric
        integer :: i
         
         h(1) = h(1) - dt/dx * u(2)*h(2)
         do i=2,n

          h(i) = h(i) - dt/dx * (u(i+1)*h(i) - u(i)*h(i-1))

         enddo

         call stream_vel (u,h,beta_fric)
end subroutine OpenAD_forward_step

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
        subroutine stream_vel (u, h, beta_fric)
!-----------------------------



        use stream_vel_variables
        use conj_grad_mod

        real(8), intent(inout), dimension(n+1) :: u
        real(8), intent(in), dimension(n) :: h
        real(8), intent(inout), dimension(n) :: beta_fric
!        real(8), intent(inout) :: fc

!        real(8), dimension(n,3) :: A
        real(8), dimension(n) :: f, utmp, b
        real(8), dimension(n+1) :: unew
        real(8) :: fend
        integer :: i,j

!$TAF INIT tape_inner = static, n_nl
        isinloop0 = 0
        isinloop1 = 1
        isinloop2 = 2
 
!        call stream_vel_init (h0, beta_0)
!        beta_fric = beta_0
!        h=h0+bb
   

        call stream_vel_taud (h, f, fend)


        u = 0.

!------ driving stress -------------

        b=0.0
        do i=1,n
         b(i) = -dx * f(i)
         if (i.lt.n) then
          b(i) = b(i) - f(i+1) * dx
         endif
        enddo
        
        b(n) = b(n) + fend

!-----------------------------------

        call phistage (u, unew, b, h, beta_fric, isinloop0)
       ! u = unew

        do i=1, n_nl

!$TAF STORE u = tape_inner
          call phistage (u, unew, b, h, beta_fric, isinloop1)
          u = unew

        enddo

        call phistage (u, unew, b, h, beta_fric, isinloop2)
        u = unew

        end subroutine stream_vel

!-----------------------------
        subroutine phi (u_i, u_ip1, b, h, beta_fric)
!-----------------------------
        use stream_vel_variables
        use conj_grad_mod
        real(8), dimension(n) :: b, h
        real(8), intent(inout), dimension(n) :: beta_fric
        real(8), intent(in), dimension(n+1) :: u_i
        real(8), intent(out), dimension(n+1) :: u_ip1
        
        real(8), dimension(n) :: nu,utmp
        real(8), dimension(n,3) :: A
        integer :: j


        call stream_vel_visc (h, u_i, nu)                 ! update viscosities
              
        call stream_assemble (nu, beta_fric, A)         ! assemble tridiag matrix
                                                         ! this represents discretization of
                                                         !  (nu^(i-1) u^(i)_x)_x - \beta^2 u^(i) = f

        utmp = 0.
        call solve (utmp, b, A)                     ! solve linear system for new u

        do j=1,n
         u_ip1(j+1) = utmp(j)                               ! effectively apply boundary condition u(1)==0
        enddo


        end subroutine phi 

!-----------------------------
        subroutine phistage (u, u_ip1, b, h, beta_fric, isinloop)
!-----------------------------
!$openad xxx template oad_template_fixedpointstage_withconvergence.f90
        use stream_vel_variables
        real(8), dimension(n) :: b, h
        real(8), intent(inout), dimension(n) :: beta_fric
        real(8), intent(in), dimension(n+1) :: u
        real(8), intent(out), dimension(n+1) :: u_ip1
!        real(8), dimension(n+1) :: u_ip1
        integer, intent(in) :: isinloop 
        integer, save :: iter=0, adj_iter=0
        logical, save :: conv_flag=.FALSE., adj_conv_flag=.FALSE.
        integer :: k
        real(8) :: normdiff, normZ, diff

        if (conv_flag .eqv. .false.) then

          iter = iter + 1

          call phi (u, u_ip1, b, h, beta_fric)

         
        endif
        normdiff=0.
        normZ=0.
        diff = 0.
        k=0
        

        end subroutine phistage 

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
         nu(i) =.5*h(i)*Aglen**(-1./nglen) * tmp ** ((1-nglen)/2./nglen)
         !print *, nu(i)
        enddo

        end subroutine stream_vel_visc
 
!-----------------------------
        subroutine stream_assemble (nu, beta_fric, A)
!-----------------------------

        use stream_vel_variables

        real(8), intent(in), dimension(n) :: nu
        real(8), intent(inout), dimension(n) :: beta_fric
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


