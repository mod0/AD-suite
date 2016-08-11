module simulation
  use parameters
  use matrix
  use finitevolume
  implicit none
contains

  !
  ! Initialize inflow and outflow.
  !
  subroutine init_flw_trnc_norm_xin_pt_out(nx, ny, nz, ndof, mu, sigma, Q)
    integer :: nx, ny, nz, ndof
    double precision, dimension( : ) :: mu, sigma
    double precision, dimension((nx * ny * nz)) :: Q

    integer :: i, j
    double precision :: x, pi, pdf, mass, arg1
    double precision, dimension(nx) :: idx
    double precision, dimension(nx) :: Q_x

    ! value of pi
    pi = 3.14159265358979323d0

    !initialize the total mass to 0
    mass = 0.0d0
    Q_x = 0.0d0

    ! Note that the portion of the  Standard Normal distribution between
    ! -3sigma/2 to 3sigma/2 is assumed to fit the 1..Nx where sigma is 1
    do i = 1, nx
      x   = (i-1.0)*2.0/(nx-1.0d0) - 1.0d0;
      pdf = 0.0d0
      do j = 1, ndof
          arg1 = -(((x-mu(j))**2/sigma(j))**2.0/2.0);
          pdf = pdf + 1.0/(sqrt(2.0*pi)*sigma(j))*exp(arg1);
      end do
      Q_x(i) = pdf;
      mass = mass + pdf;
    end do

    ! now rescale all the entities
    Q_x = Q_x/mass * ir

    ! Assign Q_x to Q
    j = 1
    do i = 1, nx* ny, ny
       Q(i) = Q_x(j)
       j = j + 1
    end do

    ! now set the output
    Q((nx * ny * nz)) = -ir
  end subroutine init_flw_trnc_norm_xin_pt_out


  !
  ! This subroutine simulates the reservoir
  ! model.
  !
  subroutine simulate_reservoir(nx, ny, nz, nd, pt, st, Q, S, P, V, Tt, Pc, oil)
    use parameters
    integer :: nx, ny, nz
    integer :: nd, pt, st
    double precision, dimension((nx * ny * nz)) :: Q

    double precision, dimension((nx * ny * nz)) :: S
    double precision, dimension(nx, ny, nz) :: P
    double precision, dimension(3, nx + 1, ny + 1, nz + 1) :: V

    double precision, dimension((nd/st) + 1) :: Tt   
    double precision, dimension(2, (nd/st) + 1) :: Pc
    double precision ::  oil

    integer :: i, j, k
    double precision :: Mw, Mo, Mt, tempoil1, tempoil2

    S = swc_                            ! initial saturation

    Pc(1, 1) = 0.0d0                    ! initial production
    Pc(2, 1) = 1.0d0

    Tt(1) = 0.0d0                       ! initial time.

    tempoil1 = 0.0d0
    tempoil2 = 0.0d0

    k = 1
    do i = 1, nd/pt
       do j = 1, pt/st
          k = k + 1

          if (j == 1) then
             call stepforward(nx, ny, nz, nd, pt, st, 1, Q, S, P, V, Mw, Mo)
          else
             call stepforward(nx, ny, nz, nd, pt, st, 0, Q, S, P, V, Mw, Mo)            
          endif


          ! update quantites
          Mt = Mw + Mo
          Tt(k) = 1.0d0 * k * St  
          Pc(1,k) = Mw/Mt
          Pc(2,k) = Mo/Mt

          call update_oil(nd, pt, st, Pc, k, tempoil1, tempoil2)
          tempoil1 = tempoil2
       end do
    end do

    oil = tempoil2
  end subroutine simulate_reservoir

  subroutine stepforward(nx, ny, nz, nd, pt, st, pressure_step, Q, S, P, V, Mw, Mo)
    integer :: nx, ny, nz
    integer :: nd, pt, st
    integer :: pressure_step
    double precision, dimension((nx * ny * nz)) :: Q
    double precision, dimension((nx * ny * nz)) :: S
    double precision, dimension(nx, ny, nz) :: P
    double precision, dimension(3, nx + 1, ny + 1, nz + 1) :: V
    double precision :: Mw, Mo

    if (pressure_step == 1) then
       ! solve pressure
       call Pres(nx, ny, nz, Q, S, P, V)    ! Pressure solver
    endif

    call NewtRaph(nx, ny, nz, nd, pt, st, Q, V, S)      ! Solve for saturation
    call RelPerm(S((nx * ny * nz)), Mw, Mo)         ! Mobilities in well-block
  end subroutine stepforward

  subroutine update_oil(nd, pt, st, Pc, k, oilin, oilout)
    integer :: k
    integer :: nd, pt, st
    double precision ::  oilin
    double precision ::  oilout
    double precision, dimension(2, (nd/st) + 1) :: Pc

    oilout = oilin +  Pc(2, k) * st                     ! Reimann sum
  end subroutine update_oil

  subroutine wrapper(nx, ny, nz, nd, ndof, pt, st, mu, sigma, Q, S, P, V, Tt, Pc, oil) 
    integer :: nx, ny, nz, ndof
    integer :: nd, pt, st
    double precision, dimension( ndof ) :: mu
    double precision, dimension( ndof ) :: sigma
    double precision, dimension((nx * ny * nz)) :: Q
    double precision, dimension((nx * ny * nz)) :: S
    double precision, dimension(nx, ny, nz) :: P
    double precision, dimension(3, nx + 1, ny + 1, nz + 1) :: V
    double precision, dimension((nd/st) + 1) :: Tt
    double precision, dimension(2, (nd/st) + 1) :: Pc
    double precision :: oil

    call init_flw_trnc_norm_xin_pt_out(nx, ny, nz, ndof, mu, sigma, Q)
    call simulate_reservoir(nx, ny, nz, nd, pt, st, Q, S, P, V, Tt, Pc, oil)
  end subroutine wrapper
end module simulation
