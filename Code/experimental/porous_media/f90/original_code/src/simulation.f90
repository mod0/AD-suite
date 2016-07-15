module simulation
  use parameters
  use matrix
  use finitevolume
  implicit none
contains

  !
  ! This routine opens the permeability and porosity used by
  ! the MATLAB program and uses it for the simulation.
  !
  subroutine read_permeability_and_porosity(nx, ny, nz, maxnx, maxny, maxnz, PERM, POR)
    integer :: i, j, k, l, m
    integer :: nx, ny, nz
    integer :: maxnx, maxny, maxnz
    double precision, dimension((nx * ny * nz)) :: POR                 ! Porosities
    double precision, dimension(3, nx, ny, nz) :: PERM  ! Permeabilities

    double precision, dimension(nx, ny, nz) :: P
    double precision, dimension(maxnx * maxny * maxnz) :: pUr
    double precision, dimension(3 * maxnx, maxny * maxnz) :: KUr
    double precision, dimension(3 * maxnx * maxny * maxnz) :: KUrl

    integer, dimension(nx * ny * nz) :: Pindices
    integer, dimension(3 * nx * ny * nz) :: Kindices

    ! initialize porosity and permeability to zero
    PERM = 0.0d0
    POR = 0.0d0

    ! read KUr
    open(1,file=permeability_file,status='old')
    read(1,*) ((KUr(i,j), j=1,maxny * maxnz), i=1,3 * maxnx)
    close(1)

    ! reshape 2 dimension to 1 dimension
    call myreshape_2_1(KUr, KUrl)

    ! select according to specified dimension
    m = 0
    do l = 1, nz
       do k = 1,ny
          do j = 1,nx
             do i = 1,3
                m = m + 1
                Kindices(m) = ((l - 1) * (maxnx * maxny * 3) &
                     + (k - 1) * (maxnx * 3) &
                     + 3 * (j-1) + i)
             end do
          end do
       end do
    end do

    ! then reshape 1 dimension to 4 dimension (hack for time being)
    call myreshape_1_4(KUrl(Kindices), PERM)

    ! read KUr
    open(1,file=porosity_file,status='old')
    read(1,*) (pUr(i), i=1,maxnx * maxny * maxnz)
    close(1)

    m = 0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             m = m + 1
             Pindices(m) = ((k - 1) * (maxnx * maxny) &
                  + (j - 1) * (maxnx) + i)
          end do
       end do
    end do

    call mymax_1_0_double(pUr(Pindices), 1.0d-3, POR)
  end subroutine read_permeability_and_porosity


  !
  ! Initialize inflow and outflow.
  !
  subroutine init_flw_trnc_norm_xin_pt_out(nx, ny, nz, mu, sigma, Q)
    integer :: nx, ny, nz
    double precision :: mu, sigma
    double precision, dimension((nx * ny * nz)) :: Q

    integer :: i, j
    double precision :: x, pi, pdf, mass
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
       ! get the real x coordinate
       x = -1.5d0 + ((i - 1) * 3.0d0)/(nx - 1)    ! Mapping x = [-1.5, 1.5] to nx dimension

       ! Now use mu and sigma to find the pdf value at x
       pdf = 1.0d0/(sigma * sqrt(2.0d0 * pi)) * exp(-(((x - mu)/sigma)**2.0d0)/2.0d0)

       ! set the value at the index equal to the pdf value at that point
       Q_x(i) = pdf

       ! increment the mass by the value of the pdf
       mass = mass + pdf

       ! index to test initialization by plot
       idx(i) = i * 1.0
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

  subroutine wrapper(nx, ny, nz, nd, pt, st, mu, sigma, Q, S, P, V, Tt, Pc, oil) 
    integer :: nx, ny, nz
    integer :: nd, pt, st
    double precision :: mu, sigma
    double precision, dimension((nx * ny * nz)) :: Q
    double precision, dimension((nx * ny * nz)) :: S
    double precision, dimension(nx, ny, nz) :: P
    double precision, dimension(3, nx + 1, ny + 1, nz + 1) :: V
    double precision, dimension((nd/st) + 1) :: Tt
    double precision, dimension(2, (nd/st) + 1) :: Pc
    double precision :: oil

    call init_flw_trnc_norm_xin_pt_out(nx, ny, nz, mu, sigma, Q)
    call simulate_reservoir(nx, ny, nz, nd, pt, st, Q, S, P, V, Tt, Pc, oil)
  end subroutine wrapper
end module simulation
