module simulation
  use parameters
  use matrix
  use finitevolume
  use netcdfwrapper
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
          arg1 = -(((x-mu(j))/sigma(j))**2.0/2.0);
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
    Q((nx * ny * nz)) = -ir/2.0
    Q(ny) = -ir/2.0
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
    double precision, dimension(4, (nd/st) + 1) :: Pc
    double precision ::  oil

    integer :: i, j, k
    double precision :: tempoil1, tempoil2
    double precision, dimension(2) :: Mw, Mo, Mt

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
          Pc(1,k) = Mw(1)/Mt(1)
          Pc(2,k) = Mo(1)/Mt(1)
          Pc(3,k) = Mw(2)/Mt(2)
          Pc(4,k) = Mo(2)/Mt(2)

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
    double precision, dimension(2) :: Mw, Mo

    if (pressure_step == 1) then
       ! solve pressure
       call Pres(nx, ny, nz, Q, S, P, V)    ! Pressure solver
    endif

    call NewtRaph(nx, ny, nz, nd, pt, st, Q, V, S)      ! Solve for saturation
    call RelPerm(S((nx * ny * nz)), Mw(1), Mo(1))       ! Mobilities in well-block
    call RelPerm(S(ny), Mw(2), Mo(2))                   ! Mobilities in well-block
  end subroutine stepforward

  subroutine update_oil(nd, pt, st, Pc, k, oilin, oilout)
    integer :: k
    integer :: nd, pt, st
    double precision ::  oilin
    double precision ::  oilout
    double precision, dimension(4, (nd/st) + 1) :: Pc

    oilout = oilin +  (-Pc(1, k) + Pc(2, k) + Pc(3, k) - Pc(4, k)) * st    ! Reimann sum
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
    double precision, dimension(4, (nd/st) + 1) :: Pc
    double precision :: oil

    call init_flw_trnc_norm_xin_pt_out(nx, ny, nz, ndof, mu, sigma, Q)
    call simulate_reservoir(nx, ny, nz, nd, pt, st, Q, S, P, V, Tt, Pc, oil)
  end subroutine wrapper

  subroutine initialize_size_param_vars(nx, ny, nz, st, pt, nd, data_directory)
    integer :: ncid
    integer :: nx, ny, nz
    integer :: st, pt, nd
    character(len = *) :: data_directory     

    ! Open parameters file
    call ncopen(trim(adjustl(data_directory))//"/parameters1.nc", &
         NC_NOWRITE, ncid)

    ! read NX, NY, and NZ
    call ncread(ncid, "NX", nx)
    call ncread(ncid, "NY", ny)
    call ncread(ncid, "NZ", nz)

    ! read St, Pt, ND
    call ncread(ncid, "St",  st)
    call ncread(ncid, "Pt",  pt)
    call ncread(ncid, "ND",  nd)
    
    call ncclose(ncid)
  end subroutine initialize_size_param_vars

  subroutine initialize_scalar_param_vars(data_directory)
    implicit none
    integer :: ncid  
    character(len = *) :: data_directory 

    ! Open parameters file
    call ncopen(trim(adjustl(data_directory))//"/parameters1.nc", &
         NC_NOWRITE, ncid)

    ! initialize all constants insize parameters file.
    ! read ir_const
    call ncread(ncid, "ir_const", ir)

    ! read hX, hY, hZ
    call ncread(ncid, "hX", hx)
    call ncread(ncid, "hY", hy)
    call ncread(ncid, "hZ", hz)

    ! real volume
    call ncread(ncid, "vol", vol)

    ! read vw, vo, swc, sor
    call ncread(ncid, "vw", vw)
    call ncread(ncid, "vo", vo)
    call ncread(ncid, "swc", swc)
    call ncread(ncid, "sor", sor)

    ! read solver_parameters
    call ncread(ncid, "solver_inner", solver_inner)    
    call ncread(ncid, "solver_outer", solver_outer)    

    ! Close parameters file
    call ncclose(ncid)
  end subroutine initialize_scalar_param_vars

  subroutine allocate_vector_param_vars(nx, ny, nz)
    integer :: nx, ny, nz

    ! Allocate space for constant permeability/porosity
    allocate(POR(nx * ny * nz), PERM(3, nx, ny, nz))    
  end subroutine allocate_vector_param_vars

  subroutine deallocate_vector_param_vars()
    ! also deallocate Porosity and Permeability
    deallocate(POR, PERM)
  end subroutine deallocate_vector_param_vars

  subroutine initialize_vector_param_vars(nx, ny, nz, data_directory)
    character(len = *) :: data_directory                            
    integer :: nx, ny, nz

    ! Now read the permeabilities and porosities 
    call read_permeability_and_porosity(nx, ny, nz, PERM, POR, data_directory)
  end subroutine initialize_vector_param_vars

  subroutine read_permeability_and_porosity(nx, ny, nz, PERM, POR, data_directory)   
    integer :: ncid                                                  
    integer :: nx, ny, nz
    character(len = *) :: data_directory                             
    
    double precision, dimension((nx * ny * nz)) :: POR       ! Porosities
    double precision, dimension(nx, ny, nz) :: POR_temp    ! Porosities_temp
    double precision, dimension(3, nx, ny, nz)  :: PERM      ! Permeabilities

    ! initialize porosity and permeability to zero
    PERM = 0.0d0
    POR = 0.0d0
    POR_temp = 0.0d0

    ! read permeability
    call ncopen(trim(adjustl(data_directory))//"/parameters2.nc", &
         NC_NOWRITE, ncid)
    call ncread(ncid, "permeability", PERM)
    call ncclose(ncid)

    ! read porosity
    call ncopen(trim(adjustl(data_directory))//"/parameters3.nc", &
         NC_NOWRITE, ncid)
    call ncread(ncid, "porosity", POR_temp)
    call ncclose(ncid)    

    POR = reshape(POR_temp, (/nx * ny * nz/))

    call mymax_1_0_double(POR, 1.0d-3, POR)
  end subroutine read_permeability_and_porosity

  subroutine allocate_vector_sim_vars(nx, ny, nz, st, pt, nd, Q, S, P, V, Tt, Pc)
    implicit none
    integer :: nx, ny, nz
    integer :: st, pt, nd
     
    double precision, dimension(:), allocatable          :: Q
    double precision, dimension(:), allocatable          :: S
    double precision, dimension(:, :, :), allocatable    :: P
    double precision, dimension(:, :, :, :), allocatable :: V
    double precision, dimension(:), allocatable          :: Tt
    double precision, dimension(:, :), allocatable       :: Pc

    ! Allocate space for Q, S, P, V
    allocate(Q(nx * ny * nz), &
         S(nx * ny * nz), &
         P(nx, ny, nz), &
         V(3, nx+1, ny+1, nz+1))

    ! Allocate space for Time and Production data 
    allocate(Tt((nd/st) + 1), & 
         Pc(4, (nd/st) + 1))
  end subroutine allocate_vector_sim_vars

  subroutine deallocate_vector_sim_vars(Q, S, P, V, Tt, Pc)
    implicit none

    double precision, dimension(:), allocatable          :: Q
    double precision, dimension(:), allocatable          :: S
    double precision, dimension(:, :, :), allocatable    :: P
    double precision, dimension(:, :, :, :), allocatable :: V
    double precision, dimension(:), allocatable          :: Tt
    double precision, dimension(:, :), allocatable       :: Pc

    ! deallocate variables
    deallocate(Q, S, P, V)
    deallocate(Tt, Pc)
  end subroutine deallocate_vector_sim_vars

  subroutine initialize_vector_sim_vars(nx, ny, nz, Q, S, P, V, Tt, Pc, data_directory)
    integer :: ncid
    integer :: nx, ny, nz
    double precision, dimension(:)          :: Q
    double precision, dimension(:)          :: S
    double precision, dimension(:, :, :)    :: P
    double precision, dimension(:, :, :, :) :: V
    double precision, dimension(:)          :: Tt
    double precision, dimension(:, :)       :: Pc
    double precision, dimension(:,:,:), allocatable :: S_temp
 
    character(len = *) :: data_directory                             

    ! initialize memory for inflow and saturation
    Q = 0.0d0

    ! read initial saturation from file
    ! Open parameters file
    allocate(S_temp(nx, ny, nz))
    call ncopen(trim(adjustl(data_directory))//"/parameters4.nc", &
         NC_NOWRITE, ncid)
    call ncread(ncid, "saturation", S_temp) 
    call ncclose(ncid)
    S = reshape(S_temp, (/ nx * ny * nz /))
    deallocate(S_temp)

    ! initialize memory for P and V
    P = 0.0d0

    ! note that V vector has an additional
    ! length in each dimension x,y,z
    V = 0.0d0

    Tt = 0.0d0               ! simulation time
    Pc = 0.0d0               ! production data
  end subroutine initialize_vector_sim_vars

  subroutine write_results(results_directory, nx, ny, nz, st, pt, nd, n_dof, mu, sigma, Tt, Pc, totaloil)
    implicit none

    integer :: ncid                                                  
    integer :: nx, ny, nz, st, pt, nd, n_dof
    double precision :: totaloil
    double precision, dimension(:)             :: mu
    double precision, dimension(:)             :: sigma
    double precision, dimension(:)             :: Tt
    double precision, dimension(:, :)          :: Pc
    character(len = *) :: results_directory
                                                                        
    call ncopen(trim(adjustl(results_directory))//"/results_eval_original_code.nc", &
         NC_WRITE, ncid)

    call ncwrite(ncid, "NX", nx)
    call ncwrite(ncid, "NY", ny)
    call ncwrite(ncid, "NZ", nz)
    call ncwrite(ncid, "MU", mu)
    call ncwrite(ncid, "SIGMA", sigma)
    call ncwrite(ncid, "TIME", Tt)
    call ncwrite(ncid, "PROD_CURVES", Pc)
    call ncwrite(ncid, "OIL", totaloil)

    ! Close the file.
    call ncclose(ncid)
  end subroutine write_results
end module simulation
