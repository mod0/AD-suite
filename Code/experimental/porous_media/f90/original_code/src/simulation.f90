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

  subroutine allocate_shared_arrays(nx, ny, nz, st, pt, nd)
    integer :: nx, ny, nz
    integer :: st, pt, nd

    ! Allocate space for constant permeability/porosity
    allocate(POR(nx * ny * nz), PERM(3, nx, ny, nz))    
  end subroutine allocate_shared_arrays

  subroutine deallocate_shared_arrays()
    ! also deallocate Porosity and Permeability
    deallocate(POR, PERM)
  end subroutine deallocate_shared_arrays

  subroutine initialize_shared_arrays(data_directory, nx, ny, nz, st, pt, nd)
    character(len = *) :: data_directory                             ! directory location of parameters 

    integer :: nx, ny, nz
    integer :: st, pt, nd

    ! Now read the permeabilities and porosities which are global
    call read_permeability_and_porosity(data_directory, nx, ny, nz, PERM, POR)
  end subroutine initialize_shared_arrays

  subroutine allocate_arrays(nx, ny, nz, st, pt, nd, Q, S, P, V, Tt, Pc)
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
         Pc(2, (nd/st) + 1))
  end subroutine allocate_arrays

  subroutine deallocate_arrays(Q, S, P, V, Tt, Pc)
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
  end subroutine deallocate_arrays

  subroutine initialize_arrays(Q, S, P, V, Tt, Pc)
    double precision, dimension(:)          :: Q
    double precision, dimension(:)          :: S
    double precision, dimension(:, :, :)    :: P
    double precision, dimension(:, :, :, :) :: V
    double precision, dimension(:)          :: Tt
    double precision, dimension(:, :)       :: Pc

    ! initialize memory for inflow and saturation
    Q = 0.0d0
    S = 0.0d0

    ! initialize memory for P and V
    P = 0.0d0

    ! note that V vector has an additional
    ! length in each dimension x,y,z
    V = 0.0d0

    Tt = 0.0d0               ! simulation time
    Pc = 0.0d0               ! production data
  end subroutine initialize_arrays

  subroutine write_results(results_directory, nx, ny, nz, st, pt, nd, n_dof, mu, sigma, Tt, Pc, totaloil)
    implicit none
    character(len = *) :: results_directory

    integer :: nx, ny, nz, st, pt, nd, n_dof
    double precision :: totaloil, totaloil_mud, totaloil_sigmad
    double precision, dimension(:)             :: mu
    double precision, dimension(:)             :: sigma
    double precision, dimension(:)             :: Tt
    double precision, dimension(:, :)          :: Pc

    ! netCDF variables
    integer :: ncid                                                     ! file handle
    integer :: var_scene_id                                             ! scene variable id
    character(len = *), parameter :: var_scene_name = "Scenario"        ! scene variable name
    integer :: var_nx_id                                                ! nx variable id
    character(len = *), parameter :: var_nx_name = "NX"                 ! nx variable name
    integer :: var_ny_id                                                ! ny variable id
    character(len = *), parameter :: var_ny_name = "NY"                 ! ny variable name
    integer :: var_nz_id                                                ! nz variable id
    character(len = *), parameter :: var_nz_name = "NZ"                 ! nz variable name
!     integer :: var_mu_id                                                ! mu variable id
!     character(len = *), parameter :: var_mu_name = "Mu"                 ! mu variable name
!     integer :: var_sigma_id                                             ! sigma variable id
!     character(len = *), parameter :: var_sigma_name = "Sigma"           ! sigma variable name
                                                                        
    integer :: dim_time_id                                              ! time dimension id
    integer :: dim_time_len                                             ! time dimension length
    character(len = *), parameter :: dim_time_name = "Time"             ! time dimension name
                                                                        
    integer :: var_time_id                                              ! time variable id
    character(len = *), parameter :: var_time_name = "Time"             ! time variable name
                                                                        
    integer :: dim_mobility_id                                          ! mobility dimension id
    integer :: dim_mobility_len                                         ! mobility dimension length
    character(len = *), parameter :: dim_mobility_name = "Mobilities"   ! mobility dimension name
                                                                        
    integer :: var_mobility_id                                          ! mobility variable id
    character(len = *), parameter :: var_mobility_name = "WaterOil"     ! mobility variable name
    integer, dimension(2) :: var_mobility_dimids                        ! id of dimensions for mobilities

    integer :: var_oil_id                                               ! oil variable id
    character(len = *), parameter :: var_oil_name = "Oil"               ! oil variable name

    ! Start writing netCDF file having all the computed values
    ! Open file
    call iserror(nf90_create(trim(adjustl(results_directory))//"results_eval_original_code.nc", &
         nf90_clobber, ncid))

    ! setup the sizes of the dimensions
    dim_time_len = (nd/st) + 1
    dim_mobility_len = 2

    ! Define all dimension
    ! Define time dimensions
    call iserror(nf90_def_dim(ncid, dim_time_name, dim_time_len, dim_time_id))
    ! Define mobility dimensions
    call iserror(nf90_def_dim(ncid, dim_mobility_name, dim_mobility_len, dim_mobility_id))

    ! Define scalar scenario input
    call iserror(nf90_def_var(ncid, var_scene_name, NF90_INT, var_scene_id))
    ! Define scalar nx input
    call iserror(nf90_def_var(ncid, var_nx_name, NF90_INT, var_nx_id))
    ! Define scalar ny input
    call iserror(nf90_def_var(ncid, var_ny_name, NF90_INT, var_ny_id))
    ! Define scalar nz input
    call iserror(nf90_def_var(ncid, var_nz_name, NF90_INT, var_nz_id))
    ! Define scalar mu input
!     call iserror(nf90_def_var(ncid, var_mu_name, NF90_DOUBLE, var_mu_id))
    ! Define scalar sigma input
!     call iserror(nf90_def_var(ncid, var_sigma_name, NF90_DOUBLE, var_sigma_id))

    ! Define the time variables. Varid is returned.
    call iserror(nf90_def_var(ncid, var_time_name, NF90_DOUBLE, dim_time_id, var_time_id))

    ! Define the netCDF variables. The dimids array is used to pass the
    ! dimids of the dimensions of the netCDF variables.
    var_mobility_dimids = (/ dim_mobility_id, dim_time_id /)
    call iserror(nf90_def_var(ncid, var_mobility_name, NF90_DOUBLE, var_mobility_dimids &
         , var_mobility_id))

    ! Define scalar oil output
    call iserror(nf90_def_var(ncid, var_oil_name, NF90_DOUBLE, var_oil_id))

    ! End define mode.
    call iserror(nf90_enddef(ncid))

    ! Write scalar inputs
    call iserror(nf90_put_var(ncid, var_scene_id, scenario_id))
    call iserror(nf90_put_var(ncid, var_nx_id, nx))
    call iserror(nf90_put_var(ncid, var_ny_id, ny))
    call iserror(nf90_put_var(ncid, var_nz_id, nz))
!     call iserror(nf90_put_var(ncid, var_mu_id, mu))
!     call iserror(nf90_put_var(ncid, var_sigma_id, sigma))

    ! Write the time data. 
    call iserror(nf90_put_var(ncid, var_time_id, Tt))

    ! Write the data. This will write our mobility and oil data.
    ! The arrays of data are the same size as
    ! the netCDF variables we have defined.
    call iserror(nf90_put_var(ncid, var_mobility_id, Pc))
    call iserror(nf90_put_var(ncid, var_oil_id, totaloil))

    ! Close the file.
    call iserror(nf90_close(ncid))
  end subroutine write_results

  !
  ! This routine opens the permeability and porosity used by
  ! the MATLAB program and uses it for the simulation.
  !
  subroutine read_permeability_and_porosity(data_directory, nx, ny, nz, PERM, POR)   
    integer :: ncid                                                  ! netcdf file handle

    character(len = *) :: data_directory                             ! directory location of parameters 
    integer :: varid                                                 ! generic variable id
    integer :: nc_chunksize                                          ! chunk size for reading

    integer :: i, j, k, l, m
    integer :: nx, ny, nz
    
    double precision, dimension((nx * ny * nz)) :: POR       ! Porosities
    double precision, dimension((nx * ny * nz)) :: POR_temp  ! Porosities_temp
    double precision, dimension(3, nx, ny, nz)  :: PERM      ! Permeabilities

    ! initialize porosity and permeability to zero
    PERM = 0.0d0
    POR = 0.0d0
    POR_temp = 0.0d0

    ! read permeability
    nc_chunksize = 4096
    call iserror(nf90_open(trim(adjustl(data_directory))//"/permeability.nc", &
         NF90_NOWRITE, ncid, nc_chunksize))
    write(*,*) "Chosen chunk size is ", nc_chunksize
    
    call iserror(nf90_inq_varid(ncid, "permeability", varid))
    call iserror(nf90_get_var(ncid, varid, PERM))
    call iserror(nf90_close(ncid))

    ! read porosity
    nc_chunksize = 4096
    call iserror(nf90_open(trim(adjustl(data_directory))//"/porosity.nc", &
         NF90_NOWRITE, ncid, nc_chunksize))
    write(*,*) "Chosen chunk size is ", nc_chunksize
    
    call iserror(nf90_inq_varid(ncid, "porosity", varid))
    call iserror(nf90_get_var(ncid, varid, POR_temp))
    call iserror(nf90_close(ncid))

    call mymax_1_0_double(POR_temp, 1.0d-3, POR)
  end subroutine read_permeability_and_porosity
end module simulation
