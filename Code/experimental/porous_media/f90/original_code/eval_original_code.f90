program runspe10
  use parameters
  use utils
  use matrix
  use finitevolume
  use simulation
  use netcdf

  implicit none

  integer :: i, j, k
  integer :: nx, ny, nz
  integer :: nd, st, pt

  ! input/intermediate variables
  double precision :: mu, sigma
  double precision, dimension(:), allocatable :: Q
  double precision, dimension(:), allocatable :: S
  double precision, dimension(:, :, :), allocatable :: P
  double precision, dimension(:, :, :, :), allocatable :: V

  ! output variables
  double precision :: totaloil
  double precision, dimension(:), allocatable :: Tt
  double precision, dimension(:, :), allocatable :: Pc

  ! read the command line argument for the location of the data files
  integer :: iargs
  character(len = 100) :: data_directory
  character(len = 100) :: results_directory

  iargs = iargc()
  if(iargs /= 2) then
     write(*,*) "Incorrect number of arguments passed. Expected data and results directory"
  endif

  call getarg(1, data_directory)
  call getarg(2, results_directory)
  
  write(*,*) "Data directory: ", data_directory
  write(*,*) "Results directory: ", results_directory


  ! Not passing shared constants present in the parameters file.
  call initialize_scenario(data_directory, nx, ny, nz, st, pt, nd, solver_inner, solver_outer)

  ! allocate arrays
  call allocate_arrays(nx, ny, nz, st, pt, nd, Q, S, P, V, Tt, Pc)
  call allocate_shared_arrays(nx, ny, nz, st, pt, nd)

  ! initialize arrays
  call initialize_arrays(Q, S, P, V, Tt, Pc)
  call initialize_shared_arrays(data_directory, nx, ny, nz, st, pt, nd)  ! will not be reinitialized
  
  ! initialize scalar inputs and outputs
  mu = 0.0d0
  sigma = 1.0d0

  ! solver verbose parameter, inner and outer iterations are read from
  ! command file.
  verbose = .false.              ! Verbose solver output  

  call wrapper(nx, ny, nz, nd, pt, st, mu, sigma, Q, S, P, V, Tt, Pc, totaloil)
  
  ! write results
  call write_results(results_directory, nx, ny, nz, mu, sigma, Tt, Pc, totaloil)

  call deallocate_arrays(Q, S, P, V, Tt, Pc)
  call deallocate_shared_arrays()
  return
contains

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

  subroutine initialize_scenario(data_directory, nx, ny, nz, st, pt, nd, solver_inner, solver_outer)
    implicit none
    
    integer :: ncid                                                  ! netcdf file handle

    character(len = *) :: data_directory                             ! directory location of parameters 
    integer :: varid                                                 ! generic variable id
    integer :: nc_chunksize                                          ! chunk size for reading

    integer :: nx, ny, nz
    integer :: st, pt, nd
    integer :: solver_inner, solver_outer  

    ! Open parameters file
    nc_chunksize = 4096
    call iserror(nf90_open(trim(adjustl(data_directory))//"/parameters.nc", &
         NF90_NOWRITE, ncid, nc_chunksize))
    write(*,*) "Chosen chunk size is ", nc_chunksize

    ! read scenario_id
    call iserror(nf90_inq_varid(ncid, "scenario_id", varid))
    call iserror(nf90_get_var(ncid, varid, scenario_id))
    write (*,*) "Scenario: ", scenario_id

    ! read NX, NY, and NZ
    call iserror(nf90_inq_varid(ncid, "NX", varid))
    call iserror(nf90_get_var(ncid, varid, nx))
    call iserror(nf90_inq_varid(ncid, "NY", varid))
    call iserror(nf90_get_var(ncid, varid, ny))
    call iserror(nf90_inq_varid(ncid, "NZ", varid))
    call iserror(nf90_get_var(ncid, varid, nz))
    write (*,*) "(NX, NY, NZ) = ", NX, NY, NZ

    ! read St, Pt, ND
    call iserror(nf90_inq_varid(ncid, "St", varid))
    call iserror(nf90_get_var(ncid, varid, st))
    call iserror(nf90_inq_varid(ncid, "Pt", varid))
    call iserror(nf90_get_var(ncid, varid, pt))
    call iserror(nf90_inq_varid(ncid, "ND", varid))
    call iserror(nf90_get_var(ncid, varid, nd))
    write (*,*) "(St, Pt, ND) = ", st, pt, nd

    ! initialize all constants insize parameters file.
    ! read ir_const
    call iserror(nf90_inq_varid(ncid, "ir_const", varid))
    call iserror(nf90_get_var(ncid, varid, ir))

    ! compute ir
    ir = ir * nx * ny * nz
    write (*,*) "(ir) = ", ir

    ! read hX, hY, hZ
    call iserror(nf90_inq_varid(ncid, "hX", varid))
    call iserror(nf90_get_var(ncid, varid, hx_))
    call iserror(nf90_inq_varid(ncid, "hY", varid))
    call iserror(nf90_get_var(ncid, varid, hy_))
    call iserror(nf90_inq_varid(ncid, "hZ", varid))
    call iserror(nf90_get_var(ncid, varid, hz_))
    write (*,*) "(hX, hY, hZ) = ", hx_, hy_, hz_

    ! compute V
    V_ = hx_ * hy_ * hz_
    write (*,*) "(V) = ", V_

    ! read vw, vo, swc, sor
    call iserror(nf90_inq_varid(ncid, "vw", varid))
    call iserror(nf90_get_var(ncid, varid, vw_))
    call iserror(nf90_inq_varid(ncid, "vo", varid))
    call iserror(nf90_get_var(ncid, varid, vo_))
    call iserror(nf90_inq_varid(ncid, "swc", varid))
    call iserror(nf90_get_var(ncid, varid, swc_))
    call iserror(nf90_inq_varid(ncid, "sor", varid))
    call iserror(nf90_get_var(ncid, varid, sor_))
    write (*,*) "(vw, vo, swc, sor) = ", vw_, vo_, swc_, sor_ 

    ! read solver_parameters
    call iserror(nf90_inq_varid(ncid, "solver_inner", varid))
    call iserror(nf90_get_var(ncid, varid, solver_inner))    
    call iserror(nf90_inq_varid(ncid, "solver_outer", varid))
    call iserror(nf90_get_var(ncid, varid, solver_outer))    

    ! Close parameters file
    call iserror(nf90_close(ncid))
  end subroutine initialize_scenario

  subroutine write_results(results_directory, nx, ny, nz, mu, sigma, Tt, Pc, totaloil)
    implicit none
    character(len = *) :: results_directory

    integer :: nx, ny, nz
    double precision :: mu, sigma, totaloil, totaloil_mud, totaloil_sigmad
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
    integer :: var_mu_id                                                ! mu variable id
    character(len = *), parameter :: var_mu_name = "Mu"                 ! mu variable name
    integer :: var_sigma_id                                             ! sigma variable id
    character(len = *), parameter :: var_sigma_name = "Sigma"           ! sigma variable name
                                                                        
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
    call iserror(nf90_def_var(ncid, var_mu_name, NF90_DOUBLE, var_mu_id))
    ! Define scalar sigma input
    call iserror(nf90_def_var(ncid, var_sigma_name, NF90_DOUBLE, var_sigma_id))

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
    call iserror(nf90_put_var(ncid, var_mu_id, mu))
    call iserror(nf90_put_var(ncid, var_sigma_id, sigma))

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

  ! !
  ! ! This routine opens the permeability and porosity used by
  ! ! the MATLAB program and uses it for the simulation.
  ! !
  ! subroutine read_permeability_and_porosity(data_directory, nx, ny, nz, PERM, POR)   
    
  !   integer :: i, j, k, l, m
  !   integer :: nx, ny, nz
  !   integer :: maxnx, maxny, maxnz
  !   parameter(maxnx = 60, maxny = 220, maxnz = 85)
  !   character(len = *) :: data_directory                             ! directory location of parameters 
  !   double precision, dimension((nx * ny * nz)) :: POR                 ! Porosities
  !   double precision, dimension(3, nx, ny, nz) :: PERM  ! Permeabilities

  !   double precision, dimension(nx, ny, nz) :: P
  !   double precision, dimension(maxnx * maxny * maxnz) :: pUr
  !   double precision, dimension(3 * maxnx, maxny * maxnz) :: KUr
  !   double precision, dimension(3 * maxnx * maxny * maxnz) :: KUrl

  !   integer, dimension(nx * ny * nz) :: Pindices
  !   integer, dimension(3 * nx * ny * nz) :: Kindices

  !   ! initialize porosity and permeability to zero
  !   PERM = 0.0d0
  !   POR = 0.0d0

  !   ! read KUr
  !   open(1,file=trim(adjustl(data_directory))//"KUr.txt",status='old')
  !   read(1,*) ((KUr(i,j), j=1,maxny * maxnz), i=1,3 * maxnx)
  !   close(1)

  !   ! reshape 2 dimension to 1 dimension
  !   call myreshape_2_1(KUr, KUrl)

  !   ! select according to specified dimension
  !   m = 0
  !   do l = 1, nz
  !      do k = 1,ny
  !         do j = 1,nx
  !            do i = 1,3
  !               m = m + 1
  !               Kindices(m) = ((l - 1) * (maxnx * maxny * 3) &
  !                    + (k - 1) * (maxnx * 3) &
  !                    + 3 * (j-1) + i)
  !            end do
  !         end do
  !      end do
  !   end do

  !   ! then reshape 1 dimension to 4 dimension (hack for time being)
  !   call myreshape_1_4(KUrl(Kindices), PERM)

  !   ! read KUr
  !   open(1,file=trim(adjustl(data_directory))//"pUr.txt",status='old')
  !   read(1,*) (pUr(i), i=1,maxnx * maxny * maxnz)

  !   close(1)

  !   m = 0
  !   do k = 1,nz
  !      do j = 1,ny
  !         do i = 1,nx
  !            m = m + 1
  !            Pindices(m) = ((k - 1) * (maxnx * maxny) &
  !                 + (j - 1) * (maxnx) + i)
  !         end do
  !      end do
  !   end do

  !   call mymax_1_0_double(pUr(Pindices), 1.0d-3, POR) 
  ! end subroutine read_permeability_and_porosity

  ! netCDF Error Check Routine
  subroutine iserror(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Error encountered. Stopped."
    end if
  end subroutine iserror

end program runspe10
