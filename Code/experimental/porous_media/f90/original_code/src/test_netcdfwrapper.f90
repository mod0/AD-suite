program test_netcdfwrapper
  use netcdf
  use netcdfwrapper
  implicit none
  integer :: ncid
  integer :: int0_write
  integer, dimension(10) :: int1_write
  integer, dimension(10,10) :: int2_write
  integer, dimension(10,10,10) :: int3_write
  integer, dimension(10,10,10,10) :: int4_write
  integer :: int0_read
  integer, dimension(10) :: int1_read
  integer, dimension(10,10) :: int2_read
  integer, dimension(10,10,10) :: int3_read
  integer, dimension(10,10,10,10) :: int4_read
  double precision :: double0_write
  double precision, dimension(10) :: double1_write
  double precision, dimension(10,10) :: double2_write
  double precision, dimension(10,10,10) :: double3_write
  double precision, dimension(10,10,10,10) :: double4_write
  double precision :: double0_read
  double precision, dimension(10) :: double1_read
  double precision, dimension(10,10) :: double2_read
  double precision, dimension(10,10,10) :: double3_read
  double precision, dimension(10,10,10,10) :: double4_read
  real*4 :: float0_write
  real*4, dimension(10) :: float1_write
  real*4, dimension(10,10) :: float2_write
  real*4, dimension(10,10,10) :: float3_write
  real*4, dimension(10,10,10,10) :: float4_write
  real*4 :: float0_read
  real*4, dimension(10) :: float1_read
  real*4, dimension(10,10) :: float2_read
  real*4, dimension(10,10,10) :: float3_read
  real*4, dimension(10,10,10,10) :: float4_read


  call ncopen("test_netcdfwrapper.nc", .true., ncid)
  call test_write_integer(ncid, int0_write, int1_write, int2_write, &
       int3_write, int4_write)
  call ncclose(ncid)

  call ncopen("test_netcdfwrapper.nc", .false., ncid)
  call test_read_integer(ncid, int0_read, int1_read, int2_read, &
       int3_read, int4_read)
  call ncclose(ncid)

  write (*,*) int0_write - int0_read, int0_write, int0_read
  write (*,*) sum(int1_write - int1_read), sum(int1_write), sum(int1_read)
  write (*,*) sum(int2_write - int2_read), sum(int2_write), sum(int2_read)
  write (*,*) sum(int3_write - int3_read), sum(int3_write), sum(int3_read)
  write (*,*) sum(int4_write - int4_read), sum(int4_write), sum(int4_read)

  call ncopen("test_netcdfwrapper.nc", .true., ncid)
  call test_write_double(ncid, double0_write, double1_write, double2_write, &
       double3_write, double4_write)
  call ncclose(ncid)

  call ncopen("test_netcdfwrapper.nc", .false., ncid)
  call test_read_double(ncid, double0_read, double1_read, double2_read, &
       double3_read, double4_read)
  call ncclose(ncid)

  write (*,*) double0_write - double0_read, double0_write, double0_read
  write (*,*) sum(double1_write - double1_read), sum(double1_write), sum(double1_read)
  write (*,*) sum(double2_write - double2_read), sum(double2_write), sum(double2_read)
  write (*,*) sum(double3_write - double3_read), sum(double3_write), sum(double3_read)
  write (*,*) sum(double4_write - double4_read), sum(double4_write), sum(double4_read)

  call ncopen("test_netcdfwrapper.nc", .true., ncid)
  call test_write_float(ncid, float0_write, float1_write, float2_write, &
       float3_write, float4_write)
  call ncclose(ncid)

  call ncopen("test_netcdfwrapper.nc", .false., ncid)
  call test_read_float(ncid, float0_read, float1_read, float2_read, &
       float3_read, float4_read)
  call ncclose(ncid)

  write (*,*) float0_write - float0_read, float0_write, float0_read
  write (*,*) sum(float1_write - float1_read), sum(float1_write), sum(float1_read)
  write (*,*) sum(float2_write - float2_read), sum(float2_write), sum(float2_read)
  write (*,*) sum(float3_write - float3_read), sum(float3_write), sum(float3_read)
  write (*,*) sum(float4_write - float4_read), sum(float4_write), sum(float4_read)
contains

  subroutine test_write_integer(ncid, int0_write, int1_write, int2_write, &
       int3_write, int4_write)
    integer :: ncid
    integer :: i, j, k, l
    integer :: int0_write
    integer, dimension(10), intent(out) :: int1_write
    integer, dimension(10,10), intent(out) :: int2_write
    integer, dimension(10,10,10), intent(out) :: int3_write
    integer, dimension(10,10,10,10), intent(out) :: int4_write

    int0_write = 1
    do i = 1,10
       int1_write(i) = i
       do j = 1,10
          int2_write(i, j) = (i+j) + (i * j)
          do k = 1,10
             int3_write(i, j, k) = (i + j + k) + (i * j * k) + (i * j + j * k + i * k)
             do l = 1,10
                int4_write(i, j, k, l) = (i + j + k + l) + (i * j * k * l) + &
                (i * j * k + i * k * l + i * j * l + k * j * l) + &
                (i * j + i * k  + i * l + j * k + j * l + k * l)
             end do
          end do
       end do
    end do
    
    call ncwrite(ncid, "int0", int0_write)
    call ncwrite(ncid, "int1", int1_write)
    call ncwrite(ncid, "int2", int2_write)
    call ncwrite(ncid, "int3", int3_write)
    call ncwrite(ncid, "int4", int4_write)
  end subroutine test_write_integer

  subroutine test_read_integer(ncid, int0_read, int1_read, int2_read, &
       int3_read, int4_read)
    integer :: ncid
    integer :: int0_read
    integer, dimension(10), intent(out) :: int1_read
    integer, dimension(10,10), intent(out) :: int2_read
    integer, dimension(10,10,10), intent(out) :: int3_read
    integer, dimension(10,10,10,10), intent(out) :: int4_read

    call ncread(ncid, "int0", int0_read)
    call ncread(ncid, "int1", int1_read)
    call ncread(ncid, "int2", int2_read)
    call ncread(ncid, "int3", int3_read)
    call ncread(ncid, "int4", int4_read)
  end subroutine test_read_integer

  subroutine test_write_double(ncid, double0_write, double1_write, double2_write, &
       double3_write, double4_write)
    integer :: ncid
    integer :: i, j, k, l
    double precision :: double0_write
    double precision, dimension(10), intent(out) :: double1_write
    double precision, dimension(10,10), intent(out) :: double2_write
    double precision, dimension(10,10,10), intent(out) :: double3_write
    double precision, dimension(10,10,10,10), intent(out) :: double4_write

    double0_write = 1.0d0
    do i = 1,10
       double1_write(i) = 1.0d0/(1.0d0*i)
       do j = 1,10
          double2_write(i, j) = 1.0d0/(1.0d0*((i+j) + (i * j)))
          do k = 1,10
             double3_write(i, j, k) = 1.0d0/(1.0d0*((i + j + k) + (i * j * k) + (i * j + j * k + i * k)))
             do l = 1,10
                double4_write(i, j, k, l) = 1.0d0/(1.0d0*((i + j + k + l) + (i * j * k * l) + &
                (i * j * k + i * k * l + i * j * l + k * j * l) + &
                (i * j + i * k  + i * l + j * k + j * l + k * l)))
             end do
          end do
       end do
    end do
    
    
    call ncwrite(ncid, "double0", double0_write)
    call ncwrite(ncid, "double1", double1_write)
    call ncwrite(ncid, "double2", double2_write)
    call ncwrite(ncid, "double3", double3_write)
    call ncwrite(ncid, "double4", double4_write)
  end subroutine test_write_double

  subroutine test_read_double(ncid, double0_read, double1_read, double2_read, &
       double3_read, double4_read)
    integer :: ncid
    double precision :: double0_read
    double precision, dimension(10), intent(out) :: double1_read
    double precision, dimension(10,10), intent(out) :: double2_read
    double precision, dimension(10,10,10), intent(out) :: double3_read
    double precision, dimension(10,10,10,10), intent(out) :: double4_read

    call ncread(ncid, "double0", double0_read)
    call ncread(ncid, "double1", double1_read)
    call ncread(ncid, "double2", double2_read)
    call ncread(ncid, "double3", double3_read)
    call ncread(ncid, "double4", double4_read)
  end subroutine test_read_double

  subroutine test_write_float(ncid, float0_write, float1_write, float2_write, &
       float3_write, float4_write)
    integer :: ncid
    integer :: i, j, k, l
    real*4 :: float0_write
    real*4, dimension(10), intent(out) :: float1_write
    real*4, dimension(10,10), intent(out) :: float2_write
    real*4, dimension(10,10,10), intent(out) :: float3_write
    real*4, dimension(10,10,10,10), intent(out) :: float4_write

    float0_write = 1.0
    do i = 1,10
       float1_write(i) = 1.0/(1.0*i)
       do j = 1,10
          float2_write(i, j) = 1.0/(1.0*((i+j) + (i * j)))
          do k = 1,10
             float3_write(i, j, k) = 1.0/(1.0*((i + j + k) + (i * j * k) + (i * j + j * k + i * k)))
             do l = 1,10
                float4_write(i, j, k, l) = 1.0/(1.0*((i + j + k + l) + (i * j * k * l) + &
                (i * j * k + i * k * l + i * j * l + k * j * l) + &
                (i * j + i * k  + i * l + j * k + j * l + k * l)))
             end do
          end do
       end do
    end do
    
    call ncwrite(ncid, "float0", float0_write)
    call ncwrite(ncid, "float1", float1_write)
    call ncwrite(ncid, "float2", float2_write)
    call ncwrite(ncid, "float3", float3_write)
    call ncwrite(ncid, "float4", float4_write)
  end subroutine test_write_float

  subroutine test_read_float(ncid, float0_read, float1_read, float2_read, &
       float3_read, float4_read)
    integer :: ncid
    real*4 :: float0_read
    real*4, dimension(10), intent(out) :: float1_read
    real*4, dimension(10,10), intent(out) :: float2_read
    real*4, dimension(10,10,10), intent(out) :: float3_read
    real*4, dimension(10,10,10,10), intent(out) :: float4_read

    call ncread(ncid, "float0", float0_read)
    call ncread(ncid, "float1", float1_read)
    call ncread(ncid, "float2", float2_read)
    call ncread(ncid, "float3", float3_read)
    call ncread(ncid, "float4", float4_read)
  end subroutine test_read_float

end program test_netcdfwrapper
