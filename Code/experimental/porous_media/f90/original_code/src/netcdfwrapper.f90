module netcdfwrapper
  use netcdf

  public :: ncread, ncwrite, ncopen, ncclose

  interface ncread
     module procedure ncread_scalar_int
     module procedure ncread_vector_int
     module procedure ncread_matrix_int
     module procedure ncread_3tensor_int
     module procedure ncread_4tensor_int
     module procedure ncread_scalar_float
     module procedure ncread_vector_float
     module procedure ncread_matrix_float
     module procedure ncread_3tensor_float
     module procedure ncread_4tensor_float
     module procedure ncread_scalar_double
     module procedure ncread_vector_double
     module procedure ncread_matrix_double
     module procedure ncread_3tensor_double
     module procedure ncread_4tensor_double
  end interface ncread

  interface ncwrite
     module procedure ncwrite_scalar_int
     module procedure ncwrite_vector_int
     module procedure ncwrite_matrix_int
     module procedure ncwrite_3tensor_int
     module procedure ncwrite_4tensor_int
     module procedure ncwrite_scalar_float
     module procedure ncwrite_vector_float
     module procedure ncwrite_matrix_float
     module procedure ncwrite_3tensor_float
     module procedure ncwrite_4tensor_float
     module procedure ncwrite_scalar_double
     module procedure ncwrite_vector_double
     module procedure ncwrite_matrix_double
     module procedure ncwrite_3tensor_double
     module procedure ncwrite_4tensor_double
  end interface ncwrite

contains

  ! ============= OPEN/CLOSE/HANDLE ERROR =============

  ! Open file in mode and return ncid
  subroutine ncopen(filename, mode, ncid)
    integer :: nc_chunksize
    parameter(nc_chunksize = 4096)
    character(len=50) :: filename
    integer :: ncid
    integer :: mode

    call iserror(nf90_open(filename, &
         mode, ncid, nc_chunksize))

    ! if the file is opened to write, end define mode
    ! each write operation will enter define mode 
    ! independently
    if (mode == NF90_CLOBBER) then
       call iserror(nf90_enddef(ncid))
    end if
  end subroutine ncopen
  
  ! close the file identified by ncid
  subroutine ncclose(ncid)
    integer :: ncid

    call iserror(nf90_close(ncid))
  end subroutine ncclose

  ! netCDF Error Check Routine
  subroutine iserror(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Error encountered. Stopped."
    end if
  end subroutine iserror

  !============== READ ROUTINES========================

  ! routine to read a scalar integer value
  subroutine ncread_scalar_int(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid, var

    
    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_scalar_int

  ! routine to read a scalar double value
  subroutine ncread_scalar_double(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    double precision :: var

    
    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_scalar_double

  ! routine to read a scalar float value
  subroutine ncread_scalar_float(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    real*4 :: var

    
    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_scalar_float

  ! routine to read a vector integer value
  subroutine ncread_vector_int(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1
    integer :: varid, vardimid1, vardim1
    integer, dimension(:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    if (size(var, 1) /= vardim1) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_vector_int


  ! routine to read a vector double value
  subroutine ncread_vector_double(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1
    integer :: varid, vardimid1, vardim1
    double precision, dimension(:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    if (size(var, 1) /= vardim1) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_vector_double


  ! routine to read a vector float value
  subroutine ncread_vector_float(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1
    integer :: varid, vardimid1, vardim1
    real*4, dimension(:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    if (size(var, 1) /= vardim1) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_vector_float

  ! routine to read a matrix of integer values
  subroutine ncread_matrix_int(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1, vardimname2
    integer :: varid, vardimid1, vardimid2, vardim1, vardim2
    integer, dimension(:,:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim2", vardimid2))
    call iserror(nf90_inquire_dimension(ncid, vardimid2, vardimname2, vardim2))

    if (size(var, 1) /= vardim1 &
         .or. size(var, 2) /= vardim2) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_matrix_int

  ! routine to read a matrix of double values
  subroutine ncread_matrix_double(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1, vardimname2
    integer :: varid, vardimid1, vardimid2, vardim1, vardim2
    double precision, dimension(:,:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim2", vardimid2))
    call iserror(nf90_inquire_dimension(ncid, vardimid2, vardimname2, vardim2))

    if (size(var, 1) /= vardim1 &
         .or. size(var, 2) /= vardim2) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_matrix_double

  ! routine to read a matrix of float values
  subroutine ncread_matrix_float(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1, vardimname2
    integer :: varid, vardimid1, vardimid2, vardim1, vardim2
    real*4, dimension(:,:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim2", vardimid2))
    call iserror(nf90_inquire_dimension(ncid, vardimid2, vardimname2, vardim2))

    if (size(var, 1) /= vardim1 &
         .or. size(var, 2) /= vardim2) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_matrix_float

  ! routine to read a 3-tensor of integer values
  subroutine ncread_3tensor_int(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1, vardimname2, vardimname3
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    integer :: vardim1, vardim2, vardim3
    integer, dimension(:,:,:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim2", vardimid2))
    call iserror(nf90_inquire_dimension(ncid, vardimid2, vardimname2, vardim2))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim3", vardimid3))
    call iserror(nf90_inquire_dimension(ncid, vardimid3, vardimname3, vardim3))

    if (size(var, 1) /= vardim1 &
         .or. size(var, 2) /= vardim2 &
         .or. size(var, 3) /= vardim3) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_3tensor_int

  ! routine to read a 3-tensor of double values
  subroutine ncread_3tensor_double(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1, vardimname2, vardimname3
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    integer :: vardim1, vardim2, vardim3
    double precision, dimension(:,:,:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim2", vardimid2))
    call iserror(nf90_inquire_dimension(ncid, vardimid2, vardimname2, vardim2))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim3", vardimid3))
    call iserror(nf90_inquire_dimension(ncid, vardimid3, vardimname3, vardim3))

    if (size(var, 1) /= vardim1 &
         .or. size(var, 2) /= vardim2 &
         .or. size(var, 3) /= vardim3) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_3tensor_double

  ! routine to read a 3-tensor of float values
  subroutine ncread_3tensor_float(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1, vardimname2, vardimname3
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    integer :: vardim1, vardim2, vardim3
    real*4, dimension(:,:,:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim2", vardimid2))
    call iserror(nf90_inquire_dimension(ncid, vardimid2, vardimname2, vardim2))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim3", vardimid3))
    call iserror(nf90_inquire_dimension(ncid, vardimid3, vardimname3, vardim3))

    if (size(var, 1) /= vardim1 &
         .or. size(var, 2) /= vardim2 &
         .or. size(var, 3) /= vardim3) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_3tensor_float

  ! routine to read a 4-tensor of integer values
  subroutine ncread_4tensor_int(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1, vardimname2, vardimname3, vardimname4
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3, vardimid4
    integer :: vardim1, vardim2, vardim3, vardim4
    integer, dimension(:,:,:,:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim2", vardimid2))
    call iserror(nf90_inquire_dimension(ncid, vardimid2, vardimname2, vardim2))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim3", vardimid3))
    call iserror(nf90_inquire_dimension(ncid, vardimid3, vardimname3, vardim3))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim4", vardimid4))
    call iserror(nf90_inquire_dimension(ncid, vardimid4, vardimname4, vardim4))

    if (size(var, 1) /= vardim1 &
         .or. size(var, 2) /= vardim2 &
         .or. size(var, 3) /= vardim3 &
         .or. size(var, 4) /= vardim4) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_4tensor_int

  ! routine to read a 4-tensor of double values
  subroutine ncread_4tensor_double(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1, vardimname2, vardimname3, vardimname4
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3, vardimid4
    integer :: vardim1, vardim2, vardim3, vardim4
    double precision, dimension(:,:,:,:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim2", vardimid2))
    call iserror(nf90_inquire_dimension(ncid, vardimid2, vardimname2, vardim2))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim3", vardimid3))
    call iserror(nf90_inquire_dimension(ncid, vardimid3, vardimname3, vardim3))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim4", vardimid4))
    call iserror(nf90_inquire_dimension(ncid, vardimid4, vardimname4, vardim4))

    if (size(var, 1) /= vardim1 &
         .or. size(var, 2) /= vardim2 &
         .or. size(var, 3) /= vardim3 &
         .or. size(var, 4) /= vardim4) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_4tensor_double

  ! routine to read a 4-tensor of float values
  subroutine ncread_4tensor_float(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    character(len=50) :: vardimname1, vardimname2, vardimname3, vardimname4
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3, vardimid4
    integer :: vardim1, vardim2, vardim3, vardim4
    real*4, dimension(:,:,:,:) :: var

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim1", vardimid1))
    call iserror(nf90_inquire_dimension(ncid, vardimid1, vardimname1, vardim1))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim2", vardimid2))
    call iserror(nf90_inquire_dimension(ncid, vardimid2, vardimname2, vardim2))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim3", vardimid3))
    call iserror(nf90_inquire_dimension(ncid, vardimid3, vardimname3, vardim3))

    call iserror(nf90_inq_dimid(ncid, trim(adjustl(varname))//"_dim4", vardimid4))
    call iserror(nf90_inquire_dimension(ncid, vardimid4, vardimname4, vardim4))

    if (size(var, 1) /= vardim1 &
         .or. size(var, 2) /= vardim2 &
         .or. size(var, 3) /= vardim3 &
         .or. size(var, 4) /= vardim4) then
       stop "Length of array mismatch."
    end if

    call iserror(nf90_inq_varid(ncid, varname, varid))
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_4tensor_float

  !=============== WRITE ROUTINES =====================

  ! routine to write a scalar integer value
  subroutine ncwrite_scalar_int(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid, var

    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_scalar_int

  ! routine to write a vector integer value
  subroutine ncwrite_vector_int(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1
    integer, dimension(:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define vector dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_vector_int

  ! routine to write a matrix of integer values
  subroutine ncwrite_matrix_int(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2
    integer, dimension(:,:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define matrix dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1, vardimid2/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_matrix_int

  ! routine to write a 3-tensor of integer values
  subroutine ncwrite_3tensor_int(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    integer, dimension(:,:,:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define 3-tensor dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1, vardimid2, vardimid3/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_3tensor_int


  ! routine to write a 4-tensor of integer values
  subroutine ncwrite_4tensor_int(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    integer, dimension(:,:,:,:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define 4-tensor dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim4", size(var, 4), vardimid4))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1, vardimid2, vardimid3, vardimid4/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_4tensor_int

  ! routine to write a scalar double value
  subroutine ncwrite_scalar_double(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    double precision :: var

    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_scalar_double

  ! routine to write a vector double value
  subroutine ncwrite_vector_double(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1
    double precision, dimension(:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define vector dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_vector_double

  ! routine to write a matrix of double values
  subroutine ncwrite_matrix_double(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2
    double precision, dimension(:,:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define matrix dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1, vardimid2/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_matrix_double

  ! routine to write a 3-tensor of double values
  subroutine ncwrite_3tensor_double(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    double precision, dimension(:,:,:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define 3-tensor dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1, vardimid2, vardimid3/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_3tensor_double


  ! routine to write a 4-tensor of double values
  subroutine ncwrite_4tensor_double(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    double precision, dimension(:,:,:,:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define 4-tensor dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim4", size(var, 4), vardimid4))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1, vardimid2, vardimid3, vardimid4/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_4tensor_double


  ! routine to write a scalar float value
  subroutine ncwrite_scalar_float(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    real*4 :: var

    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_scalar_float

  ! routine to write a vector float value
  subroutine ncwrite_vector_float(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1
    real*4, dimension(:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define vector dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_vector_float

  ! routine to write a matrix of float values
  subroutine ncwrite_matrix_float(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2
    real*4, dimension(:,:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define matrix dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1, vardimid2/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_matrix_float

  ! routine to write a 3-tensor of float values
  subroutine ncwrite_3tensor_float(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    real*4, dimension(:,:,:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define 3-tensor dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1, vardimid2, vardimid3/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_3tensor_float


  ! routine to write a 4-tensor of float values
  subroutine ncwrite_4tensor_float(ncid, varname, var)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    real*4, dimension(:,:,:,:) :: var

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    ! define 4-tensor dimensions
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))
    call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim4", size(var, 4), vardimid4))

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1, vardimid2, vardimid3, vardimid4/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_4tensor_float
end module netcdfwrapper
