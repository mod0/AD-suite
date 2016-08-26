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
    character(len=50) :: filename
    integer :: ncid
    integer :: mode
    
    ! Try chunk size = 4K
    nc_chunksize = 4096

    call iserror(nf90_open(filename, mode, ncid, nc_chunksize))

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
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    integer :: var
    integer, dimension(1) :: vararray
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(1) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 0 or 1
    if (dimcount /= 1 .or. dimcount /= 0) then
       stop "The input variable is not a scalar"
    end if

    if(dimcount == 1) then
       ! now get the dimensions
       call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

       ! check that the
       ! for each dimension, get the length of the dimension
       do i = 1, dimcount
          call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

          ! check that this length equals 1
          if (dimlens(i) /= 1) then
             stop "Length of array mismatch."
          end if
       enddo

       ! read the variable.
       call iserror(nf90_get_var(ncid, varid, vararray))

       var = vararray(1)
    else
       ! read the variable.
       call iserror(nf90_get_var(ncid, varid, var))
    end if
  end subroutine ncread_scalar_int

  ! routine to read a scalar double value
  subroutine ncread_scalar_double(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    double precision :: var
    double precision, dimension(1) :: vararray
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(1) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 0 or 1
    if (dimcount /= 1 .or. dimcount /= 0) then
       stop "The input variable is not a scalar"
    end if

    if(dimcount == 1) then
       ! now get the dimensions
       call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

       ! check that the
       ! for each dimension, get the length of the dimension
       do i = 1, dimcount
          call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

          ! check that this length equals 1
          if (dimlens(i) /= 1) then
             stop "Length of array mismatch."
          end if
       enddo

       ! read the variable.
       call iserror(nf90_get_var(ncid, varid, vararray))

       var = vararray(1)
    else
       ! read the variable.
       call iserror(nf90_get_var(ncid, varid, var))
    end if
  end subroutine ncread_scalar_double
  
  ! routine to read a scalar float value
  subroutine ncread_scalar_float(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    real*4 :: var
    real*4, dimension(1) :: vararray
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(1) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 0 or 1
    if (dimcount /= 1 .or. dimcount /= 0) then
       stop "The input variable is not a scalar"
    end if

    if(dimcount == 1) then
       ! now get the dimensions
       call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

       ! check that the
       ! for each dimension, get the length of the dimension
       do i = 1, dimcount
          call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

          ! check that this length equals 1
          if (dimlens(i) /= 1) then
             stop "Length of array mismatch."
          end if
       enddo

       ! read the variable.
       call iserror(nf90_get_var(ncid, varid, vararray))

       var = vararray(1)
    else
       ! read the variable.
       call iserror(nf90_get_var(ncid, varid, var))
    end if
  end subroutine ncread_scalar_float

  ! routine to read a vector integer value
  subroutine ncread_vector_int(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    integer, dimension(:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(1) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 1
    if (dimcount /= 1) then
       stop "The input variable is not a vector"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_vector_int

  ! routine to read a vector double value
  subroutine ncread_vector_double(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    double precision, dimension(:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(1) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 1
    if (dimcount /= 1) then
       stop "The input variable is not a vector"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_vector_double

  ! routine to read a vector float value
  subroutine ncread_vector_float(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    real*4, dimension(:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(1) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 1
    if (dimcount /= 1) then
       stop "The input variable is not a vector"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_vector_float

  ! routine to read a matrix integer value
  subroutine ncread_matrix_int(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    integer, dimension(:,:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(2) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 2
    if (dimcount /= 2) then
       stop "The input variable is not a matrix"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_matrix_int

  ! routine to read a matrix double value
  subroutine ncread_matrix_double(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    double precision, dimension(:,:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(2) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 2
    if (dimcount /= 2) then
       stop "The input variable is not a matrix"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_matrix_double

! routine to read a matrix float value
  subroutine ncread_matrix_float(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    real*4, dimension(:,:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(2) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 2
    if (dimcount /= 2) then
       stop "The input variable is not a vector"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_matrix_float

  ! routine to read a 3-tensor integer value
  subroutine ncread_3tensor_int(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    integer, dimension(:,:,:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(3) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 3
    if (dimcount /= 3) then
       stop "The input variable is not a 3-tensor"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_3tensor_int

  ! routine to read a 3-tensor double value
  subroutine ncread_3tensor_double(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    double precision, dimension(:,:,:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(3) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 3
    if (dimcount /= 3) then
       stop "The input variable is not a 3-tensor"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_3tensor_double

  ! routine to read a 3-tensor float value
  subroutine ncread_3tensor_float(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    real*4, dimension(:,:,:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(3) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 3
    if (dimcount /= 3) then
       stop "The input variable is not a 3-tensor"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_3tensor_float

  ! routine to read a 4-tensor integer value
  subroutine ncread_4tensor_int(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    integer, dimension(:,:,:,:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(4) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 4
    if (dimcount /= 4) then
       stop "The input variable is not a 4-tensor"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_4tensor_int

  ! routine to read a 4-tensor double value
  subroutine ncread_4tensor_double(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    double precision, dimension(:,:,:,:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(4) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 4
    if (dimcount /= 4) then
       stop "The input variable is not a 4-tensor"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
    call iserror(nf90_get_var(ncid, varid, var))    
  end subroutine ncread_4tensor_double

  ! routine to read a 4-tensor float value
  subroutine ncread_4tensor_float(ncid, varname, var)
    integer :: i
    integer :: ncid
    integer :: varid
    integer :: dimcount
    integer :: vartype
    real*4, dimension(:,:,:,:) :: var
    character(len=nf90_max_name) :: varname
    integer, dimension(nf90_max_var_dims) :: dimids
    integer, dimension(nf90_max_var_dims) :: dimlens
    character(len=nf90_max_name), dimension(4) :: dimnames

    ! get the variable id
    call iserror(nf90_inq_varid(ncid, varname, varid))

    ! now inquire the variable and get the number of dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount))
 
    ! check that the dimensions is equal to 4
    if (dimcount /= 4) then
       stop "The input variable is not a 4-tensor"
    end if

    ! now get the dimensions
    call iserror(nf90_inquire_variable(ncid, varid, varname, vartype, dimcount, dimids(:dimcount)))

    ! for each dimension, get the length of the dimension
    do i = 1, dimcount
       call iserror(nf90_inquire_dimension(ncid, dimids(i), dimnames(i), dimlens(i)))

       ! check that this length matches the size of the variable in that dimension
       if (size(var, i) /= dimlens(i)) then
          stop "Length of array mismatch."
       end if
    enddo

    ! read the variable.
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
  subroutine ncwrite_vector_int(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1
    integer, dimension(:) :: var
    character(len=nf90_max_name), dimension(1), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define vector dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
    else
       ! define vector dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_vector_int

  ! routine to write a matrix of integer values
  subroutine ncwrite_matrix_int(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2
    integer, dimension(:,:) :: var
    character(len=nf90_max_name), dimension(2), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define matrix dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, dimnames(2), size(var, 2), vardimid2))
    else
       ! define matrix dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, (/vardimid1, vardimid2/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_matrix_int

  ! routine to write a 3-tensor of integer values
  subroutine ncwrite_3tensor_int(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    integer, dimension(:,:,:) :: var
    character(len=nf90_max_name), dimension(3), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define 3-tensor dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, dimnames(2), size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, dimnames(3), size(var, 3), vardimid3))
    else
       ! define 3-tensor dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, &
         (/vardimid1, vardimid2, vardimid3/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_3tensor_int


  ! routine to write a 4-tensor of integer values
  subroutine ncwrite_4tensor_int(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3, vardimid4
    integer, dimension(:,:,:,:) :: var
    character(len=nf90_max_name), dimension(4), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define 4-tensor dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, dimnames(2), size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, dimnames(3), size(var, 3), vardimid3))
       call iserror(nf90_def_dim(ncid, dimnames(4), size(var, 4), vardimid4))
    else
       ! define 4-tensor dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim4", size(var, 4), vardimid4))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_int, &
         (/vardimid1, vardimid2, vardimid3, vardimid4/), varid))
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
  subroutine ncwrite_vector_double(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1
    double precision, dimension(:) :: var
    character(len=nf90_max_name), dimension(1), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define vector dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
    else
       ! define vector dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_double, (/vardimid1/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_vector_double

  ! routine to write a matrix of double values
  subroutine ncwrite_matrix_double(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2
    double precision, dimension(:,:) :: var
    character(len=nf90_max_name), dimension(2), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define matrix dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, dimnames(2), size(var, 2), vardimid2))
    else
       ! define matrix dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_double, (/vardimid1, vardimid2/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_matrix_double

  ! routine to write a 3-tensor of double values
  subroutine ncwrite_3tensor_double(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    double precision, dimension(:,:,:) :: var
    character(len=nf90_max_name), dimension(3), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define 3-tensor dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, dimnames(2), size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, dimnames(3), size(var, 3), vardimid3))
    else
       ! define 3-tensor dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_double, &
         (/vardimid1, vardimid2, vardimid3/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_3tensor_double


  ! routine to write a 4-tensor of double values
  subroutine ncwrite_4tensor_double(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3, vardimid4
    double precision, dimension(:,:,:,:) :: var
    character(len=nf90_max_name), dimension(4), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define 4-tensor dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, dimnames(2), size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, dimnames(3), size(var, 3), vardimid3))
       call iserror(nf90_def_dim(ncid, dimnames(4), size(var, 4), vardimid4))
    else
       ! define 4-tensor dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim4", size(var, 4), vardimid4))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_double, &
         (/vardimid1, vardimid2, vardimid3, vardimid4/), varid))
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
  subroutine ncwrite_vector_float(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1
    real*4, dimension(:) :: var
    character(len=nf90_max_name), dimension(1), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define vector dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
    else
       ! define vector dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_float, (/vardimid1/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_vector_float

  ! routine to write a matrix of float values
  subroutine ncwrite_matrix_float(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2
    real*4, dimension(:,:) :: var
    character(len=nf90_max_name), dimension(2), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define matrix dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, dimnames(2), size(var, 2), vardimid2))
    else
       ! define matrix dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_float, &
         (/vardimid1, vardimid2/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_matrix_float

  ! routine to write a 3-tensor of float values
  subroutine ncwrite_3tensor_float(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3
    real*4, dimension(:,:,:) :: var
    character(len=nf90_max_name), dimension(3), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define 3-tensor dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, dimnames(2), size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, dimnames(3), size(var, 3), vardimid3))
    else
       ! define 3-tensor dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_float, &
         (/vardimid1, vardimid2, vardimid3/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_3tensor_float


  ! routine to write a 4-tensor of float values
  subroutine ncwrite_4tensor_float(ncid, varname, var, dimnames)
    integer :: ncid
    character(len=31) :: varname
    integer :: varid
    integer :: vardimid1, vardimid2, vardimid3, vardimid4
    real*4, dimension(:,:,:,:) :: var
    character(len=nf90_max_name), dimension(4), optional :: dimnames

    ! enter define mode.
    call iserror(nf90_redef(ncid))

    if (present(dimnames)) then
       ! define 4-tensor dimensions
       call iserror(nf90_def_dim(ncid, dimnames(1), size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, dimnames(2), size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, dimnames(3), size(var, 3), vardimid3))
       call iserror(nf90_def_dim(ncid, dimnames(4), size(var, 4), vardimid4))
    else
       ! define 4-tensor dimensions
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim1", size(var, 1), vardimid1))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim2", size(var, 2), vardimid2))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim3", size(var, 3), vardimid3))
       call iserror(nf90_def_dim(ncid, trim(adjustl(varname))//"_dim4", size(var, 4), vardimid4))
    end if

    ! exit define mode
    call iserror(nf90_enddef(ncid))

    ! define variable
    call iserror(nf90_def_var(ncid, varname, nf90_float, &
         (/vardimid1, vardimid2, vardimid3, vardimid4/), varid))
    call iserror(nf90_put_var(ncid, varid, var))
  end subroutine ncwrite_4tensor_float
end module netcdfwrapper
