module parameters
  ! FIXED PARAMETERS
  ! grid parameters 
  integer :: Nx_, Ny_, Nz_, N_, scenario_id
  double precision :: hx_, hy_, hz_, V_, ir
  integer :: maxNx, maxNy, maxNz

  parameter(maxNx=60, maxNy=220, maxNz=85)
  parameter(scenario_id = 1, &
       Nx_ = 10, &              ! Dimension in x-direction
       Ny_ = 10, &                   ! Dimension in y-direction
       Nz_ = 2, &                    ! Dimension in z-direction
       hx_ = 20.0d0 * 0.3048d0, &    ! step size in x-direction
       hy_ = 10.0d0 * 0.3048d0, &    ! step size in y-direction
       hz_ = 2.0d0 * 0.3048d0, &     ! step size in z-direction
       N_ = Nx_ * Ny_ * Nz_, &       ! Total number of grid cells
       V_ = hx_ * hy_ * hz_, &        ! Volume of each grid cell
       ir = (795.0 * Nx_ * Ny_ * Nz_) / (maxNx * maxNy * maxNz)) ! Magic number

  ! fluid parameters
  double precision :: vw_, vo_, swc_, sor_
  parameter(vw_ = 3d-4, &      ! Viscosity of Water
            vo_ = 3d-3, &      ! Viscosity of Oil
            swc_ = 0.2d0, &    ! Saturation of water cut
            sor_ = 0.2d0)      ! Saturation of oil cut

  ! timestepping parameters
  integer :: St, Pt, ND
  parameter(St = 5,            &                  ! Max saturation time step
            Pt = 100,          &                  ! Pressure time step
            ND = 2000)                            ! Number of days in simulation

  ! filenames
  character(*), parameter :: data_directory = "../data/data_1/"
  character(*), parameter :: results_directory = "results/"
  character(*), parameter :: porosity_file = data_directory//"pUr.txt"
  character(*), parameter :: permeability_file = data_directory//"KUr.txt"
  character(*), parameter :: results_eval_original_code = &
       results_directory//"results_eval_original_code.nc"
  character(*), parameter :: results_eval_deriv_tap_1_fwd = &
       results_directory//"results_eval_deriv_tapenade_1_forward.nc"
  character(*), parameter :: results_eval_deriv_tap_1_rev = &
       results_directory//"results_eval_deriv_tapenade_1_reverse.nc"

  ! PARAMETERS READ FROM FILE
  ! porosity and permeability parameters
  double precision, dimension(N_) :: POR      ! Porosities
  double precision, dimension(3, Nx_, Ny_, Nz_) :: PERM  ! Permeabilities  

  ! PARAMETERS SET IN DRIVER
  ! linear solver parameters
  logical :: verbose
  integer :: solver_inner, solver_outer

  parameter(verbose = .false.,   &              ! Verbose solver output
       solver_inner = 64,        &              ! Number of inner iterations
       solver_outer = 100000)                   ! Number of outer iterations
end module parameters
module mathutil
implicit none
public :: dnrm2
contains

subroutine scalar_max(scalarin1, scalarin2, scalarout)
  double precision :: scalarin1, scalarin2, scalarout

  if (scalarin2 <= scalarin1) then
    scalarout = scalarin1
  else
    scalarout = scalarin2
  end if
end subroutine scalar_max

!
! Subroutine computes the 2 norm of the vector
! adapted from the original function version
! of the corresponding blas routine from
! NETLIB
!
subroutine dnrm2(v, len_v, n)
    implicit none
    integer :: i, len_v
    double precision :: n
    double precision, dimension(len_v) :: v
    double precision :: scalein, scaleout

    n = 0.0d0
    scalein = 0.0d0
    scaleout = 0.0d0

    do  i = 1, len_v
        call scalar_max(scalein, abs(v(i)), scaleout)
        scalein = scaleout
    enddo

    if (scaleout .eq. 0.0d0) then
        n = 0.0d0
    else
        do i = 1, len_v
            n = n + (v(i)/scaleout)**2
        enddo

        n = scaleout*sqrt(n)
    endif
end subroutine dnrm2
end module mathutil
module matrix
    use parameters
    implicit none

    ! export module interface
    public ::  myreshape, add_x, spmat_multiply, mymin, mymax !,disp_spmat

    ! interface disp_spmat
    !   module procedure disp_spmat2
    ! end interface disp_spmat

    ! Common interface for all reshape functions
    interface myreshape
      module procedure myreshape_1_2
      module procedure myreshape_2_1
      module procedure myreshape_1_3
      module procedure myreshape_3_1
      module procedure myreshape_1_4
      module procedure myreshape_4_1
    end interface myreshape

    ! Common interface for all add methods to spmat
    interface add_x
      module procedure addx_elem
      module procedure addx_diagonal
    end interface add_x

    ! Common interface to multiply SPMAT with other vectors, diagonal matrices.
    ! SPMAT * SPMAT is currently not implemented as arbitrary fill-ins are
    ! not implemented.
    ! SPMAT * MAT }- These can be implemented by repeatedly calling
    ! MAT * SPMAT }- the vector versions of the method.
    interface spmat_multiply
      module procedure spmat_multiply_diagonal
      module procedure spmat_multiply_vector
      module procedure scalar_multiply_spmat
    end interface

    interface mymin
      module procedure mymin_0_0_double
      module procedure mymin_1_0_double
      module procedure mymin_1_1_double
    end interface mymin

    interface mymax
      module procedure mymax_0_0_double
      module procedure mymax_1_0_double
      module procedure mymax_1_1_double
    end interface mymax
contains

! !
! ! Display the matrix entries
! !
! subroutine disp_spmat(innz, irow_index, irow_compressed, icol_index, ivalues, output)
!     implicit none
!     integer :: k, output
!
!     integer ::  innz
!     integer, dimension(7 * N_) :: irow_index
!     integer, dimension(7 * N_) :: icol_index
!     double precision, dimension(7 * N_) :: ivalues
!     integer, dimension(N_ + 1) :: irow_compressed
!
!     if(output /= 0) then
!         do k = 1, innz
!             write (*, '(a, i7, a, i7, a, a, e23.16)'), "(", irow_index(k), ",", &
!                         icol_index(k), ")", " ", ivalues(k)
!         end do
!         write (*, *) irow_compressed
!     end if
! end subroutine disp_spmat

!
! Subroutine adds x to a particular element.
! This subroutine is a bit flawed because at the time of construction of SPMAT
! , for other than the main diagonal, entries are skipped if they are 0.0 in
! the column matrix. Further the element may be non-existent. Structure of
! SPMAT has to be changed to allow arbitrary fill-ins.
!
subroutine addx_elem(innz, irow_index, irow_compressed,&
                      icol_index, ivalues, x, row, col)
    implicit none
    double precision :: x
    integer :: i, row, col

    integer :: n, innz
    integer, dimension(7 * N_) :: irow_index
    integer, dimension(7 * N_) :: icol_index
    double precision, dimension(7 * N_) :: ivalues
    integer, dimension(N_ + 1) :: irow_compressed

    do i = 1,innz
        if(irow_index(i) == row .and. icol_index(i) == col) then
            ivalues(i) = ivalues(i) + x
            exit
        end if
    end do
end subroutine addx_elem

!
! Subroutine adds x to a particular diagonal
! This subroutine is a bit flawed because at the time of construction of SPMAT
! , for other than the main diagonal, entries are skipped if they are 0.0 in
! the column matrix
!
subroutine addx_diagonal(innz, irow_index, irow_compressed,&
                          icol_index, ivalues, x, diag)
    implicit none
    integer :: i, diag
    double precision :: x

    integer :: innz
    integer, dimension(7 * N_) :: irow_index
    integer, dimension(7 * N_) :: icol_index
    double precision, dimension(7 * N_) :: ivalues
    integer, dimension(N_ + 1) :: irow_compressed

    do i = 1,innz
        if(icol_index(i) - irow_index(i) == diag) then
            ivalues(i) = ivalues(i) + x
        end if
    end do
end subroutine addx_diagonal

!
! Gets the diagonal on which the row, col lie
!
subroutine getdiag(row, col, diag)
    implicit none
    integer :: diag, row, col

    diag = col - row
end subroutine getdiag


!
! Gets the indices of the first element on
! diagonal
!
subroutine firstelm(diag, row, col)
    implicit none
    integer :: diag, row, col

    if (diag < 0) then
        row = 1 - diag
        col = 1
    else if (diag == 0) then
        row = 1
        col = 1
    else
        row = 1
        col = 1 + diag
    end if
end subroutine firstelm


!
! Gets the number of elements on the diagonal for a given output
! matrix size.
!
subroutine noelems(diag, orows, ocols, n)
    implicit none
    integer :: diag, orows, ocols, row, col, n

    ! first element along diagonal
    call firstelm(diag, row, col)

    if (diag < 0) then ! subdiagonal read from top
        ! number of elements along diagonal
        n = orows - row + 1
    else if (diag == 0) then ! main diagonal read from top
        ! number of elements along diagonal
        n = min(orows, ocols)
    else  ! super diagonal read from bottom (possible middle)
        ! number of elements along diagonal
        n = ocols - col + 1
    end if
end subroutine noelems


!
! Reshape a 2d matrix to a 1D array
!
subroutine myreshape_2_1(amatrix, bmatrix)
    implicit none
    integer :: i, j, k
    double precision, dimension(:,:) :: amatrix
    double precision, dimension(:) :: bmatrix

    k = 0

    do j = 1, size(amatrix, 2)
        do i = 1, size(amatrix, 1)
            k = k + 1
            bmatrix(k) = amatrix(i, j)
        end do
    end do
end subroutine myreshape_2_1


!
! Reshape a 1d matrix to a 2D array
!
subroutine myreshape_1_2(amatrix, bmatrix)
    implicit none
    integer :: i, j, k
    double precision, dimension(:) :: amatrix
    double precision, dimension(:,:) :: bmatrix

    k = 0

    do j = 1, size(bmatrix, 2)
        do i = 1, size(bmatrix, 1)
            k = k + 1
            bmatrix(i, j) = amatrix(k)
        end do
    end do
end subroutine myreshape_1_2

!
! Reshape a 3d matrix to a 1D array
!
subroutine myreshape_3_1(amatrix, bmatrix)
    implicit none
    integer :: i, j, k, l
    double precision, dimension(:,:,:) :: amatrix
    double precision, dimension(:) :: bmatrix

    l = 0

    do k = 1, size(amatrix, 3)
        do j = 1, size(amatrix, 2)
            do i = 1, size(amatrix, 1)
                l = l + 1
                bmatrix(l) = amatrix(i, j, k)
            end do
        end do
    end do
end subroutine myreshape_3_1


!
! Reshape a 1d matrix to a 3D array
!
subroutine myreshape_1_3(amatrix, bmatrix)
    implicit none
    integer :: i, j, k, l
    double precision, dimension(:) :: amatrix
    double precision, dimension(:,:,:) :: bmatrix

    l = 0

    do k = 1, size(bmatrix, 3)
        do j = 1, size(bmatrix, 2)
            do i = 1, size(bmatrix, 1)
                l = l + 1
                bmatrix(i, j, k) = amatrix(l)
            end do
        end do
    end do
end subroutine myreshape_1_3


!
! Reshape a 1d matrix to a 4D array
!
subroutine myreshape_4_1(amatrix, bmatrix)
    implicit none
    integer :: i, j, k, l, m
    double precision, dimension(:,:,:,:) :: amatrix
    double precision, dimension(:) :: bmatrix

    m = 0

    do l = 1, size(amatrix, 4)
        do k = 1, size(amatrix, 3)
            do j = 1, size(amatrix, 2)
                do i = 1, size(amatrix, 1)
                    m = m + 1
                    bmatrix(m) = amatrix(i, j, k, l)
                end do
            end do
        end do
    end do
end subroutine myreshape_4_1


!
! Reshape a 1d matrix to a 4D array
!
subroutine myreshape_1_4(amatrix, bmatrix)
    implicit none
    integer :: i, j, k, l, m
    double precision, dimension(:) :: amatrix
    double precision, dimension(:,:,:,:) :: bmatrix

    m = 0

    do l = 1, size(bmatrix, 4)
        do k = 1, size(bmatrix, 3)
            do j = 1, size(bmatrix, 2)
                do i = 1, size(bmatrix, 1)
                    m = m + 1
                    bmatrix(i, j, k, l) = amatrix(m)
                end do
            end do
        end do
    end do
end subroutine myreshape_1_4

!
! This routine pre-multiplies a diagonal matrix by a sparse matrix
!
subroutine spmat_multiply_diagonal(annz, arow_index, arow_compressed,&
                                    acol_index, avalues, dmatrix, &
                                    rnnz, rrow_index, rrow_compressed,&
                                    rcol_index, rvalues, order)
    implicit none
    integer :: i, alloc_err
    character(len = 3) :: order

    integer :: annz
    integer, dimension(7 * N_) :: arow_index
    integer, dimension(7 * N_) :: acol_index
    double precision, dimension(7 * N_) :: avalues
    integer, dimension(N_ + 1) :: arow_compressed

    integer :: rnnz
    integer, dimension(7 * N_) :: rrow_index
    integer, dimension(7 * N_) :: rcol_index
    double precision, dimension(7 * N_) :: rvalues
    integer, dimension(N_ + 1) :: rrow_compressed

    double precision, dimension(N_) :: dmatrix

    rnnz = annz
    rrow_compressed = arow_compressed

    if (order == "PRE") then
        do i = 1,annz
            rrow_index(i) = arow_index(i)
            rcol_index(i) = acol_index(i)
            ! take combination of columns of amatrix
            rvalues(i) = avalues(i) * dmatrix(acol_index(i))
        end do
    else if (order == "POS") then
        do i = 1,annz
            rrow_index(i) = arow_index(i)
            rcol_index(i) = acol_index(i)
            ! take combination of rows of amatrix
            rvalues(i) = avalues(i) * dmatrix(arow_index(i))
        end do
    end if
end subroutine spmat_multiply_diagonal


!
! The routine multiplies a vector by a sparse matrix (PRE/POST)
!
subroutine spmat_multiply_vector(annz, arow_index, arow_compressed, &
                                  acol_index, avalues, bvector, cvector, order)
    implicit none
    integer :: i
    character(len=3):: order

    integer ::  annz
    integer, dimension(7 * N_) :: arow_index
    integer, dimension(7 * N_) :: acol_index
    double precision, dimension(7 * N_) :: avalues
    integer, dimension(N_ + 1) :: arow_compressed

    double precision, dimension(N_) :: bvector
    double precision, dimension(N_) :: cvector

    cvector = 0.0d0

    if (order == "PRE") then
        do i = 1,annz
            ! Combination of the columns of amatrix
            cvector(arow_index(i)) = cvector(arow_index(i)) &
                                + avalues(i) * bvector(acol_index(i))
        end do
    else if (order == "POS") then
        do i = 1,annz
            ! Combination of the rows of amatrix
            cvector(acol_index(i)) = cvector(acol_index(i)) &
                                + avalues(i) * bvector(arow_index(i))
        end do
    end if
end subroutine spmat_multiply_vector

!
! This routine multiplies each element of the SPMAT
! by a scalar.
! Allows amatrix to be the same as rmatrix
!
subroutine scalar_multiply_spmat(annz, arow_index, arow_compressed, &
                                  acol_index, avalues, scalar, &
                                  rnnz, rrow_index, rrow_compressed, &
                                  rcol_index, rvalues)
    implicit none
    integer:: i, alloc_err
    double precision :: scalar

    integer :: annz
    integer, dimension(7 * N_) :: arow_index
    integer, dimension(7 * N_) :: acol_index
    double precision, dimension(7 * N_) :: avalues
    integer, dimension(N_ + 1) :: arow_compressed

    integer :: rnnz
    integer, dimension(7 * N_) :: rrow_index
    integer, dimension(7 * N_) :: rcol_index
    double precision, dimension(7 * N_) :: rvalues
    integer, dimension(N_ + 1) :: rrow_compressed

    rnnz = annz
    rrow_compressed = arow_compressed

    do i = 1,annz
        rrow_index(i) = arow_index(i)
        rcol_index(i) = acol_index(i)
        rvalues(i) = scalar * avalues(i)
    end do
end subroutine scalar_multiply_spmat

subroutine mymin_0_0_double(scalarin1, scalarin2, scalarout)
  double precision :: scalarin1, scalarin2, scalarout

  if (scalarin2 >= scalarin1) then
    scalarout = scalarin1
  else
    scalarout = scalarin2
  end if
end subroutine mymin_0_0_double

subroutine mymin_1_0_double(vectorin, scalarin, vectorout)
  integer :: i
  double precision :: scalarin
  double precision, dimension(:) :: vectorin, vectorout

  do i = 1, size(vectorin, 1)
    if (vectorin(i) >= scalarin) then
      vectorout(i) = scalarin
    else
      vectorout(i) = vectorin(i)
    end if
  end do
end subroutine mymin_1_0_double

subroutine mymin_1_1_double(vectorin1, vectorin2, vectorout)
  integer :: i
  double precision, dimension(:) :: vectorin1, vectorin2, vectorout

  do i = 1, size(vectorin1, 1)
    if (vectorin1(i) >= vectorin2(i)) then
      vectorout(i) = vectorin2(i)
    else
      vectorout(i) = vectorin1(i)
    end if
  end do
end subroutine mymin_1_1_double

subroutine mymax_0_0_double(scalarin1, scalarin2, scalarout)
  double precision :: scalarin1, scalarin2, scalarout

  if (scalarin2 <= scalarin1) then
    scalarout = scalarin1
  else
    scalarout = scalarin2
  end if
end subroutine mymax_0_0_double

subroutine mymax_1_0_double(vectorin, scalarin, vectorout)
  integer :: i
  double precision :: scalarin
  double precision, dimension(:) :: vectorin, vectorout

  do i = 1, size(vectorin, 1)
    if (vectorin(i) <= scalarin) then
      vectorout(i) = scalarin
    else
      vectorout(i) = vectorin(i)
    end if
  end do
end subroutine mymax_1_0_double

subroutine mymax_1_1_double(vectorin1, vectorin2, vectorout)
  integer :: i
  double precision, dimension(:) :: vectorin1, vectorin2, vectorout

  do i = 1, size(vectorin1, 1)
    if (vectorin1(i) <= vectorin2(i)) then
      vectorout(i) = vectorin2(i)
    else
      vectorout(i) = vectorin1(i)
    end if
  end do
end subroutine mymax_1_1_double

end module matrix
module linsolve
use parameters
use mathutil
use matrix


implicit none

! export solve
public :: solve

! interface to linear solvers
interface solve
    module procedure sparse_solve
end interface solve

contains

!
! Calls a specific solver - here the jacobi method
!
subroutine sparse_solve(annz, arow_index, arow_compressed, &
                        acol_index, avalues, b, x)
  implicit none
 
  integer :: annz
  integer, parameter :: matdim = N_
  integer, parameter :: maxlen = 7 * matdim
  integer, dimension(7 * N_) :: arow_index
  integer, dimension(7 * N_) :: acol_index
  double precision, dimension(7 * N_) :: avalues
  integer, dimension(N_ + 1) :: arow_compressed

  double precision, dimension(N_) :: b
  double precision, dimension(N_) :: x

  call sparse_dummy_method(matdim, annz, maxlen, arow_index, arow_compressed,&
                           acol_index, avalues, b, x, solver_inner, solver_outer, verbose)
end subroutine sparse_solve

subroutine sparse_dummy_method(n, annz, alen, arow_index, arow_compressed, &
                               acol_index, avalues, b, x, solver_inner, solver_outer, verbose)
  integer :: i
  logical :: verbose
  integer :: solver_inner, solver_outer
  integer :: n, annz, alen
  integer, dimension(alen) :: arow_index
  integer, dimension(alen) :: acol_index
  double precision, dimension(alen) :: avalues
  integer, dimension(n + 1) :: arow_compressed
  
  double precision, dimension(n) :: b
  double precision, dimension(n) :: x

  call spmat_multiply_vector2(annz, arow_index, arow_compressed, &
                              acol_index, avalues, b, x, "PRE")
end subroutine sparse_dummy_method
end module linsolve
module finitevolume
use parameters
use mathutil
use matrix
use linsolve
implicit none

interface RelPerm
    module procedure RelPerm_scalar
    module procedure RelPerm_vector
end interface RelPerm

contains

!
! Performs Newton Raphson to solve for saturations
!
subroutine NewtRaph(Q, V, S)
  implicit none
  
  integer :: i, j, it
  logical :: converged
  double precision :: dt, dsn

  double precision, dimension(N_) :: S
  double precision, dimension(N_) :: Q
  double precision, dimension(N_) :: S_copy
  double precision, dimension(N_) :: S_iter_copy
  double precision, dimension(N_) :: dtx
  double precision, dimension(N_) :: fi
  double precision, dimension(N_) :: fw
  double precision, dimension(N_) :: Mw
  double precision, dimension(N_) :: Mo
  double precision, dimension(N_) :: dMw
  double precision, dimension(N_) :: dMo
  double precision, dimension(N_) :: dF
  double precision, dimension(N_) :: G
  double precision, dimension(N_) :: dS
  double precision, dimension(N_) :: bfw

  double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V
  
  integer ::  annz
  integer, dimension(7 * N_) :: arow_index
  integer, dimension(7 * N_) :: acol_index
  double precision, dimension(7 * N_) :: avalues
  integer, dimension(N_ + 1) :: arow_compressed

  integer ::  bnnz
  integer, dimension(7 * N_) :: brow_index
  integer, dimension(7 * N_) :: bcol_index
  double precision, dimension(7 * N_) :: bvalues
  integer, dimension(N_ + 1) :: brow_compressed

  integer ::  dgnnz
  integer, dimension(7 * N_) :: dgrow_index
  integer, dimension(7 * N_) :: dgcol_index
  double precision, dimension(7 * N_) :: dgvalues
  integer, dimension(N_ + 1) :: dgrow_compressed

  ! not yet converged
  converged = .false.

  ! Assemble system matrix
  call GenA(V, Q,  annz, arow_index, arow_compressed, acol_index, avalues)

  ! copy S over
  S_copy = S

  ! set scaling factor
  it = 0

  do while(.not. converged)
    dt = (1.0d0 * St)/(2**it)
    dtx = dt/(V_ * POR)

    call mymax_1_0_double(Q, 0.0d0, fi)
    fi = fi * dtx

    ! Matrix-diagonal matrix product
    call spmat_multiply_diagonal( annz, arow_index, arow_compressed,&
                                  acol_index, avalues, dtx, &
                                    bnnz, brow_index, brow_compressed,&
                                  bcol_index, bvalues, "POS")
    i = 0
    
    do while (i < 2**it)
      j = 0
      i = i + 1
      dsn = 1.0d0
      S_iter_copy = S

      do while (dsn > 1.0d-3 .and. j < 10)
        call RelPerm(S, Mw, Mo, dMw, dMo)

        dF = dMw/(Mw + Mo) - (Mw/((Mw + Mo)**2) * (dMw + dMo))

        ! Matrix-diagonal matrix product
        call spmat_multiply_diagonal( bnnz, brow_index, brow_compressed, &
                                      bcol_index, bvalues, dF, &
                                        dgnnz, dgrow_index, dgrow_compressed, &
                                      dgcol_index, dgvalues, "PRE")

        call addx_diagonal( dgnnz, dgrow_index, dgrow_compressed, &
                            dgcol_index, dgvalues, -1.0d0, 0)

        fw = Mw / (Mw + Mo)

        ! Matrix-vector matrix product
        call spmat_multiply_vector( bnnz, brow_index, brow_compressed, &
                                    bcol_index, bvalues, fw, bfw, "PRE")

        G = S - S_iter_copy - bfw - fi

        call solve( dgnnz, dgrow_index, dgrow_compressed, &
                    dgcol_index, dgvalues, G, dS)

        S = S + dS

        call dnrm2(dS, N_, dsn)

        j = j + 1
      end do

      if (dsn > 1.0d-3) then
          i = 2**it               ! Breaks out of while loop.
          S = S_copy
      end if
    end do

    if (dsn < 1.0d-3) then
        converged = .true.
    else
        it = it + 1
    end if
  end do
end subroutine NewtRaph


!
! Pressure Solver
!
subroutine Pres(Q, S, P, V)
  integer i

  double precision, dimension(N_) :: S
  double precision, dimension(N_) :: Q
  double precision, dimension(3 * N_) :: M
  double precision, dimension(Nx_, Ny_, Nz_) :: P
  double precision, dimension(3, Nx_, Ny_, Nz_) :: KM
  double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

  double precision, dimension(N_) :: Mw
  double precision, dimension(N_) :: Mo

  call RelPerm(S, Mw, Mo)

  do i = 1,N_
    M(1 + (i - 1) * 3) = Mw(i) + Mo(i)
    M(2 + (i - 1) * 3) = M(1 + (i - 1) * 3)
    M(3 + (i - 1) * 3) = M(1 + (i - 1) * 3)
  end do

  call myreshape_1_4(M, KM)

  ! point-wise multiply
  KM = KM * PERM

  call TPFA(KM, Q, P, V)
end subroutine Pres


!
! Relative Permeabilities
!
subroutine RelPerm_vector(S, Mw, Mo, dMw, dMo)
    implicit none
    double precision, dimension(N_) :: S
    double precision, dimension(N_) :: Mw
    double precision, dimension(N_) :: Mo
    double precision, dimension(N_) :: S_temp
    double precision, dimension(N_), optional :: dMw
    double precision, dimension(N_), optional :: dMo

    S_temp = (S - swc_)/(1.0d0 - swc_ - sor_)   ! rescale saturation
    Mw = S_temp**2/vw_
    Mo = (1 - S_temp)**2/vo_

    if (present(dMo) .and. present(dMw)) then
        dMw = 2 * S_temp/vw_/(1 - swc_ - sor_)
        dMo = -2 * (1 - S_temp)/vo_/(1 - swc_ - sor_)
    endif
end subroutine RelPerm_vector

!
! Relative Permeabilities
!
subroutine RelPerm_scalar(S, Mw, Mo, dMw, dMo)
    implicit none
    double precision :: S, Mw, Mo, S_temp
    double precision, optional :: dMw, dMo

    S_temp = (S - swc_)/(1.0d0 - swc_ - sor_)   ! rescale saturation
    Mw = S_temp**2/vw_
    Mo = (1 - S_temp)**2/vo_

    if (present(dMo) .and. present(dMw)) then
        dMw = 2 * S_temp/vw_/(1 - swc_ - sor_)
        dMo = -2 * (1 - S_temp)/vo_/(1 - swc_ - sor_)
    endif
end subroutine RelPerm_scalar

!
! Generate A matrix
!
subroutine GenA(V, Q,  annz, arow_index, arow_compressed, acol_index, avalues)
  implicit none

  integer, dimension(7) :: idiags
  double precision, dimension(N_) :: Q
  double precision, dimension(N_, 7) :: diags ! the matrix containing the diagonal entries
  double precision, dimension(N_) :: diag_tmp
  double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V  ! V has an extra length
                                                                  ! across each x, y, z
  double precision, dimension(Nx_, Ny_, Nz_) :: VXYZ

  integer ::  annz
  integer, dimension(7 * N_) :: arow_index
  integer, dimension(7 * N_) :: acol_index
  double precision, dimension(7 * N_) :: avalues
  integer, dimension(N_ + 1) :: arow_compressed

  ! initialize diags
  diags = 0.0d0


  ! reshape arrays first
  VXYZ = V(3,1:Nx_, 1:Ny_, 2:Nz_ + 1)
  call myreshape_3_1(VXYZ, diags(:, 1)) ! z2
  VXYZ = V(2,1:Nx_, 2:Ny_ + 1, 1:Nz_)
  call myreshape_3_1(VXYZ, diags(:, 2)) ! y2
  VXYZ = V(1,2:Nx_ + 1, 1:Ny_, 1:Nz_)
  call myreshape_3_1(VXYZ, diags(:, 3)) ! x2
  VXYZ = V(1,1:Nx_, 1:Ny_, 1:Nz_)
  call myreshape_3_1(VXYZ, diags(:, 5)) ! x1
  VXYZ = V(2,1:Nx_, 1:Ny_, 1:Nz_)
  call myreshape_3_1(VXYZ, diags(:, 6)) ! y1
  VXYZ = V(3,1:Nx_, 1:Ny_, 1:Nz_)
  call myreshape_3_1(VXYZ, diags(:, 7)) ! z1


  diag_tmp = 0.0d0
  call mymax_1_0_double(diags(:,1), 0.0d0, diag_tmp)
  diags(:, 1) = diag_tmp

  diag_tmp = 0.0d0
  call mymax_1_0_double(diags(:,2), 0.0d0, diag_tmp)
  diags(:, 2) = diag_tmp

  diag_tmp = 0.0d0
  call mymax_1_0_double(diags(:,3), 0.0d0, diag_tmp)
  diags(:, 3) = diag_tmp

  diag_tmp = 0.0d0
  call mymin_1_0_double(diags(:,5), 0.0d0, diag_tmp)
  diags(:, 5) = -diag_tmp

  diag_tmp = 0.0d0
  call mymin_1_0_double(diags(:,6), 0.0d0, diag_tmp)
  diags(:, 6) = -diag_tmp

  diag_tmp = 0.0d0
  call mymin_1_0_double(diags(:,7), 0.0d0, diag_tmp)
  diags(:, 7) = -diag_tmp

  diag_tmp = 0.0d0
  call mymin_1_0_double(Q, 0.0d0, diag_tmp)
  diags(:, 4) = diag_tmp - diags(:, 5) - diags(:, 3) &
                          - diags(:, 6) - diags(:, 2) &
                          - diags(:, 7) - diags(:, 1)
  
  call spdiags_fvm_csr(diags, annz, arow_index, arow_compressed,&
                       acol_index, avalues)
end subroutine GenA

!
! Two point flux approximation.
!
subroutine TPFA(K, Q, P, V)
  implicit none
  
  integer :: i
  integer, dimension(7) :: idiags
  double precision, dimension(N_, 7) :: diags ! the matrix containing the diagonal entries

  double precision, dimension(N_) :: Q
  double precision, dimension(Nx_, Ny_, Nz_) :: P
  double precision, dimension(3, Nx_, Ny_, Nz_) :: K
  double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

  ! local variables
  double precision :: tx_, ty_, tz_
  double precision, dimension(Nx_ + 1, Ny_, Nz_) :: TX
  double precision, dimension(Nx_, Ny_ + 1, Nz_) :: TY
  double precision, dimension(Nx_, Ny_, Nz_ + 1) :: TZ
  double precision, dimension(Nx_, Ny_, Nz_) :: TXYZ

  ! solution to the linear system
  double precision, dimension(N_) :: u

  ! point-wise inverse of permeability
  double precision, dimension(3, Nx_, Ny_, Nz_) :: L

  ! sparse matrix
  integer :: annz
  integer, dimension(7 * N_) :: arow_index
  integer, dimension(7 * N_) :: acol_index
  double precision, dimension(7 * N_) :: avalues
  integer, dimension(N_ + 1) :: arow_compressed


  ! get the point-wise inverse of the permeability matrix
  L = 1.0d0/K

  tx_ = 2.0d0 * hy_ * hz_ / hx_
  ty_ = 2.0d0 * hx_ * hz_ / hy_
  tz_ = 2.0d0 * hy_ * hx_ / hz_

  TX = 0.0d0
  TY = 0.0d0
  TZ = 0.0d0

  ! Compute transmissibilities by averaging harmonically
  TX(2:Nx_,1:Ny_,1:Nz_) = tx_/(L(1, 1:Nx_ - 1, 1:Ny_, 1:Nz_) + L(1, 2:Nx_, 1:Ny_, 1:Nz_))
  TY(1:Nx_,2:Ny_,1:Nz_) = ty_/(L(2, 1:Nx_, 1:Ny_ - 1, 1:Nz_) + L(2, 1:Nx_, 2:Ny_, 1:Nz_))
  TZ(1:Nx_,1:Ny_,2:Nz_) = tz_/(L(3, 1:Nx_, 1:Ny_, 1:Nz_ - 1) + L(3, 1:Nx_, 1:Ny_, 2:Nz_))

  ! initialize diags
  diags = 0.0d0

  TXYZ = -TX(1:Nx_,1:Ny_,1:Nz_)
  call myreshape_3_1(TXYZ, diags(:, 5))          ! -x1
  TXYZ = -TY(1:Nx_,1:Ny_,1:Nz_)
  call myreshape_3_1(TXYZ, diags(:, 6))          ! -y1
  TXYZ = -TZ(1:Nx_,1:Ny_,1:Nz_)
  call myreshape_3_1(TXYZ, diags(:, 7))          ! -z1
  TXYZ = -TX(2:Nx_ + 1,1:Ny_,1:Nz_)
  call myreshape_3_1(TXYZ, diags(:, 3))      ! -x2
  TXYZ = -TY(1:Nx_,2:Ny_ + 1,1:Nz_)
  call myreshape_3_1(TXYZ, diags(:, 2))      ! -y2
  TXYZ = -TZ(1:Nx_,1:Ny_,2:Nz_ + 1)
  call myreshape_3_1(TXYZ, diags(:, 1))      ! -z2

  ! Assemble discretization matrix
  diags(:, 4) = -(diags(:,1) + diags(:,2) + diags(:,3) &
                  + diags(:,5) + diags(:,6) + diags(:,7))

  call spdiags_fvm_csr(diags, annz, arow_index, arow_compressed, &
                    acol_index, avalues)


  ! ! Increment the 1,1 element of A
!     call addx_elem(annz, arow_index, arow_compressed,&
!                     acol_index, avalues, &
!                     PERM(1,1,1,1) + PERM(2,1,1,1) + PERM(3,1,1,1), 1, 1)

  ! Fix the pressure at the inlets
  do i = 1,annz
    if(arow_index(i) < Nx_ * Ny_ .and. mod(arow_index(i), Ny_) == 1) then
        if(arow_index(i) == acol_index(i)) then
          avalues(i) = 1
        else
          avalues(i) = 0
        endif
    endif
  enddo

  ! solve the linear system
  ! Pass the rows_index, cols_index, values separately.
  call solve(annz, arow_index, arow_compressed, &
              acol_index, avalues, Q, u)

  ! reshape the solution
  call myreshape_1_3(u, P)

  ! V.x
  V(1, 2:Nx_, 1:Ny_, 1:Nz_) = (P(1:Nx_ - 1, :, :) - P(2:Nx_, :, :)) * TX(2:Nx_,:,:)
  ! V.y
  V(2, 1:Nx_, 2:Ny_, 1:Nz_) = (P(:, 1:Ny_ - 1, :) - P(:, 2:Ny_, :)) * TY(:,2:Ny_,:)
  ! V.z
  V(3, 1:Nx_, 1:Ny_, 2:Nz_) = (P(:, :, 1:Nz_ - 1) - P(:, :, 2:Nz_)) * TZ(:,:,2:Nz_)
end subroutine TPFA

!
! Creates sparse diags matrix from rectangular matrix having the diagonals
! orow_compressed is not populated.
subroutine spdiags_fvm(imatrix, onnz, orow_index, &
                       orow_compressed, ocol_index, ovalues)
  implicit none
  logical :: done
  double precision :: elm
  integer :: i, j, start_row_imatrix, end_row_imatrix, row, col

  double precision, dimension(N_, 7) :: imatrix
  integer, dimension(7) :: idiags 
  integer :: onnz
  integer, dimension(7 * N_) :: orow_index
  integer, dimension(7 * N_) :: ocol_index
  double precision, dimension(7 * N_) :: ovalues
  integer, dimension(N_ + 1) :: orow_compressed

  idiags(1) = -Nx_ * Ny_
  idiags(2) =  -Nx_
  idiags(3) = -1
  idiags(4) = 0
  idiags(5) = 1
  idiags(6) = Nx_
  idiags(7) = Nx_ * Ny_
    
  onnz = 0
  orow_compressed = 0
  
  do i = 1, 7
    if (idiags(i) > 0) then
      start_row_imatrix = idiags(i) + 1
      end_row_imatrix = N_
    else if (idiags(i) <= 0) then
      start_row_imatrix = 1
      end_row_imatrix = N_ + idiags(i)
    end if

    call firstelm(idiags(i), row, col)

    do j = start_row_imatrix,end_row_imatrix
      if (row == col .or. imatrix(j, i) /= 0) then
        onnz = onnz + 1
        orow_index(onnz) = row
        ocol_index(onnz) = col
        ovalues(onnz) = imatrix(j, i)
      end if
      row = row + 1
      col = col + 1
    enddo
  enddo
end subroutine spdiags_fvm

!
! Creates sparse diags matrix from rectangular matrix having the diagonals
!
subroutine spdiags_fvm_csr(imatrix, onnz, orow_index,&
                           orow_compressed, ocol_index, ovalues)
  implicit none
  logical :: done
  double precision :: elm
  integer :: i, j, rownnz

  double precision, dimension(N_, 7) :: imatrix
  integer, dimension(7) :: idiags
  integer, dimension(7) :: start_row_imatrix, end_row_imatrix
  integer, dimension(7) :: row_diag, col_diag      ! row, column along diagonal



  integer :: onnz
  integer, dimension(7 * N_) :: orow_index
  integer, dimension(7 * N_) :: ocol_index
  double precision, dimension(7 * N_) :: ovalues
  integer, dimension(N_ + 1) :: orow_compressed

    idiags(1) = -Nx_ * Ny_
    idiags(2) =  -Nx_
    idiags(3) = -1
    idiags(4) = 0
    idiags(5) = 1
    idiags(6) = Nx_
    idiags(7) = Nx_ * Ny_

  onnz = 0
  orow_compressed(1) = 1                      ! compressed row storage
  
  do i = 1, 7
    if (idiags(i) > 0) then
      start_row_imatrix(i) = idiags(i) + 1
      end_row_imatrix(i) = N_  
    else if (idiags(i) <= 0) then
      start_row_imatrix(i) = 1
      end_row_imatrix(i) = N_ + idiags(i)
    end if
    
    call firstelm(idiags(i), row_diag(i), col_diag(i))
  enddo 
  
  ! Do for each row in imatrix
  do i = 1, N_
    rownnz = 0                        ! count the number of nonzeros in row
    ! Do for each column in imatrix
    do j = 1, 7
    ! Need to check if that column has any entry in this row
      if(row_diag(j) <= i .and.  &    ! checks that you are past the beginning of diagonal, remains fixed
        start_row_imatrix(j) <= end_row_imatrix(j)) then ! checks that you have not exhausted the diagonal
        if(imatrix(start_row_imatrix(j), j) /= 0.0d0) then ! checks that the diagonal entry is non zero
          onnz = onnz + 1
          rownnz = rownnz + 1
          orow_index(onnz) = i
          ocol_index(onnz) = col_diag(j)
          ovalues(onnz) = imatrix(start_row_imatrix(j), j)
        endif
        start_row_imatrix(j) = start_row_imatrix(j) + 1
        col_diag(j) = col_diag(j) + 1
      endif
    enddo
    orow_compressed(i + 1) = orow_compressed(i) + rownnz
  enddo
end subroutine spdiags_fvm_csr
end module finitevolume
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
subroutine read_permeability_and_porosity(PERM, POR)
    integer :: i, j, k, l, m

    double precision, dimension(N_) :: POR                 ! Porosities
    double precision, dimension(3, Nx_, Ny_, Nz_) :: PERM  ! Permeabilities

    double precision, dimension(Nx_, Ny_, Nz_) :: P
    double precision, dimension(maxNx * maxNy * maxNz) :: pUr
    double precision, dimension(3 * maxNx, maxNy * maxNz) :: KUr
    double precision, dimension(3 * maxNx * maxNy * maxNz) :: KUrl

    integer, dimension(Nx_ * Ny_ * Nz_) :: Pindices
    integer, dimension(3 * Nx_ * Ny_ * Nz_) :: Kindices

    ! initialize porosity and permeability to zero
    PERM = 0.0d0
    POR = 0.0d0

    ! read KUr
    open(1,file=permeability_file,status='old')
    read(1,*) ((KUr(i,j), j=1,maxNy * maxNz), i=1,3 * maxNx)
    close(1)

    ! reshape 2 dimension to 1 dimension
    call myreshape_2_1(KUr, KUrl)

    ! select according to specified dimension
    m = 0
    do l = 1, Nz_
        do k = 1,Ny_
            do j = 1,Nx_
                do i = 1,3
                    m = m + 1
                    Kindices(m) = ((l - 1) * (maxNx * maxNy * 3) &
                                  + (k - 1) * (maxNx * 3) &
                                  + 3 * (j-1) + i)
                end do
            end do
        end do
    end do

    ! then reshape 1 dimension to 4 dimension (hack for time being)
    call myreshape_1_4(KUrl(Kindices), PERM)

    ! read KUr
    open(1,file=porosity_file,status='old')
    read(1,*) (pUr(i), i=1,maxNx * maxNy * maxNz)
    close(1)

    m = 0
    do k = 1,Nz_
        do j = 1,Ny_
            do i = 1,Nx_
                m = m + 1
                Pindices(m) = ((k - 1) * (maxNx * maxNy) &
                              + (j - 1) * (maxNx) + i)
            end do
        end do
    end do

    call mymax_1_0_double(pUr(Pindices), 1.0d-3, POR)
end subroutine read_permeability_and_porosity


!
! Initialize inflow and outflow.
!
subroutine init_flw_trnc_norm_xin_pt_out(mu, sigma, Q)
  double precision :: mu, sigma
  double precision, dimension(N_) :: Q

  integer :: i, j
  double precision :: x, pi, pdf, mass
  double precision, dimension(Nx_) :: idx
  double precision, dimension(Nx_) :: Q_x

  ! value of pi
  pi = 3.14159265358979323d0

  !initialize the total mass to 0
  mass = 0.0d0
  Q_x = 0.0d0

  ! Note that the portion of the  Standard Normal distribution between
  ! -3sigma/2 to 3sigma/2 is assumed to fit the 1..Nx where sigma is 1
  do i = 1, Nx_
      ! get the real x coordinate
      x = -1.5d0 + ((i - 1) * 3.0d0)/(Nx_ - 1)    ! Mapping x = [-1.5, 1.5] to Nx_ dimension

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
  do i = 1, Nx_
     Q_x(i) = Q_x(i)/mass * ir
  end do

  ! Assign Q_x to Q
  j = 1
  do i = 1, Nx_* Ny_, Ny_
    Q(i) = Q_x(j)
    j = j + 1
  end do

  ! now set the output
  Q(N_) = -ir
end subroutine init_flw_trnc_norm_xin_pt_out


!
! This subroutine simulates the reservoir
! model.
!
subroutine simulate_reservoir(Q, S, P, V, Tt, Pc, oil)
    use parameters
    
    double precision, dimension(N_) :: Q
    double precision, dimension(N_) :: S
    double precision, dimension(Nx_, Ny_, Nz_) :: P
    double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V
    
    double precision, dimension((ND/St) + 1) :: Tt   
    double precision, dimension(2, (ND/St) + 1) :: Pc
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
    do i = 1, ND/Pt
         do j = 1, Pt/St
            k = k + 1
            
            if (j == 1) then
              call stepforward(1, Q, S, P, V, Mw, Mo)
            else
              call stepforward(0, Q, S, P, V, Mw, Mo)            
            endif

            
            ! update quantites
            Mt = Mw + Mo
            Tt(k) = 1.0d0 * k * St  
            Pc(1,k) = Mw/Mt
            Pc(2,k) = Mo/Mt

            call update_oil(Pc, k, St, tempoil1, tempoil2)
            tempoil1 = tempoil2
        end do
    end do

    oil = tempoil2
end subroutine simulate_reservoir

subroutine stepforward(pressure_step, Q, S, P, V, Mw, Mo)
  integer :: pressure_step
  double precision, dimension(N_) :: Q
  double precision, dimension(N_) :: S
  double precision, dimension(Nx_, Ny_, Nz_) :: P
  double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V
  double precision :: Mw, Mo
  
  if (pressure_step == 1) then
    ! solve pressure
    call Pres(Q, S, P, V)    ! Pressure solver
  endif

  call NewtRaph(Q, V, S)      ! Solve for saturation
  call RelPerm(S(N_), Mw, Mo)         ! Mobilities in well-block
end subroutine stepforward

subroutine update_oil(Pc, k, St, oilin, oilout)
  integer :: St, k
  double precision ::  oilin
  double precision ::  oilout
  double precision, dimension(2, (ND/St) + 1) :: Pc

  oilout = oilin +  Pc(2, k) * St                     ! Reimann sum
end subroutine update_oil

subroutine wrapper(mu, sigma, Q, S, P, V, Tt, Pc, oil) 
  double precision :: mu, sigma
  double precision, dimension(N_) :: Q
  double precision, dimension(N_) :: S
  double precision, dimension(Nx_, Ny_, Nz_) :: P
  double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V
  double precision, dimension((ND/St) + 1) :: Tt
  double precision, dimension(2, (ND/St) + 1) :: Pc
  double precision :: oil

  !$openad independent(mu)
  !$openad independent(sigma)
  call init_flw_trnc_norm_xin_pt_out(mu, sigma, Q)
  call simulate_reservoir(Q, S, P, V, Tt, Pc, oil)
  !$openad dependent(oil)
end subroutine wrapper
end module simulation
