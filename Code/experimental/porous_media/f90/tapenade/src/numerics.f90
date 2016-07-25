module parameters
  ! FIXED PARAMETERS
  ! grid parameters 
  integer :: scenario_id
  double precision :: hx_, hy_, hz_, V_, ir

  ! fluid parameters
  double precision :: vw_, vo_, swc_, sor_

  ! PARAMETERS READ FROM FILE
  ! porosity and permeability parameters
  double precision, dimension(:), allocatable :: POR      ! Porosities
  double precision, dimension(:, :, :, :), allocatable :: PERM  ! Permeabilities  

  ! PARAMETERS SET IN DRIVER
  ! linear solver parameters
  logical :: verbose
  integer :: solver_inner, solver_outer
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
  end interface spmat_multiply

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
  ! subroutine disp_spmat(n, innz, irow_index, irow_compressed, icol_index, ivalues, output)
  !     implicit none
  !     integer :: k, output
  !     integer :: n
  !     integer :: innz
  !     integer, dimension(7 * n) :: irow_index
  !     integer, dimension(7 * n) :: icol_index
  !     double precision, dimension(7 * n) :: ivalues
  !     integer, dimension(n + 1) :: irow_compressed
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
  subroutine addx_elem(n, innz, irow_index, irow_compressed,&
       icol_index, ivalues, x, row, col)
    implicit none
    double precision :: x
    integer :: i, row, col

    integer :: n, innz
    integer, dimension(7 * n) :: irow_index
    integer, dimension(7 * n) :: icol_index
    double precision, dimension(7 * n) :: ivalues
    integer, dimension(n + 1) :: irow_compressed

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
  subroutine addx_diagonal(n, innz, irow_index, irow_compressed,&
       icol_index, ivalues, x, diag)
    implicit none
    integer :: i, diag
    double precision :: x

    integer :: n, innz
    integer, dimension(7 * n) :: irow_index
    integer, dimension(7 * n) :: icol_index
    double precision, dimension(7 * n) :: ivalues
    integer, dimension(n + 1) :: irow_compressed

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
  subroutine spmat_multiply_diagonal(n, annz, arow_index, arow_compressed,&
       acol_index, avalues, dmatrix, &
       rnnz, rrow_index, rrow_compressed,&
       rcol_index, rvalues, order)
    implicit none
    integer :: i, alloc_err
    character(len = 3) :: order
    integer :: n

    integer :: annz
    integer, dimension(7 * n) :: arow_index
    integer, dimension(7 * n) :: acol_index
    double precision, dimension(7 * n) :: avalues
    integer, dimension(n + 1) :: arow_compressed

    integer :: rnnz
    integer, dimension(7 * n) :: rrow_index
    integer, dimension(7 * n) :: rcol_index
    double precision, dimension(7 * n) :: rvalues
    integer, dimension(n + 1) :: rrow_compressed

    double precision, dimension(n) :: dmatrix

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
  subroutine spmat_multiply_vector(n, annz, arow_index, arow_compressed, &
       acol_index, avalues, bvector, cvector, order)
    implicit none
    integer :: i
    character(len=3):: order
    integer :: n
    integer ::  annz
    integer, dimension(7 * n) :: arow_index
    integer, dimension(7 * n) :: acol_index
    double precision, dimension(7 * n) :: avalues
    integer, dimension(n + 1) :: arow_compressed

    double precision, dimension(n) :: bvector
    double precision, dimension(n) :: cvector

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
  subroutine scalar_multiply_spmat(n, annz, arow_index, arow_compressed, &
       acol_index, avalues, scalar, &
       rnnz, rrow_index, rrow_compressed, &
       rcol_index, rvalues)
    implicit none
    integer:: i, alloc_err
    double precision :: scalar
    integer :: n

    integer :: annz
    integer, dimension(7 * n) :: arow_index
    integer, dimension(7 * n) :: acol_index
    double precision, dimension(7 * n) :: avalues
    integer, dimension(n + 1) :: arow_compressed

    integer :: rnnz
    integer, dimension(7 * n) :: rrow_index
    integer, dimension(7 * n) :: rcol_index
    double precision, dimension(7 * n) :: rvalues
    integer, dimension(n + 1) :: rrow_compressed

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
  subroutine sparse_solve(n, annz, arow_index, arow_compressed, &
       acol_index, avalues, b, x)
    implicit none

    integer :: annz
    integer :: n
    integer, dimension(7 * n) :: arow_index
    integer, dimension(7 * n) :: acol_index
    double precision, dimension(7 * n) :: avalues
    integer, dimension(n + 1) :: arow_compressed

    double precision, dimension(n) :: b
    double precision, dimension(n) :: x

    call sparse_pmgmres_method(n, annz, 7 * n, arow_index, arow_compressed,&
         acol_index, avalues, b, x, solver_inner, solver_outer, verbose)
  end subroutine sparse_solve
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
  subroutine NewtRaph(nx, ny, nz, nd, pt, st, Q, V, S)
    implicit none

    integer :: i, j, it
    logical :: converged
    double precision :: dt, dsn
    integer :: nx, ny, nz, n
    integer :: nd, pt, st
    double precision, dimension((nx * ny * nz)) :: S
    double precision, dimension((nx * ny * nz)) :: Q
    double precision, dimension((nx * ny * nz)) :: S_copy
    double precision, dimension((nx * ny * nz)) :: S_iter_copy
    double precision, dimension((nx * ny * nz)) :: dtx
    double precision, dimension((nx * ny * nz)) :: fi
    double precision, dimension((nx * ny * nz)) :: fw
    double precision, dimension((nx * ny * nz)) :: Mw
    double precision, dimension((nx * ny * nz)) :: Mo
    double precision, dimension((nx * ny * nz)) :: dMw
    double precision, dimension((nx * ny * nz)) :: dMo
    double precision, dimension((nx * ny * nz)) :: dF
    double precision, dimension((nx * ny * nz)) :: G
    double precision, dimension((nx * ny * nz)) :: dS
    double precision, dimension((nx * ny * nz)) :: bfw

    double precision, dimension(3, nx + 1, ny + 1, nz + 1) :: V

    integer ::  annz
    integer, dimension(7 * (nx * ny * nz)) :: arow_index
    integer, dimension(7 * (nx * ny * nz)) :: acol_index
    double precision, dimension(7 * (nx * ny * nz)) :: avalues
    integer, dimension((nx * ny * nz) + 1) :: arow_compressed

    integer ::  bnnz
    integer, dimension(7 * (nx * ny * nz)) :: brow_index
    integer, dimension(7 * (nx * ny * nz)) :: bcol_index
    double precision, dimension(7 * (nx * ny * nz)) :: bvalues
    integer, dimension((nx * ny * nz) + 1) :: brow_compressed

    integer ::  dgnnz
    integer, dimension(7 * (nx * ny * nz)) :: dgrow_index
    integer, dimension(7 * (nx * ny * nz)) :: dgcol_index
    double precision, dimension(7 * (nx * ny * nz)) :: dgvalues
    integer, dimension((nx * ny * nz) + 1) :: dgrow_compressed

    ! not yet converged
    converged = .false.
    n = nx * ny * nz

    ! Assemble system matrix
    call GenA(nx, ny, nz, V, Q,  annz, arow_index, arow_compressed, acol_index, avalues)

    ! copy S over
    S_copy = S

    ! set scaling factor
    it = 0

    do while(.not. converged)
       dt = (1.0d0 * st)/(2**it)
       dtx = dt/(V_ * POR)

       call mymax_1_0_double(Q, 0.0d0, fi)
       fi = fi * dtx

       ! Matrix-diagonal matrix product
       call spmat_multiply_diagonal(n, annz, arow_index, arow_compressed,&
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
             call RelPerm(nx, ny, nz, S, Mw, Mo, dMw, dMo)

             dF = dMw/(Mw + Mo) - (Mw/((Mw + Mo)**2) * (dMw + dMo))

             ! Matrix-diagonal matrix product
             call spmat_multiply_diagonal(n, bnnz, brow_index, brow_compressed, &
                  bcol_index, bvalues, dF, &
                  dgnnz, dgrow_index, dgrow_compressed, &
                  dgcol_index, dgvalues, "PRE")

             call addx_diagonal(n, dgnnz, dgrow_index, dgrow_compressed, &
                  dgcol_index, dgvalues, -1.0d0, 0)

             fw = Mw / (Mw + Mo)

             ! Matrix-vector matrix product
             call spmat_multiply_vector(n, bnnz, brow_index, brow_compressed, &
                  bcol_index, bvalues, fw, bfw, "PRE")

             G = S - S_iter_copy - bfw - fi

             call solve(n,  dgnnz, dgrow_index, dgrow_compressed, &
                  dgcol_index, dgvalues, G, dS)

             S = S + dS

             call dnrm2(dS, n, dsn)

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
  subroutine Pres(nx, ny, nz, Q, S, P, V)
    integer :: i
    integer :: nx, ny, nz
    double precision, dimension((nx * ny * nz)) :: S
    double precision, dimension((nx * ny * nz)) :: Q
    double precision, dimension(3 * (nx * ny * nz)) :: M
    double precision, dimension(nx, ny, nz) :: P
    double precision, dimension(3, nx, ny, nz) :: KM
    double precision, dimension(3, nx + 1, ny + 1, nz + 1) :: V

    double precision, dimension((nx * ny * nz)) :: Mw
    double precision, dimension((nx * ny * nz)) :: Mo

    call RelPerm(nx, ny, nz, S, Mw, Mo)

    do i = 1,(nx * ny * nz)
       M(1 + (i - 1) * 3) = Mw(i) + Mo(i)
       M(2 + (i - 1) * 3) = M(1 + (i - 1) * 3)
       M(3 + (i - 1) * 3) = M(1 + (i - 1) * 3)
    end do

    call myreshape_1_4(M, KM)

    ! point-wise multiply
    KM = KM * PERM

    call TPFA(nx, ny, nz, KM, Q, P, V)
  end subroutine Pres


  !
  ! Relative Permeabilities
  !
  subroutine RelPerm_vector(nx, ny, nz, S, Mw, Mo, dMw, dMo)
    implicit none
    integer :: nx, ny, nz
    double precision, dimension((nx * ny * nz)) :: S
    double precision, dimension((nx * ny * nz)) :: Mw
    double precision, dimension((nx * ny * nz)) :: Mo
    double precision, dimension((nx * ny * nz)) :: S_temp
    double precision, dimension((nx * ny * nz)), optional :: dMw
    double precision, dimension((nx * ny * nz)), optional :: dMo

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
  subroutine GenA(nx, ny, nz, V, Q,  annz, arow_index, arow_compressed, acol_index, avalues)
    implicit none
    integer :: nx, ny, nz
    integer, dimension(7) :: idiags
    double precision, dimension((nx * ny * nz)) :: Q
    double precision, dimension((nx * ny * nz), 7) :: diags ! the matrix containing the diagonal entries
    double precision, dimension((nx * ny * nz)) :: diag_tmp
    double precision, dimension(3, nx + 1, ny + 1, nz + 1) :: V  ! V has an extra length
    ! across each x, y, z
    double precision, dimension(nx, ny, nz) :: VXYZ

    integer ::  annz
    integer, dimension(7 * (nx * ny * nz)) :: arow_index
    integer, dimension(7 * (nx * ny * nz)) :: acol_index
    double precision, dimension(7 * (nx * ny * nz)) :: avalues
    integer, dimension((nx * ny * nz) + 1) :: arow_compressed

    ! initialize diags
    diags = 0.0d0


    ! reshape arrays first
    VXYZ = V(3,1:nx, 1:ny, 2:nz + 1)
    call myreshape_3_1(VXYZ, diags(:, 1)) ! z2
    VXYZ = V(2,1:nx, 2:ny + 1, 1:nz)
    call myreshape_3_1(VXYZ, diags(:, 2)) ! y2
    VXYZ = V(1,2:nx + 1, 1:ny, 1:nz)
    call myreshape_3_1(VXYZ, diags(:, 3)) ! x2
    VXYZ = V(1,1:nx, 1:ny, 1:nz)
    call myreshape_3_1(VXYZ, diags(:, 5)) ! x1
    VXYZ = V(2,1:nx, 1:ny, 1:nz)
    call myreshape_3_1(VXYZ, diags(:, 6)) ! y1
    VXYZ = V(3,1:nx, 1:ny, 1:nz)
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

    call spdiags_fvm_csr(nx, ny, nz, diags, annz, arow_index, arow_compressed,&
         acol_index, avalues)
  end subroutine GenA

  !
  ! Two point flux approximation.
  !
  subroutine TPFA(nx, ny, nz, K, Q, P, V)
    implicit none
    integer :: nx, ny, nz, n
    integer :: i
    integer, dimension(7) :: idiags
    double precision, dimension((nx * ny * nz), 7) :: diags ! the matrix containing the diagonal entries

    double precision, dimension((nx * ny * nz)) :: Q
    double precision, dimension(nx, ny, nz) :: P
    double precision, dimension(3, nx, ny, nz) :: K
    double precision, dimension(3, nx + 1, ny + 1, nz + 1) :: V

    ! local variables
    double precision :: tx_, ty_, tz_
    double precision, dimension(nx + 1, ny, nz) :: TX
    double precision, dimension(nx, ny + 1, nz) :: TY
    double precision, dimension(nx, ny, nz + 1) :: TZ
    double precision, dimension(nx, ny, nz) :: TXYZ

    ! solution to the linear system
    double precision, dimension((nx * ny * nz)) :: u

    ! point-wise inverse of permeability
    double precision, dimension(3, nx, ny, nz) :: L

    ! sparse matrix
    integer :: annz
    integer, dimension(7 * (nx * ny * nz)) :: arow_index
    integer, dimension(7 * (nx * ny * nz)) :: acol_index
    double precision, dimension(7 * (nx * ny * nz)) :: avalues
    integer, dimension((nx * ny * nz) + 1) :: arow_compressed

    n = nx * ny * nz

    ! get the point-wise inverse of the permeability matrix
    L = 1.0d0/K

    tx_ = 2.0d0 * hy_ * hz_ / hx_
    ty_ = 2.0d0 * hx_ * hz_ / hy_
    tz_ = 2.0d0 * hy_ * hx_ / hz_

    TX = 0.0d0
    TY = 0.0d0
    TZ = 0.0d0

    ! Compute transmissibilities by averaging harmonically
    TX(2:nx,1:ny,1:nz) = tx_/(L(1, 1:nx - 1, 1:ny, 1:nz) + L(1, 2:nx, 1:ny, 1:nz))
    TY(1:nx,2:ny,1:nz) = ty_/(L(2, 1:nx, 1:ny - 1, 1:nz) + L(2, 1:nx, 2:ny, 1:nz))
    TZ(1:nx,1:ny,2:nz) = tz_/(L(3, 1:nx, 1:ny, 1:nz - 1) + L(3, 1:nx, 1:ny, 2:nz))

    ! initialize diags
    diags = 0.0d0

    TXYZ = -TX(1:nx,1:ny,1:nz)
    call myreshape_3_1(TXYZ, diags(:, 5))          ! -x1
    TXYZ = -TY(1:nx,1:ny,1:nz)
    call myreshape_3_1(TXYZ, diags(:, 6))          ! -y1
    TXYZ = -TZ(1:nx,1:ny,1:nz)
    call myreshape_3_1(TXYZ, diags(:, 7))          ! -z1
    TXYZ = -TX(2:nx + 1,1:ny,1:nz)
    call myreshape_3_1(TXYZ, diags(:, 3))      ! -x2
    TXYZ = -TY(1:nx,2:ny + 1,1:nz)
    call myreshape_3_1(TXYZ, diags(:, 2))      ! -y2
    TXYZ = -TZ(1:nx,1:ny,2:nz + 1)
    call myreshape_3_1(TXYZ, diags(:, 1))      ! -z2

    ! Assemble discretization matrix
    diags(:, 4) = -(diags(:,1) + diags(:,2) + diags(:,3) &
         + diags(:,5) + diags(:,6) + diags(:,7))

    call spdiags_fvm_csr(nx, ny, nz, diags, annz, arow_index, arow_compressed, &
         acol_index, avalues)


    ! ! Increment the 1,1 element of A
    !     call addx_elem(annz, arow_index, arow_compressed,&
    !                     acol_index, avalues, &
    !                     PERM(1,1,1,1) + PERM(2,1,1,1) + PERM(3,1,1,1), 1, 1)

    ! Fix the pressure at the inlets
    do i = 1,annz
       if(arow_index(i) < nx * ny .and. mod(arow_index(i), ny) == 1) then
          if(arow_index(i) == acol_index(i)) then
             avalues(i) = 1
          else
             avalues(i) = 0
          endif
       endif
    enddo

    ! solve the linear system
    ! Pass the rows_index, cols_index, values separately.
    call solve(n, annz, arow_index, arow_compressed, &
         acol_index, avalues, Q, u)

    ! reshape the solution
    call myreshape_1_3(u, P)

    ! V.x
    V(1, 2:nx, 1:ny, 1:nz) = (P(1:nx - 1, :, :) - P(2:nx, :, :)) * TX(2:nx,:,:)
    ! V.y
    V(2, 1:nx, 2:ny, 1:nz) = (P(:, 1:ny - 1, :) - P(:, 2:ny, :)) * TY(:,2:ny,:)
    ! V.z
    V(3, 1:nx, 1:ny, 2:nz) = (P(:, :, 1:nz - 1) - P(:, :, 2:nz)) * TZ(:,:,2:nz)
  end subroutine TPFA

  !
  ! Creates sparse diags matrix from rectangular matrix having the diagonals
  ! orow_compressed is not populated.
  subroutine spdiags_fvm(nx, ny, nz, imatrix, onnz, orow_index, &
       orow_compressed, ocol_index, ovalues)
    implicit none
    logical :: done
    integer :: nx, ny, nz
    double precision :: elm
    integer :: i, j, start_row_imatrix, end_row_imatrix, row, col

    double precision, dimension((nx * ny * nz), 7) :: imatrix
    integer, dimension(7) :: idiags 
    integer :: onnz
    integer, dimension(7 * (nx * ny * nz)) :: orow_index
    integer, dimension(7 * (nx * ny * nz)) :: ocol_index
    double precision, dimension(7 * (nx * ny * nz)) :: ovalues
    integer, dimension((nx * ny * nz) + 1) :: orow_compressed

    idiags(1) = -nx * ny
    idiags(2) = -nx
    idiags(3) = -1
    idiags(4) =  0
    idiags(5) =  1
    idiags(6) =  nx
    idiags(7) =  nx * ny

    onnz = 0
    orow_compressed = 0

    do i = 1, 7
       if (idiags(i) > 0) then
          start_row_imatrix = idiags(i) + 1
          end_row_imatrix = (nx * ny * nz)
       else if (idiags(i) <= 0) then
          start_row_imatrix = 1
          end_row_imatrix = (nx * ny * nz) + idiags(i)
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
  subroutine spdiags_fvm_csr(nx, ny, nz, imatrix, onnz, orow_index,&
       orow_compressed, ocol_index, ovalues)
    implicit none
    integer :: nx, ny, nz
    logical :: done
    double precision :: elm
    integer :: i, j, rownnz

    double precision, dimension((nx * ny * nz), 7) :: imatrix
    integer, dimension(7) :: idiags 
    integer, dimension(7) :: start_row_imatrix, end_row_imatrix
    integer, dimension(7) :: row_diag, col_diag      ! row, column along diagonal


    integer :: onnz
    integer, dimension(7 * (nx * ny * nz)) :: orow_index
    integer, dimension(7 * (nx * ny * nz)) :: ocol_index
    double precision, dimension(7 * (nx * ny * nz)) :: ovalues
    integer, dimension((nx * ny * nz) + 1) :: orow_compressed

    idiags(1) = -nx * ny
    idiags(2) = -nx
    idiags(3) = -1
    idiags(4) =  0
    idiags(5) =  1
    idiags(6) =  nx
    idiags(7) =  nx * ny

    onnz = 0
    orow_compressed(1) = 1                      ! compressed row storage

    do i = 1, 7
       if (idiags(i) > 0) then
          start_row_imatrix(i) = idiags(i) + 1
          end_row_imatrix(i) = (nx * ny * nz)  
       else if (idiags(i) <= 0) then
          start_row_imatrix(i) = 1
          end_row_imatrix(i) = (nx * ny * nz) + idiags(i)
       end if

       call firstelm(idiags(i), row_diag(i), col_diag(i))
    enddo

    ! Do for each row in imatrix
    do i = 1, (nx * ny * nz)
       rownnz = 0                        ! count the number of nonzeros in row
       ! Do for each column in imatrix
       do j = 1, 7
          ! Need to check if that column has any entry in this row
          if(row_diag(j) <= i .and.  &    ! checks that you are past the beginning of diagonal, remains fixed
               start_row_imatrix(j) <= end_row_imatrix(j)) then ! checks that you have not exhausted the diagonal
             if(imatrix(start_row_imatrix(j), j) /= 0.0d0 .or. &
                  idiags(j) == 0) then ! checks that the diagonal entry is non zero or is the main diagonal
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
