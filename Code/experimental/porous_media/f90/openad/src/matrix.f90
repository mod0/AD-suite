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