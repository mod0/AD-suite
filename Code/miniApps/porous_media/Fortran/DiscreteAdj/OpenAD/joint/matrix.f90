module matrix
    implicit none

    ! export module interface
    public :: spdiags,  myreshape, add_x, spmat_multiply, mymin, mymax !,disp_spmat

    interface spdiags
      module procedure spdiags2
    end interface spdiags

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
      module procedure addx_elem2
      module procedure addx_diagonal2
    end interface add_x

    ! Common interface to multiply SPMAT with other vectors, diagonal matrices.
    ! SPMAT * SPMAT is currently not implemented as arbitrary fill-ins are
    ! not implemented.
    ! SPMAT * MAT }- These can be implemented by repeatedly calling
    ! MAT * SPMAT }- the vector versions of the method.
    interface spmat_multiply
      module procedure spmat_multiply_diagonal2
      module procedure spmat_multiply_vector2
      module procedure scalar_multiply_spmat2
    end interface

    interface mymin
      module procedure mymin_1_0_double
      module procedure mymin_1_1_double
    end interface mymin

    interface mymax
      module procedure mymax_1_0_double
      module procedure mymax_1_1_double
    end interface mymax
contains

!
! This routine creates a sparse matrix using the input columns
! in imatrix.
! imatrix: the columns of values to be placed in diagonal matrix
! idiags: the diagonals in which the values would go (sorted)
! orows: the rows in the output matrix (output/input)
! ocols: the columns in the output matrix (output/input)
! onnz: number of nonzeros (output)
! ocols: number of columns (output)
! orow_index: row indices (output)
! ocol_index: col indices (output)
! ovalues: values (output)
! The placement of values in the diagonals will follow MATLABs
! convention: spdiags
! The structure of the sparse matrix is column-major format
! TODO: Extend this to initialize diagonal entries no matter what.
!
subroutine spdiags2(imatrix, imatrix_rows, imatrix_cols, idiags, idiags_count, &
                    orows, ocols, onnz, olen, orow_index, ocol_index, ovalues)
    implicit none
    logical :: done
    double precision :: elm
    integer :: crow, ccol, nrow, ncol, k, alloc_err, orows, ocols
    integer :: imatrix_rows, imatrix_cols, idiags_count, onnz, olen

    integer, dimension(idiags_count) :: idiags
    double precision, dimension(imatrix_rows, imatrix_cols) :: imatrix

    integer, dimension(olen) :: orow_index
    integer, dimension(olen) :: ocol_index
    double precision, dimension(olen) :: ovalues

! Skipping validation of arguments
! *Note*: imatrix should have as many rows as min(orows, ocols)
!       the number of columns in imatrix should match length of idiags
!       and should not exceed m + n - 1 (soft requirement)
!       the values of idiags should range between -orows + 1 to ocols - 1

    call countnnz(imatrix, imatrix_rows, imatrix_cols, idiags, idiags_count, &
                  orows, ocols, onnz)

    ! initialize k and done.
    k = 0
    done = .false.

    ! initialize current row and col
    crow = 0
    ccol = 0

    do while (.not. done)
        call nextnzelm(imatrix, imatrix_rows, imatrix_cols, idiags, idiags_count, &
                       orows, ocols, crow, ccol, nrow, ncol, elm)

        if (nrow > 0 .and. ncol > 0) then
            k = k + 1
            orow_index(k) = nrow
            ocol_index(k) = ncol
            ovalues(k) = elm
            crow = nrow
            ccol = ncol
        else
            done = .true.
        end if
    end do
end subroutine spdiags2


!
! The subroutine counts the number of non-zeros in the
! matrix
!
subroutine countnnz(imatrix, imatrix_rows, imatrix_cols, idiags, idiags_count, &
                    orows, ocols, onnz)
    implicit none
    logical :: done
    double precision :: elm
    integer :: crow, ccol, nrow, ncol, orows, ocols, onnz
    integer :: imatrix_rows, imatrix_cols, idiags_count

    integer, dimension(idiags_count) :: idiags
    double precision, dimension(imatrix_rows, imatrix_cols) :: imatrix

    ! initialize nnz
    onnz = 0

    ! initialize current row and col
    crow = 0
    ccol = 0

    ! initialize done.
    done = .false.

    do while (.not. done)
        call nextnzelm(imatrix, imatrix_rows, imatrix_cols, idiags, idiags_count, &
                       orows, ocols, crow, ccol, nrow, ncol, elm)

        if (nrow > 0 .and. ncol > 0) then
            onnz = onnz + 1
            crow = nrow
            ccol = ncol
        else
            done = .true.
        end if
    end do
end subroutine countnnz

! !
! ! Display the matrix entries
! !
! subroutine disp_spmat2(irows, icols, innz, ilen, irow_index, icol_index, ivalues, output)
!     implicit none
!     integer :: k, output
!
!     integer :: irows, icols, innz, ilen
!     integer, dimension(ilen) :: irow_index
!     integer, dimension(ilen) :: icol_index
!     double precision, dimension(ilen) :: ivalues
!
!     if(output /= 0) then
!         do k = 1, innz
!             write (*, '(a, i7, a, i7, a, a, e23.16)'), "(", irow_index(k), ",", &
!                         icol_index(k), ")", " ", ivalues(k)
!         end do
!     end if
! end subroutine disp_spmat2

!
! Subroutine adds x to a particular element.
! This subroutine is a bit flawed because at the time of construction of SPMAT
! , for other than the main diagonal, entries are skipped if they are 0.0 in
! the column matrix. Further the element may be non-existent. Structure of
! SPMAT has to be changed to allow arbitrary fill-ins.
!
subroutine addx_elem2(irows, icols, innz, ilen, irow_index, icol_index, ivalues, x, row, col)
    implicit none
    double precision :: x
    integer :: i, row, col

    integer :: irows, icols, innz, ilen
    integer, dimension(ilen) :: irow_index
    integer, dimension(ilen) :: icol_index
    double precision, dimension(ilen) :: ivalues

    do i = 1,innz
        if(irow_index(i) == row .and. icol_index(i) == col) then
            ivalues(i) = ivalues(i) + x
            exit
        end if
    end do
end subroutine addx_elem2

!
! Subroutine adds x to a particular diagonal
! This subroutine is a bit flawed because at the time of construction of SPMAT
! , for other than the main diagonal, entries are skipped if they are 0.0 in
! the column matrix
!
subroutine addx_diagonal2(irows, icols, innz, ilen, irow_index, icol_index, ivalues, x, diag)
    implicit none
    integer :: i, diag
    double precision :: x

    integer :: irows, icols, innz, ilen
    integer, dimension(ilen) :: irow_index
    integer, dimension(ilen) :: icol_index
    double precision, dimension(ilen) :: ivalues

    do i = 1,innz
        if(icol_index(i) - irow_index(i) == diag) then
            ivalues(i) = ivalues(i) + x
        end if
    end do
end subroutine addx_diagonal2

!
! Returns the next non-zero elm to be filled in the matrix or used
! for counting the non-zeros. Goes column wise in output matrix,
! intersects with idiag, looks at element in imatrix to check if
! non-zero
!
subroutine nextnzelm(imatrix, imatrix_rows, imatrix_cols, idiags, idiags_count, &
                     orows, ocols, crow, ccol, nrow, ncol, nelm)
    implicit none
    double precision :: nelm
    logical :: nnzfound, diagfound

    integer :: maindiagind, maxsubd, minsupd
    integer :: orows, ocols, nrow, ncol, crow, ccol, trow, tcol
    integer :: i, j, diag, diagind, imatrix_rows, imatrix_cols, idiags_count

    integer, dimension(idiags_count) :: idiags
    double precision, dimension(imatrix_rows, imatrix_cols) :: imatrix

    ! Copy current row and column
    trow = crow
    tcol = ccol

    ! initialize nrow and ncol
    nrow = 0
    ncol = 0

    ! This method was intended to get the next element in a sparse matrix
    ! to be constructed in a column major mode using the imatrix for the
    ! diagonal entries to be placed along idiag diagonals of the resultant
    ! matrix with dimensions orows and ocols. This is extremely complicated
    ! and costly. Switching to a simpler alternative of constructing by
    ! diagonals and sorting the entries by column

    ! haven't found an element yet
    nnzfound = .false.

    ! Find the super diagonal or the max sub diagonal or min super diagonal
    maindiagind = 0
    maxsubd = -orows
    minsupd = ocols

    do i = 1, size(idiags, 1)
        if (idiags(i) == 0) then
            maindiagind = i
        else if (idiags(i) > 0 .and. minsupd > idiags(i)) then
            minsupd = idiags(i)
        else if (idiags(i) < 0 .and. maxsubd < idiags(i)) then
            maxsubd = idiags(i)
        end if
    end do

    do while (.not. nnzfound)
        ! if current element is empty, find the first non-zero element
        if (trow == 0 .and. tcol == 0) then
            ! Making changes to always include the main-diagonal elements.
            nrow = 1
            ncol = 1
            if (maindiagind > 0) then
                call getelm(imatrix, imatrix_rows, imatrix_cols, idiags, idiags_count, &
                            orows, ocols, nrow, ncol, maindiagind, nelm)
            else
                nelm = 0.0d0
            end if
        else if (trow > 0 .and. tcol > 0) then
            ! Given the current row and current column, find the next row, next column
            diagfound = .false.

            ! start with the current row and column
            i = trow
            j = tcol

            do while (.not. diagfound)
                i = i + 1

                if (i > orows) then
                    j = j + 1

                    ! No more diagonals
                    if (j > ocols) then
                        diagfound = .true.

                        nrow = 0
                        ncol = 0
                    else
                        i = 1
                    endif
                end if

                if (i <= orows .and. j <= ocols) then
                    ! get the diagonal
                    call getdiag(i, j, diag)

                    ! check if diagonal exists
                    call getdiagind(idiags, idiags_count, maindiagind, diag, diagind)

                    ! if diagonal exists or is the main diagonal
                    if (diagind > 0 .and.  diagind <= size(idiags, 1)) then
                        diagfound = .true.

                        nrow = i
                        ncol = j
                    else if (i == j) then ! if diagonal does not exist but is MD
                        diagfound = .true.

                        nrow = i
                        ncol = j
                    end if
                end if
            end do

            ! Now get the element if nrow and ncol are > 0
            if (nrow > 0 .and. ncol > 0) then
                if(diagind > 0 .and. diagind <= size(idiags, 1)) then
                    call getelm(imatrix, imatrix_rows, imatrix_cols, idiags, idiags_count, &
                                orows, ocols, nrow, ncol, maindiagind, nelm)
                else
                    nelm = 0.0d0
                end if
            end if
            ! if there are no more non-zero elements return nrow ncol as is.
        end if

        if (nelm /= 0.0d0) then
            nnzfound = .true.
        else if (nelm == 0.0d0 .and. nrow == ncol) then
            nnzfound = .true.
        else
            trow = nrow
            tcol = ncol
        end if
    end do
end subroutine nextnzelm


!
! Given an imatrix, idiag and the output rows and columns and
! a given row, column and diagonal index get the element from
! the imatrix according to the rules followed by MATLAB
! Note should be called only when a diagonal exists for the element row, col
!
subroutine getelm(imatrix, imatrix_rows, imatrix_cols, idiags, idiags_count, &
                  orows, ocols, row, col, maindiagind, elm)
    implicit none
    double precision :: elm

    integer :: i, maindiagind, diag, diagind, n, idiags_count
    integer :: orows, ocols, row, col, srow, scol, imatrix_rows, imatrix_cols

    integer, dimension(idiags_count) :: idiags
    double precision, dimension(imatrix_rows, imatrix_cols) :: imatrix

    ! get the diagonal on which the requested entry lies
    diag = col - row

    if (diag == 0) then
        ! This has to be come here only if the maindiagonal existed in the first
        ! place.
        if(maindiagind > 0 .and. maindiagind <= size(idiags, 1)) then
            diagind = maindiagind
            i = row ! equivalently column
            elm = imatrix(i, diagind)
        else
            stop "Cannot return element from non-existing main diagonal."
        end if
    else if (diag < 0) then
        ! Get the index of the diagonal.
        call getdiagind(idiags, idiags_count, maindiagind, diag, diagind)

        ! ensure that the diagonal exists in the list.
        if (diagind <= 0 .and. diagind > size(idiags, 1)) then
                stop "Diagonal not found in the list of diagonals."
        end if

        ! Get the coordinates of the first element along the diagonal to est. ref.
        call firstelm(diag, srow, scol)

        ! follow MATLAB conventions
        if (orows >= ocols) then ! subdiagonal read from top
            ! The index along the diagonal from the top
            i = row - srow + 1

            ! read from the top
            elm = imatrix(i, diagind)
        else if (orows < ocols) then ! subdiagonal read from bottom
            ! The index along the diagonal from the top
            i = row - srow + 1

            ! get number of elements along the diagonal.
            call noelems(diag, orows, ocols, n)

            ! read from the bottom
            elm = imatrix(min(orows, ocols) - n + i, diagind)
        end if
    else if (diag > 0) then
        ! Get the index of the diagonal.
        call getdiagind(idiags, idiags_count, maindiagind, diag, diagind)

        ! ensure that the diagonal exists in the list.
        if (diagind <= 0 .and. diagind > size(idiags, 1)) then
            stop "Diagonal not found in the list of diagonals."
        end if

        ! Get the coordinates of the first element along the diagonal to est. ref.
        call firstelm(diag, srow, scol)

        ! follow MATLAB conventions
        if (orows < ocols) then ! superdiagonal read from top
            ! The index along the diagonal from the top
            i = row - srow + 1

            ! read from the top
            elm = imatrix(i, diagind)
        else if (orows >= ocols) then ! superdiagonal read from bottom
            ! The index along the diagonal from the top
            i = row - srow + 1

            ! get number of elements along the diagonal.
            call noelems(diag, orows, ocols, n)

            ! read from the bottom
            elm = imatrix(min(orows, ocols) - n + i, diagind)
        end if
    end if
end subroutine getelm


!
! This subroutine returns the index of the given diagonal
! by searching relative to the main diagonal in the idiags
! list.
!
subroutine getdiagind(idiags, idiags_count, maindiagind, diag, diagind)
    implicit none
    integer :: i, maindiagind, diag, diagind, idiags_count

    integer, dimension(idiags_count) :: idiags

    diagind = 0

    if (diag < 0) then
        i = 1

        ! go over the list of idiags forwards and find the needed diagonal.
        do while (i <= size(idiags, 1))
            if(idiags(i) /= diag ) then
                i = i + 1
            else
                exit
            end if
        end do
    else if (diag > 0) then
        i = size(idiags, 1)

        ! go over the list of idiags backwards and find the needed diagonal.
        do while (i > 0)
            if(idiags(i) /= diag) then
                i = i - 1
            else
                exit
            end if
        end do
    else if (diag == 0) then
        i = maindiagind
    end if

    if (i > 0 .and. i <= size(idiags, 1)) then
        diagind = i
    end if
end subroutine


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
end subroutine


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
end subroutine

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
end subroutine


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
end subroutine


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
end subroutine


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
end subroutine

!
! This routine pre-multiplies a diagonal matrix by a sparse matrix
!
subroutine spmat_multiply_diagonal2(arows, acols, annz, alen, arow_index, acol_index, avalues, dmatrix, &
                                    rrows, rcols, rnnz, rlen, rrow_index, rcol_index, rvalues, order)
    implicit none
    integer :: i, alloc_err
    character(len = 3) :: order

    integer :: arows, acols, annz, alen
    integer, dimension(alen) :: arow_index
    integer, dimension(alen) :: acol_index
    double precision, dimension(alen) :: avalues

    integer :: rrows, rcols, rnnz, rlen
    integer, dimension(rlen) :: rrow_index
    integer, dimension(rlen) :: rcol_index
    double precision, dimension(rlen) :: rvalues

    double precision, dimension(arows) :: dmatrix

    rnnz = annz
    rrows = arows
    rcols = acols

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
end subroutine spmat_multiply_diagonal2


!
! The routine multiplies a vector by a sparse matrix (PRE/POST)
!
subroutine spmat_multiply_vector2(arows, acols, annz, alen, arow_index, acol_index, &
                                  avalues, bvector, cvector, order)
    implicit none
    integer :: i
    character(len=3):: order

    integer :: arows, acols, annz, alen
    integer, dimension(alen) :: arow_index
    integer, dimension(alen) :: acol_index
    double precision, dimension(alen) :: avalues

    double precision, dimension(arows) :: bvector
    double precision, dimension(arows) :: cvector

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
end subroutine spmat_multiply_vector2

!
! This routine multiplies each element of the SPMAT
! by a scalar.
! Allows amatrix to be the same as rmatrix
!
subroutine scalar_multiply_spmat2(arows, acols, annz, alen, arow_index, acol_index, avalues, scalar, &
                                  rrows, rcols, rnnz, rlen, rrow_index, rcol_index, rvalues)
    implicit none
    integer:: i, alloc_err
    double precision :: scalar

    integer :: arows, acols, annz, alen
    integer, dimension(alen) :: arow_index
    integer, dimension(alen) :: acol_index
    double precision, dimension(alen) :: avalues

    integer :: rrows, rcols, rnnz, rlen
    integer, dimension(rlen) :: rrow_index
    integer, dimension(rlen) :: rcol_index
    double precision, dimension(rlen) :: rvalues

    rnnz = annz
    rrows = arows
    rcols = acols

    do i = 1,annz
        rrow_index(i) = arow_index(i)
        rcol_index(i) = acol_index(i)
        rvalues(i) = scalar * avalues(i)
    end do
end subroutine scalar_multiply_spmat2


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
