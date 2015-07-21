module sparse_matrix
    implicit none
    !
    ! The sparse matrix format
    ! rows: The number of rows in the full matrix
    ! cols: The number of columns in the full matrix
    ! row_index: The row index of the element
    ! col_index: The column index of the element
    ! values: The actual values at row, column
    ! nnz: The number of non-zeros in the matrix
    !
    type spmat
        integer :: rows, columns, nnz
        integer, dimension(:), pointer :: row_index
        integer, dimension(:), pointer :: col_index
        double precision, dimension(:), pointer :: values
    end type

contains


!
! This routine creates a sparse matrix using the input columns
! in imatrix.
! imatrix: the columns of values to be placed in diagonal matrix
! idiags: the diagonals in which the values would go (sorted)
! orows: the rows in the output matrix
! ocols: the columns in the output matrix
! omatrix: the output sparse matrix in spmat format.
! The placement of values in the diagonals will follow MATLABs
! convention: spdiags
! The structure of the sparse matrix is column-major format
!
subroutine spdiags(imatrix, idiags, orows, ocols, omatrix)
    implicit none
    integer :: orows, ocols, i, j, k, n, d, row, col, nnz, alloc_err
    double precision, dimension(:,:) :: imatrix
    integer, dimension(:) :: idiags
    type(spmat) :: omatrix

! Skipping validation of arguments
! *Note*: imatrix should have as many rows as min(orows, ocols)
!       the number of columns in imatrix should match length of idiags
!       and should not exceed m + n - 1 (soft requirement)
!       the values of idiags should range between -orows + 1 to ocols - 1

! Allocate space with min(orows, ocols) rows and size(idiags,1) as columns
! Ensure omatrix fields are unallocated
    if (allocated(omatrix%row_index) .or. &
       allocated(omatrix%col_index) .or. &
       allocated(omatrix%values)) then
        stop "The output matrix has already been allocated space in the heap."
    end if

    call countnnz(imatrix, idiags, orows, ocols, nnz)

    allocate(omatrix%row_index(nnz), omatrix%col_index(nnz), &
             omatrix%values(nnz), stat = alloc_err)

    if (alloc_err /= 0) then
        stop "Could not allocate memory for the sparse matrix."
    end if

    ! Set the number of non-zeros, rows and columns
    omatrix%nnz = nnz
    omatrix%rows = orows
    omatrix%columns = ocols

    k = 0

    do i = 1, size(idiags, 1)
        d = idiags(i)
! if orows == ocols || orows > ocols, elements of superdiagonal filled from lower
! part of corresponding column in imatrix and from upper part for subdiagonals
        if (orows == ocols .or. orows > ocols) then
            if (d < 0) then ! subdiagonal read from top
                ! first element along diagonal
                row = 1 - d
                col = 1

                ! number of elements along diagonal
                n = orows - row + 1

                do j = 1, n
                    if (imatrix(j, i) /= 0.0d0) then
                        k = k + 1

                        omatrix%row_index(k) = row + j - 1
                        omatrix%col_index(k) = col + j - 1
                        omatrix%values(k) = imatrix(j, i)
                    end if
                end do
            else if (d == 0) then ! main diagonal read from top
                ! first element along diagonal
                row = 1
                col = 1

                ! number of elements along diagonal
                n = min(orows, ocols)

                do j = 1, n
                    if (imatrix(j, i) /= 0.0d0) then
                        k = k + 1

                        omatrix%row_index(k) = row + j - 1
                        omatrix%col_index(k) = col + j - 1
                        omatrix%values(k) = imatrix(j, i)
                    end if
                end do
            else  ! super diagonal read from bottom (possible middle)
                ! first element along diagonal
                row = 1
                col = 1 + d

                ! number of elements along diagonal
                n = ocols - col + 1

                do j = 1, n
                    if (imatrix(orows - j + 1, i) /= 0.0d0) then
                        k = k + 1

                        omatrix%row_index(k) = row + n - j
                        omatrix%col_index(k) = col + n - j
                        omatrix%values(k) = imatrix(j, i)
                    end if
                end do
            end if
        else
    ! if orows < ocols, elements of superdiagonal filled from top part of corresponding
    ! column in imatrix and from lower part for subdiagonals
            if (d < 0) then ! subdiagonal read from bottom (possible middle)
                ! first element along diagonal
                row = 1 - d
                col = 1

                ! number of elements along diagonal
                n = orows - row + 1

                do j = 1, n
                    if (imatrix(orows - j + 1, i) /= 0.0d0) then
                        k = k + 1

                        omatrix%row_index(k) = row + n - j
                        omatrix%col_index(k) = col + n - j
                        omatrix%values(k) = imatrix(j, i)
                    end if
                end do
            else if (d == 0) then ! super diagonal read from top
                ! first element along diagonal
                row = 1
                col = 1

                ! number of elements along diagonal
                n = min(orows, ocols)

                do j = 1, n
                    if (imatrix(j, i) /= 0.0d0) then
                        k = k + 1

                        omatrix%row_index(k) = row + j - 1
                        omatrix%col_index(k) = col + j - 1
                        omatrix%values(k) = imatrix(j, i)
                    end if
                end do
            else
                ! first element along diagonal
                row = 1
                col = 1 + d

                ! number of elements along diagonal
                n = ocols - col + 1

                do j = 1, n
                    if (imatrix(j, i) /= 0.0d0) then
                        k = k + 1

                        omatrix%row_index(k) = row + j - 1
                        omatrix%col_index(k) = col + j - 1
                        omatrix%values(k) = imatrix(j, i)
                    end if
                end do
            end if
        end if
    end do

    ! now sort the entries in imatrix by columns

end subroutine


!
! The subroutine counts the number of non-zeros in the
! matrix
!
subroutine countnnz(imatrix, idiags, orows, ocols, nnz)
    implicit none
    integer :: orows, ocols, i, j, n, d, row, col, nnz
    double precision, dimension(:,:) :: imatrix
    integer, dimension(:) :: idiags

    ! initialize nnz
    nnz = 0

    do i = 1, size(idiags, 1)
        d = idiags(i)
! if orows == ocols || orows > ocols, elements of superdiagonal filled from lower
! part of corresponding column in imatrix and from upper part for subdiagonals
        if (orows == ocols .or. orows > ocols) then
            if (d < 0) then ! subdiagonal read from top
                ! first element along diagonal
                row = 1 - d
                col = 1

                ! number of elements along diagonal
                n = orows - row + 1

                do j = 1, n
                    if (imatrix(j, i) /= 0.0d0) then
                        nnz = nnz + 1
                    end if
                end do
            else if (d == 0) then ! main diagonal read from top
                ! first element along diagonal
                row = 1
                col = 1

                ! number of elements along diagonal
                n = min(orows, ocols)

                do j = 1, n
                    if (imatrix(j, i) /= 0.0d0) then
                        nnz = nnz + 1
                    end if
                end do
            else  ! super diagonal read from bottom (possible middle)
                ! first element along diagonal
                row = 1
                col = 1 + d

                ! number of elements along diagonal
                n = ocols - col + 1

                do j = 1, n
                    if (imatrix(orows - j + 1, i) /= 0.0d0) then
                        nnz = nnz + 1
                    end if
                end do
            end if
        else
    ! if orows < ocols, elements of superdiagonal filled from top part of corresponding
    ! column in imatrix and from lower part for subdiagonals
            if (d < 0) then ! subdiagonal read from bottom (possible middle)
                ! first element along diagonal
                row = 1 - d
                col = 1

                ! number of elements along diagonal
                n = orows - row + 1

                do j = 1, n
                    if (imatrix(orows - j + 1, i) /= 0.0d0) then
                        nnz = nnz + 1
                    end if
                end do
            else if (d == 0) then ! super diagonal read from top
                ! first element along diagonal
                row = 1
                col = 1

                ! number of elements along diagonal
                n = min(orows, ocols)

                do j = 1, n
                    if (imatrix(j, i) /= 0.0d0) then
                        nnz = nnz + 1
                    end if
                end do
            else
                ! first element along diagonal
                row = 1
                col = 1 + d

                ! number of elements along diagonal
                n = ocols - col + 1

                do j = 1, n
                    if (imatrix(j, i) /= 0.0d0) then
                        nnz = nnz + 1
                    end if
                end do
            end if
        end if
    end do
end subroutine countnnz


!
! Display the matrix entries
!
subroutine disp_spmat(imatrix)
    implicit none
    integer :: k
    type(spmat) :: imatrix

    do k = 1, imatrix%nnz
        write (*, '(a, i7, a, i7, a, a, f)'), "(", imatrix%row_index(k), ",", &
                    imatrix%col_index(k), ")", " ", imatrix%values(k)
    end do
end subroutine

!
! Deallocate the entries in the sparse matrix
!
subroutine free_spmat(imatrix)
    implicit none
    integer :: dealloc_err
    type(spmat) :: imatrix


    if (allocated(imatrix%row_index)) then
        deallocate(imatrix%row_index, stat = dealloc_err)

        if (dealloc_err /= 0) then
            stop "Could not deallocate memory for the sparse matrix."
        end if
    end if

    if (allocated(imatrix%col_index)) then
        deallocate(imatrix%col_index, stat = dealloc_err)

        if (dealloc_err /= 0) then
            stop "Could not deallocate memory for the sparse matrix."
        end if
    end if

    if (allocated(imatrix%values)) then
        deallocate(imatrix%values, stat = dealloc_err)

        if (dealloc_err /= 0) then
            stop "Could not deallocate memory for the sparse matrix."
        end if
    end if
end subroutine

!
! Common matrix operations need to be implemented.
! i)  inc/dec the elements of the matrix by some quantity
! ii) Add/Sub two matrices (account for fill in? Heck no.)
! iii) Multiply two matrices - 2 sparse matrices (fill in? Heck no.)
! iv) Sparse linear solve.
!

subroutine addx(imatrix, x)
    implicit none
    integer :: i
    type(spmat) :: imatrix
    double precision :: x

    do i = 1,imatrix%nnz
        imatrix%values(i) = imatrix%values(i) + x
    end do
end subroutine


!
! Returns the next non-zero elm to be filled in the matrix or used
! for counting the non-zeros. Goes column wise in output matrix,
! intersects with idiag, looks at element in imatrix to check if
! non-zero
!
subroutine nextnzelm(imatrix, idiags, orows, ocols, crow, ccol, nrow, ncol, nelm)
    implicit none
    logical :: nnzfound, diagfound
    integer :: maindiagind, maxsubd, minsupd
    integer :: orows, ocols, nrow, ncol, crow, ccol, trow, tcol
    integer ::  i, j, k, l, diag, diagind
    integer, dimension(:) :: idiags
    double precision :: nelm
    double precision, dimension(:,:) :: imatrix

    ! initialize nrow and ncol
    nrow = 0
    ncol = 0

    ! Copy current row and column
    trow = crow
    tcol = ccol

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
            if (maindiagind > 0) then
                call firstelm(0, nrow, ncol)
                call getelm(imatrix, idiags, orows, ocols, nrow, ncol, &
                            maindiagind, nelm)
            else if (maxsubd > -orows) then
                call firstelm(maxsubd, nrow, ncol)
                call getelm(imatrix, idiags, orows, ocols, nrow, ncol, &
                            maindiagind, nelm)
            else if (minsupd < ocols) then
                call firstelm(minsupd, nrow, ncol)
                call getelm(imatrix, idiags, orows, ocols, nrow, ncol, &
                            maindiagind, nelm)
            else
                stop "Incorrect diagonal values."
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
                    call getdiagind(idiags, maindiagind, diag, diagind)

                    if (diagind > 0 .and.  diagind <= size(idiags, 1)) then
                        diagfound = .true.

                        nrow = i
                        ncol = j
                    end if
                end if
            end do

            ! Now get the element if nrow and ncol are > 0
            if (nrow > 0 .and. ncol > 0) then
                call getelm(imatrix, idiags, orows, ocols, nrow, ncol, &
                            maindiagind, nelm)
            end if
            ! if there are no more non-zero elements return nrow ncol as is.
        end if

        if (nelm /= 0.0d0) then
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
subroutine getelm(imatrix, idiags, orows, ocols, row, col, maindiagind, elm)
    implicit none
    integer :: orows, ocols, row, col, srow, scol
    integer :: i, maindiagind, diag, diagind, n
    double precision :: elm
    integer, dimension(:) :: idiags
    double precision, dimension(:,:) :: imatrix

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
        call getdiagind(idiags, maindiagind, diag, diagind)

        ! ensure that the diagonal exists in the list.
        if (diagind > 0 .and. diagind <= size(idiags, 1)) then
                stop "Diagonal not found in the list of diagonals."
        end if

        ! Get the coordinates of the first element along the diagonal to est. ref.
        call firstelm(diag, srow, scol)

        ! follow MATLAB conventions
        if (orows > ocols) then ! subdiagonal read from top
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
        call getdiagind(idiags, maindiagind, diag, diagind)

        ! ensure that the diagonal exists in the list.
        if (diagind > 0 .and. diagind <= size(idiags, 1)) then
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
        else if (orows > ocols) then ! superdiagonal read from bottom
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
subroutine getdiagind(idiags, maindiagind, diag, diagind)
    integer :: i, maindiagind, diag, diagind
    integer, dimension(:) :: idiags

    diagind = 0

    if (diag < 0) then
        i = 1

        ! go over the list of idiags forwards and find the needed diagonal.
        do while (idiags(i) /= diag .and. i <= size(idiags, 1))
            i = i + 1
        end do
    else if (diag > 0) then
        i = size(idiags, 1)

        ! go over the list of idiags backwards and find the needed diagonal.
        do while (idiags(i) /= diag .and. i > 0)
            i = i - 1
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

end module sparse_matrix
