module matrix
    implicit none

    ! export module interface
    public :: zeros, ones, free_mat, pinverse, myreshape, free_mat, add_x

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

    ! Common interface for all the zeros functions
    interface zeros
        module procedure zeros1
        module procedure zeros2
        module procedure zeros3
        module procedure zeros4
    end interface zeros

    ! Common interface for all the ones functions
    interface ones
        module procedure ones1
        module procedure ones2
        module procedure ones3
        module procedure ones4
    end interface ones

    ! Common interface for all the point wise inverse functions
    interface pinverse
        module procedure pinverse1
        module procedure pinverse2
        module procedure pinverse3
        module procedure pinverse4
    end interface pinverse

    ! Common interface for all reshape functions
    interface myreshape
        module procedure myreshape_1_2
        module procedure myreshape_2_1
        module procedure myreshape_1_3
        module procedure myreshape_3_1
        module procedure myreshape_1_4
        module procedure myreshape_4_1
    end interface myreshape

    ! Common interface for all free functions
    interface free_mat
        module procedure free_mat1
        module procedure free_mat2
        module procedure free_mat3
        module procedure free_mat4
        module procedure free_spmat
    end interface free_mat

    ! Common interface for all add methods to spmat
    interface add_x
        module procedure addx_elem
        module procedure addx_diagonal
    end interface add_x

    ! Common interface to multiply SPMAT with other vectors, diagonal matrices.
    ! SPMAT * SPMAT is currently not implemented as arbitrary fill-ins are
    ! not implemented.
    interface spmat_multiply
        module procedure spmat_multiply_diagonal
        module procedure diagonal_multiply_spmat
    end interface

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
! TODO: Extend this to initialize diagonal entries no matter what.
!
subroutine spdiags(imatrix, idiags, orows, ocols, omatrix)
    implicit none
    logical :: done
    integer :: orows, ocols, row, col, nnz, k, alloc_err
    double precision :: elm
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
    if (associated(omatrix%row_index) .or. &
       associated(omatrix%col_index) .or. &
       associated(omatrix%values)) then
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

    ! initialize k and done.
    k = 0
    done = .false.

    ! initialize current row and col
    row = 0
    col = 0

    do while (.not. done)
        call nextnzelm(imatrix, idiags, orows, ocols, row, col, row, col, elm)

        if (row > 0 .and. col > 0) then
            k = k + 1
            omatrix%row_index(k) = row
            omatrix%col_index(k) = col
            omatrix%values(k) = elm
        else
            done = .true.
        end if
    end do
end subroutine


!
! The subroutine counts the number of non-zeros in the
! matrix
!
subroutine countnnz(imatrix, idiags, orows, ocols, nnz)
    implicit none
    logical :: done
    integer :: orows, ocols, row, col, nnz
    integer, dimension(:) :: idiags
    double precision :: elm
    double precision, dimension(:,:) :: imatrix

    ! initialize nnz
    nnz = 0

    ! initialize current row and col
    row = 0
    col = 0

    ! initialize done.
    done = .false.

    do while (.not. done)
        call nextnzelm(imatrix, idiags, orows, ocols, row, col, row, col, elm)

        if (row > 0 .and. col > 0) then
            nnz = nnz + 1
        else
            done = .true.
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
        write (*, '(a, i7, a, i7, a, a, f10.5)'), "(", imatrix%row_index(k), ",", &
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


    if (associated(imatrix%row_index)) then
        deallocate(imatrix%row_index, stat = dealloc_err)

        if (dealloc_err /= 0) then
            stop "Could not deallocate memory for the sparse matrix."
        end if
    end if

    if (associated(imatrix%col_index)) then
        deallocate(imatrix%col_index, stat = dealloc_err)

        if (dealloc_err /= 0) then
            stop "Could not deallocate memory for the sparse matrix."
        end if
    end if

    if (associated(imatrix%values)) then
        deallocate(imatrix%values, stat = dealloc_err)

        if (dealloc_err /= 0) then
            stop "Could not deallocate memory for the sparse matrix."
        end if
    end if
end subroutine


!
! Subroutine adds x to a particular element.
! This subroutine is a bit flawed because at the time of construction of SPMAT
! , for other than the main diagonal, entries are skipped if they are 0.0 in
! the column matrix. Further the element may be non-existent. Structure of
! SPMAT has to be changed to allow arbitrary fill-ins.
!
subroutine addx_elem(imatrix, x, row, col)
    implicit none
    integer :: i, row, col
    type(spmat) :: imatrix
    double precision :: x

    do i = 1,imatrix%nnz
        if(imatrix%row_index(i) == row .and. imatrix%col_index(i) == col) then
            imatrix%values(i) = imatrix%values(i) + x
            exit
        end if
    end do
end subroutine


!
! Subroutine adds x to a particular diagonal
! This subroutine is a bit flawed because at the time of construction of SPMAT
! , for other than the main diagonal, entries are skipped if they are 0.0 in
! the column matrix
!
subroutine addx_diagonal(imatrix, x, diag)
    implicit none
    integer :: i, diag
    type(spmat) :: imatrix
    double precision :: x

    do i = 1,imatrix%nnz
        if(imatrix%col_index(i) - imatrix%row_index(i) == diag) then
            imatrix%values(i) = imatrix%values(i) + x
        end if
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
    integer ::  i, j, diag, diagind
    integer, dimension(:) :: idiags
    double precision :: nelm
    double precision, dimension(:,:) :: imatrix

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
                call getelm(imatrix, idiags, orows, ocols, nrow, ncol, &
                            maindiagind, nelm)
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
                    call getdiagind(idiags, maindiagind, diag, diagind)

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
                    call getelm(imatrix, idiags, orows, ocols, nrow, ncol, &
                            maindiagind, nelm)
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
        call getdiagind(idiags, maindiagind, diag, diagind)

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
subroutine getdiagind(idiags, maindiagind, diag, diagind)
    implicit none
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


!
! Generates Zeros matrix
!
subroutine zeros1(rows, Z)
    implicit none
    integer :: rows, alloc_error
    double precision, dimension(:), pointer :: Z

    if (associated(Z)) then
        stop "The zeros matrix already has been assigned space in the heap."
    endif

    allocate(Z(rows), stat=alloc_error)

    if (alloc_error /= 0) then
        stop "Could not allocate memory for the zeros matrix."
    end if

    Z = 0.0d0
end subroutine


!
! Generates Zeros matrix
!
subroutine zeros2(rows, cols, Z)
    implicit none
    integer :: rows, cols, alloc_error
    double precision, dimension(:, :), pointer :: Z

    if (associated(Z)) then
        stop "The zeros matrix already has been assigned space in the heap."
    endif

    allocate(Z(rows, cols), stat=alloc_error)

    if (alloc_error /= 0) then
        stop "Could not allocate memory for the zeros matrix."
    end if

    Z = 0.0d0
end subroutine

!
! Generates Zeros matrix
!
subroutine zeros3(rows, cols, stacks, Z)
    implicit none
    integer :: rows, cols, stacks, alloc_error
    double precision, dimension(:, :, :), pointer :: Z

    if (associated(Z)) then
        stop "The zeros matrix already has been assigned space in the heap."
    endif

    allocate(Z(rows, cols, stacks), stat=alloc_error)

    if (alloc_error /= 0) then
        stop "Could not allocate memory for the zeros matrix."
    end if

    Z = 0.0d0
end subroutine


!
! Generates Zeros matrix
!
subroutine zeros4(rows, cols, stacks, boxes, Z)
    implicit none
    integer :: rows, cols, stacks, boxes, alloc_error
    double precision, dimension(:, :, :, :), pointer :: Z

    if (associated(Z)) then
        stop "The zeros matrix already has been assigned space in the heap."
    endif

    allocate(Z(rows, cols, stacks, boxes), stat=alloc_error)

    if (alloc_error /= 0) then
        stop "Could not allocate memory for the zeros matrix."
    end if

    Z = 0.0d0
end subroutine


!
! Generates Ones matrix
!
subroutine ones1(rows, O)
    implicit none
    integer :: rows, alloc_error
    double precision, dimension(:), pointer :: O

    if (associated(O)) then
        stop "The ones matrix already has been assigned space in the heap."
    endif

    allocate(O(rows), stat=alloc_error)

    if (alloc_error /= 0) then
        stop "Could not allocate memory for the ones matrix."
    end if

    O = 1.0d0
end subroutine


!
! Generates Ones matrix
!
subroutine ones2(rows, cols, O)
    implicit none
    integer :: rows, cols, alloc_error
    double precision, dimension(:, :), pointer :: O

    if (associated(O)) then
        stop "The ones matrix already has been assigned space in the heap."
    endif

    allocate(O(rows, cols), stat=alloc_error)

    if (alloc_error /= 0) then
        stop "Could not allocate memory for the ones matrix."
    end if

    O = 1.0d0
end subroutine

!
! Generates Ones matrix
!
subroutine ones3(rows, cols, stacks, O)
    implicit none
    integer :: rows, cols, stacks, alloc_error
    double precision, dimension(:, :, :), pointer :: O

    if (associated(O)) then
        stop "The ones matrix already has been assigned space in the heap."
    endif

    allocate(O(rows, cols, stacks), stat=alloc_error)

    if (alloc_error /= 0) then
        stop "Could not allocate memory for the ones matrix."
    end if

    O = 1.0d0
end subroutine


!
! Generates Ones matrix
!
subroutine ones4(rows, cols, stacks, boxes, O)
    implicit none
    integer :: rows, cols, stacks, boxes, alloc_error
    double precision, dimension(:, :, :, :), pointer :: O

    if (associated(O)) then
        stop "The ones matrix already has been assigned space in the heap."
    endif

    allocate(O(rows, cols, stacks, boxes), stat=alloc_error)

    if (alloc_error /= 0) then
        stop "Could not allocate memory for the ones matrix."
    end if

    O = 1.0d0
end subroutine

!
! Matrix pointwise inverse
!
subroutine pinverse1(amatrix, bmatrix)
    implicit none
    double precision, dimension(:), pointer :: amatrix
    double precision, dimension(:), pointer :: bmatrix


    bmatrix = 1.0d0 / amatrix
end subroutine pinverse1

!
! Matrix pointwise inverse
!
subroutine pinverse2(amatrix, bmatrix)
    implicit none
    double precision, dimension(:,:), pointer :: amatrix
    double precision, dimension(:,:), pointer :: bmatrix


    bmatrix = 1.0d0 / amatrix
end subroutine pinverse2

!
! Matrix pointwise inverse
!
subroutine pinverse3(amatrix, bmatrix)
    implicit none
    double precision, dimension(:,:,:), pointer :: amatrix
    double precision, dimension(:,:,:), pointer :: bmatrix


    bmatrix = 1.0d0 / amatrix
end subroutine pinverse3

!
! Matrix pointwise inverse
!
subroutine pinverse4(amatrix, bmatrix)
    implicit none
    double precision, dimension(:,:,:,:), pointer :: amatrix
    double precision, dimension(:,:,:,:), pointer :: bmatrix


    bmatrix = 1.0d0 / amatrix
end subroutine pinverse4

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
! Check that the matrix is associated and free it.
!
subroutine free_mat1(amatrix)
    implicit none
    integer :: dealloc_err
    double precision, dimension(:), pointer :: amatrix

     if (associated(amatrix)) then
        deallocate(amatrix, stat = dealloc_err)

        if (dealloc_err /= 0) then
           stop "Could not deallocate memory for 1D matrix"
        end if
    end if
end subroutine free_mat1

!
! Check that the matrix is associated and free it.
!
subroutine free_mat2(amatrix)
    implicit none
    integer :: dealloc_err
    double precision, dimension(:,:), pointer :: amatrix

     if (associated(amatrix)) then
        deallocate(amatrix, stat = dealloc_err)

        if (dealloc_err /= 0) then
           stop "Could not deallocate memory for 2D matrix"
        end if
    end if
end subroutine free_mat2

!
! Check that the matrix is associated and free it.
!
subroutine free_mat3(amatrix)
    implicit none
    integer :: dealloc_err
    double precision, dimension(:,:,:), pointer :: amatrix

     if (associated(amatrix)) then
        deallocate(amatrix, stat = dealloc_err)

        if (dealloc_err /= 0) then
           stop "Could not deallocate memory for 3D matrix"
        end if
    end if
end subroutine free_mat3

!
! Check that the matrix is associated and free it.
!
subroutine free_mat4(amatrix)
    implicit none
    integer :: dealloc_err
    double precision, dimension(:,:,:,:), pointer :: amatrix

     if (associated(amatrix)) then
        deallocate(amatrix, stat = dealloc_err)

        if (dealloc_err /= 0) then
           stop "Could not deallocate memory for 4D matrix"
        end if
    end if
end subroutine free_mat4


!
! This routine pre-multiplies a diagonal matrix by a sparse matrix
!
subroutine spmat_multiply_diagonal(amatrix, dmatrix, rmatrix)
    implicit none
    integer :: i, alloc_err
    type(spmat) :: amatrix, rmatrix
    double precision, dimension(:) :: dmatrix

    ! Ensure rmatrix fields are unallocated
    if (associated(rmatrix%row_index) .or. &
        associated(rmatrix%col_index) .or. &
        associated(rmatrix%values)) then
        stop "The output matrix has already been allocated space in the heap."
    end if

    allocate(rmatrix%row_index(amatrix%nnz), rmatrix%col_index(amatrix%nnz), &
             rmatrix%values(amatrix%nnz), stat = alloc_err)

    if (alloc_err /= 0) then
        stop "Could not allocate memory for the sparse matrix."
    end if

    rmatrix%nnz = amatrix%nnz
    rmatrix%rows = amatrix%rows
    rmatrix%columns = amatrix%columns

    do i = 1,amatrix%nnz
        rmatrix%row_index(i) = amatrix%row_index(i)
        rmatrix%col_index(i) = amatrix%col_index(i)
        ! take combination of columns of amatrix
        rmatrix%values(i) = amatrix%values(i) * dmatrix(amatrix%col_index(i))
    end do
end subroutine spmat_multiply_diagonal


!
! This routine pre-multiplies a sparse matrix by a diagonal matrix
!
subroutine diagonal_multiply_spmat(dmatrix, amatrix, rmatrix)
    implicit none
    integer :: i, alloc_err
    type(spmat) :: amatrix, rmatrix
    double precision, dimension(:) :: dmatrix

    ! Ensure rmatrix fields are unallocated
    if (associated(rmatrix%row_index) .or. &
        associated(rmatrix%col_index) .or. &
        associated(rmatrix%values)) then
        stop "The output matrix has already been allocated space in the heap."
    end if

    allocate(rmatrix%row_index(amatrix%nnz), rmatrix%col_index(amatrix%nnz), &
             rmatrix%values(amatrix%nnz), stat = alloc_err)

    if (alloc_err /= 0) then
        stop "Could not allocate memory for the sparse matrix."
    end if

    rmatrix%nnz = amatrix%nnz
    rmatrix%rows = amatrix%rows
    rmatrix%columns = amatrix%columns

    do i = 1,amatrix%nnz
        rmatrix%row_index(i) = amatrix%row_index(i)
        rmatrix%col_index(i) = amatrix%col_index(i)
        ! take combination of columns of amatrix
        rmatrix%values(i) = amatrix%values(i) * dmatrix(amatrix%row_index(i))
    end do
end subroutine diagonal_multiply_spmat


end module matrix
