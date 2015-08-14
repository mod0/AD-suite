module fluid
    double precision :: vw_, vo_, swc_, sor_

    parameter(vw_ = 3d-4, &      ! Viscosity
              vo_ = 3d-3, &      ! Viscosity
              swc_ = 0.2d0, &    ! Saturation
              sor_ = 0.2d0)      ! Saturation

end module fluid
module grid
    integer :: Nx_, Ny_, Nz_, N_
    double precision :: hx_, hy_, hz_, V_
    integer :: maxNx, maxNy, maxNz
    parameter(maxNx=60, maxNy=220, maxNz=85)

    parameter(Nx_ = 4, &                   ! Dimension in x-direction
              Ny_ = 4, &                  ! Dimension in y-direction
              Nz_ = 2, &                    ! Dimension in z-direction
              hx_ = 20.0d0 * 0.3048d0, &    ! step size in x-direction
              hy_ = 10.0d0 * 0.3048d0, &    ! step size in y-direction
              hz_ = 2.0d0 * 0.3048d0, &     ! step size in z-direction
              N_ = Nx_ * Ny_ * Nz_, &       ! Total number of grid cells
              V_ = hx_ * hy_ * hz_)         ! Volume of each grid cell

    double precision, dimension(N_) :: Por_      ! Porosities
    double precision, dimension(3, Nx_, Ny_, Nz_) :: K_  ! Permeabilities 
    
    integer, dimension(7) :: diags_
    parameter(diags_ = (/ -Nx_ * Ny_, -Nx_, -1, 0, 1, Nx_, Nx_ * Ny_ /))
end module grid
module mathutil
implicit none
contains
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
    double precision :: scale

    n = 0.0d0
    scale = 0.0d0

    do  i = 1, len_v
        scale = max(scale, abs(v(i)))
    enddo

    if (scale .eq. 0.0d0) then
        n = 0.0d0
    else
        do i = 1, len_v
            n = n + (v(i)/scale)**2
        enddo

        n = scale*sqrt(n)
    endif
end subroutine dnrm2
end module mathutil
module matrix
    implicit none

    ! export module interface
    public :: zeros, ones, free_mat, pinverse, myreshape, add_x, spmat_multiply

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
    ! SPMAT * MAT }- These can be implemented by repeatedly calling
    ! MAT * SPMAT }- the vector versions of the method.
    interface spmat_multiply
        module procedure spmat_multiply_diagonal
        module procedure spmat_multiply_vector
        module procedure scalar_multiply_spmat
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

    alloc_err = 0

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
subroutine disp_spmat(imatrix, output)
    implicit none
    integer :: k, output
    type(spmat) :: imatrix

    if(output /= 0) then
        do k = 1, imatrix%nnz
            write (*, '(a, i7, a, i7, a, a, e23.16)'), "(", imatrix%row_index(k), ",", &
                        imatrix%col_index(k), ")", " ", imatrix%values(k)
        end do
    end if
end subroutine

!
! Deallocate the entries in the sparse matrix
!
subroutine free_spmat(imatrix)
    implicit none
    integer :: dealloc_err
    type(spmat) :: imatrix

    dealloc_err = 0

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
! Generates Zeros matrix
!
subroutine zeros1(rows, Z)
    implicit none
    integer :: rows, alloc_err
    double precision, dimension(:), pointer :: Z

    if (associated(Z)) then
        stop "The zeros matrix already has been assigned space in the heap."
    endif
    
    alloc_err = 0

    allocate(Z(rows), stat=alloc_err)

    if (alloc_err /= 0) then
        stop "Could not allocate memory for the zeros matrix."
    end if

    Z = 0.0d0
end subroutine


!
! Generates Zeros matrix
!
subroutine zeros2(rows, cols, Z)
    implicit none
    integer :: rows, cols, alloc_err
    double precision, dimension(:, :), pointer :: Z

    if (associated(Z)) then
        stop "The zeros matrix already has been assigned space in the heap."
    endif

    alloc_err = 0
    
    allocate(Z(rows, cols), stat=alloc_err)

    if (alloc_err /= 0) then
        stop "Could not allocate memory for the zeros matrix."
    end if

    Z = 0.0d0
end subroutine

!
! Generates Zeros matrix
!
subroutine zeros3(rows, cols, stacks, Z)
    implicit none
    integer :: rows, cols, stacks, alloc_err
    double precision, dimension(:, :, :), pointer :: Z

    if (associated(Z)) then
        stop "The zeros matrix already has been assigned space in the heap."
    endif

    alloc_err = 0
    
    allocate(Z(rows, cols, stacks), stat=alloc_err)

    if (alloc_err /= 0) then
        stop "Could not allocate memory for the zeros matrix."
    end if

    Z = 0.0d0
end subroutine


!
! Generates Zeros matrix
!
subroutine zeros4(rows, cols, stacks, boxes, Z)
    implicit none
    integer :: rows, cols, stacks, boxes, alloc_err
    double precision, dimension(:, :, :, :), pointer :: Z

    if (associated(Z)) then
        stop "The zeros matrix already has been assigned space in the heap."
    endif

    alloc_err = 0
   
    allocate(Z(rows, cols, stacks, boxes), stat=alloc_err)

    if (alloc_err /= 0) then
        stop "Could not allocate memory for the zeros matrix."
    end if

    Z = 0.0d0
end subroutine


!
! Generates Ones matrix
!
subroutine ones1(rows, O)
    implicit none
    integer :: rows, alloc_err
    double precision, dimension(:), pointer :: O

    if (associated(O)) then
        stop "The ones matrix already has been assigned space in the heap."
    endif

    alloc_err = 0
   
    allocate(O(rows), stat=alloc_err)

    if (alloc_err /= 0) then
        stop "Could not allocate memory for the ones matrix."
    end if

    O = 1.0d0
end subroutine


!
! Generates Ones matrix
!
subroutine ones2(rows, cols, O)
    implicit none
    integer :: rows, cols, alloc_err
    double precision, dimension(:, :), pointer :: O

    if (associated(O)) then
        stop "The ones matrix already has been assigned space in the heap."
    endif

    alloc_err = 0
   
    allocate(O(rows, cols), stat=alloc_err)

    if (alloc_err /= 0) then
        stop "Could not allocate memory for the ones matrix."
    end if

    O = 1.0d0
end subroutine

!
! Generates Ones matrix
!
subroutine ones3(rows, cols, stacks, O)
    implicit none
    integer :: rows, cols, stacks, alloc_err
    double precision, dimension(:, :, :), pointer :: O

    if (associated(O)) then
        stop "The ones matrix already has been assigned space in the heap."
    endif

    alloc_err = 0
   
    allocate(O(rows, cols, stacks), stat=alloc_err)

    if (alloc_err /= 0) then
        stop "Could not allocate memory for the ones matrix."
    end if

    O = 1.0d0
end subroutine


!
! Generates Ones matrix
!
subroutine ones4(rows, cols, stacks, boxes, O)
    implicit none
    integer :: rows, cols, stacks, boxes, alloc_err
    double precision, dimension(:, :, :, :), pointer :: O

    if (associated(O)) then
        stop "The ones matrix already has been assigned space in the heap."
    endif

    alloc_err = 0
   
    allocate(O(rows, cols, stacks, boxes), stat=alloc_err)

    if (alloc_err /= 0) then
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

    dealloc_err = 0
    
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

    dealloc_err = 0
    
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

    dealloc_err = 0
    
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

    dealloc_err = 0
    
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
subroutine spmat_multiply_diagonal(amatrix, dmatrix, rmatrix, order)
    implicit none
    integer :: i, alloc_err
    type(spmat) :: amatrix, rmatrix
    double precision, dimension(:) :: dmatrix
    character(len = 3) :: order

    ! Ensure rmatrix fields are unallocated
    if (.not.(associated(rmatrix%row_index) .or. &
              associated(rmatrix%col_index) .or. &
              associated(rmatrix%values))) then
    
        alloc_err = 0
    
        allocate(rmatrix%row_index(amatrix%nnz), rmatrix%col_index(amatrix%nnz), &
                 rmatrix%values(amatrix%nnz), stat = alloc_err)

        if (alloc_err /= 0) then
            stop "Could not allocate memory for the sparse matrix."
        end if

        rmatrix%nnz = amatrix%nnz
        rmatrix%rows = amatrix%rows
        rmatrix%columns = amatrix%columns
    end if

    if (order == "PRE") then
        do i = 1,amatrix%nnz
            rmatrix%row_index(i) = amatrix%row_index(i)
            rmatrix%col_index(i) = amatrix%col_index(i)
            ! take combination of columns of amatrix
            rmatrix%values(i) = amatrix%values(i) * dmatrix(amatrix%col_index(i))
        end do
    else if (order == "POS") then
        do i = 1,amatrix%nnz
            rmatrix%row_index(i) = amatrix%row_index(i)
            rmatrix%col_index(i) = amatrix%col_index(i)
            ! take combination of rows of amatrix
            rmatrix%values(i) = amatrix%values(i) * dmatrix(amatrix%row_index(i))
        end do
    end if
end subroutine spmat_multiply_diagonal



!
! The routine multiplies a vector by a sparse matrix (PRE/POST)
!
subroutine spmat_multiply_vector(amatrix, bvector, cvector, order)
    implicit none
    integer :: i
    type(spmat) :: amatrix
    double precision, dimension(:) :: bvector, cvector
    character(len=3):: order

    cvector = 0.0d0

    if (order == "PRE") then
        do i = 1,amatrix%nnz
            ! Combination of the columns of amatrix
            cvector(amatrix%row_index(i)) = cvector(amatrix%row_index(i)) &
                                + amatrix%values(i) * bvector(amatrix%col_index(i))
        end do
    else if (order == "POS") then
        do i = 1,amatrix%nnz
            ! Combination of the rows of amatrix
            cvector(amatrix%col_index(i)) = cvector(amatrix%col_index(i)) &
                                + amatrix%values(i) * bvector(amatrix%row_index(i))
        end do
    end if
end subroutine spmat_multiply_vector


!
! This routine multiplies each element of the SPMAT
! by a scalar.
! Allows amatrix to be the same as rmatrix
!
subroutine scalar_multiply_spmat(amatrix, scalar, rmatrix)
    implicit none
    integer:: i, alloc_err
    double precision :: scalar
    type(spmat) :: amatrix, rmatrix

    ! Ensure rmatrix fields are unallocated
    if (.not.(associated(rmatrix%row_index) .or. &
              associated(rmatrix%col_index) .or. &
              associated(rmatrix%values))) then

        alloc_err = 0
              
        allocate(rmatrix%row_index(amatrix%nnz), rmatrix%col_index(amatrix%nnz), &
             rmatrix%values(amatrix%nnz), stat = alloc_err)

        if (alloc_err /= 0) then
            stop "Could not allocate memory for the sparse matrix."
        end if

        rmatrix%nnz = amatrix%nnz
        rmatrix%rows = amatrix%rows
        rmatrix%columns = amatrix%columns
    end if

    do i = 1,amatrix%nnz
        rmatrix%row_index(i) = amatrix%row_index(i)
        rmatrix%col_index(i) = amatrix%col_index(i)
        rmatrix%values(i) = scalar * amatrix%values(i)
    end do
end subroutine

end module matrix
module linsolve
use matrix
use mathutil

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
subroutine sparse_solve(A, b, x)
    implicit none
    integer :: i
    type(spmat) :: A
    double precision, dimension(:) :: b
    double precision, dimension(:) :: x
    if(1 .ne. 0) then
        x = 0.0d0
    end if
    
    !call sparse_mgmres_method(A, b, x)
       
    ! dummy relationship
    do i = 1, A%nnz
       x(i) = A%values(i) * b(A%row_index(i))
    end do
end subroutine sparse_solve

!
! The subroutine assumes the diagonals are non-zero
!
subroutine sparse_jacobi_method(A, b, x)
    implicit none
    integer :: i, j
    type(spmat) :: A
    double precision :: nrm
    double precision, dimension(A%rows) :: b, x, x_old, main_diag


    do i = 1, 100000      ! Max iterations

        x_old = x
        ! start x with right hand side
        x = b

!$omp parallel
        do j = 1, A%nnz
            ! Not the main diagonal
            if (A%row_index(j) /= A%col_index(j)) then
                ! update the new iterate at the corresponding row
                x(A%row_index(j)) = x(A%row_index(j)) &
                                    - A%values(j) * x_old(A%col_index(j))
            else
                ! place the entry in the corresponding diagonal.
                main_diag(A%row_index(j)) = A%values(j)
            end if
        end do
!$omp end parallel

        ! divide by the diagonal entries
        x = x/main_diag

        ! get the norm
        call dnrm2((x - x_old), A%rows, nrm)

        ! Exit when converged.
        if (nrm < 1.0d-8) then
            !print *, "Converged in norm"
            exit
        end if
    end do

    if (i > 5000) then
        print *, "Norm: ", nrm
    end if
end subroutine sparse_jacobi_method


!
! The subroutine assumes the diagonals are non-zero
!
subroutine sparse_gauss_seidel_method(A, b, x)
    implicit none
    integer :: i, j, k
    type(spmat) :: A
    double precision :: nrm
    double precision, dimension(A%rows) :: b, x, x_old, main_diag

    do k = 1, 100000      ! Max iterations
        ! start with old x
        x_old = x

        ! start x with right hand side
        x = b

        do j = 1, A%nnz
            ! Not the main diagonal
            ! This is same as Jacobi (upper triangle part of A)
            if (A%row_index(j) < A%col_index(j)) then
                ! update the new iterate at the corresponding row
                x(A%row_index(j)) = x(A%row_index(j)) &
                                    - A%values(j) * x_old(A%col_index(j))
            else
                ! place the entry in the corresponding diagonal.
                main_diag(A%row_index(j)) = A%values(j)
            end if
        end do

        do i = 1,A%rows                 ! An artificial bound for the row
            do j = 1,A%nnz
                if(A%row_index(j) == i .and. A%col_index(j) < i) then ! 1st part not required.
                    ! Use the new solution upto row < i
                    x(A%row_index(j)) = x(A%row_index(j)) &
                                    - A%values(j) * x(A%col_index(j))
                else if(A%col_index(j) >= i) then
                    ! exit when the column is more than or equal i
                    exit
                end if
            end do

            ! Divide the ith component of soln by the ith component of maindiag
            x(i) = x(i)/main_diag(i)
        end do

        ! get the norm
        call dnrm2((x - x_old), A%rows, nrm)

        ! Exit when converged.
        if (nrm < 1.0d-8) then
!            print *, "Iterations: ", k
!            print *, "Converged in norm"
            exit
        end if
    end do

    if (k > 5000) then
        print *, "Norm: ", nrm
    end if
end subroutine sparse_gauss_seidel_method

!
! A wrapper for mgmres
!
subroutine sparse_mgmres_method(A, b, x)
    !use mgmres
    implicit none
    type(spmat) :: A
    integer :: itr_max, mr
    double precision :: tol_abs, tol_rel
    double precision, dimension(:) :: b, x


    tol_abs = 1.0d-6
    tol_rel = 1.0d-6
    itr_max = 5
    mr = (A%rows)

    !call mgmres_st ( A%rows, A%nnz, A%row_index, A%col_index, A%values, x, b, &
    !                itr_max, mr, tol_abs, tol_rel )

end subroutine sparse_mgmres_method


end module linsolve
module fvm
use grid
use fluid
use matrix
use linsolve
use mathutil

implicit none

interface RelPerm
    module procedure RelPerm_scalar
    module procedure RelPerm_vector
end interface RelPerm

integer :: output
parameter(output = 0)

contains

!
! Performs Newton Raphson to solve for saturations
!
subroutine NewtRaph(S, V, Q, St)
    !use print_active
    implicit none
    double precision :: St                             ! Maximum saturation Time Step
    double precision, dimension(N_) :: S, Q
    double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

    integer :: i, j, it
    type(spmat) :: A, B, dG
    logical :: converged
    double precision :: dt, dsn
    double precision, dimension(N_) :: S_copy, S_iter_copy, dtx, fi, fw, Mw, &
                                       Mo, dMw, dMo, dF, G, dS, bfw

    ! nullify the fields of A, B, dG
    nullify(A%row_index, A%col_index, A%values, &
            B%row_index, B%col_index, B%values, &
            dG%row_index, dG%col_index, dG%values)

    ! not yet converged
    converged = .false.

    ! Assemble system matrix
    call GenA(V, Q, A)

    ! copy S over
    S_copy = S

    ! set scaling factor
    it = 0

    ! print A
    !call disp_spmat(A, output)

    do while(.not. converged)
        dt = St/(2**it)
        dtx = dt/(V_ * Por_)
        fi = max(Q, 0.0d0) * dtx

        !call print_array(dtx,1,0,output)

        call spmat_multiply(A, dtx, B, "POS")

        !call disp_spmat(B, output)

        i = 0

        ! This loop is very badly implemented in the original MATLAB code
        do while (i < 2**it)
            j = 0
            i = i + 1
            dsn = 1.0d0
            S_iter_copy = S

            do while (dsn > 1.0d-3 .and. j < 10)
                call RelPerm(S, Mw, Mo, dMw, dMo)
                dF = dMw/(Mw + Mo) - (Mw/((Mw + Mo)**2) * (dMw + dMo))

                !call print_array(dF,1,0,output)

                call spmat_multiply(B, dF, dG, "PRE")

                !call disp_spmat(dG, output)

                call addx_diagonal(dG, -1.0d0, 0)

                !call disp_spmat(dG, output)

                fw = Mw / (Mw + Mo)
                call spmat_multiply(B, fw, bfw, "PRE")

                !call print_array(bfw,1,0,output)

                G = S - S_iter_copy - bfw - fi

                !call print_array(G,1,0,output)

                call solve(dG, G, dS)
                S = S + dS

                !call print_array(dS,1,0,output)

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

    call free_mat(A)
    call free_mat(B)
    call free_mat(dG)
end subroutine NewtRaph

!
! Pressure Solver
!
subroutine Pres(S, Q, P, V)
    !use print_active
    double precision, dimension(N_) :: S, Q
    double precision, dimension(Nx_, Ny_, Nz_) :: P
    double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

    double precision, dimension(3 * N_), target :: M
    double precision, dimension(3, Nx_, Ny_, Nz_) :: KM
    
    double precision, dimension(:), pointer :: M1
    double precision, dimension(:), pointer :: M2
    double precision, dimension(:), pointer :: M3
    
    M1 => M(1 : 3 * N_ : 3)
    M2 => M(2 : 3 * N_ : 3)
    M3 => M(3 : 3 * N_ : 3)
    
    call RelPerm(S, M1, M2)

    !call print_array(M, 1, 0,output)

    M1 = M1 + M2
    M2 = M1
    M3 = M1

    !call print_array(M, 1, 0,output)

    call myreshape(M, KM)

    !call print_array(KM, 1,0,1,0,1,0,1,0,output)

    !call print_array(K_, 1,0,1,0,1,0,1,0,output)

    ! point-wise multiply
    KM = KM * K_

    !call print_array(KM, 1,0,1,0,1,0,1,0,output)
    !call print_array(Q, 1,0,output)

    call tpfa(KM, Q, P, V)

    !call print_array(V, 1,0,1,0,1,0,1,0,output)
    !call print_array(P, 1,0,1,0,1,0,output)
end subroutine


!
! Relative Permeabilities
!
subroutine RelPerm_vector(S, Mw, Mo, dMw, dMo)
    !use print_active
    implicit none
    double precision, dimension(N_) :: S, Mw, Mo
    double precision, dimension(N_), optional :: dMw, dMo

    double precision, dimension(size(S, 1)) :: S_

    !call print_array(S, 1, 0,output)

    S_ = (S - swc_)/(1.0d0 - swc_ - sor_)   ! rescale saturation

    !call print_array(S_, 1, 0,output)

    Mw = S_**2/vw_
    Mo = (1 - S_)**2/vo_

    !call print_array(Mw, 1, 0,output)
    !call print_array(Mo, 1, 0,output)

    if (present(dMo) .and. present(dMw)) then
        dMw = 2 * S_/vw_/(1 - swc_ - sor_)
        dMo = -2 * (1 - S_)/vo_/(1 - swc_ - sor_)
    endif
end subroutine RelPerm_vector

!
! Relative Permeabilities
!
subroutine RelPerm_scalar(S, Mw, Mo, dMw, dMo)
    implicit none
    double precision :: S, Mw, Mo
    double precision, optional :: dMw, dMo

    double precision :: S_

    S_ = (S - swc_)/(1.0d0 - swc_ - sor_)   ! rescale saturation
    Mw = S_**2/vw_
    Mo = (1 - S_)**2/vo_

    if (present(dMo) .and. present(dMw)) then
        dMw = 2 * S_/vw_/(1 - swc_ - sor_)
        dMo = -2 * (1 - S_)/vo_/(1 - swc_ - sor_)
    endif
end subroutine RelPerm_scalar

!
! Generate A matrix
!
subroutine GenA(V, Q, A)
    !use print_active
    implicit none

    type(spmat) :: A
    double precision, dimension(N_) :: Q
    double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V  ! V has an extra length
                                                                    ! across each x, y, z

    ! the matrix containing the diagonal entries
    double precision, dimension(N_, 7) :: diags

    !call print_array(V, 1,0,1,0,1,0,1,0,output)

    ! reshape arrays first
    call myreshape(V(3,1:Nx_, 1:Ny_, 2:Nz_ + 1), diags(:, 1)) ! z2
    call myreshape(V(2,1:Nx_, 2:Ny_ + 1, 1:Nz_), diags(:, 2)) ! y2
    call myreshape(V(1,2:Nx_ + 1, 1:Ny_, 1:Nz_), diags(:, 3)) ! x2
    call myreshape(V(1,1:Nx_, 1:Ny_, 1:Nz_), diags(:, 5)) ! x1
    call myreshape(V(2,1:Nx_, 1:Ny_, 1:Nz_), diags(:, 6)) ! y1
    call myreshape(V(3,1:Nx_, 1:Ny_, 1:Nz_), diags(:, 7)) ! z1

    !call print_array(diags,1,0,1,0,output)

    diags(:, 1) = max(diags(:,1), 0.0d0)
    diags(:, 2) = max(diags(:,2), 0.0d0)
    diags(:, 3) = max(diags(:,3), 0.0d0)
    diags(:, 5) = -min(diags(:,5), 0.0d0)
    diags(:, 6) = -min(diags(:,6), 0.0d0)
    diags(:, 7) = -min(diags(:,7), 0.0d0)
    diags(:, 4) = min(Q, 0.0d0) - diags(:, 5) - diags(:, 3) &
                                - diags(:, 6) - diags(:, 2) &
                                - diags(:, 7) - diags(:, 1)

    !call print_array(diags,1,0,1,0,output)

    ! This can be sped up by passing 3 arrays having rows, cols and diagind
    ! this can be done because diag positions are fixed.
    call spdiags(diags, diags_, N_, N_, A)
end subroutine GenA


!
! Two point flux approximation.
!
subroutine tpfa(K, Q, P, V)
    !use print_active
    implicit none
    double precision, dimension(N_) :: Q
    double precision, dimension(Nx_, Ny_, Nz_) :: P
    double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V
    double precision, dimension(3, Nx_, Ny_, Nz_) :: K

    ! local variables
    double precision :: tx_, ty_, tz_
    double precision, dimension(Nx_ + 1, Ny_, Nz_) :: TX
    double precision, dimension(Nx_, Ny_ + 1, Nz_) :: TY
    double precision, dimension(Nx_, Ny_, Nz_ + 1) :: TZ

    ! the matrix containing the diagonal entries
    double precision, dimension(N_, 7) :: diags

    ! solution to the linear system
    double precision, dimension(N_) :: u

    ! point-wise inverse of permeability
    double precision, dimension(3,Nx_,Ny_,Nz_) :: L

    ! sparse matrix
    type(spmat) :: A
    nullify(A%row_index)
    nullify(A%col_index)
    nullify(A%values)

    ! get the point-wise inverse of the permeability matrix
    L = 1.0d0/K

    !call print_array(K,1,0,1,0,1,0,1,0,output)
    !call print_array(L,1,0,1,0,1,0,1,0,output)

    tx_ = 2.0d0 * hy_ * hz_ / hx_
    ty_ = 2.0d0 * hx_ * hz_ / hy_
    tz_ = 2.0d0 * hy_ * hx_ / hz_

    TX = 0.0d0
    TY = 0.0d0
    TZ = 0.0d0

    ! Compute transmissibilities by averaging harmonically
    TX(2:Nx_,:,:) = tx_/(L(1, 1:Nx_ - 1, :, :) + L(1, 2:Nx_, :, :))
    TY(:,2:Ny_,:) = ty_/(L(2, :, 1:Ny_ - 1, :) + L(2, :, 2:Ny_, :))
    TZ(:,:,2:Nz_) = tz_/(L(3, :, :, 1:Nz_ - 1) + L(3, :, :, 2:Nz_))

    !call print_array(TX,1,0,1,0,1,0,output)
    !call print_array(TY,1,0,1,0,1,0,output)
    !call print_array(TZ,1,0,1,0,1,0,output)

    call myreshape(-TX(1:Nx_,:,:), diags(:, 5))          ! -x1
    call myreshape(-TY(:,1:Ny_,:), diags(:, 6))          ! -y1
    call myreshape(-TZ(:,:,1:Nz_), diags(:, 7))          ! -z1
    call myreshape(-TX(2:Nx_ + 1,:,:), diags(:, 3))      ! -x2
    call myreshape(-TY(:,2:Ny_ + 1,:), diags(:, 2))      ! -y2
    call myreshape(-TZ(:,:,2:Nz_ + 1), diags(:, 1))      ! -z2

    !call print_array(diags,1,0,1,0,output)

    ! Assemble discretization matrix
    diags(:, 4) = -(diags(:,1) + diags(:,2) + diags(:,3) &
                    + diags(:,5) + diags(:,6) + diags(:,7))

    !call print_array(diags,1,0,1,0,output)

    call spdiags(diags, diags_, N_, N_, A)

    !call disp_spmat(A, output)

    ! Increment the 1,1 element of A
    call addx_elem(A, K_(1,1,1,1) + K_(2,1,1,1) + K_(3,1,1,1), 1, 1)

    !call print_array(K_,1,0,1,0,1,0,1,0,output)
    !call disp_spmat(A, output)


    !call print_array(Q,1,0,output)



    ! solve the linear system
    call solve(A, Q, u )

    !call print_array(u,1,0,output)


    !call print_array(P,1,0,1,0,1,0,output)

    ! reshape the solution
    call myreshape(u, P)

    !call print_array(P,1,0,1,0,1,0,output)

    !call print_array(V,1,0,1,0,1,0,1,0,output)

    ! V.x
    V(1, 2:Nx_, 1:Ny_, 1:Nz_) = (P(1:Nx_ - 1, :, :) - P(2:Nx_, :, :)) * TX(2:Nx_,:,:)
    ! V.y
    V(2, 1:Nx_, 2:Ny_, 1:Nz_) = (P(:, 1:Ny_ - 1, :) - P(:, 2:Ny_, :)) * TY(:,2:Ny_,:)
    ! V.z
    V(3, 1:Nx_, 1:Ny_, 2:Nz_) = (P(:, :, 1:Nz_ - 1) - P(:, :, 2:Nz_)) * TZ(:,:,2:Nz_)

    !call print_array(V,1,0,1,0,1,0,1,0,output)

    ! free matrices
    call free_mat(A)
end subroutine tpfa

end module fvm
program runspe10
use fvm
use grid
use fluid
use matrix
!use gnufor2
!use print_active

implicit none

integer :: St, Pt, ND
parameter(St = 5,&       ! Max saturation time step
          Pt = 100,&     ! Pressure time step
          ND = 2000)    ! Number of days in simulation
          
double precision :: ir, Mw, Mo, Mt, mu, sigma, oil
double precision, dimension((ND/St) + 1) :: Tt
double precision, dimension(N_) :: S
double precision, dimension(N_) :: Q
double precision, dimension(2, (ND/St) + 1) :: Pc
double precision, dimension(Nx_, Ny_, Nz_) :: P
double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

! setup memory for P and V
P = 0.0d0

! note that V vector has an additional
! length in each dimension x,y,z
V = 0.0d0

Tt = 0.0d0               ! simulation time
Pc = 0.0d0               ! production data

! setup memory for inflow and saturation
Q = 0.0d0
S = 0.0d0

! Initialize permeability and porosity for testing
K_ = 0.0d0
Por_ = 0.0d0

! Now read the permeabilities and porosities
call inputKP()                      

! compute ir
ir = (795.0 * Nx_ * Ny_ * Nz_) / (maxNx * maxNy * maxNz)

!Q(1:N_:Nx_*Ny_) = ir
!Q(Nx_*Ny_:N_:Nx_*Ny_) = -ir
!Q(1) = ir
!Q(N_) = -ir

!$openad independent(mu)
!$openad independent(sigma)
!$openad independent(sigma)

! initialize mu, sigma, oil
mu = 0.0d0
sigma = 1.0d0
oil = 0.0d0

call run_simulation(ir, mu, sigma, Q, oil)


!! Set filename to collect results.
!open(unit = 1, file = 'Pc1', &
!        form = 'unformatted', access = 'stream')
!
!! write the production curve 1
!write(unit=1) Pc(1, :)
!
!! close files
!close(unit=1)
!
!! Set filename to collect results.
!open(unit = 1, file = 'Pc2', &
!        form = 'unformatted', access = 'stream')
!
!! write the production curve 2
!write(unit=1) Pc(2, :)
!
!! close files
!close(unit=1)

!---------------------------------------------------------------------------
! Call GNUPLOT through the interface module.
! Uncomment these plot calls after verifying you have GNUPlot installed.
!---------------------------------------------------------------------------
! Plot the final solution for the forward trajectory
!call plot(Tt, Pc(1,:), Tt, Pc(2,:), terminal='png', filename='WaterOilCut.png')


!call print_array(Pc, 1,0,1,0,output)

contains


!
! This routine opens the permeability and porosity used by
! the MATLAB program and uses it for the simulation.
!
subroutine inputKP()
    integer :: i, j, k, l, m

    double precision, dimension(maxNx * maxNy * maxNz) :: pUr
    double precision, dimension(3 * maxNx, maxNy * maxNz) :: KUr
    double precision, dimension(3 * maxNx * maxNy * maxNz), target :: KUrl
    double precision, dimension(:), pointer :: KUrlP
    
    integer, dimension(Nx_ * Ny_ * Nz_) :: Pindices
    integer, dimension(3 * Nx_ * Ny_ * Nz_) :: Kindices

    ! read KUr
    open(1,file='KUr.txt',status='old')
    read(1,*) ((KUr(i,j), j=1,maxNy * maxNz), i=1,3 * maxNx)
    close(1)

    ! reshape 2 dimension to 1 dimension
    call myreshape(KUr, KUrl)
    ! then reshape 1 dimension to 4 dimension (hack for time being)
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

    KUrlP => KUrl(Kindices)
    
    call myreshape(KUrlP, K_)

    ! read KUr
    open(1,file='pUr.txt',status='old')
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

    Por_ = max(pUr(Pindices), 1.0d-3)
end subroutine inputKP


!
! Initialize inflow and outflow.
!
subroutine inflow_truncated_normal_x_outflow_point(Q, ir, mu, sigma)
!    use gnufor2
    integer :: i
!    double precision, dimension(Nx_) :: idx
    double precision :: x, pi, pdf, mu, sigma, ir, mass
    double precision, dimension(N_), target :: Q
    double precision, dimension(:), pointer :: Q_x

    ! set Q_x to be only the values along the first plane and x direction
    Q_x => Q(1:Nx_ * Ny_:Ny_)

    ! value of pi
    pi = 3.14159265358979323d0

    !initialize the total mass to 0
    mass = 0.0d0

    ! Note that the portion of the  Standard Normal distribution between
    ! -3sigma/2 to 3sigma/2 is assumed to fit the 1..Nx
    do i = 1, Nx_
        ! get the real x coordinate
        x = -1.5d0 + ((i - 1) * 3.0d0)/(Nx_ - 1)

        ! Now use mu and sigma to find the pdf value at x
        pdf = 1.0d0/(sigma * sqrt(2.0d0 * pi)) * exp(-(((x - mu)/sigma)**2.0d0)/2.0d0)

        ! set the value at the index equal to the pdf value at that point
        Q_x(i) = pdf

        ! increment the mass by the value of the pdf
        mass = mass + pdf

!        ! index to test initialization by plot
!        idx(i) = i * 1.0
    end do

    ! now rescale all the entities
    Q_x = Q_x/mass * ir

    ! now set the output
    Q(N_) = -ir

!    !---------------------------------------------------------------------------
!    ! Call GNUPLOT through the interface module.
!    ! Uncomment these plot calls after verifying you have GNUPlot installed.
!    !---------------------------------------------------------------------------
!    ! Plot the Q and check if it is correct.
!    call plot(idx, Q_x, terminal='png', filename='inflow.png')
end subroutine inflow_truncated_normal_x_outflow_point


!
! This subroutine simulates the reservoir
! model.
!
subroutine simulate_reservoir(Q, oil)
    integer :: i, j, k
    double precision :: oil, oil_temp
    double precision, dimension(N_) :: Q

    S = swc_                            ! initial saturation

    oil_temp = 0.0d0
    
    Pc(1, 1) = 0.0d0                    ! initial production
    Pc(2, 1) = 1.0d0
    Tt(1) = 0.0d0                       ! initial time.

    k = 1
    do i = 1, ND/Pt
        call Pres(S, Q, P, V)                       ! Pressure solver

        do j = 1, Pt/St
            k = k + 1
            call NewtRaph(S, V, Q, St * 1.0d0)      ! Solve for saturation
            call RelPerm(S(N_), Mw, Mo)             ! Mobilities in well-block

            Mt = Mw + Mo

            Tt(k) = 1.0d0 * k * St
            Pc(1,k) = Mw/Mt
            Pc(2,k) = Mo/Mt

            oil_temp = oil_temp + Pc(2, k) * St     ! Reimann sum
        end do
    end do
    
    oil = oil_temp
end subroutine simulate_reservoir


!
! Top level method to be differentiated
!
subroutine run_simulation(ir, mu, sigma, Q, oil)
    double precision :: ir, mu, sigma, oil
    double precision, dimension(N_) :: Q
    
    call inflow_truncated_normal_x_outflow_point(Q, ir, mu, sigma)
    call simulate_reservoir(Q, oil)
end subroutine run_simulation

end program runspe10
