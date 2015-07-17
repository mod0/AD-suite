module sparse_matrix

contains


!
! The sparse diagonal matrix format
! rows: The number of rows in the full matrix
! cols: The number of columns in the full matrix
! diagonals: The order in which the diagonals are stored in values.
! values: The actual values in the diagonals with zeros filled in
!         like spdiags does in MATLAB
!
type spdiag
    integer :: rows, cols
    integer, dimension(:), allocatable :: diagonals
    double precision, dimension(:,:), allocatable :: values
end type


!
! This routine creates a sparse matrix using the input columns
! in imatrix.
! imatrix: the columns of values to be placed in diagonal matrix
! idiags: the diagonals in which the values would go (sorted)
! irows: the rows in the output matrix
! icols: the columns in the output matrix
! omatrixsp: the output sparse matrix in spdiag format.
! The placement of values in the diagonals will follow MATLABs
! convention: spdiags
!
subroutine spdiags(imatrix, idiags, irows, icols, omatrixsp)
    integer :: irows, icols
    double precision, dimension(:,:) :: imatrix
    integer, dimension(:) :: idiags
    type(spdiag) :: omatrix

! Skipping validation of arguments
! Note: imatrix should have as many rows as min(irows, icols)
!       the number of columns in imatrix should match length of idiags
!       and should not exceed m + n - 1 (soft requirement)
!       the values of idiags should range between -irows + 1 to icols - 1

! if irows == icols || irows > icols, elements of superdiagonal filled from lower
! part of corresponding column in imatrix and from upper part for subdiagonals
! if irows < icols, elements of superdiagonal filled from top part of corresponding
! column in imatrix and from lower part for subdiagonals

! Note: Zeros will be filled in appropriately where elements do not belong to the
!       matrix. This should be taken into consideration whenever the matrix is
!       accessed.

end subroutine

!
! Common matrix operations need to be implemented.
! i)  inc/dec the elements of the matrix by some quantity
! ii) Add/Sub two matrices (account for fill in?)
! iii) Multiply two matrices - 2 sparse matrices (fill in?)
! iv) Sparse linear solve.
!



end module sparse_matrix
