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
    integer :: debug

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
    debug = 0

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
