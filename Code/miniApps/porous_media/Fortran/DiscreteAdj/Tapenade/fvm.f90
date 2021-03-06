module fvm
use grid
use fluid
use matrix
use linsolve
use mathutil
!use print_active
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
subroutine NewtRaph(S, V, Q, St, solver_inner, solver_outer, verbose)
  !use print_active
  implicit none
  logical :: verbose
  integer :: solver_inner, solver_outer
  
  integer :: i, j, it
  integer :: St                             ! Maximum saturation Time Step
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
      call spmat_multiply_diagonal2( annz, arow_index, arow_compressed,&
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
              call spmat_multiply_diagonal2( bnnz, brow_index, brow_compressed, &
                                            bcol_index, bvalues, dF, &
                                             dgnnz, dgrow_index, dgrow_compressed, &
                                            dgcol_index, dgvalues, "PRE")

              call addx_diagonal2( dgnnz, dgrow_index, dgrow_compressed, &
                                  dgcol_index, dgvalues, -1.0d0, 0)

              fw = Mw / (Mw + Mo)

              ! Matrix-vector matrix product
              call spmat_multiply_vector2( bnnz, brow_index, brow_compressed, &
                                          bcol_index, bvalues, fw, bfw, "PRE")

              G = S - S_iter_copy - bfw - fi

              call solve( dgnnz, dgrow_index, dgrow_compressed, &
                         dgcol_index, dgvalues, G, dS, solver_inner, solver_outer, verbose)

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
subroutine Pres(S, Q, P, V, solver_inner, solver_outer, verbose)
    !use print_active
    logical :: verbose
    integer :: solver_inner, solver_outer
    
    integer i

    double precision, dimension(N_) :: S
    double precision, dimension(N_) :: Q
    double precision, dimension(3 * N_) :: M
    double precision, dimension(Nx_, Ny_, Nz_) :: P
    double precision, dimension(3, Nx_, Ny_, Nz_) :: KM
    double precision, dimension(3, Nx_ + 1, Ny_ + 1, Nz_ + 1) :: V

    double precision, dimension(N_) :: Mw
    double precision, dimension(N_) :: Mo

    !call RelPerm(S, M(1 : 3 * N_ : 3), M(2 : 3 * N_ : 3))

    call RelPerm(S, Mw, Mo)

    do i = 1,N_
      M(1 + (i - 1) * 3) = Mw(i) + Mo(i)
      M(2 + (i - 1) * 3) = M(1 + (i - 1) * 3)
      M(3 + (i - 1) * 3) = M(1 + (i - 1) * 3)
    end do

    call myreshape_1_4(M, KM)

    ! point-wise multiply
    KM = KM * PERM

    call tpfa(KM, Q, P, V, solver_inner, solver_outer, verbose)
end subroutine


!
! Relative Permeabilities
!
subroutine RelPerm_vector(S, Mw, Mo, dMw, dMo)
    !use print_active
    implicit none
    double precision, dimension(N_) :: S
    double precision, dimension(N_) :: Mw
    double precision, dimension(N_) :: Mo
    double precision, dimension(N_) :: S_
    double precision, dimension(N_), optional :: dMw
    double precision, dimension(N_), optional :: dMo

    S_ = (S - swc_)/(1.0d0 - swc_ - sor_)   ! rescale saturation

    Mw = S_**2/vw_
    Mo = (1 - S_)**2/vo_

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
    double precision :: S, Mw, Mo, S_
    double precision, optional :: dMw, dMo

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
subroutine GenA(V, Q,  annz, arow_index, arow_compressed, acol_index, avalues)
  !use print_active
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

  ! This can be sped up by passing 3 arrays having rows, cols and diagind
  ! this can be done because diag positions are fixed.
  !TODO: Have a variant of spdiags which will instead of writing it in
  !the spmat type, it will write it in separate arrays sent to it.
  
  call spdiags_fvm_csr(diags, annz, arow_index, arow_compressed,&
                    acol_index, avalues)
end subroutine GenA

!
! Two point flux approximation.
!
subroutine tpfa(K, Q, P, V, solver_inner, solver_outer, verbose)
    !use print_active
    implicit none
    logical :: verbose
    integer :: solver_inner, solver_outer
    
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


    !TODO: Have a variant of spdiags which will instead of writing it in
    !the spmat type, it will write it in separate arrays sent to it.

    call spdiags_fvm_csr(diags, annz, arow_index, arow_compressed, &
                     acol_index, avalues)


    ! ! Increment the 1,1 element of A
!     call addx_elem2(annz, arow_index, arow_compressed,&
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
               acol_index, avalues, Q, u, solver_inner, solver_outer, verbose)

    ! reshape the solution
    call myreshape_1_3(u, P)

    ! V.x
    V(1, 2:Nx_, 1:Ny_, 1:Nz_) = (P(1:Nx_ - 1, :, :) - P(2:Nx_, :, :)) * TX(2:Nx_,:,:)
    ! V.y
    V(2, 1:Nx_, 2:Ny_, 1:Nz_) = (P(:, 1:Ny_ - 1, :) - P(:, 2:Ny_, :)) * TY(:,2:Ny_,:)
    ! V.z
    V(3, 1:Nx_, 1:Ny_, 2:Nz_) = (P(:, :, 1:Nz_ - 1) - P(:, :, 2:Nz_)) * TZ(:,:,2:Nz_)
end subroutine tpfa

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
    integer, parameter, dimension(7) :: idiags = (/-Nx_ * Ny_, -Nx_, -1, &
                                                   0, 1, Nx_, Nx_ * Ny_/)
                                                   
    integer :: onnz
    integer, dimension(7 * N_) :: orow_index
    integer, dimension(7 * N_) :: ocol_index
    double precision, dimension(7 * N_) :: ovalues
    integer, dimension(N_ + 1) :: orow_compressed


    
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
    integer, parameter, dimension(7) :: idiags = (/-Nx_ * Ny_, -Nx_, -1, &
                                                   0, 1, Nx_, Nx_ * Ny_/)
    integer, dimension(7) :: start_row_imatrix, end_row_imatrix
    integer, dimension(7) :: row_diag, col_diag      ! row, column along diagonal

                                                   
    integer :: onnz
    integer, dimension(7 * N_) :: orow_index
    integer, dimension(7 * N_) :: ocol_index
    double precision, dimension(7 * N_) :: ovalues
    integer, dimension(N_ + 1) :: orow_compressed



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

end module fvm