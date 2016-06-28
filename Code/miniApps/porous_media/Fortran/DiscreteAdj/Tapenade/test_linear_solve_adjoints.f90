program test_linsolve_adj

  integer :: n, annz, alen
  integer, dimension(5) :: arow_compressed
  integer, dimension(6) :: arow_index, acol_index
  double precision, dimension(6) :: avalues
  double precision, dimension(6) :: avaluesd
  double precision, dimension(6) :: avaluesb
  
  double precision, dimension(4) :: b1
  double precision, dimension(4) :: bd1
  double precision, dimension(4) :: bb1
  
  double precision, dimension(4) :: x1
  double precision, dimension(4) :: xd1
  double precision, dimension(4,4) :: xDD1
  double precision, dimension(4) :: xb1
  double precision, dimension(4,4) :: xBB1
  

  n = 4
  annz = 6
  alen = 6
  
  arow_index = (/1,2,3,3,4,4/)
  acol_index = (/1,2,3,4,2,4/)
  arow_compressed = (/1,2,3,5,7/)
  avalues = (/9.474578019835103e-01, 9.630885392869130e-01, 5.468057187389680e-01, &
              2.372835797715215e-01, 5.785250610234389e-01, 5.211358308040015e-01/)
  
  b1 = (/2.315943867085238e-01, 4.888977439201669e-01, 6.240600881736895e-01, &
       6.791355408657477e-01/)
  
  avaluesd = 0.0d0

  bd1 = (/1,0,0,0/)
  CALL SPARSE_PMGMRES_METHOD_D(n, annz, alen, arow_index, arow_compressed, &
                               acol_index, avalues, avaluesd, b1, bd1, x1, xd1,&
                               4, 4, .false.)
  xDD1(:, 1) = xd1
  write (*,*) "x1: ", x1
  write (*,*) "xd1: ", xd1
  
  bd1 = (/0,1,0,0/)
  CALL SPARSE_PMGMRES_METHOD_D(n, annz, alen, arow_index, arow_compressed, &
                               acol_index, avalues, avaluesd, b1, bd1, x1, xd1, &
                               4, 4, .false.)
  xDD1(:, 2) = xd1
  write (*,*) "x1: ", x1
  write (*,*) "xd1: ", xd1
  
  bd1 = (/0,0,1,0/)
  CALL SPARSE_PMGMRES_METHOD_D(n, annz, alen, arow_index, arow_compressed, &
                               acol_index, avalues, avaluesd, b1, bd1, x1, xd1, &
                               4, 4, .false.)

  xDD1(:, 3) = xd1
  write (*,*) "x1: ", x1
  write (*,*) "xd1: ", xd1
  
  bd1 = (/0,0,0,1/)
  CALL SPARSE_PMGMRES_METHOD_D(n, annz, alen, arow_index, arow_compressed, &
                               acol_index, avalues, avaluesd, b1, bd1, x1, xd1, &
                               4, 4, .false.)
                               
  xDD1(:, 4) = xd1
  write (*,*) "x1: ", x1
  write (*,*) "xd1: ", xd1
  
  xb1 = (/1.0d0, 0.0d0, 0.0d0, 0.0d0/)
  bb1 = 0.0d0
  avaluesb = 0.0d0
  CALL SPARSE_PMGMRES_METHOD_B(n, annz, alen, arow_index,   &
                          &   arow_compressed, acol_index, avalues,&
                          &   avaluesb, b1, bb1, x1, xb1,              &
                          4, 4, .false.)
  xBB1(1, :) = bb1
  write (*,*) "x1: ", x1
  write (*,*) "bb1: ", bb1

                          
  xb1 = (/0.0d0, 1.0d0, 0.0d0, 0.0d0/)
  bb1 = 0.0d0
  avaluesb = 0.0d0
  CALL SPARSE_PMGMRES_METHOD_B(n, annz, alen, arow_index,   &
                          &   arow_compressed, acol_index, avalues,&
                          &   avaluesb, b1, bb1, x1, xb1,              &
                          4, 4, .false.)                      
  xBB1(2, :) = bb1
  write (*,*) "x1: ", x1
  write (*,*) "bb1: ", bb1
  
  xb1 = (/0.0d0, 0.0d0, 1.0d0, 0.0d0/)
  bb1 = 0.0d0
  avaluesb = 0.0d0
  CALL SPARSE_PMGMRES_METHOD_B(n, annz, alen, arow_index,   &
                          &   arow_compressed, acol_index, avalues,&
                          &   avaluesb, b1, bb1, x1, xb1,              &
                          4, 4, .false.)
  xBB1(3, :) = bb1                          
  write (*,*) "x1: ", x1
  write (*,*) "bb1: ", bb1
  
  xb1 = (/0.0d0, 0.0d0, 0.0d0, 1.0d0/)
  bb1 = 0.0d0
  avaluesb = 0.0d0
  CALL SPARSE_PMGMRES_METHOD_B(n, annz, alen, arow_index,   &
                          &   arow_compressed, acol_index, avalues,&
                          &   avaluesb, b1, bb1, x1, xb1,              &
                          4, 4, .false.)
  xBB1(4, :) = bb1
  write (*,*) "x1: ", x1
  write (*,*) "bb1: ", bb1  
  
  
  write(*,*) "Adjoint Jacobian"
  write (*, *) xBB1(1, :)
  write (*, *) xBB1(2, :)
  write (*, *) xBB1(3, :)
  write (*, *) xBB1(4, :)

  write(*,*) "Tangent Linear Jacobian"
  write (*, *) xDD1(1, :)
  write (*, *) xDD1(2, :)
  write (*, *) xDD1(3, :)
  write (*, *) xDD1(4, :)
  

  
!   SUBROUTINE SPARSE_PMGMRES_METHOD(n, annz, alen, arow_index, arow_compressed, &
!                                   acol_index, avalues, b, x, solver_inner, &
!                                   solver_outer, verbose)
  
end program test_linsolve_adj