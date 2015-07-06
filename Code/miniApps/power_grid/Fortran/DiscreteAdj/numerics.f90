module blas_and_lapack
contains

!
!  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”
!  or “3-clause license”)
!  Please read attached file License.txt
!

    double precision function dnrm2(n,x,incx)
      integer :: n,incx
      double precision :: x(n)
!     **********
!
!     Function dnrm2
!
!     Given a vector x of length n, this function calculates the
!     Euclidean norm of x with stride incx.
!
!     The function statement is
!
!       double precision function dnrm2(n,x,incx)
!
!     where
!
!       n is a positive integer input variable.
!
!       x is an input array of length n.
!
!       incx is a positive integer variable that specifies the
!         stride of the vector.
!
!     Subprograms called
!
!       FORTRAN-supplied ... abs, max, sqrt
!
!     MINPACK-2 Project. February 1991.
!     Argonne National Laboratory.
!     Brett M. Averick.
!
!     **********
      integer :: i
      double precision :: scale

      dnrm2 = 0.0d0
      scale = 0.0d0

      do  i = 1, n, incx
         scale = max(scale, abs(x(i)))
      enddo

      if (scale .eq. 0.0d0) return

      do i = 1, n, incx
         dnrm2 = dnrm2 + (x(i)/scale)**2
      enddo

      dnrm2 = scale*sqrt(dnrm2)

      return
    end function dnrm2

!====================== The end of dnrm2 ===============================

      subroutine daxpy(n,da,dx,incx,dy,incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end subroutine daxpy

!====================== The end of daxpy ===============================

      subroutine dcopy(n,dx,incx,dy,incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end subroutine dcopy

!====================== The end of dcopy ===============================

      double precision function ddot(n,dx,incx,dy,incy)
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
            dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end function ddot

!====================== The end of ddot ================================

      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end subroutine  dscal

!====================== The end of dscal ===============================

      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  Purpose
!  =======
!
!     DASUM takes the sum of the absolute values.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS,MOD
!     ..
      DASUM = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
!        code for increment equal to 1
!
!
!        clean-up loop
!
         M = MOD(N,6)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DABS(DX(I))
            END DO
            IF (N.LT.6) THEN
               DASUM = DTEMP
               RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,6
            DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) + &
                   DABS(DX(I+2)) + DABS(DX(I+3)) + &
                   DABS(DX(I+4)) + DABS(DX(I+5))
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DTEMP = DTEMP + DABS(DX(I))
         END DO
      END IF
      DASUM = DTEMP
      RETURN
      END FUNCTION DASUM
!====================== The end of dasum ===============================

      INTEGER FUNCTION IDAMAX(N,DX,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  Purpose
!  =======
!
!     IDAMAX finds the index of element having max. absolute value.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS
!     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
!
!        code for increment equal to 1
!
         DMAX = DABS(DX(1))
         DO I = 2,N
            IF (DABS(DX(I)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(I))
            END IF
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DABS(DX(IX)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END FUNCTION IDAMAX
!====================== The end of idamax ===============================

!
!  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”
!  or “3-clause license”)
!  Please read attached file License.txt
!
      subroutine dpofa(a,lda,n,info)
      integer lda,n,info
      double precision a(lda,*)
!
!     dpofa factors a double precision symmetric positive definite
!     matrix.
!
!     dpofa is usually called by dpoco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
!
!     on entry
!
!        a       double precision(lda, n)
!                the symmetric matrix to be factored.  only the
!                diagonal and upper triangle are used.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix  r  so that  a = trans(r)*r
!                where  trans(r)  is the transpose.
!                the strict lower triangle is unaltered.
!                if  info .ne. 0 , the factorization is not complete.
!
!        info    integer
!                = 0  for normal return.
!                = k  signals an error condition.  the leading minor
!                     of order  k  is not positive definite.
!
!     linpack.  this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas ddot
!     fortran sqrt
!
!     internal variables
!
      double precision t
      double precision s
      integer j,jm1,k
!     begin block with ...exits to 40
!
!
         do 30 j = 1, n
            info = j
            s = 0.0d0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
   10       continue
   20       continue
            s = a(j,j) - s
!     ......exit
            if (s .le. 0.0d0) go to 40
            a(j,j) = sqrt(s)
   30    continue
         info = 0
   40 continue
      return
      end subroutine dpofa

!====================== The end of dpofa ===============================

      subroutine dtrsl(t,ldt,n,b,job,info)
      integer ldt,n,job,info
      double precision t(ldt,*),b(*)
!
!
!     dtrsl solves systems of the form
!
!                   t * x = b
!     or
!                   trans(t) * x = b
!
!     where t is a triangular matrix of order n. here trans(t)
!     denotes the transpose of the matrix t.
!
!     on entry
!
!         t         double precision(ldt,n)
!                   t contains the matrix of the system. the zero
!                   elements of the matrix are not referenced, and
!                   the corresponding elements of the array can be
!                   used to store other information.
!
!         ldt       integer
!                   ldt is the leading dimension of the array t.
!
!         n         integer
!                   n is the order of the system.
!
!         b         double precision(n).
!                   b contains the right hand side of the system.
!
!         job       integer
!                   job specifies what kind of system is to be solved.
!                   if job is
!
!                        00   solve t*x=b, t lower triangular,
!                        01   solve t*x=b, t upper triangular,
!                        10   solve trans(t)*x=b, t lower triangular,
!                        11   solve trans(t)*x=b, t upper triangular.
!
!     on return
!
!         b         b contains the solution, if info .eq. 0.
!                   otherwise b is unaltered.
!
!         info      integer
!                   info contains zero if the system is nonsingular.
!                   otherwise info contains the index of
!                   the first zero diagonal element of t.
!
!     linpack. this version dated 08/14/78 .
!     g. w. stewart, university of maryland, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,ddot
!     fortran mod
!
!     internal variables
!
      double precision temp
      integer case,j,jj
!
!     begin block permitting ...exits to 150
!
!        check for zero diagonal elements.
!
         do 10 info = 1, n
!     ......exit
            if (t(info,info) .eq. 0.0d0) go to 150
   10    continue
         info = 0
!
!        determine the task and go to it.
!
         case = 1
         if (mod(job,10) .ne. 0) case = 2
         if (mod(job,100)/10 .ne. 0) case = case + 2
         go to (20,50,80,110), case
!
!        solve t*x=b for t lower triangular
!
   20    continue
            b(1) = b(1)/t(1,1)
            if (n .lt. 2) go to 40
            do 30 j = 2, n
               temp = -b(j-1)
               call daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
               b(j) = b(j)/t(j,j)
   30       continue
   40       continue
         go to 140
!
!        solve t*x=b for t upper triangular.
!
   50    continue
            b(n) = b(n)/t(n,n)
            if (n .lt. 2) go to 70
            do 60 jj = 2, n
               j = n - jj + 1
               temp = -b(j+1)
               call daxpy(j,temp,t(1,j+1),1,b(1),1)
               b(j) = b(j)/t(j,j)
   60       continue
   70       continue
         go to 140
!
!        solve trans(t)*x=b for t lower triangular.
!
   80    continue
            b(n) = b(n)/t(n,n)
            if (n .lt. 2) go to 100
            do 90 jj = 2, n
               j = n - jj + 1
               b(j) = b(j) - ddot(jj-1,t(j+1,j),1,b(j+1),1)
               b(j) = b(j)/t(j,j)
   90       continue
  100       continue
         go to 140
!
!        solve trans(t)*x=b for t upper triangular.
!
  110    continue
            b(1) = b(1)/t(1,1)
            if (n .lt. 2) go to 130
            do 120 j = 2, n
               b(j) = b(j) - ddot(j-1,t(1,j),1,b(1),1)
               b(j) = b(j)/t(j,j)
  120       continue
  130       continue
  140    continue
  150 continue
      return
      end subroutine dtrsl

!====================== The end of dtrsl ===============================

      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
!
!     dgesl solves the double precision system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgeco or dgefa.
!
!     on entry
!
!        a       double precision(lda, n)
!                the output from dgeco or dgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from dgeco or dgefa.
!
!        b       double precision(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgeco has set rcond .gt. 0.0
!        or dgefa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call dgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,ddot
!
!     internal variables
!
      double precision t
      integer k,kb,l,nm1
!
      nm1 = n - 1
      if (job .ne. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve  l*y = b
!
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
!
!        now solve  u*x = y
!
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
!
!        now solve trans(l)*x = y
!
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end subroutine dgesl
!====================== The end of dgesl ===============================

      subroutine dgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(1)
      double precision a(lda,1),z(1)
      double precision rcond
!
!     dgeco factors a double precision matrix by gaussian elimination
!     and estimates the condition of the matrix.
!
!     if  rcond  is not needed, dgefa is slightly faster.
!     to solve  a*x = b , follow dgeco by dgesl.
!     to compute  inverse(a)*c , follow dgeco by dgesl.
!     to compute  determinant(a) , follow dgeco by dgedi.
!     to compute  inverse(a) , follow dgeco by dgedi.
!
!     on entry
!
!        a       double precision(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        rcond   double precision
!                an estimate of the reciprocal condition of  a .
!                for the system  a*x = b , relative perturbations
!                in  a  and  b  of size  epsilon  may cause
!                relative perturbations in  x  of size  epsilon/rcond .
!                if  rcond  is so small that the logical expression
!                           1.0 + rcond .eq. 1.0
!                is true, then  a  may be singular to working
!                precision.  in particular,  rcond  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        z       double precision(n)
!                a work vector whose contents are usually unimportant.
!                if  a  is close to a singular matrix, then  z  is
!                an approximate null vector in the sense that
!                norm(a*z) = rcond*norm(a)*norm(z) .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     linpack dgefa
!     blas daxpy,ddot,dscal,dasum
!     fortran dabs,dmax1,dsign
!
!     internal variables
!
      double precision ek,t,wk,wkm
      double precision anorm,s,sm,ynorm
      integer info,j,k,kb,kp1,l

!
!
!     compute 1-norm of a
!
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
!
!     factor
!
      call dgefa(a,lda,n,ipvt,info)
!
!     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
!     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
!     trans(a)  is the transpose of a .  the components of  e  are
!     chosen to cause maximum local growth in the elements of w  where
!     trans(u)*w = e .  the vectors are frequently rescaled to avoid
!     overflow.
!
!     solve trans(u)*w = e
!
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
!
!     solve trans(l)*y = w
!
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
!
      ynorm = 1.0d0
!
!     solve l*v = y
!
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
!
!     solve  u*z = v
!
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
!     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
!
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end subroutine dgeco
!====================== The end of dgeco ===============================

      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
!
!     dgefa factors a double precision matrix by gaussian elimination.
!
!     dgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
!
!     on entry
!
!        a       double precision(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgesl or dgedi will divide by zero
!                     if called.  use  rcond  in dgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,idamax
!
!     internal variables
!
      double precision t
      integer j,k,kp1,l,nm1

!
!
!     gaussian elimination with partial pivoting
!
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
!
!        find l = pivot index
!
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
         if (a(l,k) .eq. 0.0d0) go to 40
!
!           interchange if necessary
!
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
!
!           compute multipliers
!
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
!
!           row elimination with column indexing
!
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end subroutine dgefa
!====================== The end of dgefa ===============================
end module blas_and_lapack


module numerics
    use blas_and_lapack
    implicit none

    ! Constants
    double precision :: pi, pmax, omegaB, omegaS, H, D, phiS, beta, c, t0, &
                        tend, tf, tcl, dt, perturb

    ! Shared
    double precision, save :: pm, quad_f, quad_b

    ! Set the parameter values
    parameter(                          &
        pi = 3.141592653589793d0,       &
        pmax = 2.0877d0,                &
        omegaB = 120.0d0 * pi,          &
        omegaS = 1d0,                   &
        H = 5.0d0,                      &
        D = 5.0d0,                      &
        phiS = 1.0d0,                   &
        beta = 2.0d0,                   &
        c = 10000.0d0,                  &
        t0 = 0.0d0,                     &
        tend = 10.0d0,                  &
        tf = 0.095d0,                   &  ! trying to evade floating pt issues
        tcl = 0.205d0,                  &  ! trying to evade floating pt issues
        dt = 0.01d0,                    &
        perturb = 1.0d0)

contains

    !
    ! Wrapper method which obtains the cost function and gradient
    ! used by the upper level optimization routine.
    !
    subroutine get_cost_function_and_gradient(pm, f, g, tlen)
        integer :: tlen, i
        double precision :: pm
        double precision, intent(out) :: f
        double precision, dimension(:), intent(out) :: g

        ! -------------------------------------------------------------------------
        ! declare variables for the integration scheme - use fixed time stepping
        ! may be use the SAVE attribute to avoid repeated allocation/deallocation
        ! -------------------------------------------------------------------------
        double precision, dimension(2) :: x
        double precision, dimension(2) :: lm
        double precision, dimension(tlen) :: tout_f
        double precision, dimension(tlen) :: tout_b
        double precision, dimension(tlen,2) :: yout_f
        double precision, dimension(tlen,2) :: yout_b

#ifdef USE_OPENAD
!$openad independent(pm)
!$openad dependent(f)
#endif
        ! Initialize the time array.
        tout_f(1) = t0
        do i = 2, tlen - 1
            tout_f(i) = tout_f(i - 1) + dt
        enddo
        tout_f(tlen) = tend

        ! Reverse the forward time array for adjoint use.
        tout_b = tout_f(tlen:1:-1)

        ! print the current parameter value
        ! print *, "pm is ", pm

        ! Initialize the x value
        x(1) = asin(perturb*pm/pmax)
        x(2) = perturb*1.0d0

        ! Compute the function value.
        call power_grid_cost_function(x, tout_f, yout_f, f)

#if !defined(USE_OPENAD)
        ! call the function that computes the sensitivites using
        ! continuous adjoints.

        ! Initialize the lambda value for the adjoint integration
        lm = (/0.0d0, 0.0d0/)

        ! Compute the gradient value.
        call power_grid_cost_gradient(lm, tout_b, yout_b, tout_f, yout_f, g)
#endif
    end subroutine get_cost_function_and_gradient


    !
    ! The subroutine computes power grid cost function
    ! This subroutine should call the time integrator
    ! in the forward mode along with the mode being
    ! set to "FWD".
    !
    subroutine power_grid_cost_function(x, tout_f, yout_f, f)
        implicit none

        double precision, dimension(2), intent(in) :: x
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:,:), intent(out) :: yout_f
        double precision, intent(out) :: f

        call crank_nicolson(x, tout_f, yout_f, "FWD")

        f = -pm + quad_f
        !print *, "The cost function is ", f
    end subroutine power_grid_cost_function

#if !defined(USE_OPENAD)
    !
    ! The subroutine computes the sensitivity of cost
    ! function with respect to the parameter p
    ! This subroutine should call the time integrator
    ! in the adjoint mode along with the mode being
    ! set to "ADJ".
    !
    subroutine power_grid_cost_gradient(x, tout_b, yout_b, tout_f, yout_f, g)
        implicit none
        ! The gradient can be computed using continuous adjoints
        ! or by discrete adjoints

        double precision, dimension(2), intent(in) :: x
        double precision, dimension(:), intent(in) :: tout_b
        double precision, dimension(:,:), intent(out) :: yout_b
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:,:), intent(in) :: yout_f
        double precision, dimension(:), intent(out) :: g

        call crank_nicolson(x, tout_b, yout_b, "ADJ", tout_f, yout_f)

        g(1) = -1.0d0 + quad_b &
                + 1.0d0/sqrt(pmax**2 - pm**2) * yout_b(size(yout_b, 1), 1)

        !print *, "The gradient is ", g(1)
    end subroutine power_grid_cost_gradient
#endif

    !
    ! Perform time integration by Crank Nicolson scheme.
    !
    subroutine crank_nicolson(x0, tout, yout, mode, tout_traj, yout_traj)
        implicit none

        logical :: converged
        integer :: max_conv_iter, iter, idx
        double precision :: dh, t, rcond
        double precision, dimension(:), intent(in) :: x0
        integer, dimension(2) :: ipvt
        double precision, dimension(2) :: x_n_1, x_n, dx_n, &
                                          f_n_1, f_n, wz
        double precision, dimension(2, 2) :: J_n, I_n
        double precision, dimension(:), intent(in) :: tout
        double precision, dimension(:, :), intent(out) :: yout
        double precision, dimension(:), intent(in), optional :: tout_traj
        double precision, dimension(:, :), intent(in), optional :: yout_traj
        character(len = 3) :: mode

        ! set the maximum number of iterations (newton) before convergence.
        max_conv_iter = 10

        ! get the first time step
        t = tout(1)

        ! set the solution at the first time step as starting value
        yout(1, :) = x0

#if !defined(USE_OPENAD)
        if (mode == "FWD") then
            ! get the step size depending on whether we are going forward or
            ! backward. (FWD)
            dh = dt

            ! Initialize the forward quadrature value
            call forward_ode_quad(t, x0, "INIT")
        elseif (mode == "ADJ") then
            ! Check whether optional arguments are present
            if(.not.(present(tout_traj) .and. present(yout_traj))) then
                stop "The solution trajectories are not present for &
                &the adjoint mode"
            endif

            ! get the step size depending on whether we are going forward or
            ! backward. (FWD)
            dh = -dt

            ! Initialize the backward quadrature value
            call adjoint_ode_quad(t, x0, "INIT")
        endif
#else
        ! get the step size depending on whether we are going forward or
        ! backward. (FWD)
        dh = dt

        ! Initialize the forward quadrature value
        call forward_ode_quad(t, x0, "INIT")
#endif
        ! read the first time instant x value
        x_n_1 = x0

        ! start newton iteration with previous time instant value
        ! this assignment was initially inside the newton iteration
        ! loop. It has however been pulled out as it is duplicating
        ! the work done inside the norm check if-block
        x_n = x_n_1

        ! have an identity matrix handy
        call identity_matrix(I_n)

        ! compute all the other x values based on
        ! solving the nonlinear system (semi-implicit)
        ! using newton iteration.
        do idx = 2, size(tout)
#if !defined(USE_OPENAD)
            ! compute the f value at the previous time instant
            ! and initial previous x value.
            if (mode == "FWD") then
                call forward_ode_rhs(t, x_n_1, f_n_1)
            elseif (mode == "ADJ") then
                call adjoint_ode_rhs(t, x_n_1, tout_traj, yout_traj, f_n_1)
            endif
#else
            ! compute the f value at the previous time instant
            ! and initial previous x value.
            call forward_ode_rhs(t, x_n_1, f_n_1)
#endif
            ! get the new t value
            t = tout(idx)

            ! set converged to false
            converged = .false.

            ! try to converge within max_conv_iter
            do iter = 1, max_conv_iter
                if(.not. converged) then
#if !defined(USE_OPENAD)
                    ! get the function value and the jacobian
                    ! at the unconverged x and new t
                    if (mode == "FWD") then
                        call forward_ode_rhs(t, x_n, f_n)
                        call forward_ode_jac(t, x_n, J_n)
                    elseif (mode == "ADJ") then
                        call adjoint_ode_rhs(t, x_n, tout_traj, yout_traj, f_n)
                        ! Note that the jacobian in the adjoint mode is
                        ! independent of the current iterate. However,
                        ! to keep things symmetric between fwd and adj,
                        ! we are recomputing the jacobian each time.
                        ! This is wasteful, yes. But its a 2x2 system
                        ! in this particular case.
                        call adjoint_ode_jac(t, tout_traj, yout_traj, J_n)
                    endif
#else
                    ! get the function value and the jacobian
                    ! at the unconverged x and new t
                    call forward_ode_rhs(t, x_n, f_n)
                    call forward_ode_jac(t, x_n, J_n)
#endif

                    ! construct the left hand side matrix
                    J_n =  (I_n - dh/2 * J_n)

                    ! construct the right hand side vector
                    dx_n = x_n - x_n_1 - dh/2 * (f_n + f_n_1)

                    ! solve the system for dx_n
                    ! first factorize the matrix
                    call dgeco(J_n, 2, 2, ipvt, rcond, wz)

                    ! check the conditioning of the matrix
                    if (1.0d0 + rcond .eq. 1.0d0) then
                        !stop "Matrix is badly conditioned."
                        !print *,"Matrix is badly conditioned."
                    endif

                    ! solve the system
                    call dgesl(J_n, 2, 2, ipvt, dx_n, 0)

                    ! update the iterate
                    x_n = x_n - dx_n

                    ! check whether the right hand side vector with
                    ! new iterate has converged
                    ! Note that moving the if block down may change
                    ! the optimal slightly.
                    if (dnrm2(2,dx_n, 1) < 1.0d-8) then
                        ! print *, "Converged at t = ", t
                        x_n_1 = x_n
                        f_n_1 = f_n
                        converged = .true.
                    endif
                else
                    ! do nothing.
                endif
            enddo

            ! update the new value in the solution trajectory
            yout(idx, :) = x_n


            ! Perform quadrature along
#if !defined(USE_OPENAD)
            if (mode == "FWD") then
                ! Iterate the forward quadrature value
                call forward_ode_quad(t, x_n, "ITER")
            elseif (mode == "ADJ") then
                ! Iterate the backward quadrature value
                call adjoint_ode_quad(t, x_n, "ITER")
            endif
#else
            call forward_ode_quad(t, x_n, "ITER")
#endif
        enddo


#if !defined(USE_OPENAD)
        ! Finalize the quadrature term.
        if (mode == "FWD") then
            ! Finalize the forward quadrature value
            call forward_ode_quad(t, x_n, "DONE")
        elseif (mode == "ADJ") then
            ! Finalize the backward quadrature value
            call adjoint_ode_quad(t, x_n, "DONE")
        endif
#else
        call forward_ode_quad(t, x_n, "DONE")
#endif
    end subroutine crank_nicolson

    !
    ! Subroutine initializes the values of I to
    ! be a square/non-square (meaningless) identity
    ! matrix.
    !
    subroutine identity_matrix(I_n)
        integer :: idx, jdx
        double precision, dimension(:,:), intent(out) :: I_n

        do idx = 1, size(I_n,1)
            do jdx = 1, size(I_n, 2)
                if (idx == jdx) then
                    I_n(idx,jdx) = 1
                else
                    I_n(idx,jdx) = 0
                endif
            enddo
        enddo
    end subroutine identity_matrix

    !
    ! The RHS function of the forward ODE
    !
    subroutine forward_ode_rhs(t, x, y)
        implicit none

        double precision, intent(in) :: t
        double precision, dimension(2), intent(in) :: x
        double precision, dimension(2), intent(out) :: y
        double precision :: phi
        double precision :: omega
        double precision :: pmax_

        phi = x(1)
        omega = x(2)

        if (t > tf .AND. t <= tcl) then
            pmax_ = 0.0d0
        else
            pmax_ = pmax
        endif

        y(1) = omegaB * (omega - omegaS)
        y(2) = omegaS * (pm - pmax_ * sin(phi) - D * (omega - omegaS))/(2*H)
    end subroutine forward_ode_rhs

    !
    ! The JAC function of the forward ODE
    !
    subroutine forward_ode_jac(t, x, J)
        implicit none

        double precision, intent(in) :: t
        double precision, dimension(2), intent(in) :: x
        double precision, dimension(2,2), intent(out) :: J
        double precision :: phi
        double precision :: omega
        double precision :: pmax_

        phi = x(1)
        omega = x(2)

        if (t > tf .AND. t <= tcl) then
            pmax_ = 0.0d0
        else
            pmax_ = pmax
        endif

        J(1, 1) = 0.0d0
        J(1, 2) = omegaB
        J(2, 1) = -omegaS * pmax_ * cos(phi) / (2.0d0 * H)
        J(2, 2) = -D * omegaS / (2.0d0 * H)
    end subroutine forward_ode_jac

    !
    ! The subroutine to compute the quadrature in FWD mode
    !
    subroutine forward_ode_quad(t, x, state)
        implicit none

        character(len = 4), intent(in) :: state
        double precision, intent(in) :: t
        double precision, dimension(2), intent(in) :: x
        double precision :: phi, omega, new_h
        double precision,save :: prev_h

        phi = x(1)
        omega = x(2)

        if (state == "ITER") then
            ! perform trapezoidal integration using previous height
            ! and new height.
            new_h =  c * max(0.0d0, (phi - phiS))**beta
            quad_f = quad_f + 0.5d0 * dt * (prev_h + &
                       new_h)
            ! update the prev height.
            prev_h = new_h
        elseif (state == "INIT") then
            quad_f = 0.0d0
            prev_h = c * max(0.0d0, (phi - phiS))**beta
        elseif (state == "DONE") then
            !print *, "The forward quadrature value is ", quad_f
        endif
    end subroutine forward_ode_quad


#if !defined(USE_OPENAD)
    !
    ! The RHS function of the adjoint ODE
    ! Uses the trajectory of the forward solution to
    ! compute adjoint RHS function
    !
    subroutine adjoint_ode_rhs(t, x, tout_f, yout_f, y)
        implicit none

        integer :: idx
        double precision :: phi
        double precision, intent(in) :: t
        double precision, dimension(2,2) :: J
        double precision, dimension(2) :: forcing
        double precision, dimension(2), intent(in) :: x
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:, :), intent(in) :: yout_f
        double precision, dimension(2), intent(out) :: y

        ! find the index of the current t
        idx = findindex(t, tout_f)

        ! get the jacobian of the forward ode
        call forward_ode_jac(t, yout_f(idx, :), J)

        ! get phi at the given time instant
        phi = yout_f(idx, 1)

        ! compute the forcing term
        forcing = (/ c * beta * max(0.0d0, &
                    (phi - phiS))**(beta - 1.0d0), 0.0d0 /)

        ! now compute y
        y = matmul(transpose(-J), x) - forcing
    end subroutine

    !
    ! The jacobian of the adjoint ODE
    ! Uses the trajectory of the forward solution to
    ! compute the jacobian at a given `t`
    !
    subroutine adjoint_ode_jac(t, tout_f, yout_f, J)
        implicit none

        integer :: idx
        double precision, intent(in) :: t
        double precision, dimension(2,2), intent(out) :: J
        double precision, dimension(:), intent(in) :: tout_f
        double precision, dimension(:, :), intent(in) :: yout_f

        ! find the index of the current t
        idx = findindex(t, tout_f)

        ! get the jacobian of the forward ode
        call forward_ode_jac(t, yout_f(idx, :), J)

        J = transpose(-J)
    end subroutine

    !
    ! The subroutine to compute the quadrature in ADJ mode
    ! Uses the trajectory of the forward solution to compute
    ! the quadrature term
    !
    subroutine adjoint_ode_quad(t, x, state)
        implicit none

        character(len = 4), intent(in) :: state
        double precision, intent(in) :: t
        double precision, dimension(2), intent(in) :: x
        double precision :: phi, omega, new_h
        double precision,save :: prev_h

        phi = x(1)
        omega = x(2)

        if (state == "ITER") then
            ! perform trapezoidal integration using previous height
            ! and new height.
            new_h = dot_product((/0.0d0, omegaS/(2.0d0*H)/), x)
            quad_b = quad_b + 0.5d0 * dt * (prev_h + &
                       new_h)
            ! update the previous height.
            prev_h = new_h
        elseif (state == "INIT") then
            quad_b = 0.0d0
            prev_h = dot_product((/0.0d0, omegaS/(2.0d0*H)/), x)
        elseif (state == "DONE") then
            !print *, "The adjoint quadrature value is ", quad_b
        endif
    end subroutine adjoint_ode_quad

    !
    ! Function returns the first index when traversing
    ! from lowest index (1) to highest index (length)
    ! whose value matches the desired value
    !
    integer function findindex(t, array)
        implicit none

        double precision, intent(in) :: t
        double precision, dimension(:), intent(in) :: array

        integer :: idx, default_index

        default_index = 0
        do idx = 1, size(array, 1)
            if(array(idx) == t) then
                default_index = idx
                exit
            endif
        enddo

        findindex = default_index
    end function findindex
#endif

end module numerics
