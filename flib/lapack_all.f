c
c   This file contains LAPACK routines used in quantum-espresso
c   that are part of ATLAS - from www.netlib.org
c   These are:
* [S,D,C,Z]GESV
* [S,D,C,Z]GETRF
* [S,D,C,Z]GETRS
* [S,D,C,Z]GETRI
* [S,D,C,Z]TRTRI
* [S,D,C,Z]POSV
* [S,D,C,Z]POTRF
* [S,D,C,Z]POTRS
* [S,D,C,Z]POTRI
* [S,D,C,Z]LAUUM 
c
      SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGESV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor A as
*     A = P * L * U,
*  where P is a permutation matrix, L is unit lower triangular, and U is
*  upper triangular.  The factored form of A is then used to solve the
*  system of equations A * X = B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N coefficient matrix A.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS matrix of right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, so the solution could not be computed.
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL           DGETRF, DGETRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESV ', -INFO )
         RETURN
      END IF
*
*     Compute the LU factorization of A.
*
      CALL DGETRF( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL DGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
     $                INFO )
      END IF
      RETURN
*
*     End of DGESV
*
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, JP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRI computes the inverse of a matrix using the LU factorization
*  computed by DGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by DGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,
     $                   NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRI', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DGETRI
*
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A' * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END
 
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
*  Further Details
*  ===============
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
      SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOTF2 computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U' * U ,  if UPLO = 'U', or
*     A = L  * L',  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n by n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U'*U  or A = L*L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization A = U'*U.
*
         DO 10 J = 1, N
*
*           Compute U(J,J) and test for non-positive-definiteness.
*
            AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of row J.
*
            IF( J.LT.N ) THEN
               CALL DGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
*
*        Compute the Cholesky factorization A = L*L'.
*
         DO 20 J = 1, N
*
*           Compute L(J,J) and test for non-positive-definiteness.
*
            AJJ = A( J, J ) - DDOT( J-1, A( J, 1 ), LDA, A( J, 1 ),
     $            LDA )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of column J.
*
            IF( J.LT.N ) THEN
               CALL DGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ),
     $                     LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = J
*
   40 CONTINUE
      RETURN
*
*     End of DPOTF2
*
      END
      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOTRF computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U**T * U,  if UPLO = 'U', or
*     A = L  * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the block version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the factorization could not be
*                completed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DPOTF2, DSYRK, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         CALL DPOTF2( UPLO, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code.
*
         IF( UPPER ) THEN
*
*           Compute the Cholesky factorization A = U'*U.
*
            DO 10 J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Upper', 'Transpose', JB, J-1, -ONE,
     $                     A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block row.
*
                  CALL DGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1,
     $                        J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ),
     $                        LDA, ONE, A( J, J+JB ), LDA )
                  CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit',
     $                        JB, N-J-JB+1, ONE, A( J, J ), LDA,
     $                        A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
*
         ELSE
*
*           Compute the Cholesky factorization A = L*L'.
*
            DO 20 J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Lower', 'No transpose', JB, J-1, -ONE,
     $                     A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block column.
*
                  CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,
     $                        J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ),
     $                        LDA, ONE, A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'Transpose', 'Non-unit',
     $                        N-J-JB+1, JB, ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = INFO + J - 1
*
   40 CONTINUE
      RETURN
*
*     End of DPOTRF
*
      END
      SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOTRS solves a system of linear equations A*X = B with a symmetric
*  positive definite matrix A using the Cholesky factorization
*  A = U**T*U or A = L*L**T computed by DPOTRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The triangular factor U or L from the Cholesky factorization
*          A = U**T*U or A = L*L**T, as computed by DPOTRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B where A = U'*U.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A*X = B where A = L*L'.
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
      END IF
*
      RETURN
*
*     End of DPOTRS
*
      END

      SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTI2 computes the inverse of a real upper or lower triangular
*  matrix.
*
*  This is the Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading n by n upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DTRMV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DTRTI2
*
      END
      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTRI computes the inverse of a real upper or lower triangular
*  matrix A.
*
*  This is the Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*               matrix is singular and its inverse can not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DTRTRI
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
      SUBROUTINE ZDROT( N, CX, INCX, CY, INCY, C, S )
*
*     applies a plane rotation, where the cos and sin (c and s) are real
*     and the vectors cx and cy are complex.
*     jack dongarra, linpack, 3/11/78.
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   C, S
*     ..
*     .. Array Arguments ..
      COMPLEX*16         CX( * ), CY( * )
*
* =====================================================================
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, IY
      COMPLEX*16         CTEMP
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 )
     $   GO TO 20
*
*        code for unequal increments or equal increments not equal
*          to 1
*
      IX = 1
      IY = 1
      IF( INCX.LT.0 )
     $   IX = ( -N+1 )*INCX + 1
      IF( INCY.LT.0 )
     $   IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         CTEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - S*CX( IX )
         CX( IX ) = CTEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
   20 CONTINUE
      DO 30 I = 1, N
         CTEMP = C*CX( I ) + S*CY( I )
         CY( I ) = C*CY( I ) - S*CX( I )
         CX( I ) = CTEMP
   30 CONTINUE
      RETURN
      END
      SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            J, JP
*     ..
*     .. External Functions ..
      INTEGER            IZAMAX
      EXTERNAL           IZAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M )
     $         CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),
     $                  LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of ZGETF2
*
      END
      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGETF2, ZLASWP, ZTRSM
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL ZGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL ZGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZGETRF
*
      END

      SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRI computes the inverse of a matrix using the LU factorization
*  computed by ZGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by ZGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,
     $                   NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGEMV, ZSWAP, ZTRSM, ZTRTRI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NB = ILAENV( 1, 'ZGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRI', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from ZTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL ZTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL ZGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL ZGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL ZTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL ZSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZGETRI
*
      END
      SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRS solves a system of linear equations
*     A * X = B,  A**T * X = B,  or  A**H * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by ZGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by ZGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLASWP, ZTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL ZTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A**T * X = B  or A**H * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL ZTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL ZTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,
     $               LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of ZGETRS
*
      END


      SUBROUTINE ZPOTF2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZPOTF2 computes the Cholesky factorization of a complex Hermitian
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U' * U ,  if UPLO = 'U', or
*     A = L  * L',  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          n by n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U'*U  or A = L*L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZGEMV, ZLACGV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPOTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization A = U'*U.
*
         DO 10 J = 1, N
*
*           Compute U(J,J) and test for non-positive-definiteness.
*
            AJJ = DBLE( A( J, J ) ) - ZDOTC( J-1, A( 1, J ), 1,
     $            A( 1, J ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of row J.
*
            IF( J.LT.N ) THEN
               CALL ZLACGV( J-1, A( 1, J ), 1 )
               CALL ZGEMV( 'Transpose', J-1, N-J, -CONE, A( 1, J+1 ),
     $                     LDA, A( 1, J ), 1, CONE, A( J, J+1 ), LDA )
               CALL ZLACGV( J-1, A( 1, J ), 1 )
               CALL ZDSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
*
*        Compute the Cholesky factorization A = L*L'.
*
         DO 20 J = 1, N
*
*           Compute L(J,J) and test for non-positive-definiteness.
*
            AJJ = DBLE( A( J, J ) ) - ZDOTC( J-1, A( J, 1 ), LDA,
     $            A( J, 1 ), LDA )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of column J.
*
            IF( J.LT.N ) THEN
               CALL ZLACGV( J-1, A( J, 1 ), LDA )
               CALL ZGEMV( 'No transpose', N-J, J-1, -CONE, A( J+1, 1 ),
     $                     LDA, A( J, 1 ), LDA, CONE, A( J+1, J ), 1 )
               CALL ZLACGV( J-1, A( J, 1 ), LDA )
               CALL ZDSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = J
*
   40 CONTINUE
      RETURN
*
*     End of ZPOTF2
*
      END
      SUBROUTINE ZPOTRF( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZPOTRF computes the Cholesky factorization of a complex Hermitian
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U**H * U,  if UPLO = 'U', or
*     A = L  * L**H,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the block version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**H*U or A = L*L**H.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the factorization could not be
*                completed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      COMPLEX*16         CONE
      PARAMETER          ( ONE = 1.0D+0, CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZHERK, ZPOTF2, ZTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPOTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         CALL ZPOTF2( UPLO, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code.
*
         IF( UPPER ) THEN
*
*           Compute the Cholesky factorization A = U'*U.
*
            DO 10 J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL ZHERK( 'Upper', 'Conjugate transpose', JB, J-1,
     $                     -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL ZPOTF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block row.
*
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', JB,
     $                        N-J-JB+1, J-1, -CONE, A( 1, J ), LDA,
     $                        A( 1, J+JB ), LDA, CONE, A( J, J+JB ),
     $                        LDA )
                  CALL ZTRSM( 'Left', 'Upper', 'Conjugate transpose',
     $                        'Non-unit', JB, N-J-JB+1, CONE, A( J, J ),
     $                        LDA, A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
*
         ELSE
*
*           Compute the Cholesky factorization A = L*L'.
*
            DO 20 J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL ZHERK( 'Lower', 'No transpose', JB, J-1, -ONE,
     $                     A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL ZPOTF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block column.
*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose',
     $                        N-J-JB+1, JB, J-1, -CONE, A( J+JB, 1 ),
     $                        LDA, A( J, 1 ), LDA, CONE, A( J+JB, J ),
     $                        LDA )
                  CALL ZTRSM( 'Right', 'Lower', 'Conjugate transpose',
     $                        'Non-unit', N-J-JB+1, JB, CONE, A( J, J ),
     $                        LDA, A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = INFO + J - 1
*
   40 CONTINUE
      RETURN
*
*     End of ZPOTRF
*
      END

      SUBROUTINE ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZTRTI2 computes the inverse of a complex upper or lower triangular
*  matrix.
*
*  This is the Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading n by n upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      COMPLEX*16         AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZSCAL, ZTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL ZTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL ZSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL ZTRMV( 'Lower', 'No transpose', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL ZSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTRTI2
*
      END

      SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZTRTRI computes the inverse of a complex upper or lower triangular
*  matrix A.
*
*  This is the Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*               matrix is singular and its inverse can not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZTRMM, ZTRSM, ZTRTI2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL ZTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL ZTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              Compute inverse of current diagonal block
*
               CALL ZTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL ZTRMM( 'Left', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL ZTRSM( 'Right', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL ZTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZTRTRI
*
      END
      SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,
     $                   LIWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSTEDC computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the divide and conquer method.
*  The eigenvectors of a full or band real symmetric matrix can also be
*  found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this
*  matrix to tridiagonal form.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.  See DLAED3 for details.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'I':  Compute eigenvectors of tridiagonal matrix also.
*          = 'V':  Compute eigenvectors of original dense symmetric
*                  matrix also.  On entry, Z contains the orthogonal
*                  matrix used to reduce the original matrix to
*                  tridiagonal form.
*
*  N       (input) INTEGER
*          The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the subdiagonal elements of the tridiagonal matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
*          On entry, if COMPZ = 'V', then Z contains the orthogonal
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original symmetric matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If  COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1.
*          If eigenvectors are desired, then LDZ >= max(1,N).
*
*  WORK    (workspace/output) DOUBLE PRECISION array,
*                                         dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.
*          If COMPZ = 'V' and N > 1 then LWORK must be at least
*                         ( 1 + 3*N + 2*N*lg N + 3*N**2 ),
*                         where lg( N ) = smallest integer k such
*                         that 2**k >= N.
*          If COMPZ = 'I' and N > 1 then LWORK must be at least
*                         ( 1 + 4*N + N**2 ).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.
*          If COMPZ = 'V' and N > 1 then LIWORK must be at least
*                         ( 6 + 6*N + 5*N*lg N ).
*          If COMPZ = 'I' and N > 1 then LIWORK must be at least
*                         ( 3 + 5*N ).
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  The algorithm failed to compute an eigenvalue while
*                working on the submatrix lying in rows and columns
*                INFO/(N+1) through mod(INFO,N+1).
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*  Modified by Francoise Tisseur, University of Tennessee.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            DTRTRW, END, I, ICOMPZ, II, J, K, LGN, LIWMIN,
     $                   LWMIN, M, SMLSIZ, START, STOREZ
      DOUBLE PRECISION   EPS, ORGNRM, P, TINY
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLAED0, DLASCL, DLASET, DLASRT,
     $                   DSTEQR, DSTERF, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, MAX, MOD, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( N.LE.1 .OR. ICOMPZ.LE.0 ) THEN
         LIWMIN = 1
         LWMIN = 1
      ELSE
         LGN = INT( LOG( DBLE( N ) ) / LOG( TWO ) )
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( ICOMPZ.EQ.1 ) THEN
            LWMIN = 1 + 3*N + 2*N*LGN + 3*N**2
            LIWMIN = 6 + 6*N + 5*N*LGN
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            LWMIN = 1 + 4*N + N**2
            LIWMIN = 3 + 5*N
         END IF
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -8
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.NE.0 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
      SMLSIZ = ILAENV( 9, 'DSTEDC', ' ', 0, 0, 0, 0 )
*
*     If the following conditional clause is removed, then the routine
*     will use the Divide and Conquer routine to compute only the
*     eigenvalues, which requires (3N + 3N**2) real workspace and
*     (2 + 5N + 2N lg(N)) integer workspace.
*     Since on many architectures DSTERF is much faster than any other
*     algorithm for finding eigenvalues only, it is used here
*     as the default.
*
*     If COMPZ = 'N', use DSTERF to compute the eigenvalues.
*
      IF( ICOMPZ.EQ.0 ) THEN
         CALL DSTERF( N, D, E, INFO )
         RETURN
      END IF
*
*     If N is smaller than the minimum divide size (SMLSIZ+1), then
*     solve the problem with another solver.
*
      IF( N.LE.SMLSIZ ) THEN
         IF( ICOMPZ.EQ.0 ) THEN
            CALL DSTERF( N, D, E, INFO )
            RETURN
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            CALL DSTEQR( 'I', N, D, E, Z, LDZ, WORK, INFO )
            RETURN
         ELSE
            CALL DSTEQR( 'V', N, D, E, Z, LDZ, WORK, INFO )
            RETURN
         END IF
      END IF
*
*     If COMPZ = 'V', the Z matrix must be stored elsewhere for later
*     use.
*
      IF( ICOMPZ.EQ.1 ) THEN
         STOREZ = 1 + N*N
      ELSE
         STOREZ = 1
      END IF
*
      IF( ICOMPZ.EQ.2 ) THEN
         CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
      END IF
*
*     Scale.
*
      ORGNRM = DLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO )
     $   RETURN
*
      EPS = DLAMCH( 'Epsilon' )
*
      START = 1
*
*     while ( START <= N )
*
   10 CONTINUE
      IF( START.LE.N ) THEN
*
*     Let END be the position of the next subdiagonal entry such that
*     E( END ) <= TINY or END = N if no such subdiagonal exists.  The
*     matrix identified by the elements between START and END
*     constitutes an independent sub-problem.
*
         END = START
   20    CONTINUE
         IF( END.LT.N ) THEN
            TINY = EPS*SQRT( ABS( D( END ) ) )*SQRT( ABS( D( END+1 ) ) )
            IF( ABS( E( END ) ).GT.TINY ) THEN
               END = END + 1
               GO TO 20
            END IF
         END IF
*
*        (Sub) Problem determined.  Compute its size and solve it.
*
         M = END - START + 1
         IF( M.EQ.1 ) THEN
            START = END + 1
            GO TO 10
         END IF
         IF( M.GT.SMLSIZ ) THEN
            INFO = SMLSIZ
*
*           Scale.
*
            ORGNRM = DLANST( 'M', M, D( START ), E( START ) )
            CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M,
     $                   INFO )
            CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ),
     $                   M-1, INFO )
*
            IF( ICOMPZ.EQ.1 ) THEN
               DTRTRW = 1
            ELSE
               DTRTRW = START
            END IF
            CALL DLAED0( ICOMPZ, N, M, D( START ), E( START ),
     $                   Z( DTRTRW, START ), LDZ, WORK( 1 ), N,
     $                   WORK( STOREZ ), IWORK, INFO )
            IF( INFO.NE.0 ) THEN
               INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) +
     $                MOD( INFO, ( M+1 ) ) + START - 1
               RETURN
            END IF
*
*           Scale back.
*
            CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M,
     $                   INFO )
*
         ELSE
            IF( ICOMPZ.EQ.1 ) THEN
*
*     Since QR won't update a Z matrix which is larger than the
*     length of D, we must solve the sub-problem in a workspace and
*     then multiply back into Z.
*
               CALL DSTEQR( 'I', M, D( START ), E( START ), WORK, M,
     $                      WORK( M*M+1 ), INFO )
               CALL DLACPY( 'A', N, M, Z( 1, START ), LDZ,
     $                      WORK( STOREZ ), N )
               CALL DGEMM( 'N', 'N', N, M, M, ONE, WORK( STOREZ ), LDZ,
     $                     WORK, M, ZERO, Z( 1, START ), LDZ )
            ELSE IF( ICOMPZ.EQ.2 ) THEN
               CALL DSTEQR( 'I', M, D( START ), E( START ),
     $                      Z( START, START ), LDZ, WORK, INFO )
            ELSE
               CALL DSTERF( M, D( START ), E( START ), INFO )
            END IF
            IF( INFO.NE.0 ) THEN
               INFO = START*( N+1 ) + END
               RETURN
            END IF
         END IF
*
         START = END + 1
         GO TO 10
      END IF
*
*     endwhile
*
*     If the problem split any number of times, then the eigenvalues
*     will not be properly ordered.  Here we permute the eigenvalues
*     (and the associated eigenvectors) into ascending order.
*
      IF( M.NE.N ) THEN
         IF( ICOMPZ.EQ.0 ) THEN
*
*        Use Quick Sort
*
            CALL DLASRT( 'I', N, D, INFO )
*
         ELSE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
            DO 40 II = 2, N
               I = II - 1
               K = I
               P = D( I )
               DO 30 J = II, N
                  IF( D( J ).LT.P ) THEN
                     K = J
                     P = D( J )
                  END IF
   30          CONTINUE
               IF( K.NE.I ) THEN
                  D( K ) = D( I )
                  D( I ) = P
                  CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
               END IF
   40       CONTINUE
         END IF
      END IF
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of DSTEDC
*
      END
      SUBROUTINE DLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS,
     $                   WORK, IWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAED0 computes all eigenvalues and corresponding eigenvectors of a
*  symmetric tridiagonal matrix using the divide and conquer method.
*
*  Arguments
*  =========
*
*  ICOMPQ  (input) INTEGER
*          = 0:  Compute eigenvalues only.
*          = 1:  Compute eigenvectors of original dense symmetric matrix
*                also.  On entry, Q contains the orthogonal matrix used
*                to reduce the original matrix to tridiagonal form.
*          = 2:  Compute eigenvalues and eigenvectors of tridiagonal
*                matrix.
*
*  QSIZ   (input) INTEGER
*         The dimension of the orthogonal matrix used to reduce
*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D      (input/output) DOUBLE PRECISION array, dimension (N)
*         On entry, the main diagonal of the tridiagonal matrix.
*         On exit, its eigenvalues.
*
*  E      (input) DOUBLE PRECISION array, dimension (N-1)
*         The off-diagonal elements of the tridiagonal matrix.
*         On exit, E has been destroyed.
*
*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*         On entry, Q must contain an N-by-N orthogonal matrix.
*         If ICOMPQ = 0    Q is not referenced.
*         If ICOMPQ = 1    On entry, Q is a subset of the columns of the
*                          orthogonal matrix used to reduce the full
*                          matrix to tridiagonal form corresponding to
*                          the subset of the full matrix which is being
*                          decomposed at this time.
*         If ICOMPQ = 2    On entry, Q will be the identity matrix.
*                          On exit, Q contains the eigenvectors of the
*                          tridiagonal matrix.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  If eigenvectors are
*         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.
*
*  QSTORE (workspace) DOUBLE PRECISION array, dimension (LDQS, N)
*         Referenced only when ICOMPQ = 1.  Used to store parts of
*         the eigenvector matrix when the updating matrix multiplies
*         take place.
*
*  LDQS   (input) INTEGER
*         The leading dimension of the array QSTORE.  If ICOMPQ = 1,
*         then  LDQS >= max(1,N).  In any case,  LDQS >= 1.
*
*  WORK   (workspace) DOUBLE PRECISION array,
*         If ICOMPQ = 0 or 1, the dimension of WORK must be at least
*                     1 + 3*N + 2*N*lg N + 2*N**2
*                     ( lg( N ) = smallest integer k
*                                 such that 2^k >= N )
*         If ICOMPQ = 2, the dimension of WORK must be at least
*                     4*N + N**2.
*
*  IWORK  (workspace) INTEGER array,
*         If ICOMPQ = 0 or 1, the dimension of IWORK must be at least
*                        6 + 6*N + 5*N*lg N.
*                        ( lg( N ) = smallest integer k
*                                    such that 2^k >= N )
*         If ICOMPQ = 2, the dimension of IWORK must be at least
*                        3 + 5*N.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  The algorithm failed to compute an eigenvalue while
*                working on the submatrix lying in rows and columns
*                INFO/(N+1) through mod(INFO,N+1).
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CURLVL, CURPRB, CURR, I, IGIVCL, IGIVNM,
     $                   IGIVPT, INDXQ, IPERM, IPRMPT, IQ, IQPTR, IWREM,
     $                   J, K, LGN, MATSIZ, MSD2, SMLSIZ, SMM1, SPM1,
     $                   SPM2, SUBMAT, SUBPBS, TLVLS
      DOUBLE PRECISION   TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLAED1, DLAED7, DSTEQR,
     $                   XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.2 ) THEN
         INFO = -1
      ELSE IF( ( ICOMPQ.EQ.1 ) .AND. ( QSIZ.LT.MAX( 0, N ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDQS.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED0', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      SMLSIZ = ILAENV( 9, 'DLAED0', ' ', 0, 0, 0, 0 )
*
*     Determine the size and placement of the submatrices, and save in
*     the leading elements of IWORK.
*
      IWORK( 1 ) = N
      SUBPBS = 1
      TLVLS = 0
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.SMLSIZ ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         TLVLS = TLVLS + 1
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
   30 CONTINUE
*
*     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
*     using rank-1 modifications (cuts).
*
      SPM1 = SUBPBS - 1
      DO 40 I = 1, SPM1
         SUBMAT = IWORK( I ) + 1
         SMM1 = SUBMAT - 1
         D( SMM1 ) = D( SMM1 ) - ABS( E( SMM1 ) )
         D( SUBMAT ) = D( SUBMAT ) - ABS( E( SMM1 ) )
   40 CONTINUE
*
      INDXQ = 4*N + 3
      IF( ICOMPQ.NE.2 ) THEN
*
*        Set up workspaces for eigenvalues only/accumulate new vectors
*        routine
*
         TEMP = LOG( DBLE( N ) ) / LOG( TWO )
         LGN = INT( TEMP )
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IPRMPT = INDXQ + N + 1
         IPERM = IPRMPT + N*LGN
         IQPTR = IPERM + N*LGN
         IGIVPT = IQPTR + N + 2
         IGIVCL = IGIVPT + N*LGN
*
         IGIVNM = 1
         IQ = IGIVNM + 2*N*LGN
         IWREM = IQ + N**2 + 1
*
*        Initialize pointers
*
         DO 50 I = 0, SUBPBS
            IWORK( IPRMPT+I ) = 1
            IWORK( IGIVPT+I ) = 1
   50    CONTINUE
         IWORK( IQPTR ) = 1
      END IF
*
*     Solve each submatrix eigenproblem at the bottom of the divide and
*     conquer tree.
*
      CURR = 0
      DO 70 I = 0, SPM1
         IF( I.EQ.0 ) THEN
            SUBMAT = 1
            MATSIZ = IWORK( 1 )
         ELSE
            SUBMAT = IWORK( I ) + 1
            MATSIZ = IWORK( I+1 ) - IWORK( I )
         END IF
         IF( ICOMPQ.EQ.2 ) THEN
            CALL DSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ),
     $                   Q( SUBMAT, SUBMAT ), LDQ, WORK, INFO )
            IF( INFO.NE.0 )
     $         GO TO 130
         ELSE
            CALL DSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ),
     $                   WORK( IQ-1+IWORK( IQPTR+CURR ) ), MATSIZ, WORK,
     $                   INFO )
            IF( INFO.NE.0 )
     $         GO TO 130
            IF( ICOMPQ.EQ.1 ) THEN
               CALL DGEMM( 'N', 'N', QSIZ, MATSIZ, MATSIZ, ONE,
     $                     Q( 1, SUBMAT ), LDQ, WORK( IQ-1+IWORK( IQPTR+
     $                     CURR ) ), MATSIZ, ZERO, QSTORE( 1, SUBMAT ),
     $                     LDQS )
            END IF
            IWORK( IQPTR+CURR+1 ) = IWORK( IQPTR+CURR ) + MATSIZ**2
            CURR = CURR + 1
         END IF
         K = 1
         DO 60 J = SUBMAT, IWORK( I+1 )
            IWORK( INDXQ+J ) = K
            K = K + 1
   60    CONTINUE
   70 CONTINUE
*
*     Successively merge eigensystems of adjacent submatrices
*     into eigensystem for the corresponding larger matrix.
*
*     while ( SUBPBS > 1 )
*
      CURLVL = 1
   80 CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         SPM2 = SUBPBS - 2
         DO 90 I = 0, SPM2, 2
            IF( I.EQ.0 ) THEN
               SUBMAT = 1
               MATSIZ = IWORK( 2 )
               MSD2 = IWORK( 1 )
               CURPRB = 0
            ELSE
               SUBMAT = IWORK( I ) + 1
               MATSIZ = IWORK( I+2 ) - IWORK( I )
               MSD2 = MATSIZ / 2
               CURPRB = CURPRB + 1
            END IF
*
*     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
*     into an eigensystem of size MATSIZ.
*     DLAED1 is used only for the full eigensystem of a tridiagonal
*     matrix.
*     DLAED7 handles the cases in which eigenvalues only or eigenvalues
*     and eigenvectors of a full symmetric matrix (which was reduced to
*     tridiagonal form) are desired.
*
            IF( ICOMPQ.EQ.2 ) THEN
               CALL DLAED1( MATSIZ, D( SUBMAT ), Q( SUBMAT, SUBMAT ),
     $                      LDQ, IWORK( INDXQ+SUBMAT ),
     $                      E( SUBMAT+MSD2-1 ), MSD2, WORK,
     $                      IWORK( SUBPBS+1 ), INFO )
            ELSE
               CALL DLAED7( ICOMPQ, MATSIZ, QSIZ, TLVLS, CURLVL, CURPRB,
     $                      D( SUBMAT ), QSTORE( 1, SUBMAT ), LDQS,
     $                      IWORK( INDXQ+SUBMAT ), E( SUBMAT+MSD2-1 ),
     $                      MSD2, WORK( IQ ), IWORK( IQPTR ),
     $                      IWORK( IPRMPT ), IWORK( IPERM ),
     $                      IWORK( IGIVPT ), IWORK( IGIVCL ),
     $                      WORK( IGIVNM ), WORK( IWREM ),
     $                      IWORK( SUBPBS+1 ), INFO )
            END IF
            IF( INFO.NE.0 )
     $         GO TO 130
            IWORK( I / 2+1 ) = IWORK( I+2 )
   90    CONTINUE
         SUBPBS = SUBPBS / 2
         CURLVL = CURLVL + 1
         GO TO 80
      END IF
*
*     end while
*
*     Re-merge the eigenvalues/vectors which were deflated at the final
*     merge step.
*
      IF( ICOMPQ.EQ.1 ) THEN
         DO 100 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
            CALL DCOPY( QSIZ, QSTORE( 1, J ), 1, Q( 1, I ), 1 )
  100    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
      ELSE IF( ICOMPQ.EQ.2 ) THEN
         DO 110 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
            CALL DCOPY( N, Q( 1, J ), 1, WORK( N*I+1 ), 1 )
  110    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
         CALL DLACPY( 'A', N, N, WORK( N+1 ), N, Q, LDQ )
      ELSE
         DO 120 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
  120    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
      END IF
      GO TO 140
*
  130 CONTINUE
      INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
*
  140 CONTINUE
      RETURN
*
*     End of DLAED0
*
      END
      SUBROUTINE DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     December 23, 1999
*
*     .. Scalar Arguments ..
      INTEGER            I, INFO, N
      DOUBLE PRECISION   DLAM, RHO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DELTA( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  This subroutine computes the I-th updated eigenvalue of a symmetric
*  rank-one modification to a diagonal matrix whose elements are
*  given in the array d, and that
*
*             D(i) < D(j)  for  i < j
*
*  and that RHO > 0.  This is arranged by the calling routine, and is
*  no loss in generality.  The rank-one modified system is thus
*
*             diag( D )  +  RHO *  Z * Z_transpose.
*
*  where we assume the Euclidean norm of Z is 1.
*
*  The method consists of approximating the rational functions in the
*  secular equation by simpler interpolating rational functions.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The length of all arrays.
*
*  I      (input) INTEGER
*         The index of the eigenvalue to be computed.  1 <= I <= N.
*
*  D      (input) DOUBLE PRECISION array, dimension (N)
*         The original eigenvalues.  It is assumed that they are in
*         order, D(I) < D(J)  for I < J.
*
*  Z      (input) DOUBLE PRECISION array, dimension (N)
*         The components of the updating vector.
*
*  DELTA  (output) DOUBLE PRECISION array, dimension (N)
*         If N .ne. 1, DELTA contains (D(j) - lambda_I) in its  j-th
*         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA
*         contains the information necessary to construct the
*         eigenvectors.
*
*  RHO    (input) DOUBLE PRECISION
*         The scalar in the symmetric updating formula.
*
*  DLAM   (output) DOUBLE PRECISION
*         The computed lambda_I, the I-th updated eigenvalue.
*
*  INFO   (output) INTEGER
*         = 0:  successful exit
*         > 0:  if INFO = 1, the updating process failed.
*
*  Internal Parameters
*  ===================
*
*  Logical variable ORGATI (origin-at-i?) is used for distinguishing
*  whether D(i) or D(i+1) is treated as the origin.
*
*            ORGATI = .true.    origin at i
*            ORGATI = .false.   origin at i+1
*
*   Logical variable SWTCH3 (switch-for-3-poles?) is for noting
*   if we are working with THREE poles!
*
*   MAXIT is the maximum number of iterations allowed for each
*   eigenvalue.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ren-Cang Li, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0,
     $                   TEN = 10.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ORGATI, SWTCH, SWTCH3
      INTEGER            II, IIM1, IIP1, IP1, ITER, J, NITER
      DOUBLE PRECISION   A, B, C, DEL, DLTLB, DLTUB, DPHI, DPSI, DW,
     $                   EPS, ERRETM, ETA, MIDPT, PHI, PREW, PSI,
     $                   RHOINV, TAU, TEMP, TEMP1, W
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   ZZ( 3 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAED5, DLAED6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Since this routine is called in an inner loop, we do no argument
*     checking.
*
*     Quick return for N=1 and 2.
*
      INFO = 0
      IF( N.EQ.1 ) THEN
*
*         Presumably, I=1 upon entry
*
         DLAM = D( 1 ) + RHO*Z( 1 )*Z( 1 )
         DELTA( 1 ) = ONE
         RETURN
      END IF
      IF( N.EQ.2 ) THEN
         CALL DLAED5( I, D, Z, DELTA, RHO, DLAM )
         RETURN
      END IF
*
*     Compute machine epsilon
*
      EPS = DLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO
*
*     The case I = N
*
      IF( I.EQ.N ) THEN
*
*        Initialize some basic variables
*
         II = N - 1
         NITER = 1
*
*        Calculate initial guess
*
         MIDPT = RHO / TWO
*
*        If ||Z||_2 is not one, then TEMP should be set to
*        RHO * ||Z||_2^2 / TWO
*
         DO 10 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
   10    CONTINUE
*
         PSI = ZERO
         DO 20 J = 1, N - 2
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
   20    CONTINUE
*
         C = RHOINV + PSI
         W = C + Z( II )*Z( II ) / DELTA( II ) +
     $       Z( N )*Z( N ) / DELTA( N )
*
         IF( W.LE.ZERO ) THEN
            TEMP = Z( N-1 )*Z( N-1 ) / ( D( N )-D( N-1 )+RHO ) +
     $             Z( N )*Z( N ) / RHO
            IF( C.LE.TEMP ) THEN
               TAU = RHO
            ELSE
               DEL = D( N ) - D( N-1 )
               A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
               B = Z( N )*Z( N )*DEL
               IF( A.LT.ZERO ) THEN
                  TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
               ELSE
                  TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
               END IF
            END IF
*
*           It can be proved that
*               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO
*
            DLTLB = MIDPT
            DLTUB = RHO
         ELSE
            DEL = D( N ) - D( N-1 )
            A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            ELSE
               TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            END IF
*
*           It can be proved that
*               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2
*
            DLTLB = ZERO
            DLTUB = MIDPT
         END IF
*
         DO 30 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - TAU
   30    CONTINUE
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 40 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   40    CONTINUE
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +
     $            ABS( TAU )*( DPSI+DPHI )
*
         W = RHOINV + PHI + PSI
*
*        Test for convergence
*
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            DLAM = D( I ) + TAU
            GO TO 250
         END IF
*
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
*
*        Calculate the new step
*
         NITER = NITER + 1
         C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
         A = ( DELTA( N-1 )+DELTA( N ) )*W -
     $       DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
         B = DELTA( N-1 )*DELTA( N )*W
         IF( C.LT.ZERO )
     $      C = ABS( C )
         IF( C.EQ.ZERO ) THEN
*          ETA = B/A
*           ETA = RHO - TAU
            ETA = DLTUB - TAU
         ELSE IF( A.GE.ZERO ) THEN
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
*
*        Note, eta should be positive if w is negative, and
*        eta should be negative otherwise. However,
*        if for some reason caused by roundoff, eta*w > 0,
*        we simply use one Newton step instead. This way
*        will guarantee eta*w < 0.
*
         IF( W*ETA.GT.ZERO )
     $      ETA = -W / ( DPSI+DPHI )
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
         DO 50 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
   50    CONTINUE
*
         TAU = TAU + ETA
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 60 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   60    CONTINUE
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +
     $            ABS( TAU )*( DPSI+DPHI )
*
         W = RHOINV + PHI + PSI
*
*        Main loop to update the values of the array   DELTA
*
         ITER = NITER + 1
*
         DO 90 NITER = ITER, MAXIT
*
*           Test for convergence
*
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               DLAM = D( I ) + TAU
               GO TO 250
            END IF
*
            IF( W.LE.ZERO ) THEN
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            END IF
*
*           Calculate the new step
*
            C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
            A = ( DELTA( N-1 )+DELTA( N ) )*W -
     $          DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
            B = DELTA( N-1 )*DELTA( N )*W
            IF( A.GE.ZERO ) THEN
               ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
*
*           Note, eta should be positive if w is negative, and
*           eta should be negative otherwise. However,
*           if for some reason caused by roundoff, eta*w > 0,
*           we simply use one Newton step instead. This way
*           will guarantee eta*w < 0.
*
            IF( W*ETA.GT.ZERO )
     $         ETA = -W / ( DPSI+DPHI )
            TEMP = TAU + ETA
            IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               END IF
            END IF
            DO 70 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
   70       CONTINUE
*
            TAU = TAU + ETA
*
*           Evaluate PSI and the derivative DPSI
*
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 80 J = 1, II
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
   80       CONTINUE
            ERRETM = ABS( ERRETM )
*
*           Evaluate PHI and the derivative DPHI
*
            TEMP = Z( N ) / DELTA( N )
            PHI = Z( N )*TEMP
            DPHI = TEMP*TEMP
            ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +
     $               ABS( TAU )*( DPSI+DPHI )
*
            W = RHOINV + PHI + PSI
   90    CONTINUE
*
*        Return with INFO = 1, NITER = MAXIT and not converged
*
         INFO = 1
         DLAM = D( I ) + TAU
         GO TO 250
*
*        End for the case I = N
*
      ELSE
*
*        The case for I < N
*
         NITER = 1
         IP1 = I + 1
*
*        Calculate initial guess
*
         DEL = D( IP1 ) - D( I )
         MIDPT = DEL / TWO
         DO 100 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
  100    CONTINUE
*
         PSI = ZERO
         DO 110 J = 1, I - 1
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
  110    CONTINUE
*
         PHI = ZERO
         DO 120 J = N, I + 2, -1
            PHI = PHI + Z( J )*Z( J ) / DELTA( J )
  120    CONTINUE
         C = RHOINV + PSI + PHI
         W = C + Z( I )*Z( I ) / DELTA( I ) +
     $       Z( IP1 )*Z( IP1 ) / DELTA( IP1 )
*
         IF( W.GT.ZERO ) THEN
*
*           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
*
*           We choose d(i) as origin.
*
            ORGATI = .TRUE.
            A = C*DEL + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
            B = Z( I )*Z( I )*DEL
            IF( A.GT.ZERO ) THEN
               TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            ELSE
               TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            END IF
            DLTLB = ZERO
            DLTUB = MIDPT
         ELSE
*
*           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
*
*           We choose d(i+1) as origin.
*
            ORGATI = .FALSE.
            A = C*DEL - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
            B = Z( IP1 )*Z( IP1 )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
            ELSE
               TAU = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
            END IF
            DLTLB = -MIDPT
            DLTUB = ZERO
         END IF
*
         IF( ORGATI ) THEN
            DO 130 J = 1, N
               DELTA( J ) = ( D( J )-D( I ) ) - TAU
  130       CONTINUE
         ELSE
            DO 140 J = 1, N
               DELTA( J ) = ( D( J )-D( IP1 ) ) - TAU
  140       CONTINUE
         END IF
         IF( ORGATI ) THEN
            II = I
         ELSE
            II = I + 1
         END IF
         IIM1 = II - 1
         IIP1 = II + 1
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 150 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  150    CONTINUE
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         DPHI = ZERO
         PHI = ZERO
         DO 160 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  160    CONTINUE
*
         W = RHOINV + PHI + PSI
*
*        W is the value of the secular function with
*        its ii-th element removed.
*
         SWTCH3 = .FALSE.
         IF( ORGATI ) THEN
            IF( W.LT.ZERO )
     $         SWTCH3 = .TRUE.
         ELSE
            IF( W.GT.ZERO )
     $         SWTCH3 = .TRUE.
         END IF
         IF( II.EQ.1 .OR. II.EQ.N )
     $      SWTCH3 = .FALSE.
*
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = W + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +
     $            THREE*ABS( TEMP ) + ABS( TAU )*DW
*
*        Test for convergence
*
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            IF( ORGATI ) THEN
               DLAM = D( I ) + TAU
            ELSE
               DLAM = D( IP1 ) + TAU
            END IF
            GO TO 250
         END IF
*
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
*
*        Calculate the new step
*
         NITER = NITER + 1
         IF( .NOT.SWTCH3 ) THEN
            IF( ORGATI ) THEN
               C = W - DELTA( IP1 )*DW - ( D( I )-D( IP1 ) )*
     $             ( Z( I ) / DELTA( I ) )**2
            ELSE
               C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )*
     $             ( Z( IP1 ) / DELTA( IP1 ) )**2
            END IF
            A = ( DELTA( I )+DELTA( IP1 ) )*W -
     $          DELTA( I )*DELTA( IP1 )*DW
            B = DELTA( I )*DELTA( IP1 )*W
            IF( C.EQ.ZERO ) THEN
               IF( A.EQ.ZERO ) THEN
                  IF( ORGATI ) THEN
                     A = Z( I )*Z( I ) + DELTA( IP1 )*DELTA( IP1 )*
     $                   ( DPSI+DPHI )
                  ELSE
                     A = Z( IP1 )*Z( IP1 ) + DELTA( I )*DELTA( I )*
     $                   ( DPSI+DPHI )
                  END IF
               END IF
               ETA = B / A
            ELSE IF( A.LE.ZERO ) THEN
               ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
         ELSE
*
*           Interpolation using THREE most relevant poles
*
            TEMP = RHOINV + PSI + PHI
            IF( ORGATI ) THEN
               TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) -
     $             ( D( IIM1 )-D( IIP1 ) )*TEMP1
               ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
               ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*
     $                   ( ( DPSI-TEMP1 )+DPHI )
            ELSE
               TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) -
     $             ( D( IIP1 )-D( IIM1 ) )*TEMP1
               ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*
     $                   ( DPSI+( DPHI-TEMP1 ) )
               ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
            END IF
            ZZ( 2 ) = Z( II )*Z( II )
            CALL DLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA,
     $                   INFO )
            IF( INFO.NE.0 )
     $         GO TO 250
         END IF
*
*        Note, eta should be positive if w is negative, and
*        eta should be negative otherwise. However,
*        if for some reason caused by roundoff, eta*w > 0,
*        we simply use one Newton step instead. This way
*        will guarantee eta*w < 0.
*
         IF( W*ETA.GE.ZERO )
     $      ETA = -W / DW
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
*
         PREW = W
*
  170    CONTINUE
         DO 180 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
  180    CONTINUE
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 190 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  190    CONTINUE
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         DPHI = ZERO
         PHI = ZERO
         DO 200 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  200    CONTINUE
*
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +
     $            THREE*ABS( TEMP ) + ABS( TAU+ETA )*DW
*
         SWTCH = .FALSE.
         IF( ORGATI ) THEN
            IF( -W.GT.ABS( PREW ) / TEN )
     $         SWTCH = .TRUE.
         ELSE
            IF( W.GT.ABS( PREW ) / TEN )
     $         SWTCH = .TRUE.
         END IF
*
         TAU = TAU + ETA
*
*        Main loop to update the values of the array   DELTA
*
         ITER = NITER + 1
*
         DO 240 NITER = ITER, MAXIT
*
*           Test for convergence
*
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               IF( ORGATI ) THEN
                  DLAM = D( I ) + TAU
               ELSE
                  DLAM = D( IP1 ) + TAU
               END IF
               GO TO 250
            END IF
*
            IF( W.LE.ZERO ) THEN
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            END IF
*
*           Calculate the new step
*
            IF( .NOT.SWTCH3 ) THEN
               IF( .NOT.SWTCH ) THEN
                  IF( ORGATI ) THEN
                     C = W - DELTA( IP1 )*DW -
     $                   ( D( I )-D( IP1 ) )*( Z( I ) / DELTA( I ) )**2
                  ELSE
                     C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )*
     $                   ( Z( IP1 ) / DELTA( IP1 ) )**2
                  END IF
               ELSE
                  TEMP = Z( II ) / DELTA( II )
                  IF( ORGATI ) THEN
                     DPSI = DPSI + TEMP*TEMP
                  ELSE
                     DPHI = DPHI + TEMP*TEMP
                  END IF
                  C = W - DELTA( I )*DPSI - DELTA( IP1 )*DPHI
               END IF
               A = ( DELTA( I )+DELTA( IP1 ) )*W -
     $             DELTA( I )*DELTA( IP1 )*DW
               B = DELTA( I )*DELTA( IP1 )*W
               IF( C.EQ.ZERO ) THEN
                  IF( A.EQ.ZERO ) THEN
                     IF( .NOT.SWTCH ) THEN
                        IF( ORGATI ) THEN
                           A = Z( I )*Z( I ) + DELTA( IP1 )*
     $                         DELTA( IP1 )*( DPSI+DPHI )
                        ELSE
                           A = Z( IP1 )*Z( IP1 ) +
     $                         DELTA( I )*DELTA( I )*( DPSI+DPHI )
                        END IF
                     ELSE
                        A = DELTA( I )*DELTA( I )*DPSI +
     $                      DELTA( IP1 )*DELTA( IP1 )*DPHI
                     END IF
                  END IF
                  ETA = B / A
               ELSE IF( A.LE.ZERO ) THEN
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               ELSE
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               END IF
            ELSE
*
*              Interpolation using THREE most relevant poles
*
               TEMP = RHOINV + PSI + PHI
               IF( SWTCH ) THEN
                  C = TEMP - DELTA( IIM1 )*DPSI - DELTA( IIP1 )*DPHI
                  ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*DPSI
                  ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*DPHI
               ELSE
                  IF( ORGATI ) THEN
                     TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) -
     $                   ( D( IIM1 )-D( IIP1 ) )*TEMP1
                     ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                     ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*
     $                         ( ( DPSI-TEMP1 )+DPHI )
                  ELSE
                     TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) -
     $                   ( D( IIP1 )-D( IIM1 ) )*TEMP1
                     ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*
     $                         ( DPSI+( DPHI-TEMP1 ) )
                     ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
                  END IF
               END IF
               CALL DLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA,
     $                      INFO )
               IF( INFO.NE.0 )
     $            GO TO 250
            END IF
*
*           Note, eta should be positive if w is negative, and
*           eta should be negative otherwise. However,
*           if for some reason caused by roundoff, eta*w > 0,
*           we simply use one Newton step instead. This way
*           will guarantee eta*w < 0.
*
            IF( W*ETA.GE.ZERO )
     $         ETA = -W / DW
            TEMP = TAU + ETA
            IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               END IF
            END IF
*
            DO 210 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
  210       CONTINUE
*
            TAU = TAU + ETA
            PREW = W
*
*           Evaluate PSI and the derivative DPSI
*
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 220 J = 1, IIM1
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
  220       CONTINUE
            ERRETM = ABS( ERRETM )
*
*           Evaluate PHI and the derivative DPHI
*
            DPHI = ZERO
            PHI = ZERO
            DO 230 J = N, IIP1, -1
               TEMP = Z( J ) / DELTA( J )
               PHI = PHI + Z( J )*TEMP
               DPHI = DPHI + TEMP*TEMP
               ERRETM = ERRETM + PHI
  230       CONTINUE
*
            TEMP = Z( II ) / DELTA( II )
            DW = DPSI + DPHI + TEMP*TEMP
            TEMP = Z( II )*TEMP
            W = RHOINV + PHI + PSI + TEMP
            ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +
     $               THREE*ABS( TEMP ) + ABS( TAU )*DW
            IF( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN )
     $         SWTCH = .NOT.SWTCH
*
  240    CONTINUE
*
*        Return with INFO = 1, NITER = MAXIT and not converged
*
         INFO = 1
         IF( ORGATI ) THEN
            DLAM = D( I ) + TAU
         ELSE
            DLAM = D( IP1 ) + TAU
         END IF
*
      END IF
*
  250 CONTINUE
*
      RETURN
*
*     End of DLAED4
*
      END
      SUBROUTINE DLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, INFO, LDQ, N
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            INDXQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), Q( LDQ, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAED1 computes the updated eigensystem of a diagonal
*  matrix after modification by a rank-one symmetric matrix.  This
*  routine is used only for the eigenproblem which requires all
*  eigenvalues and eigenvectors of a tridiagonal matrix.  DLAED7 handles
*  the case in which eigenvalues only or eigenvalues and eigenvectors
*  of a full symmetric matrix (which was reduced to tridiagonal form)
*  are desired.
*
*    T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out)
*
*     where Z = Q'u, u is a vector of length N with ones in the
*     CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
*
*     The eigenvectors of the original matrix are stored in Q, and the
*     eigenvalues are in D.  The algorithm consists of three stages:
*
*        The first stage consists of deflating the size of the problem
*        when there are multiple eigenvalues or if there is a zero in
*        the Z vector.  For each such occurence the dimension of the
*        secular equation problem is reduced by one.  This stage is
*        performed by the routine DLAED2.
*
*        The second stage consists of calculating the updated
*        eigenvalues. This is done by finding the roots of the secular
*        equation via the routine DLAED4 (as called by DLAED3).
*        This routine also calculates the eigenvectors of the current
*        problem.
*
*        The final stage consists of computing the updated eigenvectors
*        directly using the updated eigenvalues.  The eigenvectors for
*        the current problem are multiplied with the eigenvectors from
*        the overall problem.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D      (input/output) DOUBLE PRECISION array, dimension (N)
*         On entry, the eigenvalues of the rank-1-perturbed matrix.
*         On exit, the eigenvalues of the repaired matrix.
*
*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*         On entry, the eigenvectors of the rank-1-perturbed matrix.
*         On exit, the eigenvectors of the repaired tridiagonal matrix.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  INDXQ  (input/output) INTEGER array, dimension (N)
*         On entry, the permutation which separately sorts the two
*         subproblems in D into ascending order.
*         On exit, the permutation which will reintegrate the
*         subproblems back into sorted order,
*         i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.
*
*  RHO    (input) DOUBLE PRECISION
*         The subdiagonal entry used to create the rank-1 modification.
*
*  CUTPNT (input) INTEGER
*         The location of the last eigenvalue in the leading sub-matrix.
*         min(1,N) <= CUTPNT <= N/2.
*
*  WORK   (workspace) DOUBLE PRECISION array, dimension (4*N + N**2)
*
*  IWORK  (workspace) INTEGER array, dimension (4*N)
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, an eigenvalue did not converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*  Modified by Francoise Tisseur, University of Tennessee.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            COLTYP, I, IDLMDA, INDX, INDXC, INDXP, IQ2, IS,
     $                   IW, IZ, K, N1, N2, ZPP1
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAED2, DLAED3, DLAMRG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( MIN( 1, N / 2 ).GT.CUTPNT .OR. ( N / 2 ).LT.CUTPNT ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED1', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     The following values are integer pointers which indicate
*     the portion of the workspace
*     used by a particular array in DLAED2 and DLAED3.
*
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ2 = IW + N
*
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
*
*
*     Form the z-vector which consists of the last row of Q_1 and the
*     first row of Q_2.
*
      CALL DCOPY( CUTPNT, Q( CUTPNT, 1 ), LDQ, WORK( IZ ), 1 )
      ZPP1 = CUTPNT + 1
      CALL DCOPY( N-CUTPNT, Q( ZPP1, ZPP1 ), LDQ, WORK( IZ+CUTPNT ), 1 )
*
*     Deflate eigenvalues.
*
      CALL DLAED2( K, N, CUTPNT, D, Q, LDQ, INDXQ, RHO, WORK( IZ ),
     $             WORK( IDLMDA ), WORK( IW ), WORK( IQ2 ),
     $             IWORK( INDX ), IWORK( INDXC ), IWORK( INDXP ),
     $             IWORK( COLTYP ), INFO )
*
      IF( INFO.NE.0 )
     $   GO TO 20
*
*     Solve Secular Equation.
*
      IF( K.NE.0 ) THEN
         IS = ( IWORK( COLTYP )+IWORK( COLTYP+1 ) )*CUTPNT +
     $        ( IWORK( COLTYP+1 )+IWORK( COLTYP+2 ) )*( N-CUTPNT ) + IQ2
         CALL DLAED3( K, N, CUTPNT, D, Q, LDQ, RHO, WORK( IDLMDA ),
     $                WORK( IQ2 ), IWORK( INDXC ), IWORK( COLTYP ),
     $                WORK( IW ), WORK( IS ), INFO )
         IF( INFO.NE.0 )
     $      GO TO 20
*
*     Prepare the INDXQ sorting permutation.
*
         N1 = K
         N2 = N - K
         CALL DLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         DO 10 I = 1, N
            INDXQ( I ) = I
   10    CONTINUE
      END IF
*
   20 CONTINUE
      RETURN
*
*     End of DLAED1
*
      END
      SUBROUTINE DLAED7( ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,
     $                   LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR,
     $                   PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N,
     $                   QSIZ, TLVLS
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),
     $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
      DOUBLE PRECISION   D( * ), GIVNUM( 2, * ), Q( LDQ, * ),
     $                   QSTORE( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAED7 computes the updated eigensystem of a diagonal
*  matrix after modification by a rank-one symmetric matrix. This
*  routine is used only for the eigenproblem which requires all
*  eigenvalues and optionally eigenvectors of a dense symmetric matrix
*  that has been reduced to tridiagonal form.  DLAED1 handles
*  the case in which all eigenvalues and eigenvectors of a symmetric
*  tridiagonal matrix are desired.
*
*    T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out)
*
*     where Z = Q'u, u is a vector of length N with ones in the
*     CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
*
*     The eigenvectors of the original matrix are stored in Q, and the
*     eigenvalues are in D.  The algorithm consists of three stages:
*
*        The first stage consists of deflating the size of the problem
*        when there are multiple eigenvalues or if there is a zero in
*        the Z vector.  For each such occurence the dimension of the
*        secular equation problem is reduced by one.  This stage is
*        performed by the routine DLAED8.
*
*        The second stage consists of calculating the updated
*        eigenvalues. This is done by finding the roots of the secular
*        equation via the routine DLAED4 (as called by DLAED9).
*        This routine also calculates the eigenvectors of the current
*        problem.
*
*        The final stage consists of computing the updated eigenvectors
*        directly using the updated eigenvalues.  The eigenvectors for
*        the current problem are multiplied with the eigenvectors from
*        the overall problem.
*
*  Arguments
*  =========
*
*  ICOMPQ  (input) INTEGER
*          = 0:  Compute eigenvalues only.
*          = 1:  Compute eigenvectors of original dense symmetric matrix
*                also.  On entry, Q contains the orthogonal matrix used
*                to reduce the original matrix to tridiagonal form.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  QSIZ   (input) INTEGER
*         The dimension of the orthogonal matrix used to reduce
*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
*
*  TLVLS  (input) INTEGER
*         The total number of merging levels in the overall divide and
*         conquer tree.
*
*  CURLVL (input) INTEGER
*         The current level in the overall merge routine,
*         0 <= CURLVL <= TLVLS.
*
*  CURPBM (input) INTEGER
*         The current problem in the current level in the overall
*         merge routine (counting from upper left to lower right).
*
*  D      (input/output) DOUBLE PRECISION array, dimension (N)
*         On entry, the eigenvalues of the rank-1-perturbed matrix.
*         On exit, the eigenvalues of the repaired matrix.
*
*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*         On entry, the eigenvectors of the rank-1-perturbed matrix.
*         On exit, the eigenvectors of the repaired tridiagonal matrix.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  INDXQ  (output) INTEGER array, dimension (N)
*         The permutation which will reintegrate the subproblem just
*         solved back into sorted order, i.e., D( INDXQ( I = 1, N ) )
*         will be in ascending order.
*
*  RHO    (input) DOUBLE PRECISION
*         The subdiagonal element used to create the rank-1
*         modification.
*
*  CUTPNT (input) INTEGER
*         Contains the location of the last eigenvalue in the leading
*         sub-matrix.  min(1,N) <= CUTPNT <= N.
*
*  QSTORE (input/output) DOUBLE PRECISION array, dimension (N**2+1)
*         Stores eigenvectors of submatrices encountered during
*         divide and conquer, packed together. QPTR points to
*         beginning of the submatrices.
*
*  QPTR   (input/output) INTEGER array, dimension (N+2)
*         List of indices pointing to beginning of submatrices stored
*         in QSTORE. The submatrices are numbered starting at the
*         bottom left of the divide and conquer tree, from left to
*         right and bottom to top.
*
*  PRMPTR (input) INTEGER array, dimension (N lg N)
*         Contains a list of pointers which indicate where in PERM a
*         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
*         indicates the size of the permutation and also the size of
*         the full, non-deflated problem.
*
*  PERM   (input) INTEGER array, dimension (N lg N)
*         Contains the permutations (from deflation and sorting) to be
*         applied to each eigenblock.
*
*  GIVPTR (input) INTEGER array, dimension (N lg N)
*         Contains a list of pointers which indicate where in GIVCOL a
*         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
*         indicates the number of Givens rotations.
*
*  GIVCOL (input) INTEGER array, dimension (2, N lg N)
*         Each pair of numbers indicates a pair of columns to take place
*         in a Givens rotation.
*
*  GIVNUM (input) DOUBLE PRECISION array, dimension (2, N lg N)
*         Each number indicates the S value to be used in the
*         corresponding Givens rotation.
*
*  WORK   (workspace) DOUBLE PRECISION array, dimension (3*N+QSIZ*N)
*
*  IWORK  (workspace) INTEGER array, dimension (4*N)
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, an eigenvalue did not converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            COLTYP, CURR, I, IDLMDA, INDX, INDXC, INDXP,
     $                   IQ2, IS, IW, IZ, K, LDQ2, N1, N2, PTR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLAED8, DLAED9, DLAEDA, DLAMRG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPQ.EQ.1 .AND. QSIZ.LT.N ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( MIN( 1, N ).GT.CUTPNT .OR. N.LT.CUTPNT ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED7', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     The following values are for bookkeeping purposes only.  They are
*     integer pointers which indicate the portion of the workspace
*     used by a particular array in DLAED8 and DLAED9.
*
      IF( ICOMPQ.EQ.1 ) THEN
         LDQ2 = QSIZ
      ELSE
         LDQ2 = N
      END IF
*
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ2 = IW + N
      IS = IQ2 + N*LDQ2
*
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
*
*     Form the z-vector which consists of the last row of Q_1 and the
*     first row of Q_2.
*
      PTR = 1 + 2**TLVLS
      DO 10 I = 1, CURLVL - 1
         PTR = PTR + 2**( TLVLS-I )
   10 CONTINUE
      CURR = PTR + CURPBM
      CALL DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,
     $             GIVCOL, GIVNUM, QSTORE, QPTR, WORK( IZ ),
     $             WORK( IZ+N ), INFO )
*
*     When solving the final problem, we no longer need the stored data,
*     so we will overwrite the data from this level onto the previously
*     used storage space.
*
      IF( CURLVL.EQ.TLVLS ) THEN
         QPTR( CURR ) = 1
         PRMPTR( CURR ) = 1
         GIVPTR( CURR ) = 1
      END IF
*
*     Sort and Deflate eigenvalues.
*
      CALL DLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO, CUTPNT,
     $             WORK( IZ ), WORK( IDLMDA ), WORK( IQ2 ), LDQ2,
     $             WORK( IW ), PERM( PRMPTR( CURR ) ), GIVPTR( CURR+1 ),
     $             GIVCOL( 1, GIVPTR( CURR ) ),
     $             GIVNUM( 1, GIVPTR( CURR ) ), IWORK( INDXP ),
     $             IWORK( INDX ), INFO )
      PRMPTR( CURR+1 ) = PRMPTR( CURR ) + N
      GIVPTR( CURR+1 ) = GIVPTR( CURR+1 ) + GIVPTR( CURR )
*
*     Solve Secular Equation.
*
      IF( K.NE.0 ) THEN
         CALL DLAED9( K, 1, K, N, D, WORK( IS ), K, RHO, WORK( IDLMDA ),
     $                WORK( IW ), QSTORE( QPTR( CURR ) ), K, INFO )
         IF( INFO.NE.0 )
     $      GO TO 30
         IF( ICOMPQ.EQ.1 ) THEN
            CALL DGEMM( 'N', 'N', QSIZ, K, K, ONE, WORK( IQ2 ), LDQ2,
     $                  QSTORE( QPTR( CURR ) ), K, ZERO, Q, LDQ )
         END IF
         QPTR( CURR+1 ) = QPTR( CURR ) + K**2
*
*     Prepare the INDXQ sorting permutation.
*
         N1 = K
         N2 = N - K
         CALL DLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         QPTR( CURR+1 ) = QPTR( CURR )
         DO 20 I = 1, N
            INDXQ( I ) = I
   20    CONTINUE
      END IF
*
   30 CONTINUE
      RETURN
*
*     End of DLAED7
*
      END
      SUBROUTINE DLAED5( I, D, Z, DELTA, RHO, DLAM )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            I
      DOUBLE PRECISION   DLAM, RHO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( 2 ), DELTA( 2 ), Z( 2 )
*     ..
*
*  Purpose
*  =======
*
*  This subroutine computes the I-th eigenvalue of a symmetric rank-one
*  modification of a 2-by-2 diagonal matrix
*
*             diag( D )  +  RHO *  Z * transpose(Z) .
*
*  The diagonal elements in the array D are assumed to satisfy
*
*             D(i) < D(j)  for  i < j .
*
*  We also assume RHO > 0 and that the Euclidean norm of the vector
*  Z is one.
*
*  Arguments
*  =========
*
*  I      (input) INTEGER
*         The index of the eigenvalue to be computed.  I = 1 or I = 2.
*
*  D      (input) DOUBLE PRECISION array, dimension (2)
*         The original eigenvalues.  We assume D(1) < D(2).
*
*  Z      (input) DOUBLE PRECISION array, dimension (2)
*         The components of the updating vector.
*
*  DELTA  (output) DOUBLE PRECISION array, dimension (2)
*         The vector DELTA contains the information necessary
*         to construct the eigenvectors.
*
*  RHO    (input) DOUBLE PRECISION
*         The scalar in the symmetric updating formula.
*
*  DLAM   (output) DOUBLE PRECISION
*         The computed lambda_I, the I-th updated eigenvalue.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ren-Cang Li, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, FOUR
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   FOUR = 4.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   B, C, DEL, TAU, TEMP, W
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      DEL = D( 2 ) - D( 1 )
      IF( I.EQ.1 ) THEN
         W = ONE + TWO*RHO*( Z( 2 )*Z( 2 )-Z( 1 )*Z( 1 ) ) / DEL
         IF( W.GT.ZERO ) THEN
            B = DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 1 )*Z( 1 )*DEL
*
*           B > ZERO, always
*
            TAU = TWO*C / ( B+SQRT( ABS( B*B-FOUR*C ) ) )
            DLAM = D( 1 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / TAU
            DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
         ELSE
            B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 2 )*Z( 2 )*DEL
            IF( B.GT.ZERO ) THEN
               TAU = -TWO*C / ( B+SQRT( B*B+FOUR*C ) )
            ELSE
               TAU = ( B-SQRT( B*B+FOUR*C ) ) / TWO
            END IF
            DLAM = D( 2 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
            DELTA( 2 ) = -Z( 2 ) / TAU
         END IF
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      ELSE
*
*     Now I=2
*
         B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
         C = RHO*Z( 2 )*Z( 2 )*DEL
         IF( B.GT.ZERO ) THEN
            TAU = ( B+SQRT( B*B+FOUR*C ) ) / TWO
         ELSE
            TAU = TWO*C / ( -B+SQRT( B*B+FOUR*C ) )
         END IF
         DLAM = D( 2 ) + TAU
         DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
         DELTA( 2 ) = -Z( 2 ) / TAU
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      END IF
      RETURN
*
*     End OF DLAED5
*
      END
      SUBROUTINE DLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      LOGICAL            ORGATI
      INTEGER            INFO, KNITER
      DOUBLE PRECISION   FINIT, RHO, TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( 3 ), Z( 3 )
*     ..
*
*  Purpose
*  =======
*
*  DLAED6 computes the positive or negative root (closest to the origin)
*  of
*                   z(1)        z(2)        z(3)
*  f(x) =   rho + --------- + ---------- + ---------
*                  d(1)-x      d(2)-x      d(3)-x
*
*  It is assumed that
*
*        if ORGATI = .true. the root is between d(2) and d(3);
*        otherwise it is between d(1) and d(2)
*
*  This routine will be called by DLAED4 when necessary. In most cases,
*  the root sought is the smallest in magnitude, though it might not be
*  in some extremely rare situations.
*
*  Arguments
*  =========
*
*  KNITER       (input) INTEGER
*               Refer to DLAED4 for its significance.
*
*  ORGATI       (input) LOGICAL
*               If ORGATI is true, the needed root is between d(2) and
*               d(3); otherwise it is between d(1) and d(2).  See
*               DLAED4 for further details.
*
*  RHO          (input) DOUBLE PRECISION
*               Refer to the equation f(x) above.
*
*  D            (input) DOUBLE PRECISION array, dimension (3)
*               D satisfies d(1) < d(2) < d(3).
*
*  Z            (input) DOUBLE PRECISION array, dimension (3)
*               Each of the elements in z must be positive.
*
*  FINIT        (input) DOUBLE PRECISION
*               The value of f at 0. It is more accurate than the one
*               evaluated inside this routine (if someone wants to do
*               so).
*
*  TAU          (output) DOUBLE PRECISION
*               The root of the equation f(x).
*
*  INFO         (output) INTEGER
*               = 0: successful exit
*               > 0: if INFO = 1, failure to converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ren-Cang Li, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 20 )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DSCALE( 3 ), ZSCALE( 3 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, SCALE
      INTEGER            I, ITER, NITER
      DOUBLE PRECISION   A, B, BASE, C, DDF, DF, EPS, ERRETM, ETA, F,
     $                   FC, SCLFAC, SCLINV, SMALL1, SMALL2, SMINV1,
     $                   SMINV2, TEMP, TEMP1, TEMP2, TEMP3, TEMP4
*     ..
*     .. Save statement ..
      SAVE               FIRST, SMALL1, SMINV1, SMALL2, SMINV2, EPS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
      NITER = 1
      TAU = ZERO
      IF( KNITER.EQ.2 ) THEN
         IF( ORGATI ) THEN
            TEMP = ( D( 3 )-D( 2 ) ) / TWO
            C = RHO + Z( 1 ) / ( ( D( 1 )-D( 2 ) )-TEMP )
            A = C*( D( 2 )+D( 3 ) ) + Z( 2 ) + Z( 3 )
            B = C*D( 2 )*D( 3 ) + Z( 2 )*D( 3 ) + Z( 3 )*D( 2 )
         ELSE
            TEMP = ( D( 1 )-D( 2 ) ) / TWO
            C = RHO + Z( 3 ) / ( ( D( 3 )-D( 2 ) )-TEMP )
            A = C*( D( 1 )+D( 2 ) ) + Z( 1 ) + Z( 2 )
            B = C*D( 1 )*D( 2 ) + Z( 1 )*D( 2 ) + Z( 2 )*D( 1 )
         END IF
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         IF( C.EQ.ZERO ) THEN
            TAU = B / A
         ELSE IF( A.LE.ZERO ) THEN
            TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
         TEMP = RHO + Z( 1 ) / ( D( 1 )-TAU ) +
     $          Z( 2 ) / ( D( 2 )-TAU ) + Z( 3 ) / ( D( 3 )-TAU )
         IF( ABS( FINIT ).LE.ABS( TEMP ) )
     $      TAU = ZERO
      END IF
*
*     On first call to routine, get machine parameters for
*     possible scaling to avoid overflow
*
      IF( FIRST ) THEN
         EPS = DLAMCH( 'Epsilon' )
         BASE = DLAMCH( 'Base' )
         SMALL1 = BASE**( INT( LOG( DLAMCH( 'SafMin' ) ) / LOG( BASE ) /
     $            THREE ) )
         SMINV1 = ONE / SMALL1
         SMALL2 = SMALL1*SMALL1
         SMINV2 = SMINV1*SMINV1
         FIRST = .FALSE.
      END IF
*
*     Determine if scaling of inputs necessary to avoid overflow
*     when computing 1/TEMP**3
*
      IF( ORGATI ) THEN
         TEMP = MIN( ABS( D( 2 )-TAU ), ABS( D( 3 )-TAU ) )
      ELSE
         TEMP = MIN( ABS( D( 1 )-TAU ), ABS( D( 2 )-TAU ) )
      END IF
      SCALE = .FALSE.
      IF( TEMP.LE.SMALL1 ) THEN
         SCALE = .TRUE.
         IF( TEMP.LE.SMALL2 ) THEN
*
*        Scale up by power of radix nearest 1/SAFMIN**(2/3)
*
            SCLFAC = SMINV2
            SCLINV = SMALL2
         ELSE
*
*        Scale up by power of radix nearest 1/SAFMIN**(1/3)
*
            SCLFAC = SMINV1
            SCLINV = SMALL1
         END IF
*
*        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
*
         DO 10 I = 1, 3
            DSCALE( I ) = D( I )*SCLFAC
            ZSCALE( I ) = Z( I )*SCLFAC
   10    CONTINUE
         TAU = TAU*SCLFAC
      ELSE
*
*        Copy D and Z to DSCALE and ZSCALE
*
         DO 20 I = 1, 3
            DSCALE( I ) = D( I )
            ZSCALE( I ) = Z( I )
   20    CONTINUE
      END IF
*
      FC = ZERO
      DF = ZERO
      DDF = ZERO
      DO 30 I = 1, 3
         TEMP = ONE / ( DSCALE( I )-TAU )
         TEMP1 = ZSCALE( I )*TEMP
         TEMP2 = TEMP1*TEMP
         TEMP3 = TEMP2*TEMP
         FC = FC + TEMP1 / DSCALE( I )
         DF = DF + TEMP2
         DDF = DDF + TEMP3
   30 CONTINUE
      F = FINIT + TAU*FC
*
      IF( ABS( F ).LE.ZERO )
     $   GO TO 60
*
*        Iteration begins
*
*     It is not hard to see that
*
*           1) Iterations will go up monotonically
*              if FINIT < 0;
*
*           2) Iterations will go down monotonically
*              if FINIT > 0.
*
      ITER = NITER + 1
*
      DO 50 NITER = ITER, MAXIT
*
         IF( ORGATI ) THEN
            TEMP1 = DSCALE( 2 ) - TAU
            TEMP2 = DSCALE( 3 ) - TAU
         ELSE
            TEMP1 = DSCALE( 1 ) - TAU
            TEMP2 = DSCALE( 2 ) - TAU
         END IF
         A = ( TEMP1+TEMP2 )*F - TEMP1*TEMP2*DF
         B = TEMP1*TEMP2*F
         C = F - ( TEMP1+TEMP2 )*DF + TEMP1*TEMP2*DDF
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         IF( C.EQ.ZERO ) THEN
            ETA = B / A
         ELSE IF( A.LE.ZERO ) THEN
            ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
         IF( F*ETA.GE.ZERO ) THEN
            ETA = -F / DF
         END IF
*
         TEMP = ETA + TAU
         IF( ORGATI ) THEN
            IF( ETA.GT.ZERO .AND. TEMP.GE.DSCALE( 3 ) )
     $         ETA = ( DSCALE( 3 )-TAU ) / TWO
            IF( ETA.LT.ZERO .AND. TEMP.LE.DSCALE( 2 ) )
     $         ETA = ( DSCALE( 2 )-TAU ) / TWO
         ELSE
            IF( ETA.GT.ZERO .AND. TEMP.GE.DSCALE( 2 ) )
     $         ETA = ( DSCALE( 2 )-TAU ) / TWO
            IF( ETA.LT.ZERO .AND. TEMP.LE.DSCALE( 1 ) )
     $         ETA = ( DSCALE( 1 )-TAU ) / TWO
         END IF
         TAU = TAU + ETA
*
         FC = ZERO
         ERRETM = ZERO
         DF = ZERO
         DDF = ZERO
         DO 40 I = 1, 3
            TEMP = ONE / ( DSCALE( I )-TAU )
            TEMP1 = ZSCALE( I )*TEMP
            TEMP2 = TEMP1*TEMP
            TEMP3 = TEMP2*TEMP
            TEMP4 = TEMP1 / DSCALE( I )
            FC = FC + TEMP4
            ERRETM = ERRETM + ABS( TEMP4 )
            DF = DF + TEMP2
            DDF = DDF + TEMP3
   40    CONTINUE
         F = FINIT + TAU*FC
         ERRETM = EIGHT*( ABS( FINIT )+ABS( TAU )*ERRETM ) +
     $            ABS( TAU )*DF
         IF( ABS( F ).LE.EPS*ERRETM )
     $      GO TO 60
   50 CONTINUE
      INFO = 1
   60 CONTINUE
*
*     Undo scaling
*
      IF( SCALE )
     $   TAU = TAU*SCLINV
      RETURN
*
*     End of DLAED6
*
      END
      SUBROUTINE DLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMDA, W,
     $                   Q2, INDX, INDXC, INDXP, COLTYP, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDQ, N, N1
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ),
     $                   INDXQ( * )
      DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),
     $                   W( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAED2 merges the two sets of eigenvalues together into a single
*  sorted set.  Then it tries to deflate the size of the problem.
*  There are two ways in which deflation can occur:  when two or more
*  eigenvalues are close together or if there is a tiny entry in the
*  Z vector.  For each such occurrence the order of the related secular
*  equation problem is reduced by one.
*
*  Arguments
*  =========
*
*  K      (output) INTEGER
*         The number of non-deflated eigenvalues, and the order of the
*         related secular equation. 0 <= K <=N.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  N1     (input) INTEGER
*         The location of the last eigenvalue in the leading sub-matrix.
*         min(1,N) <= N1 <= N/2.
*
*  D      (input/output) DOUBLE PRECISION array, dimension (N)
*         On entry, D contains the eigenvalues of the two submatrices to
*         be combined.
*         On exit, D contains the trailing (N-K) updated eigenvalues
*         (those which were deflated) sorted into increasing order.
*
*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*         On entry, Q contains the eigenvectors of two submatrices in
*         the two square blocks with corners at (1,1), (N1,N1)
*         and (N1+1, N1+1), (N,N).
*         On exit, Q contains the trailing (N-K) updated eigenvectors
*         (those which were deflated) in its last N-K columns.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  INDXQ  (input/output) INTEGER array, dimension (N)
*         The permutation which separately sorts the two sub-problems
*         in D into ascending order.  Note that elements in the second
*         half of this permutation must first have N1 added to their
*         values. Destroyed on exit.
*
*  RHO    (input/output) DOUBLE PRECISION
*         On entry, the off-diagonal element associated with the rank-1
*         cut which originally split the two submatrices which are now
*         being recombined.
*         On exit, RHO has been modified to the value required by
*         DLAED3.
*
*  Z      (input) DOUBLE PRECISION array, dimension (N)
*         On entry, Z contains the updating vector (the last
*         row of the first sub-eigenvector matrix and the first row of
*         the second sub-eigenvector matrix).
*         On exit, the contents of Z have been destroyed by the updating
*         process.
*
*  DLAMDA (output) DOUBLE PRECISION array, dimension (N)
*         A copy of the first K eigenvalues which will be used by
*         DLAED3 to form the secular equation.
*
*  W      (output) DOUBLE PRECISION array, dimension (N)
*         The first k values of the final deflation-altered z-vector
*         which will be passed to DLAED3.
*
*  Q2     (output) DOUBLE PRECISION array, dimension (N1**2+(N-N1)**2)
*         A copy of the first K eigenvectors which will be used by
*         DLAED3 in a matrix multiply (DGEMM) to solve for the new
*         eigenvectors.
*
*  INDX   (workspace) INTEGER array, dimension (N)
*         The permutation used to sort the contents of DLAMDA into
*         ascending order.
*
*  INDXC  (output) INTEGER array, dimension (N)
*         The permutation used to arrange the columns of the deflated
*         Q matrix into three groups:  the first group contains non-zero
*         elements only at and above N1, the second contains
*         non-zero elements only below N1, and the third is dense.
*
*  INDXP  (workspace) INTEGER array, dimension (N)
*         The permutation used to place deflated values of D at the end
*         of the array.  INDXP(1:K) points to the nondeflated D-values
*         and INDXP(K+1:N) points to the deflated eigenvalues.
*
*  COLTYP (workspace/output) INTEGER array, dimension (N)
*         During execution, a label which will indicate which of the
*         following types a column in the Q2 matrix is:
*         1 : non-zero in the upper half only;
*         2 : dense;
*         3 : non-zero in the lower half only;
*         4 : deflated.
*         On exit, COLTYP(i) is the number of columns of type i,
*         for i=1 to 4 only.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*  Modified by Francoise Tisseur, University of Tennessee.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0,
     $                   TWO = 2.0D0, EIGHT = 8.0D0 )
*     ..
*     .. Local Arrays ..
      INTEGER            CTOT( 4 ), PSM( 4 )
*     ..
*     .. Local Scalars ..
      INTEGER            CT, I, IMAX, IQ1, IQ2, J, JMAX, JS, K2, N1P1,
     $                   N2, NJ, PJ
      DOUBLE PRECISION   C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           IDAMAX, DLAMCH, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DLAMRG, DROT, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( MIN( 1, ( N / 2 ) ).GT.N1 .OR. ( N / 2 ).LT.N1 ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
*     Normalize z so that norm(z) = 1.  Since z is the concatenation of
*     two normalized vectors, norm2(z) = sqrt(2).
*
      T = ONE / SQRT( TWO )
      CALL DSCAL( N, T, Z, 1 )
*
*     RHO = ABS( norm(z)**2 * RHO )
*
      RHO = ABS( TWO*RHO )
*
*     Sort the eigenvalues into increasing order
*
      DO 10 I = N1P1, N
         INDXQ( I ) = INDXQ( I ) + N1
   10 CONTINUE
*
*     re-integrate the deflated parts from the last pass
*
      DO 20 I = 1, N
         DLAMDA( I ) = D( INDXQ( I ) )
   20 CONTINUE
      CALL DLAMRG( N1, N2, DLAMDA, 1, 1, INDXC )
      DO 30 I = 1, N
         INDX( I ) = INDXQ( INDXC( I ) )
   30 CONTINUE
*
*     Calculate the allowable deflation tolerance
*
      IMAX = IDAMAX( N, Z, 1 )
      JMAX = IDAMAX( N, D, 1 )
      EPS = DLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX ) ) )
*
*     If the rank-1 modifier is small enough, no more needs to be done
*     except to reorganize Q so that its columns correspond with the
*     elements in D.
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         IQ2 = 1
         DO 40 J = 1, N
            I = INDX( J )
            CALL DCOPY( N, Q( 1, I ), 1, Q2( IQ2 ), 1 )
            DLAMDA( J ) = D( I )
            IQ2 = IQ2 + N
   40    CONTINUE
         CALL DLACPY( 'A', N, N, Q2, N, Q, LDQ )
         CALL DCOPY( N, DLAMDA, 1, D, 1 )
         GO TO 190
      END IF
*
*     If there are multiple eigenvalues then the problem deflates.  Here
*     the number of equal eigenvalues are found.  As each equal
*     eigenvalue is found, an elementary reflector is computed to rotate
*     the corresponding eigensubspace so that the corresponding
*     components of Z are zero in this new basis.
*
      DO 50 I = 1, N1
         COLTYP( I ) = 1
   50 CONTINUE
      DO 60 I = N1P1, N
         COLTYP( I ) = 3
   60 CONTINUE
*
*
      K = 0
      K2 = N + 1
      DO 70 J = 1, N
         NJ = INDX( J )
         IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
*
*           Deflate due to small z component.
*
            K2 = K2 - 1
            COLTYP( NJ ) = 4
            INDXP( K2 ) = NJ
            IF( J.EQ.N )
     $         GO TO 100
         ELSE
            PJ = NJ
            GO TO 80
         END IF
   70 CONTINUE
   80 CONTINUE
      J = J + 1
      NJ = INDX( J )
      IF( J.GT.N )
     $   GO TO 100
      IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
*
*        Deflate due to small z component.
*
         K2 = K2 - 1
         COLTYP( NJ ) = 4
         INDXP( K2 ) = NJ
      ELSE
*
*        Check if eigenvalues are close enough to allow deflation.
*
         S = Z( PJ )
         C = Z( NJ )
*
*        Find sqrt(a**2+b**2) without overflow or
*        destructive underflow.
*
         TAU = DLAPY2( C, S )
         T = D( NJ ) - D( PJ )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
*
*           Deflation is possible.
*
            Z( NJ ) = TAU
            Z( PJ ) = ZERO
            IF( COLTYP( NJ ).NE.COLTYP( PJ ) )
     $         COLTYP( NJ ) = 2
            COLTYP( PJ ) = 4
            CALL DROT( N, Q( 1, PJ ), 1, Q( 1, NJ ), 1, C, S )
            T = D( PJ )*C**2 + D( NJ )*S**2
            D( NJ ) = D( PJ )*S**2 + D( NJ )*C**2
            D( PJ ) = T
            K2 = K2 - 1
            I = 1
   90       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( PJ ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = PJ
                  I = I + 1
                  GO TO 90
               ELSE
                  INDXP( K2+I-1 ) = PJ
               END IF
            ELSE
               INDXP( K2+I-1 ) = PJ
            END IF
            PJ = NJ
         ELSE
            K = K + 1
            DLAMDA( K ) = D( PJ )
            W( K ) = Z( PJ )
            INDXP( K ) = PJ
            PJ = NJ
         END IF
      END IF
      GO TO 80
  100 CONTINUE
*
*     Record the last eigenvalue.
*
      K = K + 1
      DLAMDA( K ) = D( PJ )
      W( K ) = Z( PJ )
      INDXP( K ) = PJ
*
*     Count up the total number of the various types of columns, then
*     form a permutation which positions the four column types into
*     four uniform groups (although one or more of these groups may be
*     empty).
*
      DO 110 J = 1, 4
         CTOT( J ) = 0
  110 CONTINUE
      DO 120 J = 1, N
         CT = COLTYP( J )
         CTOT( CT ) = CTOT( CT ) + 1
  120 CONTINUE
*
*     PSM(*) = Position in SubMatrix (of types 1 through 4)
*
      PSM( 1 ) = 1
      PSM( 2 ) = 1 + CTOT( 1 )
      PSM( 3 ) = PSM( 2 ) + CTOT( 2 )
      PSM( 4 ) = PSM( 3 ) + CTOT( 3 )
      K = N - CTOT( 4 )
*
*     Fill out the INDXC array so that the permutation which it induces
*     will place all type-1 columns first, all type-2 columns next,
*     then all type-3's, and finally all type-4's.
*
      DO 130 J = 1, N
         JS = INDXP( J )
         CT = COLTYP( JS )
         INDX( PSM( CT ) ) = JS
         INDXC( PSM( CT ) ) = J
         PSM( CT ) = PSM( CT ) + 1
  130 CONTINUE
*
*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
*     and Q2 respectively.  The eigenvalues/vectors which were not
*     deflated go into the first K slots of DLAMDA and Q2 respectively,
*     while those which were deflated go into the last N - K slots.
*
      I = 1
      IQ1 = 1
      IQ2 = 1 + ( CTOT( 1 )+CTOT( 2 ) )*N1
      DO 140 J = 1, CTOT( 1 )
         JS = INDX( I )
         CALL DCOPY( N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 )
         Z( I ) = D( JS )
         I = I + 1
         IQ1 = IQ1 + N1
  140 CONTINUE
*
      DO 150 J = 1, CTOT( 2 )
         JS = INDX( I )
         CALL DCOPY( N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 )
         CALL DCOPY( N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 )
         Z( I ) = D( JS )
         I = I + 1
         IQ1 = IQ1 + N1
         IQ2 = IQ2 + N2
  150 CONTINUE
*
      DO 160 J = 1, CTOT( 3 )
         JS = INDX( I )
         CALL DCOPY( N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 )
         Z( I ) = D( JS )
         I = I + 1
         IQ2 = IQ2 + N2
  160 CONTINUE
*
      IQ1 = IQ2
      DO 170 J = 1, CTOT( 4 )
         JS = INDX( I )
         CALL DCOPY( N, Q( 1, JS ), 1, Q2( IQ2 ), 1 )
         IQ2 = IQ2 + N
         Z( I ) = D( JS )
         I = I + 1
  170 CONTINUE
*
*     The deflated eigenvalues and their corresponding vectors go back
*     into the last N - K slots of D and Q respectively.
*
      CALL DLACPY( 'A', N, CTOT( 4 ), Q2( IQ1 ), N, Q( 1, K+1 ), LDQ )
      CALL DCOPY( N-K, Z( K+1 ), 1, D( K+1 ), 1 )
*
*     Copy CTOT into COLTYP for referencing in DLAED3.
*
      DO 180 J = 1, 4
         COLTYP( J ) = CTOT( J )
  180 CONTINUE
*
  190 CONTINUE
      RETURN
*
*     End of DLAED2
*
      END
      SUBROUTINE DLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX,
     $                   CTOT, W, S, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDQ, N, N1
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            CTOT( * ), INDX( * )
      DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),
     $                   S( * ), W( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAED3 finds the roots of the secular equation, as defined by the
*  values in D, W, and RHO, between 1 and K.  It makes the
*  appropriate calls to DLAED4 and then updates the eigenvectors by
*  multiplying the matrix of eigenvectors of the pair of eigensystems
*  being combined by the matrix of eigenvectors of the K-by-K system
*  which is solved here.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  K       (input) INTEGER
*          The number of terms in the rational function to be solved by
*          DLAED4.  K >= 0.
*
*  N       (input) INTEGER
*          The number of rows and columns in the Q matrix.
*          N >= K (deflation may result in N>K).
*
*  N1      (input) INTEGER
*          The location of the last eigenvalue in the leading submatrix.
*          min(1,N) <= N1 <= N/2.
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          D(I) contains the updated eigenvalues for
*          1 <= I <= K.
*
*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
*          Initially the first K columns are used as workspace.
*          On output the columns 1 to K contain
*          the updated eigenvectors.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  RHO     (input) DOUBLE PRECISION
*          The value of the parameter in the rank one update equation.
*          RHO >= 0 required.
*
*  DLAMDA  (input/output) DOUBLE PRECISION array, dimension (K)
*          The first K elements of this array contain the old roots
*          of the deflated updating problem.  These are the poles
*          of the secular equation. May be changed on output by
*          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,
*          Cray-2, or Cray C-90, as described above.
*
*  Q2      (input) DOUBLE PRECISION array, dimension (LDQ2, N)
*          The first K columns of this matrix contain the non-deflated
*          eigenvectors for the split problem.
*
*  INDX    (input) INTEGER array, dimension (N)
*          The permutation used to arrange the columns of the deflated
*          Q matrix into three groups (see DLAED2).
*          The rows of the eigenvectors found by DLAED4 must be likewise
*          permuted before the matrix multiply can take place.
*
*  CTOT    (input) INTEGER array, dimension (4)
*          A count of the total number of the various types of columns
*          in Q, as described in INDX.  The fourth column type is any
*          column which has been deflated.
*
*  W       (input/output) DOUBLE PRECISION array, dimension (K)
*          The first K elements of this array contain the components
*          of the deflation-adjusted updating vector. Destroyed on
*          output.
*
*  S       (workspace) DOUBLE PRECISION array, dimension (N1 + 1)*K
*          Will contain the eigenvectors of the repaired matrix which
*          will be multiplied by the previously accumulated eigenvectors
*          to update the system.
*
*  LDS     (input) INTEGER
*          The leading dimension of S.  LDS >= max(1,K).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, an eigenvalue did not converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*  Modified by Francoise Tisseur, University of Tennessee.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, IQ2, J, N12, N2, N23
      DOUBLE PRECISION   TEMP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLAED4, DLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.K ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED3', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( K.EQ.0 )
     $   RETURN
*
*     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
*     be computed with high relative accuracy (barring over/underflow).
*     This is a problem on machines without a guard digit in
*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
*     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
*     which on any of these machines zeros out the bottommost
*     bit of DLAMDA(I) if it is 1; this makes the subsequent
*     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
*     occurs. On binary machines with a guard digit (almost all
*     machines) it does not change DLAMDA(I) at all. On hexadecimal
*     and decimal machines with a guard digit, it slightly
*     changes the bottommost bits of DLAMDA(I). It does not account
*     for hexadecimal or decimal machines without guard digits
*     (we know of none). We use a subroutine call to compute
*     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
*     this code.
*
      DO 10 I = 1, K
         DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
   10 CONTINUE
*
      DO 20 J = 1, K
         CALL DLAED4( K, J, DLAMDA, W, Q( 1, J ), RHO, D( J ), INFO )
*
*        If the zero finder fails, the computation is terminated.
*
         IF( INFO.NE.0 )
     $      GO TO 120
   20 CONTINUE
*
      IF( K.EQ.1 )
     $   GO TO 110
      IF( K.EQ.2 ) THEN
         DO 30 J = 1, K
            W( 1 ) = Q( 1, J )
            W( 2 ) = Q( 2, J )
            II = INDX( 1 )
            Q( 1, J ) = W( II )
            II = INDX( 2 )
            Q( 2, J ) = W( II )
   30    CONTINUE
         GO TO 110
      END IF
*
*     Compute updated W.
*
      CALL DCOPY( K, W, 1, S, 1 )
*
*     Initialize W(I) = Q(I,I)
*
      CALL DCOPY( K, Q, LDQ+1, W, 1 )
      DO 60 J = 1, K
         DO 40 I = 1, J - 1
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   40    CONTINUE
         DO 50 I = J + 1, K
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   50    CONTINUE
   60 CONTINUE
      DO 70 I = 1, K
         W( I ) = SIGN( SQRT( -W( I ) ), S( I ) )
   70 CONTINUE
*
*     Compute eigenvectors of the modified rank-1 modification.
*
      DO 100 J = 1, K
         DO 80 I = 1, K
            S( I ) = W( I ) / Q( I, J )
   80    CONTINUE
         TEMP = DNRM2( K, S, 1 )
         DO 90 I = 1, K
            II = INDX( I )
            Q( I, J ) = S( II ) / TEMP
   90    CONTINUE
  100 CONTINUE
*
*     Compute the updated eigenvectors.
*
  110 CONTINUE
*
      N2 = N - N1
      N12 = CTOT( 1 ) + CTOT( 2 )
      N23 = CTOT( 2 ) + CTOT( 3 )
*
      CALL DLACPY( 'A', N23, K, Q( CTOT( 1 )+1, 1 ), LDQ, S, N23 )
      IQ2 = N1*N12 + 1
      IF( N23.NE.0 ) THEN
         CALL DGEMM( 'N', 'N', N2, K, N23, ONE, Q2( IQ2 ), N2, S, N23,
     $               ZERO, Q( N1+1, 1 ), LDQ )
      ELSE
         CALL DLASET( 'A', N2, K, ZERO, ZERO, Q( N1+1, 1 ), LDQ )
      END IF
*
      CALL DLACPY( 'A', N12, K, Q, LDQ, S, N12 )
      IF( N12.NE.0 ) THEN
         CALL DGEMM( 'N', 'N', N1, K, N12, ONE, Q2, N1, S, N12, ZERO, Q,
     $               LDQ )
      ELSE
         CALL DLASET( 'A', N1, K, ZERO, ZERO, Q( 1, 1 ), LDQ )
      END IF
*
*
  120 CONTINUE
      RETURN
*
*     End of DLAED3
*
      END
      SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            DTRD1, DTRD2, N1, N2
*     ..
*     .. Array Arguments ..
      INTEGER            INDEX( * )
      DOUBLE PRECISION   A( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAMRG will create a permutation list which will merge the elements
*  of A (which is composed of two independently sorted sets) into a
*  single set which is sorted in ascending order.
*
*  Arguments
*  =========
*
*  N1     (input) INTEGER
*  N2     (input) INTEGER
*         These arguements contain the respective lengths of the two
*         sorted lists to be merged.
*
*  A      (input) DOUBLE PRECISION array, dimension (N1+N2)
*         The first N1 elements of A contain a list of numbers which
*         are sorted in either ascending or descending order.  Likewise
*         for the final N2 elements.
*
*  DTRD1  (input) INTEGER
*  DTRD2  (input) INTEGER
*         These are the strides to be taken through the array A.
*         Allowable strides are 1 and -1.  They indicate whether a
*         subset of A is sorted in ascending (DTRDx = 1) or descending
*         (DTRDx = -1) order.
*
*  INDEX  (output) INTEGER array, dimension (N1+N2)
*         On exit this array will contain a permutation such that
*         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be
*         sorted in ascending order.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IND1, IND2, N1SV, N2SV
*     ..
*     .. Executable Statements ..
*
      N1SV = N1
      N2SV = N2
      IF( DTRD1.GT.0 ) THEN
         IND1 = 1
      ELSE
         IND1 = N1
      END IF
      IF( DTRD2.GT.0 ) THEN
         IND2 = 1 + N1
      ELSE
         IND2 = N1 + N2
      END IF
      I = 1
*     while ( (N1SV > 0) & (N2SV > 0) )
   10 CONTINUE
      IF( N1SV.GT.0 .AND. N2SV.GT.0 ) THEN
         IF( A( IND1 ).LE.A( IND2 ) ) THEN
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + DTRD1
            N1SV = N1SV - 1
         ELSE
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + DTRD2
            N2SV = N2SV - 1
         END IF
         GO TO 10
      END IF
*     end while
      IF( N1SV.EQ.0 ) THEN
         DO 20 N1SV = 1, N2SV
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + DTRD2
   20    CONTINUE
      ELSE
*     N2SV .EQ. 0
         DO 30 N2SV = 1, N1SV
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + DTRD1
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of DLAMRG
*
      END
      SUBROUTINE DLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO,
     $                   CUTPNT, Z, DLAMDA, Q2, LDQ2, W, PERM, GIVPTR,
     $                   GIVCOL, GIVNUM, INDXP, INDX, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, GIVPTR, ICOMPQ, INFO, K, LDQ, LDQ2, N,
     $                   QSIZ
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
     $                   INDXQ( * ), PERM( * )
      DOUBLE PRECISION   D( * ), DLAMDA( * ), GIVNUM( 2, * ),
     $                   Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAED8 merges the two sets of eigenvalues together into a single
*  sorted set.  Then it tries to deflate the size of the problem.
*  There are two ways in which deflation can occur:  when two or more
*  eigenvalues are close together or if there is a tiny element in the
*  Z vector.  For each such occurrence the order of the related secular
*  equation problem is reduced by one.
*
*  Arguments
*  =========
*
*  ICOMPQ  (input) INTEGER
*          = 0:  Compute eigenvalues only.
*          = 1:  Compute eigenvectors of original dense symmetric matrix
*                also.  On entry, Q contains the orthogonal matrix used
*                to reduce the original matrix to tridiagonal form.
*
*  K      (output) INTEGER
*         The number of non-deflated eigenvalues, and the order of the
*         related secular equation.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  QSIZ   (input) INTEGER
*         The dimension of the orthogonal matrix used to reduce
*         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
*
*  D      (input/output) DOUBLE PRECISION array, dimension (N)
*         On entry, the eigenvalues of the two submatrices to be
*         combined.  On exit, the trailing (N-K) updated eigenvalues
*         (those which were deflated) sorted into increasing order.
*
*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*         If ICOMPQ = 0, Q is not referenced.  Otherwise,
*         on entry, Q contains the eigenvectors of the partially solved
*         system which has been previously updated in matrix
*         multiplies with other partially solved eigensystems.
*         On exit, Q contains the trailing (N-K) updated eigenvectors
*         (those which were deflated) in its last N-K columns.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  INDXQ  (input) INTEGER array, dimension (N)
*         The permutation which separately sorts the two sub-problems
*         in D into ascending order.  Note that elements in the second
*         half of this permutation must first have CUTPNT added to
*         their values in order to be accurate.
*
*  RHO    (input/output) DOUBLE PRECISION
*         On entry, the off-diagonal element associated with the rank-1
*         cut which originally split the two submatrices which are now
*         being recombined.
*         On exit, RHO has been modified to the value required by
*         DLAED3.
*
*  CUTPNT (input) INTEGER
*         The location of the last eigenvalue in the leading
*         sub-matrix.  min(1,N) <= CUTPNT <= N.
*
*  Z      (input) DOUBLE PRECISION array, dimension (N)
*         On entry, Z contains the updating vector (the last row of
*         the first sub-eigenvector matrix and the first row of the
*         second sub-eigenvector matrix).
*         On exit, the contents of Z are destroyed by the updating
*         process.
*
*  DLAMDA (output) DOUBLE PRECISION array, dimension (N)
*         A copy of the first K eigenvalues which will be used by
*         DLAED3 to form the secular equation.
*
*  Q2     (output) DOUBLE PRECISION array, dimension (LDQ2,N)
*         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,
*         a copy of the first K eigenvectors which will be used by
*         DLAED7 in a matrix multiply (DGEMM) to update the new
*         eigenvectors.
*
*  LDQ2   (input) INTEGER
*         The leading dimension of the array Q2.  LDQ2 >= max(1,N).
*
*  W      (output) DOUBLE PRECISION array, dimension (N)
*         The first k values of the final deflation-altered z-vector and
*         will be passed to DLAED3.
*
*  PERM   (output) INTEGER array, dimension (N)
*         The permutations (from deflation and sorting) to be applied
*         to each eigenblock.
*
*  GIVPTR (output) INTEGER
*         The number of Givens rotations which took place in this
*         subproblem.
*
*  GIVCOL (output) INTEGER array, dimension (2, N)
*         Each pair of numbers indicates a pair of columns to take place
*         in a Givens rotation.
*
*  GIVNUM (output) DOUBLE PRECISION array, dimension (2, N)
*         Each number indicates the S value to be used in the
*         corresponding Givens rotation.
*
*  INDXP  (workspace) INTEGER array, dimension (N)
*         The permutation used to place deflated values of D at the end
*         of the array.  INDXP(1:K) points to the nondeflated D-values
*         and INDXP(K+1:N) points to the deflated eigenvalues.
*
*  INDX   (workspace) INTEGER array, dimension (N)
*         The permutation used to sort the contents of D into ascending
*         order.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0,
     $                   TWO = 2.0D0, EIGHT = 8.0D0 )
*     ..
*     .. Local Scalars ..
*
      INTEGER            I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2
      DOUBLE PRECISION   C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           IDAMAX, DLAMCH, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DLAMRG, DROT, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ICOMPQ.EQ.1 .AND. QSIZ.LT.N ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) THEN
         INFO = -10
      ELSE IF( LDQ2.LT.MAX( 1, N ) ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED8', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      N1 = CUTPNT
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
*     Normalize z so that norm(z) = 1
*
      T = ONE / SQRT( TWO )
      DO 10 J = 1, N
         INDX( J ) = J
   10 CONTINUE
      CALL DSCAL( N, T, Z, 1 )
      RHO = ABS( TWO*RHO )
*
*     Sort the eigenvalues into increasing order
*
      DO 20 I = CUTPNT + 1, N
         INDXQ( I ) = INDXQ( I ) + CUTPNT
   20 CONTINUE
      DO 30 I = 1, N
         DLAMDA( I ) = D( INDXQ( I ) )
         W( I ) = Z( INDXQ( I ) )
   30 CONTINUE
      I = 1
      J = CUTPNT + 1
      CALL DLAMRG( N1, N2, DLAMDA, 1, 1, INDX )
      DO 40 I = 1, N
         D( I ) = DLAMDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
   40 CONTINUE
*
*     Calculate the allowable deflation tolerence
*
      IMAX = IDAMAX( N, Z, 1 )
      JMAX = IDAMAX( N, D, 1 )
      EPS = DLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*ABS( D( JMAX ) )
*
*     If the rank-1 modifier is small enough, no more needs to be done
*     except to reorganize Q so that its columns correspond with the
*     elements in D.
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         IF( ICOMPQ.EQ.0 ) THEN
            DO 50 J = 1, N
               PERM( J ) = INDXQ( INDX( J ) )
   50       CONTINUE
         ELSE
            DO 60 J = 1, N
               PERM( J ) = INDXQ( INDX( J ) )
               CALL DCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
   60       CONTINUE
            CALL DLACPY( 'A', QSIZ, N, Q2( 1, 1 ), LDQ2, Q( 1, 1 ),
     $                   LDQ )
         END IF
         RETURN
      END IF
*
*     If there are multiple eigenvalues then the problem deflates.  Here
*     the number of equal eigenvalues are found.  As each equal
*     eigenvalue is found, an elementary reflector is computed to rotate
*     the corresponding eigensubspace so that the corresponding
*     components of Z are zero in this new basis.
*
      K = 0
      GIVPTR = 0
      K2 = N + 1
      DO 70 J = 1, N
         IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*           Deflate due to small z component.
*
            K2 = K2 - 1
            INDXP( K2 ) = J
            IF( J.EQ.N )
     $         GO TO 110
         ELSE
            JLAM = J
            GO TO 80
         END IF
   70 CONTINUE
   80 CONTINUE
      J = J + 1
      IF( J.GT.N )
     $   GO TO 100
      IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
*
*        Deflate due to small z component.
*
         K2 = K2 - 1
         INDXP( K2 ) = J
      ELSE
*
*        Check if eigenvalues are close enough to allow deflation.
*
         S = Z( JLAM )
         C = Z( J )
*
*        Find sqrt(a**2+b**2) without overflow or
*        destructive underflow.
*
         TAU = DLAPY2( C, S )
         T = D( J ) - D( JLAM )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
*
*           Deflation is possible.
*
            Z( J ) = TAU
            Z( JLAM ) = ZERO
*
*           Record the appropriate Givens rotation
*
            GIVPTR = GIVPTR + 1
            GIVCOL( 1, GIVPTR ) = INDXQ( INDX( JLAM ) )
            GIVCOL( 2, GIVPTR ) = INDXQ( INDX( J ) )
            GIVNUM( 1, GIVPTR ) = C
            GIVNUM( 2, GIVPTR ) = S
            IF( ICOMPQ.EQ.1 ) THEN
               CALL DROT( QSIZ, Q( 1, INDXQ( INDX( JLAM ) ) ), 1,
     $                    Q( 1, INDXQ( INDX( J ) ) ), 1, C, S )
            END IF
            T = D( JLAM )*C*C + D( J )*S*S
            D( J ) = D( JLAM )*S*S + D( J )*C*C
            D( JLAM ) = T
            K2 = K2 - 1
            I = 1
   90       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( JLAM ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = JLAM
                  I = I + 1
                  GO TO 90
               ELSE
                  INDXP( K2+I-1 ) = JLAM
               END IF
            ELSE
               INDXP( K2+I-1 ) = JLAM
            END IF
            JLAM = J
         ELSE
            K = K + 1
            W( K ) = Z( JLAM )
            DLAMDA( K ) = D( JLAM )
            INDXP( K ) = JLAM
            JLAM = J
         END IF
      END IF
      GO TO 80
  100 CONTINUE
*
*     Record the last eigenvalue.
*
      K = K + 1
      W( K ) = Z( JLAM )
      DLAMDA( K ) = D( JLAM )
      INDXP( K ) = JLAM
*
  110 CONTINUE
*
*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
*     and Q2 respectively.  The eigenvalues/vectors which were not
*     deflated go into the first K slots of DLAMDA and Q2 respectively,
*     while those which were deflated go into the last N - K slots.
*
      IF( ICOMPQ.EQ.0 ) THEN
         DO 120 J = 1, N
            JP = INDXP( J )
            DLAMDA( J ) = D( JP )
            PERM( J ) = INDXQ( INDX( JP ) )
  120    CONTINUE
      ELSE
         DO 130 J = 1, N
            JP = INDXP( J )
            DLAMDA( J ) = D( JP )
            PERM( J ) = INDXQ( INDX( JP ) )
            CALL DCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  130    CONTINUE
      END IF
*
*     The deflated eigenvalues and their corresponding vectors go back
*     into the last N - K slots of D and Q respectively.
*
      IF( K.LT.N ) THEN
         IF( ICOMPQ.EQ.0 ) THEN
            CALL DCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
         ELSE
            CALL DCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
            CALL DLACPY( 'A', QSIZ, N-K, Q2( 1, K+1 ), LDQ2,
     $                   Q( 1, K+1 ), LDQ )
         END IF
      END IF
*
      RETURN
*
*     End of DLAED8
*
      END
      SUBROUTINE DLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W,
     $                   S, LDS, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, KSTART, KSTOP, LDQ, LDS, N
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), S( LDS, * ),
     $                   W( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAED9 finds the roots of the secular equation, as defined by the
*  values in D, Z, and RHO, between KSTART and KSTOP.  It makes the
*  appropriate calls to DLAED4 and then stores the new matrix of
*  eigenvectors for use in calculating the next level of Z vectors.
*
*  Arguments
*  =========
*
*  K       (input) INTEGER
*          The number of terms in the rational function to be solved by
*          DLAED4.  K >= 0.
*
*  KSTART  (input) INTEGER
*  KSTOP   (input) INTEGER
*          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP
*          are to be computed.  1 <= KSTART <= KSTOP <= K.
*
*  N       (input) INTEGER
*          The number of rows and columns in the Q matrix.
*          N >= K (delation may result in N > K).
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          D(I) contains the updated eigenvalues
*          for KSTART <= I <= KSTOP.
*
*  Q       (workspace) DOUBLE PRECISION array, dimension (LDQ,N)
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= max( 1, N ).
*
*  RHO     (input) DOUBLE PRECISION
*          The value of the parameter in the rank one update equation.
*          RHO >= 0 required.
*
*  DLAMDA  (input) DOUBLE PRECISION array, dimension (K)
*          The first K elements of this array contain the old roots
*          of the deflated updating problem.  These are the poles
*          of the secular equation.
*
*  W       (input) DOUBLE PRECISION array, dimension (K)
*          The first K elements of this array contain the components
*          of the deflation-adjusted updating vector.
*
*  S       (output) DOUBLE PRECISION array, dimension (LDS, K)
*          Will contain the eigenvectors of the repaired matrix which
*          will be stored for subsequent Z vector calculation and
*          multiplied by the previously accumulated eigenvectors
*          to update the system.
*
*  LDS     (input) INTEGER
*          The leading dimension of S.  LDS >= max( 1, K ).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, an eigenvalue did not converge
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   TEMP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAED4, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( KSTART.LT.1 .OR. KSTART.GT.MAX( 1, K ) ) THEN
         INFO = -2
      ELSE IF( MAX( 1, KSTOP ).LT.KSTART .OR. KSTOP.GT.MAX( 1, K ) )
     $          THEN
         INFO = -3
      ELSE IF( N.LT.K ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDS.LT.MAX( 1, K ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED9', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( K.EQ.0 )
     $   RETURN
*
*     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
*     be computed with high relative accuracy (barring over/underflow).
*     This is a problem on machines without a guard digit in
*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
*     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
*     which on any of these machines zeros out the bottommost
*     bit of DLAMDA(I) if it is 1; this makes the subsequent
*     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
*     occurs. On binary machines with a guard digit (almost all
*     machines) it does not change DLAMDA(I) at all. On hexadecimal
*     and decimal machines with a guard digit, it slightly
*     changes the bottommost bits of DLAMDA(I). It does not account
*     for hexadecimal or decimal machines without guard digits
*     (we know of none). We use a subroutine call to compute
*     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
*     this code.
*
      DO 10 I = 1, N
         DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
   10 CONTINUE
*
      DO 20 J = KSTART, KSTOP
         CALL DLAED4( K, J, DLAMDA, W, Q( 1, J ), RHO, D( J ), INFO )
*
*        If the zero finder fails, the computation is terminated.
*
         IF( INFO.NE.0 )
     $      GO TO 120
   20 CONTINUE
*
      IF( K.EQ.1 .OR. K.EQ.2 ) THEN
         DO 40 I = 1, K
            DO 30 J = 1, K
               S( J, I ) = Q( J, I )
   30       CONTINUE
   40    CONTINUE
         GO TO 120
      END IF
*
*     Compute updated W.
*
      CALL DCOPY( K, W, 1, S, 1 )
*
*     Initialize W(I) = Q(I,I)
*
      CALL DCOPY( K, Q, LDQ+1, W, 1 )
      DO 70 J = 1, K
         DO 50 I = 1, J - 1
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   50    CONTINUE
         DO 60 I = J + 1, K
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   60    CONTINUE
   70 CONTINUE
      DO 80 I = 1, K
         W( I ) = SIGN( SQRT( -W( I ) ), S( I, 1 ) )
   80 CONTINUE
*
*     Compute eigenvectors of the modified rank-1 modification.
*
      DO 110 J = 1, K
         DO 90 I = 1, K
            Q( I, J ) = W( I ) / Q( I, J )
   90    CONTINUE
         TEMP = DNRM2( K, Q( 1, J ), 1 )
         DO 100 I = 1, K
            S( I, J ) = Q( I, J ) / TEMP
  100    CONTINUE
  110 CONTINUE
*
  120 CONTINUE
      RETURN
*
*     End of DLAED9
*
      END
      SUBROUTINE DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,
     $                   GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            CURLVL, CURPBM, INFO, N, TLVLS
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), GIVPTR( * ), PERM( * ),
     $                   PRMPTR( * ), QPTR( * )
      DOUBLE PRECISION   GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAEDA computes the Z vector corresponding to the merge step in the
*  CURLVLth step of the merge process with TLVLS steps for the CURPBMth
*  problem.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  TLVLS  (input) INTEGER
*         The total number of merging levels in the overall divide and
*         conquer tree.
*
*  CURLVL (input) INTEGER
*         The current level in the overall merge routine,
*         0 <= curlvl <= tlvls.
*
*  CURPBM (input) INTEGER
*         The current problem in the current level in the overall
*         merge routine (counting from upper left to lower right).
*
*  PRMPTR (input) INTEGER array, dimension (N lg N)
*         Contains a list of pointers which indicate where in PERM a
*         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
*         indicates the size of the permutation and incidentally the
*         size of the full, non-deflated problem.
*
*  PERM   (input) INTEGER array, dimension (N lg N)
*         Contains the permutations (from deflation and sorting) to be
*         applied to each eigenblock.
*
*  GIVPTR (input) INTEGER array, dimension (N lg N)
*         Contains a list of pointers which indicate where in GIVCOL a
*         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
*         indicates the number of Givens rotations.
*
*  GIVCOL (input) INTEGER array, dimension (2, N lg N)
*         Each pair of numbers indicates a pair of columns to take place
*         in a Givens rotation.
*
*  GIVNUM (input) DOUBLE PRECISION array, dimension (2, N lg N)
*         Each number indicates the S value to be used in the
*         corresponding Givens rotation.
*
*  Q      (input) DOUBLE PRECISION array, dimension (N**2)
*         Contains the square eigenblocks from previous levels, the
*         starting positions for blocks are given by QPTR.
*
*  QPTR   (input) INTEGER array, dimension (N+2)
*         Contains a list of pointers which indicate where in Q an
*         eigenblock is stored.  SQRT( QPTR(i+1) - QPTR(i) ) indicates
*         the size of the block.
*
*  Z      (output) DOUBLE PRECISION array, dimension (N)
*         On output this vector contains the updating vector (the last
*         row of the first sub-eigenvector matrix and the first row of
*         the second sub-eigenvector matrix).
*
*  ZTEMP  (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            BSIZ1, BSIZ2, CURR, I, K, MID, PSIZ1, PSIZ2,
     $                   PTR, ZPTR1
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV, DROT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -1
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAEDA', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine location of first number in second half.
*
      MID = N / 2 + 1
*
*     Gather last/first rows of appropriate eigenblocks into center of Z
*
      PTR = 1
*
*     Determine location of lowest level subproblem in the full storage
*     scheme
*
      CURR = PTR + CURPBM*2**CURLVL + 2**( CURLVL-1 ) - 1
*
*     Determine size of these matrices.  We add HALF to the value of
*     the SQRT in case the machine underestimates one of these square
*     roots.
*
      BSIZ1 = INT( HALF+SQRT( DBLE( QPTR( CURR+1 )-QPTR( CURR ) ) ) )
      BSIZ2 = INT( HALF+SQRT( DBLE( QPTR( CURR+2 )-QPTR( CURR+1 ) ) ) )
      DO 10 K = 1, MID - BSIZ1 - 1
         Z( K ) = ZERO
   10 CONTINUE
      CALL DCOPY( BSIZ1, Q( QPTR( CURR )+BSIZ1-1 ), BSIZ1,
     $            Z( MID-BSIZ1 ), 1 )
      CALL DCOPY( BSIZ2, Q( QPTR( CURR+1 ) ), BSIZ2, Z( MID ), 1 )
      DO 20 K = MID + BSIZ2, N
         Z( K ) = ZERO
   20 CONTINUE
*
*     Loop thru remaining levels 1 -> CURLVL applying the Givens
*     rotations and permutation and then multiplying the center matrices
*     against the current Z.
*
      PTR = 2**TLVLS + 1
      DO 70 K = 1, CURLVL - 1
         CURR = PTR + CURPBM*2**( CURLVL-K ) + 2**( CURLVL-K-1 ) - 1
         PSIZ1 = PRMPTR( CURR+1 ) - PRMPTR( CURR )
         PSIZ2 = PRMPTR( CURR+2 ) - PRMPTR( CURR+1 )
         ZPTR1 = MID - PSIZ1
*
*       Apply Givens at CURR and CURR+1
*
         DO 30 I = GIVPTR( CURR ), GIVPTR( CURR+1 ) - 1
            CALL DROT( 1, Z( ZPTR1+GIVCOL( 1, I )-1 ), 1,
     $                 Z( ZPTR1+GIVCOL( 2, I )-1 ), 1, GIVNUM( 1, I ),
     $                 GIVNUM( 2, I ) )
   30    CONTINUE
         DO 40 I = GIVPTR( CURR+1 ), GIVPTR( CURR+2 ) - 1
            CALL DROT( 1, Z( MID-1+GIVCOL( 1, I ) ), 1,
     $                 Z( MID-1+GIVCOL( 2, I ) ), 1, GIVNUM( 1, I ),
     $                 GIVNUM( 2, I ) )
   40    CONTINUE
         PSIZ1 = PRMPTR( CURR+1 ) - PRMPTR( CURR )
         PSIZ2 = PRMPTR( CURR+2 ) - PRMPTR( CURR+1 )
         DO 50 I = 0, PSIZ1 - 1
            ZTEMP( I+1 ) = Z( ZPTR1+PERM( PRMPTR( CURR )+I )-1 )
   50    CONTINUE
         DO 60 I = 0, PSIZ2 - 1
            ZTEMP( PSIZ1+I+1 ) = Z( MID+PERM( PRMPTR( CURR+1 )+I )-1 )
   60    CONTINUE
*
*        Multiply Blocks at CURR and CURR+1
*
*        Determine size of these matrices.  We add HALF to the value of
*        the SQRT in case the machine underestimates one of these
*        square roots.
*
         BSIZ1 = INT( HALF+SQRT( DBLE( QPTR( CURR+1 )-QPTR( CURR ) ) ) )
         BSIZ2 = INT( HALF+SQRT( DBLE( QPTR( CURR+2 )-QPTR( CURR+
     $           1 ) ) ) )
         IF( BSIZ1.GT.0 ) THEN
            CALL DGEMV( 'T', BSIZ1, BSIZ1, ONE, Q( QPTR( CURR ) ),
     $                  BSIZ1, ZTEMP( 1 ), 1, ZERO, Z( ZPTR1 ), 1 )
         END IF
         CALL DCOPY( PSIZ1-BSIZ1, ZTEMP( BSIZ1+1 ), 1, Z( ZPTR1+BSIZ1 ),
     $               1 )
         IF( BSIZ2.GT.0 ) THEN
            CALL DGEMV( 'T', BSIZ2, BSIZ2, ONE, Q( QPTR( CURR+1 ) ),
     $                  BSIZ2, ZTEMP( PSIZ1+1 ), 1, ZERO, Z( MID ), 1 )
         END IF
         CALL DCOPY( PSIZ2-BSIZ2, ZTEMP( PSIZ1+BSIZ2+1 ), 1,
     $               Z( MID+BSIZ2 ), 1 )
*
         PTR = PTR + 2**( TLVLS-K )
   70 CONTINUE
*
      RETURN
*
*     End of DLAEDA
*
      END
      SUBROUTINE dgeev( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
      implicit none
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   wi( * ), work( * ), wr( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
      CHARACTER          SIDE
      INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, ITAU, IWRK, K,
     $                   lwork_trevc, maxwrk, minwrk, nout
      DOUBLE PRECISION   ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM,
     $                   sn
*     ..
*     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgebak, dgebal, dgehrd, dhseqr, dlabad, dlacpy,
     $                   dlartg, dlascl, dorghr, drot, dscal, dtrevc3,
     $                   xerbla
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE, DLAPY2, DNRM2
      EXTERNAL           lsame, idamax, ilaenv, dlamch, dlange, dlapy2,
     $                   dnrm2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, sqrt
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      lquery = ( lwork.EQ.-1 )
      wantvl = lsame( jobvl, 'V' )
      wantvr = lsame( jobvr, 'V' )
      IF( ( .NOT.wantvl ) .AND. ( .NOT.lsame( jobvl, 'N' ) ) ) THEN
         info = -1
      ELSE IF( ( .NOT.wantvr ) .AND. ( .NOT.lsame( jobvr, 'N' ) ) ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldvl.LT.1 .OR. ( wantvl .AND. ldvl.LT.n ) ) THEN
         info = -9
      ELSE IF( ldvr.LT.1 .OR. ( wantvr .AND. ldvr.LT.n ) ) THEN
         info = -11
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by DHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.)
*
      IF( info.EQ.0 ) THEN
         IF( n.EQ.0 ) THEN
            minwrk = 1
            maxwrk = 1
         ELSE
            maxwrk = 2*n + n*ilaenv( 1, 'DGEHRD', ' ', n, 1, n, 0 )
            IF( wantvl ) THEN
               minwrk = 4*n
               maxwrk = max( maxwrk, 2*n + ( n - 1 )*ilaenv( 1,
     $                       'DORGHR', ' ', n, 1, n, -1 ) )
               CALL dhseqr( 'S', 'V', n, 1, n, a, lda, wr, wi, vl, ldvl,
     $                      work, -1, info )
               hswork = int( work(1) )
               maxwrk = max( maxwrk, n + 1, n + hswork )
               CALL dtrevc3( 'L', 'B', SELECT, n, a, lda,
     $                       vl, ldvl, vr, ldvr, n, nout,
     $                       work, -1, ierr )
               lwork_trevc = int( work(1) )
               maxwrk = max( maxwrk, n + lwork_trevc )
               maxwrk = max( maxwrk, 4*n )
            ELSE IF( wantvr ) THEN
               minwrk = 4*n
               maxwrk = max( maxwrk, 2*n + ( n - 1 )*ilaenv( 1,
     $                       'DORGHR', ' ', n, 1, n, -1 ) )
               CALL dhseqr( 'S', 'V', n, 1, n, a, lda, wr, wi, vr, ldvr,
     $                      work, -1, info )
               hswork = int( work(1) )
               maxwrk = max( maxwrk, n + 1, n + hswork )
               CALL dtrevc3( 'R', 'B', SELECT, n, a, lda,
     $                       vl, ldvl, vr, ldvr, n, nout,
     $                       work, -1, ierr )
               lwork_trevc = int( work(1) )
               maxwrk = max( maxwrk, n + lwork_trevc )
               maxwrk = max( maxwrk, 4*n )
            ELSE
               minwrk = 3*n
               CALL dhseqr( 'E', 'N', n, 1, n, a, lda, wr, wi, vr, ldvr,
     $                      work, -1, info )
               hswork = int( work(1) )
               maxwrk = max( maxwrk, n + 1, n + hswork )
            END IF
            maxwrk = max( maxwrk, minwrk )
         END IF
         work( 1 ) = maxwrk
*
         IF( lwork.LT.minwrk .AND. .NOT.lquery ) THEN
            info = -13
         END IF
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEEV ', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      eps = dlamch( 'P' )
      smlnum = dlamch( 'S' )
      bignum = one / smlnum
      CALL dlabad( smlnum, bignum )
      smlnum = sqrt( smlnum ) / eps
      bignum = one / smlnum
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      anrm = dlange( 'M', n, n, a, lda, dum )
      scalea = .false.
      IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
         scalea = .true.
         cscale = smlnum
      ELSE IF( anrm.GT.bignum ) THEN
         scalea = .true.
         cscale = bignum
      END IF
      IF( scalea )
     $   CALL dlascl( 'G', 0, 0, anrm, cscale, n, n, a, lda, ierr )
*
*     Balance the matrix
*     (Workspace: need N)
*
      ibal = 1
      CALL dgebal( 'B', n, a, lda, ilo, ihi, work( ibal ), ierr )
*
*     Reduce to upper Hessenberg form
*     (Workspace: need 3*N, prefer 2*N+N*NB)
*
      itau = ibal + n
      iwrk = itau + n
      CALL dgehrd( n, ilo, ihi, a, lda, work( itau ), work( iwrk ),
     $             lwork-iwrk+1, ierr )
*
      IF( wantvl ) THEN
*
*        Want left eigenvectors
*        Copy Householder vectors to VL
*
         side = 'L'
         CALL dlacpy( 'L', n, n, a, lda, vl, ldvl )
*
*        Generate orthogonal matrix in VL
*        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*
         CALL dorghr( n, ilo, ihi, vl, ldvl, work( itau ), work( iwrk ),
     $                lwork-iwrk+1, ierr )
*
*        Perform QR iteration, accumulating Schur vectors in VL
*        (Workspace: need N+1, prefer N+HSWORK (see comments) )
*
         iwrk = itau
         CALL dhseqr( 'S', 'V', n, ilo, ihi, a, lda, wr, wi, vl, ldvl,
     $                work( iwrk ), lwork-iwrk+1, info )
*
         IF( wantvr ) THEN
*
*           Want left and right eigenvectors
*           Copy Schur vectors to VR
*
            side = 'B'
            CALL dlacpy( 'F', n, n, vl, ldvl, vr, ldvr )
         END IF
*
      ELSE IF( wantvr ) THEN
*
*        Want right eigenvectors
*        Copy Householder vectors to VR
*
         side = 'R'
         CALL dlacpy( 'L', n, n, a, lda, vr, ldvr )
*
*        Generate orthogonal matrix in VR
*        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*
         CALL dorghr( n, ilo, ihi, vr, ldvr, work( itau ), work( iwrk ),
     $                lwork-iwrk+1, ierr )
*
*        Perform QR iteration, accumulating Schur vectors in VR
*        (Workspace: need N+1, prefer N+HSWORK (see comments) )
*
         iwrk = itau
         CALL dhseqr( 'S', 'V', n, ilo, ihi, a, lda, wr, wi, vr, ldvr,
     $                work( iwrk ), lwork-iwrk+1, info )
*
      ELSE
*
*        Compute eigenvalues only
*        (Workspace: need N+1, prefer N+HSWORK (see comments) )
*
         iwrk = itau
         CALL dhseqr( 'E', 'N', n, ilo, ihi, a, lda, wr, wi, vr, ldvr,
     $                work( iwrk ), lwork-iwrk+1, info )
      END IF
*
*     If INFO .NE. 0 from DHSEQR, then quit
*
      IF( info.NE.0 )
     $   GO TO 50
*
      IF( wantvl .OR. wantvr ) THEN
*
*        Compute left and/or right eigenvectors
*        (Workspace: need 4*N, prefer N + N + 2*N*NB)
*
         CALL dtrevc3( side, 'B', SELECT, n, a, lda, vl, ldvl, vr, ldvr,
     $                 n, nout, work( iwrk ), lwork-iwrk+1, ierr )
      END IF
*
      IF( wantvl ) THEN
*
*        Undo balancing of left eigenvectors
*        (Workspace: need N)
*
         CALL dgebak( 'B', 'L', n, ilo, ihi, work( ibal ), n, vl, ldvl,
     $                ierr )
*
*        Normalize left eigenvectors and make largest component real
*
         DO 20 i = 1, n
            IF( wi( i ).EQ.zero ) THEN
               scl = one / dnrm2( n, vl( 1, i ), 1 )
               CALL dscal( n, scl, vl( 1, i ), 1 )
            ELSE IF( wi( i ).GT.zero ) THEN
               scl = one / dlapy2( dnrm2( n, vl( 1, i ), 1 ),
     $               dnrm2( n, vl( 1, i+1 ), 1 ) )
               CALL dscal( n, scl, vl( 1, i ), 1 )
               CALL dscal( n, scl, vl( 1, i+1 ), 1 )
               DO 10 k = 1, n
                  work( iwrk+k-1 ) = vl( k, i )**2 + vl( k, i+1 )**2
   10          CONTINUE
               k = idamax( n, work( iwrk ), 1 )
               CALL dlartg( vl( k, i ), vl( k, i+1 ), cs, sn, r )
               CALL drot( n, vl( 1, i ), 1, vl( 1, i+1 ), 1, cs, sn )
               vl( k, i+1 ) = zero
            END IF
   20    CONTINUE
      END IF
*
      IF( wantvr ) THEN
*
*        Undo balancing of right eigenvectors
*        (Workspace: need N)
*
         CALL dgebak( 'B', 'R', n, ilo, ihi, work( ibal ), n, vr, ldvr,
     $                ierr )
*
*        Normalize right eigenvectors and make largest component real
*
         DO 40 i = 1, n
            IF( wi( i ).EQ.zero ) THEN
               scl = one / dnrm2( n, vr( 1, i ), 1 )
               CALL dscal( n, scl, vr( 1, i ), 1 )
            ELSE IF( wi( i ).GT.zero ) THEN
               scl = one / dlapy2( dnrm2( n, vr( 1, i ), 1 ),
     $               dnrm2( n, vr( 1, i+1 ), 1 ) )
               CALL dscal( n, scl, vr( 1, i ), 1 )
               CALL dscal( n, scl, vr( 1, i+1 ), 1 )
               DO 30 k = 1, n
                  work( iwrk+k-1 ) = vr( k, i )**2 + vr( k, i+1 )**2
   30          CONTINUE
               k = idamax( n, work( iwrk ), 1 )
               CALL dlartg( vr( k, i ), vr( k, i+1 ), cs, sn, r )
               CALL drot( n, vr( 1, i ), 1, vr( 1, i+1 ), 1, cs, sn )
               vr( k, i+1 ) = zero
            END IF
   40    CONTINUE
      END IF
*
*     Undo scaling if necessary
*
   50 CONTINUE
      IF( scalea ) THEN
         CALL dlascl( 'G', 0, 0, cscale, anrm, n-info, 1, wr( info+1 ),
     $                max( n-info, 1 ), ierr )
         CALL dlascl( 'G', 0, 0, cscale, anrm, n-info, 1, wi( info+1 ),
     $                max( n-info, 1 ), ierr )
         IF( info.GT.0 ) THEN
            CALL dlascl( 'G', 0, 0, cscale, anrm, ilo-1, 1, wr, n,
     $                   ierr )
            CALL dlascl( 'G', 0, 0, cscale, anrm, ilo-1, 1, wi, n,
     $                   ierr )
         END IF
      END IF
*
      work( 1 ) = maxwrk
      RETURN
*
*     End of DGEEV
*
      END
      SUBROUTINE zgeev( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
     $                  WORK, LWORK, RWORK, INFO )
      implicit none
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   w( * ), work( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
      CHARACTER          SIDE
      INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, IRWORK, ITAU,
     $                   iwrk, k, lwork_trevc, maxwrk, minwrk, nout
      DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, SCL, SMLNUM
      COMPLEX*16         TMP
*     ..
*     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, xerbla, zdscal, zgebak, zgebal, zgehrd,
     $                   zhseqr, zlacpy, zlascl, zscal, ztrevc3, zunghr
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DZNRM2, ZLANGE
      EXTERNAL           lsame, idamax, ilaenv, dlamch, dznrm2, zlange
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dble, dcmplx, conjg, aimag, max, sqrt
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      lquery = ( lwork.EQ.-1 )
      wantvl = lsame( jobvl, 'V' )
      wantvr = lsame( jobvr, 'V' )
      IF( ( .NOT.wantvl ) .AND. ( .NOT.lsame( jobvl, 'N' ) ) ) THEN
         info = -1
      ELSE IF( ( .NOT.wantvr ) .AND. ( .NOT.lsame( jobvr, 'N' ) ) ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldvl.LT.1 .OR. ( wantvl .AND. ldvl.LT.n ) ) THEN
         info = -8
      ELSE IF( ldvr.LT.1 .OR. ( wantvr .AND. ldvr.LT.n ) ) THEN
         info = -10
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       CWorkspace refers to complex workspace, and RWorkspace to real
*       workspace. NB refers to the optimal block size for the
*       immediately following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by ZHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.)
*
      IF( info.EQ.0 ) THEN
         IF( n.EQ.0 ) THEN
            minwrk = 1
            maxwrk = 1
         ELSE
            maxwrk = n + n*ilaenv( 1, 'ZGEHRD', ' ', n, 1, n, 0 )
            minwrk = 2*n
            IF( wantvl ) THEN
               maxwrk = max( maxwrk, n + ( n - 1 )*ilaenv( 1, 'ZUNGHR',
     $                       ' ', n, 1, n, -1 ) )
               CALL ztrevc3( 'L', 'B', SELECT, n, a, lda,
     $                       vl, ldvl, vr, ldvr,
     $                       n, nout, work, -1, rwork, -1, ierr )
               lwork_trevc = int( work(1) )
               maxwrk = max( maxwrk, n + lwork_trevc )
               CALL zhseqr( 'S', 'V', n, 1, n, a, lda, w, vl, ldvl,
     $                      work, -1, info )
            ELSE IF( wantvr ) THEN
               maxwrk = max( maxwrk, n + ( n - 1 )*ilaenv( 1, 'ZUNGHR',
     $                       ' ', n, 1, n, -1 ) )
               CALL ztrevc3( 'R', 'B', SELECT, n, a, lda,
     $                       vl, ldvl, vr, ldvr,
     $                       n, nout, work, -1, rwork, -1, ierr )
               lwork_trevc = int( work(1) )
               maxwrk = max( maxwrk, n + lwork_trevc )
               CALL zhseqr( 'S', 'V', n, 1, n, a, lda, w, vr, ldvr,
     $                      work, -1, info )
            ELSE
               CALL zhseqr( 'E', 'N', n, 1, n, a, lda, w, vr, ldvr,
     $                      work, -1, info )
            END IF
            hswork = int( work(1) )
            maxwrk = max( maxwrk, hswork, minwrk )
         END IF
         work( 1 ) = maxwrk
*
         IF( lwork.LT.minwrk .AND. .NOT.lquery ) THEN
            info = -12
         END IF
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEEV ', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      eps = dlamch( 'P' )
      smlnum = dlamch( 'S' )
      bignum = one / smlnum
      CALL dlabad( smlnum, bignum )
      smlnum = sqrt( smlnum ) / eps
      bignum = one / smlnum
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      anrm = zlange( 'M', n, n, a, lda, dum )
      scalea = .false.
      IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
         scalea = .true.
         cscale = smlnum
      ELSE IF( anrm.GT.bignum ) THEN
         scalea = .true.
         cscale = bignum
      END IF
      IF( scalea )
     $   CALL zlascl( 'G', 0, 0, anrm, cscale, n, n, a, lda, ierr )
*
*     Balance the matrix
*     (CWorkspace: none)
*     (RWorkspace: need N)
*
      ibal = 1
      CALL zgebal( 'B', n, a, lda, ilo, ihi, rwork( ibal ), ierr )
*
*     Reduce to upper Hessenberg form
*     (CWorkspace: need 2*N, prefer N+N*NB)
*     (RWorkspace: none)
*
      itau = 1
      iwrk = itau + n
      CALL zgehrd( n, ilo, ihi, a, lda, work( itau ), work( iwrk ),
     $             lwork-iwrk+1, ierr )
*
      IF( wantvl ) THEN
*
*        Want left eigenvectors
*        Copy Householder vectors to VL
*
         side = 'L'
         CALL zlacpy( 'L', n, n, a, lda, vl, ldvl )
*
*        Generate unitary matrix in VL
*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
*        (RWorkspace: none)
*
         CALL zunghr( n, ilo, ihi, vl, ldvl, work( itau ), work( iwrk ),
     $                lwork-iwrk+1, ierr )
*
*        Perform QR iteration, accumulating Schur vectors in VL
*        (CWorkspace: need 1, prefer HSWORK (see comments) )
*        (RWorkspace: none)
*
         iwrk = itau
         CALL zhseqr( 'S', 'V', n, ilo, ihi, a, lda, w, vl, ldvl,
     $                work( iwrk ), lwork-iwrk+1, info )
*
         IF( wantvr ) THEN
*
*           Want left and right eigenvectors
*           Copy Schur vectors to VR
*
            side = 'B'
            CALL zlacpy( 'F', n, n, vl, ldvl, vr, ldvr )
         END IF
*
      ELSE IF( wantvr ) THEN
*
*        Want right eigenvectors
*        Copy Householder vectors to VR
*
         side = 'R'
         CALL zlacpy( 'L', n, n, a, lda, vr, ldvr )
*
*        Generate unitary matrix in VR
*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
*        (RWorkspace: none)
*
         CALL zunghr( n, ilo, ihi, vr, ldvr, work( itau ), work( iwrk ),
     $                lwork-iwrk+1, ierr )
*
*        Perform QR iteration, accumulating Schur vectors in VR
*        (CWorkspace: need 1, prefer HSWORK (see comments) )
*        (RWorkspace: none)
*
         iwrk = itau
         CALL zhseqr( 'S', 'V', n, ilo, ihi, a, lda, w, vr, ldvr,
     $                work( iwrk ), lwork-iwrk+1, info )
*
      ELSE
*
*        Compute eigenvalues only
*        (CWorkspace: need 1, prefer HSWORK (see comments) )
*        (RWorkspace: none)
*
         iwrk = itau
         CALL zhseqr( 'E', 'N', n, ilo, ihi, a, lda, w, vr, ldvr,
     $                work( iwrk ), lwork-iwrk+1, info )
      END IF
*
*     If INFO .NE. 0 from ZHSEQR, then quit
*
      IF( info.NE.0 )
     $   GO TO 50
*
      IF( wantvl .OR. wantvr ) THEN
*
*        Compute left and/or right eigenvectors
*        (CWorkspace: need 2*N, prefer N + 2*N*NB)
*        (RWorkspace: need 2*N)
*
         irwork = ibal + n
         CALL ztrevc3( side, 'B', SELECT, n, a, lda, vl, ldvl, vr, ldvr,
     $                 n, nout, work( iwrk ), lwork-iwrk+1,
     $                 rwork( irwork ), n, ierr )
      END IF
*
      IF( wantvl ) THEN
*
*        Undo balancing of left eigenvectors
*        (CWorkspace: none)
*        (RWorkspace: need N)
*
         CALL zgebak( 'B', 'L', n, ilo, ihi, rwork( ibal ), n, vl, ldvl,
     $                ierr )
*
*        Normalize left eigenvectors and make largest component real
*
         DO 20 i = 1, n
            scl = one / dznrm2( n, vl( 1, i ), 1 )
            CALL zdscal( n, scl, vl( 1, i ), 1 )
            DO 10 k = 1, n
               rwork( irwork+k-1 ) = dble( vl( k, i ) )**2 +
     $                               aimag( vl( k, i ) )**2
   10       CONTINUE
            k = idamax( n, rwork( irwork ), 1 )
            tmp = conjg( vl( k, i ) ) / sqrt( rwork( irwork+k-1 ) )
            CALL zscal( n, tmp, vl( 1, i ), 1 )
            vl( k, i ) = dcmplx( dble( vl( k, i ) ), zero )
   20    CONTINUE
      END IF
*
      IF( wantvr ) THEN
*
*        Undo balancing of right eigenvectors
*        (CWorkspace: none)
*        (RWorkspace: need N)
*
         CALL zgebak( 'B', 'R', n, ilo, ihi, rwork( ibal ), n, vr, ldvr,
     $                ierr )
*
*        Normalize right eigenvectors and make largest component real
*
         DO 40 i = 1, n
            scl = one / dznrm2( n, vr( 1, i ), 1 )
            CALL zdscal( n, scl, vr( 1, i ), 1 )
            DO 30 k = 1, n
               rwork( irwork+k-1 ) = dble( vr( k, i ) )**2 +
     $                               aimag( vr( k, i ) )**2
   30       CONTINUE
            k = idamax( n, rwork( irwork ), 1 )
            tmp = conjg( vr( k, i ) ) / sqrt( rwork( irwork+k-1 ) )
            CALL zscal( n, tmp, vr( 1, i ), 1 )
            vr( k, i ) = dcmplx( dble( vr( k, i ) ), zero )
   40    CONTINUE
      END IF
*
*     Undo scaling if necessary
*
   50 CONTINUE
      IF( scalea ) THEN
         CALL zlascl( 'G', 0, 0, cscale, anrm, n-info, 1, w( info+1 ),
     $                max( n-info, 1 ), ierr )
         IF( info.GT.0 ) THEN
            CALL zlascl( 'G', 0, 0, cscale, anrm, ilo-1, 1, w, n, ierr )
         END IF
      END IF
*
      work( 1 ) = maxwrk
      RETURN
*
*     End of ZGEEV
*
      END
      SUBROUTINE ztrevc3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                    LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO)
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, LWORK, LRWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   work( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      COMPLEX*16         CZERO, CONE
      parameter( czero = ( 0.0d+0, 0.0d+0 ),
     $                     cone  = ( 1.0d+0, 0.0d+0 ) )
      INTEGER            NBMIN, NBMAX
      parameter( nbmin = 8, nbmax = 128 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, LQUERY, OVER, RIGHTV, SOMEV
      INTEGER            I, II, IS, J, K, KI, IV, MAXWRK, NB
      DOUBLE PRECISION   OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
      COMPLEX*16         CDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV, IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      EXTERNAL           lsame, ilaenv, izamax, dlamch, dzasum
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zcopy, zdscal, zgemv, zlatrs,
     $                   zgemm, dlabad, zlaset, zlacpy
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, conjg, dimag, max
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      bothv  = lsame( side, 'B' )
      rightv = lsame( side, 'R' ) .OR. bothv
      leftv  = lsame( side, 'L' ) .OR. bothv
*
      allv  = lsame( howmny, 'A' )
      over  = lsame( howmny, 'B' )
      somev = lsame( howmny, 'S' )
*
*     Set M to the number of columns required to store the selected
*     eigenvectors.
*
      IF( somev ) THEN
         m = 0
         DO 10 j = 1, n
            IF( SELECT( j ) )
     $         m = m + 1
   10    CONTINUE
      ELSE
         m = n
      END IF
*
      info = 0
      nb = ilaenv( 1, 'ZTREVC', side // howmny, n, -1, -1, -1 )
      maxwrk = n + 2*n*nb
      work(1) = maxwrk
      rwork(1) = n
      lquery = ( lwork.EQ.-1 .OR. lrwork.EQ.-1 )
      IF( .NOT.rightv .AND. .NOT.leftv ) THEN
         info = -1
      ELSE IF( .NOT.allv .AND. .NOT.over .AND. .NOT.somev ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( ldt.LT.max( 1, n ) ) THEN
         info = -6
      ELSE IF( ldvl.LT.1 .OR. ( leftv .AND. ldvl.LT.n ) ) THEN
         info = -8
      ELSE IF( ldvr.LT.1 .OR. ( rightv .AND. ldvr.LT.n ) ) THEN
         info = -10
      ELSE IF( mm.LT.m ) THEN
         info = -11
      ELSE IF( lwork.LT.max( 1, 2*n ) .AND. .NOT.lquery ) THEN
         info = -14
      ELSE IF ( lrwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -16
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZTREVC3', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Use blocked version of back-transformation if sufficient workspace.
*     Zero-out the workspace to avoid potential NaN propagation.
*
      IF( over .AND. lwork .GE. n + 2*n*nbmin ) THEN
         nb = (lwork - n) / (2*n)
         nb = min( nb, nbmax )
         CALL zlaset( 'F', n, 1+2*nb, czero, czero, work, n )
      ELSE
         nb = 1
      END IF
*
*     Set the constants to control overflow.
*
      unfl = dlamch( 'Safe minimum' )
      ovfl = one / unfl
      CALL dlabad( unfl, ovfl )
      ulp = dlamch( 'Precision' )
      smlnum = unfl*( n / ulp )
*
*     Store the diagonal elements of T in working array WORK.
*
      DO 20 i = 1, n
         work( i ) = t( i, i )
   20 CONTINUE
*
*     Compute 1-norm of each column of strictly upper triangular
*     part of T to control overflow in triangular solver.
*
      rwork( 1 ) = zero
      DO 30 j = 2, n
         rwork( j ) = dzasum( j-1, t( 1, j ), 1 )
   30 CONTINUE
*
      IF( rightv ) THEN
*
*        ============================================================
*        Compute right eigenvectors.
*
*        IV is index of column in current block.
*        Non-blocked version always uses IV=NB=1;
*        blocked     version starts with IV=NB, goes down to 1.
*        (Note the "0-th" column is used to store the original diagonal.)
         iv = nb
         is = m
         DO 80 ki = n, 1, -1
            IF( somev ) THEN
               IF( .NOT.SELECT( ki ) )
     $            GO TO 80
            END IF
            smin = max( ulp*( cabs1( t( ki, ki ) ) ), smlnum )
*
*           --------------------------------------------------------
*           Complex right eigenvector
*
            work( ki + iv*n ) = cone
*
*           Form right-hand side.
*
            DO 40 k = 1, ki - 1
               work( k + iv*n ) = -t( k, ki )
   40       CONTINUE
*
*           Solve upper triangular system:
*           [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE*WORK.
*
            DO 50 k = 1, ki - 1
               t( k, k ) = t( k, k ) - t( ki, ki )
               IF( cabs1( t( k, k ) ).LT.smin )
     $            t( k, k ) = smin
   50       CONTINUE
*
            IF( ki.GT.1 ) THEN
               CALL zlatrs( 'Upper', 'No transpose', 'Non-unit', 'Y',
     $                      ki-1, t, ldt, work( 1 + iv*n ), scale,
     $                      rwork, info )
               work( ki + iv*n ) = scale
            END IF
*
*           Copy the vector x or Q*x to VR and normalize.
*
            IF( .NOT.over ) THEN
*              ------------------------------
*              no back-transform: copy x to VR and normalize.
               CALL zcopy( ki, work( 1 + iv*n ), 1, vr( 1, is ), 1 )
*
               ii = izamax( ki, vr( 1, is ), 1 )
               remax = one / cabs1( vr( ii, is ) )
               CALL zdscal( ki, remax, vr( 1, is ), 1 )
*
               DO 60 k = ki + 1, n
                  vr( k, is ) = czero
   60          CONTINUE
*
            ELSE IF( nb.EQ.1 ) THEN
*              ------------------------------
*              version 1: back-transform each vector with GEMV, Q*x.
               IF( ki.GT.1 )
     $            CALL zgemv( 'N', n, ki-1, cone, vr, ldvr,
     $                        work( 1 + iv*n ), 1, dcmplx( scale ),
     $                        vr( 1, ki ), 1 )
*
               ii = izamax( n, vr( 1, ki ), 1 )
               remax = one / cabs1( vr( ii, ki ) )
               CALL zdscal( n, remax, vr( 1, ki ), 1 )
*
            ELSE
*              ------------------------------
*              version 2: back-transform block of vectors with GEMM
*              zero out below vector
               DO k = ki + 1, n
                  work( k + iv*n ) = czero
               END DO
*
*              Columns IV:NB of work are valid vectors.
*              When the number of vectors stored reaches NB,
*              or if this was last vector, do the GEMM
               IF( (iv.EQ.1) .OR. (ki.EQ.1) ) THEN
                  CALL zgemm( 'N', 'N', n, nb-iv+1, ki+nb-iv, cone,
     $                        vr, ldvr,
     $                        work( 1 + (iv)*n    ), n,
     $                        czero,
     $                        work( 1 + (nb+iv)*n ), n )
*                 normalize vectors
                  DO k = iv, nb
                     ii = izamax( n, work( 1 + (nb+k)*n ), 1 )
                     remax = one / cabs1( work( ii + (nb+k)*n ) )
                     CALL zdscal( n, remax, work( 1 + (nb+k)*n ), 1 )
                  END DO
                  CALL zlacpy( 'F', n, nb-iv+1,
     $                         work( 1 + (nb+iv)*n ), n,
     $                         vr( 1, ki ), ldvr )
                  iv = nb
               ELSE
                  iv = iv - 1
               END IF
            END IF
*
*           Restore the original diagonal elements of T.
*
            DO 70 k = 1, ki - 1
               t( k, k ) = work( k )
   70       CONTINUE
*
            is = is - 1
   80    CONTINUE
      END IF
*
      IF( leftv ) THEN
*
*        ============================================================
*        Compute left eigenvectors.
*
*        IV is index of column in current block.
*        Non-blocked version always uses IV=1;
*        blocked     version starts with IV=1, goes up to NB.
*        (Note the "0-th" column is used to store the original diagonal.)
         iv = 1
         is = 1
         DO 130 ki = 1, n
*
            IF( somev ) THEN
               IF( .NOT.SELECT( ki ) )
     $            GO TO 130
            END IF
            smin = max( ulp*( cabs1( t( ki, ki ) ) ), smlnum )
*
*           --------------------------------------------------------
*           Complex left eigenvector
*
            work( ki + iv*n ) = cone
*
*           Form right-hand side.
*
            DO 90 k = ki + 1, n
               work( k + iv*n ) = -conjg( t( ki, k ) )
   90       CONTINUE
*
*           Solve conjugate-transposed triangular system:
*           [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE*WORK.
*
            DO 100 k = ki + 1, n
               t( k, k ) = t( k, k ) - t( ki, ki )
               IF( cabs1( t( k, k ) ).LT.smin )
     $            t( k, k ) = smin
  100       CONTINUE
*
            IF( ki.LT.n ) THEN
               CALL zlatrs( 'Upper', 'Conjugate transpose', 'Non-unit',
     $                      'Y', n-ki, t( ki+1, ki+1 ), ldt,
     $                      work( ki+1 + iv*n ), scale, rwork, info )
               work( ki + iv*n ) = scale
            END IF
*
*           Copy the vector x or Q*x to VL and normalize.
*
            IF( .NOT.over ) THEN
*              ------------------------------
*              no back-transform: copy x to VL and normalize.
               CALL zcopy( n-ki+1, work( ki + iv*n ), 1, vl(ki,is), 1 )
*
               ii = izamax( n-ki+1, vl( ki, is ), 1 ) + ki - 1
               remax = one / cabs1( vl( ii, is ) )
               CALL zdscal( n-ki+1, remax, vl( ki, is ), 1 )
*
               DO 110 k = 1, ki - 1
                  vl( k, is ) = czero
  110          CONTINUE
*
            ELSE IF( nb.EQ.1 ) THEN
*              ------------------------------
*              version 1: back-transform each vector with GEMV, Q*x.
               IF( ki.LT.n )
     $            CALL zgemv( 'N', n, n-ki, cone, vl( 1, ki+1 ), ldvl,
     $                        work( ki+1 + iv*n ), 1, dcmplx( scale ),
     $                        vl( 1, ki ), 1 )
*
               ii = izamax( n, vl( 1, ki ), 1 )
               remax = one / cabs1( vl( ii, ki ) )
               CALL zdscal( n, remax, vl( 1, ki ), 1 )
*
            ELSE
*              ------------------------------
*              version 2: back-transform block of vectors with GEMM
*              zero out above vector
*              could go from KI-NV+1 to KI-1
               DO k = 1, ki - 1
                  work( k + iv*n ) = czero
               END DO
*
*              Columns 1:IV of work are valid vectors.
*              When the number of vectors stored reaches NB,
*              or if this was last vector, do the GEMM
               IF( (iv.EQ.nb) .OR. (ki.EQ.n) ) THEN
                  CALL zgemm( 'N', 'N', n, iv, n-ki+iv, cone,
     $                        vl( 1, ki-iv+1 ), ldvl,
     $                        work( ki-iv+1 + (1)*n ), n,
     $                        czero,
     $                        work( 1 + (nb+1)*n ), n )
*                 normalize vectors
                  DO k = 1, iv
                     ii = izamax( n, work( 1 + (nb+k)*n ), 1 )
                     remax = one / cabs1( work( ii + (nb+k)*n ) )
                     CALL zdscal( n, remax, work( 1 + (nb+k)*n ), 1 )
                  END DO
                  CALL zlacpy( 'F', n, iv,
     $                         work( 1 + (nb+1)*n ), n,
     $                         vl( 1, ki-iv+1 ), ldvl )
                  iv = 1
               ELSE
                  iv = iv + 1
               END IF
            END IF
*
*           Restore the original diagonal elements of T.
*
            DO 120 k = ki + 1, n
               t( k, k ) = work( k )
  120       CONTINUE
*
            is = is + 1
  130    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTREVC3
*
      END
      SUBROUTINE dgehrd( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      parameter( nbmax = 64, ldt = nbmax+1,
     $                     tsize = ldt*nbmax )
      DOUBLE PRECISION  ZERO, ONE
      parameter( zero = 0.0d+0,
     $                     one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWT, J, LDWORK, LWKOPT, NB,
     $                   NBMIN, NH, NX
      DOUBLE PRECISION  EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           daxpy, dgehd2, dgemm, dlahr2, dlarfb, dtrmm,
     $                   xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      lquery = ( lwork.EQ.-1 )
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
*
      IF( info.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         nb = min( nbmax, ilaenv( 1, 'DGEHRD', ' ', n, ilo, ihi, -1 ) )
         lwkopt = n*nb + tsize
         work( 1 ) = lwkopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEHRD', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
*
      DO 10 i = 1, ilo - 1
         tau( i ) = zero
   10 CONTINUE
      DO 20 i = max( 1, ihi ), n - 1
         tau( i ) = zero
   20 CONTINUE
*
*     Quick return if possible
*
      nh = ihi - ilo + 1
      IF( nh.LE.1 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
*     Determine the block size
*
      nb = min( nbmax, ilaenv( 1, 'DGEHRD', ' ', n, ilo, ihi, -1 ) )
      nbmin = 2
      IF( nb.GT.1 .AND. nb.LT.nh ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code)
*
         nx = max( nb, ilaenv( 3, 'DGEHRD', ' ', n, ilo, ihi, -1 ) )
         IF( nx.LT.nh ) THEN
*
*           Determine if workspace is large enough for blocked code
*
            IF( lwork.LT.n*nb+tsize ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code
*
               nbmin = max( 2, ilaenv( 2, 'DGEHRD', ' ', n, ilo, ihi,
     $                 -1 ) )
               IF( lwork.GE.(n*nbmin + tsize) ) THEN
                  nb = (lwork-tsize) / n
               ELSE
                  nb = 1
               END IF
            END IF
         END IF
      END IF
      ldwork = n
*
      IF( nb.LT.nbmin .OR. nb.GE.nh ) THEN
*
*        Use unblocked code below
*
         i = ilo
*
      ELSE
*
*        Use blocked code
*
         iwt = 1 + n*nb
         DO 40 i = ilo, ihi - 1 - nx, nb
            ib = min( nb, ihi-i )
*
*           Reduce columns i:i+ib-1 to Hessenberg form, returning the
*           matrices V and T of the block reflector H = I - V*T*V**T
*           which performs the reduction, and also the matrix Y = A*V*T
*
            CALL dlahr2( ihi, i, ib, a( 1, i ), lda, tau( i ),
     $                   work( iwt ), ldt, work, ldwork )
*
*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
*           right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set
*           to 1
*
            ei = a( i+ib, i+ib-1 )
            a( i+ib, i+ib-1 ) = one
            CALL dgemm( 'No transpose', 'Transpose',
     $                  ihi, ihi-i-ib+1,
     $                  ib, -one, work, ldwork, a( i+ib, i ), lda, one,
     $                  a( 1, i+ib ), lda )
            a( i+ib, i+ib-1 ) = ei
*
*           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
*           right
*
            CALL dtrmm( 'Right', 'Lower', 'Transpose',
     $                  'Unit', i, ib-1,
     $                  one, a( i+1, i ), lda, work, ldwork )
            DO 30 j = 0, ib-2
               CALL daxpy( i, -one, work( ldwork*j+1 ), 1,
     $                     a( 1, i+j+1 ), 1 )
   30       CONTINUE
*
*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
*           left
*
            CALL dlarfb( 'Left', 'Transpose', 'Forward',
     $                   'Columnwise',
     $                   ihi-i, n-i-ib+1, ib, a( i+1, i ), lda,
     $                   work( iwt ), ldt, a( i+1, i+ib ), lda,
     $                   work, ldwork )
   40    CONTINUE
      END IF
*
*     Use unblocked code to reduce the rest of the matrix
*
      CALL dgehd2( n, i, ihi, a, lda, tau, work, iinfo )
      work( 1 ) = lwkopt
*
      RETURN
*
*     End of DGEHRD
*
      END
      SUBROUTINE dorghr( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LWKOPT, NB, NH
*     ..
*     .. External Subroutines ..
      EXTERNAL           dorgqr, xerbla
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      nh = ihi - ilo
      lquery = ( lwork.EQ.-1 )
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, nh ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
*
      IF( info.EQ.0 ) THEN
         nb = ilaenv( 1, 'DORGQR', ' ', nh, nh, nh, -1 )
         lwkopt = max( 1, nh )*nb
         work( 1 ) = lwkopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORGHR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
*     Shift the vectors which define the elementary reflectors one
*     column to the right, and set the first ilo and the last n-ihi
*     rows and columns to those of the unit matrix
*
      DO 40 j = ihi, ilo + 1, -1
         DO 10 i = 1, j - 1
            a( i, j ) = zero
   10    CONTINUE
         DO 20 i = j + 1, ihi
            a( i, j ) = a( i, j-1 )
   20    CONTINUE
         DO 30 i = ihi + 1, n
            a( i, j ) = zero
   30    CONTINUE
   40 CONTINUE
      DO 60 j = 1, ilo
         DO 50 i = 1, n
            a( i, j ) = zero
   50    CONTINUE
         a( j, j ) = one
   60 CONTINUE
      DO 80 j = ihi + 1, n
         DO 70 i = 1, n
            a( i, j ) = zero
   70    CONTINUE
         a( j, j ) = one
   80 CONTINUE
*
      IF( nh.GT.0 ) THEN
*
*        Generate Q(ilo+1:ihi,ilo+1:ihi)
*
         CALL dorgqr( nh, nh, nh, a( ilo+1, ilo+1 ), lda, tau( ilo ),
     $                work, lwork, iinfo )
      END IF
      work( 1 ) = lwkopt
      RETURN
*
*     End of DORGHR
*
      END
      SUBROUTINE dhseqr( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z,
     $                   LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
      CHARACTER          COMPZ, JOB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
     $                   z( ldz, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    DLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 15 )
*
*     ==== NL allocates some local workspace to help small matrices
*     .    through a rare DLAHQR failure.  NL > NTINY = 15 is
*     .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
*     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
*     .    allows up to six simultaneous shifts and a 16-by-16
*     .    deflation window.  ====
      INTEGER            NL
      parameter( nl = 49 )
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   HL( NL, NL ), WORKL( NL )
*     ..
*     .. Local Scalars ..
      INTEGER            I, KBOT, NMIN
      LOGICAL            INITZ, LQUERY, WANTT, WANTZ
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      LOGICAL            LSAME
      EXTERNAL           ilaenv, lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlacpy, dlahqr, dlaqr0, dlaset, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dble, max, min
*     ..
*     .. Executable Statements ..
*
*     ==== Decode and check the input parameters. ====
*
      wantt = lsame( job, 'S' )
      initz = lsame( compz, 'I' )
      wantz = initz .OR. lsame( compz, 'V' )
      work( 1 ) = dble( max( 1, n ) )
      lquery = lwork.EQ.-1
*
      info = 0
      IF( .NOT.lsame( job, 'E' ) .AND. .NOT.wantt ) THEN
         info = -1
      ELSE IF( .NOT.lsame( compz, 'N' ) .AND. .NOT.wantz ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -5
      ELSE IF( ldh.LT.max( 1, n ) ) THEN
         info = -7
      ELSE IF( ldz.LT.1 .OR. ( wantz .AND. ldz.LT.max( 1, n ) ) ) THEN
         info = -11
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -13
      END IF
*
      IF( info.NE.0 ) THEN
*
*        ==== Quick return in case of invalid argument. ====
*
         CALL xerbla( 'DHSEQR', -info )
         RETURN
*
      ELSE IF( n.EQ.0 ) THEN
*
*        ==== Quick return in case N = 0; nothing to do. ====
*
         RETURN
*
      ELSE IF( lquery ) THEN
*
*        ==== Quick return in case of a workspace query ====
*
         CALL dlaqr0( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, ilo,
     $                ihi, z, ldz, work, lwork, info )
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
         work( 1 ) = max( dble( max( 1, n ) ), work( 1 ) )
         RETURN
*
      ELSE
*
*        ==== copy eigenvalues isolated by DGEBAL ====
*
         DO 10 i = 1, ilo - 1
            wr( i ) = h( i, i )
            wi( i ) = zero
   10    CONTINUE
         DO 20 i = ihi + 1, n
            wr( i ) = h( i, i )
            wi( i ) = zero
   20    CONTINUE
*
*        ==== Initialize Z, if requested ====
*
         IF( initz )
     $      CALL dlaset( 'A', n, n, zero, one, z, ldz )
*
*        ==== Quick return if possible ====
*
         IF( ilo.EQ.ihi ) THEN
            wr( ilo ) = h( ilo, ilo )
            wi( ilo ) = zero
            RETURN
         END IF
*
*        ==== DLAHQR/DLAQR0 crossover point ====
*
         nmin = ilaenv( 12, 'DHSEQR', job( : 1 ) // compz( : 1 ), n,
     $          ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
*
*        ==== DLAQR0 for big matrices; DLAHQR for small ones ====
*
         IF( n.GT.nmin ) THEN
            CALL dlaqr0( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, ilo,
     $                   ihi, z, ldz, work, lwork, info )
         ELSE
*
*           ==== Small matrix ====
*
            CALL dlahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, ilo,
     $                   ihi, z, ldz, info )
*
            IF( info.GT.0 ) THEN
*
*              ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
*              .    when DLAHQR fails. ====
*
               kbot = info
*
               IF( n.GE.nl ) THEN
*
*                 ==== Larger matrices have enough subdiagonal scratch
*                 .    space to call DLAQR0 directly. ====
*
                  CALL dlaqr0( wantt, wantz, n, ilo, kbot, h, ldh, wr,
     $                         wi, ilo, ihi, z, ldz, work, lwork, info )
*
               ELSE
*
*                 ==== Tiny matrices don't have enough subdiagonal
*                 .    scratch space to benefit from DLAQR0.  Hence,
*                 .    tiny matrices must be copied into a larger
*                 .    array before calling DLAQR0. ====
*
                  CALL dlacpy( 'A', n, n, h, ldh, hl, nl )
                  hl( n+1, n ) = zero
                  CALL dlaset( 'A', nl, nl-n, zero, zero, hl( 1, n+1 ),
     $                         nl )
                  CALL dlaqr0( wantt, wantz, nl, ilo, kbot, hl, nl, wr,
     $                         wi, ilo, ihi, z, ldz, workl, nl, info )
                  IF( wantt .OR. info.NE.0 )
     $               CALL dlacpy( 'A', n, n, hl, nl, h, ldh )
               END IF
            END IF
         END IF
*
*        ==== Clear out the trash, if necessary. ====
*
         IF( ( wantt .OR. info.NE.0 ) .AND. n.GT.2 )
     $      CALL dlaset( 'L', n-2, n-2, zero, zero, h( 3, 1 ), ldh )
*
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
*
         work( 1 ) = max( dble( max( 1, n ) ), work( 1 ) )
      END IF
*
*     ==== End of DHSEQR ====
*
      END
      SUBROUTINE dtrevc3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL,
     $                    VR, LDVR, MM, M, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   work( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      INTEGER            NBMIN, NBMAX
      parameter( nbmin = 8, nbmax = 128 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, LQUERY, OVER, PAIR,
     $                   rightv, somev
      INTEGER            I, IERR, II, IP, IS, J, J1, J2, JNXT, K, KI,
     $                   iv, maxwrk, nb, ki2
      DOUBLE PRECISION   BETA, BIGNUM, EMAX, OVFL, REC, REMAX, SCALE,
     $                   smin, smlnum, ulp, unfl, vcrit, vmax, wi, wr,
     $                   xnorm
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DDOT, DLAMCH
      EXTERNAL           lsame, idamax, ilaenv, ddot, dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           daxpy, dcopy, dgemv, dlaln2, dscal, xerbla,
     $                   dgemm, dlaset, dlabad, dlacpy
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, sqrt
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   X( 2, 2 )
      INTEGER            ISCOMPLEX( NBMAX )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      bothv  = lsame( side, 'B' )
      rightv = lsame( side, 'R' ) .OR. bothv
      leftv  = lsame( side, 'L' ) .OR. bothv
*
      allv  = lsame( howmny, 'A' )
      over  = lsame( howmny, 'B' )
      somev = lsame( howmny, 'S' )
*
      info = 0
      nb = ilaenv( 1, 'DTREVC', side // howmny, n, -1, -1, -1 )
      maxwrk = n + 2*n*nb
      work(1) = maxwrk
      lquery = ( lwork.EQ.-1 )
      IF( .NOT.rightv .AND. .NOT.leftv ) THEN
         info = -1
      ELSE IF( .NOT.allv .AND. .NOT.over .AND. .NOT.somev ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( ldt.LT.max( 1, n ) ) THEN
         info = -6
      ELSE IF( ldvl.LT.1 .OR. ( leftv .AND. ldvl.LT.n ) ) THEN
         info = -8
      ELSE IF( ldvr.LT.1 .OR. ( rightv .AND. ldvr.LT.n ) ) THEN
         info = -10
      ELSE IF( lwork.LT.max( 1, 3*n ) .AND. .NOT.lquery ) THEN
         info = -14
      ELSE
*
*        Set M to the number of columns required to store the selected
*        eigenvectors, standardize the array SELECT if necessary, and
*        test MM.
*
         IF( somev ) THEN
            m = 0
            pair = .false.
            DO 10 j = 1, n
               IF( pair ) THEN
                  pair = .false.
                  SELECT( j ) = .false.
               ELSE
                  IF( j.LT.n ) THEN
                     IF( t( j+1, j ).EQ.zero ) THEN
                        IF( SELECT( j ) )
     $                     m = m + 1
                     ELSE
                        pair = .true.
                        IF( SELECT( j ) .OR. SELECT( j+1 ) ) THEN
                           SELECT( j ) = .true.
                           m = m + 2
                        END IF
                     END IF
                  ELSE
                     IF( SELECT( n ) )
     $                  m = m + 1
                  END IF
               END IF
   10       CONTINUE
         ELSE
            m = n
         END IF
*
         IF( mm.LT.m ) THEN
            info = -11
         END IF
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DTREVC3', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Use blocked version of back-transformation if sufficient workspace.
*     Zero-out the workspace to avoid potential NaN propagation.
*
      IF( over .AND. lwork .GE. n + 2*n*nbmin ) THEN
         nb = (lwork - n) / (2*n)
         nb = min( nb, nbmax )
         CALL dlaset( 'F', n, 1+2*nb, zero, zero, work, n )
      ELSE
         nb = 1
      END IF
*
*     Set the constants to control overflow.
*
      unfl = dlamch( 'Safe minimum' )
      ovfl = one / unfl
      CALL dlabad( unfl, ovfl )
      ulp = dlamch( 'Precision' )
      smlnum = unfl*( n / ulp )
      bignum = ( one-ulp ) / smlnum
*
*     Compute 1-norm of each column of strictly upper triangular
*     part of T to control overflow in triangular solver.
*
      work( 1 ) = zero
      DO 30 j = 2, n
         work( j ) = zero
         DO 20 i = 1, j - 1
            work( j ) = work( j ) + abs( t( i, j ) )
   20    CONTINUE
   30 CONTINUE
*
*     Index IP is used to specify the real or complex eigenvalue:
*       IP = 0, real eigenvalue,
*            1, first  of conjugate complex pair: (wr,wi)
*           -1, second of conjugate complex pair: (wr,wi)
*       ISCOMPLEX array stores IP for each column in current block.
*
      IF( rightv ) THEN
*
*        ============================================================
*        Compute right eigenvectors.
*
*        IV is index of column in current block.
*        For complex right vector, uses IV-1 for real part and IV for complex part.
*        Non-blocked version always uses IV=2;
*        blocked     version starts with IV=NB, goes down to 1 or 2.
*        (Note the "0-th" column is used for 1-norms computed above.)
         iv = 2
         IF( nb.GT.2 ) THEN
            iv = nb
         END IF
 
         ip = 0
         is = m
         DO 140 ki = n, 1, -1
            IF( ip.EQ.-1 ) THEN
*              previous iteration (ki+1) was second of conjugate pair,
*              so this ki is first of conjugate pair; skip to end of loop
               ip = 1
               GO TO 140
            ELSE IF( ki.EQ.1 ) THEN
*              last column, so this ki must be real eigenvalue
               ip = 0
            ELSE IF( t( ki, ki-1 ).EQ.zero ) THEN
*              zero on sub-diagonal, so this ki is real eigenvalue
               ip = 0
            ELSE
*              non-zero on sub-diagonal, so this ki is second of conjugate pair
               ip = -1
            END IF
 
            IF( somev ) THEN
               IF( ip.EQ.0 ) THEN
                  IF( .NOT.SELECT( ki ) )
     $               GO TO 140
               ELSE
                  IF( .NOT.SELECT( ki-1 ) )
     $               GO TO 140
               END IF
            END IF
*
*           Compute the KI-th eigenvalue (WR,WI).
*
            wr = t( ki, ki )
            wi = zero
            IF( ip.NE.0 )
     $         wi = sqrt( abs( t( ki, ki-1 ) ) )*
     $              sqrt( abs( t( ki-1, ki ) ) )
            smin = max( ulp*( abs( wr )+abs( wi ) ), smlnum )
*
            IF( ip.EQ.0 ) THEN
*
*              --------------------------------------------------------
*              Real right eigenvector
*
               work( ki + iv*n ) = one
*
*              Form right-hand side.
*
               DO 50 k = 1, ki - 1
                  work( k + iv*n ) = -t( k, ki )
   50          CONTINUE
*
*              Solve upper quasi-triangular system:
*              [ T(1:KI-1,1:KI-1) - WR ]*X = SCALE*WORK.
*
               jnxt = ki - 1
               DO 60 j = ki - 1, 1, -1
                  IF( j.GT.jnxt )
     $               GO TO 60
                  j1 = j
                  j2 = j
                  jnxt = j - 1
                  IF( j.GT.1 ) THEN
                     IF( t( j, j-1 ).NE.zero ) THEN
                        j1   = j - 1
                        jnxt = j - 2
                     END IF
                  END IF
*
                  IF( j1.EQ.j2 ) THEN
*
*                    1-by-1 diagonal block
*
                     CALL dlaln2( .false., 1, 1, smin, one, t( j, j ),
     $                            ldt, one, one, work( j+iv*n ), n, wr,
     $                            zero, x, 2, scale, xnorm, ierr )
*
*                    Scale X(1,1) to avoid overflow when updating
*                    the right-hand side.
*
                     IF( xnorm.GT.one ) THEN
                        IF( work( j ).GT.bignum / xnorm ) THEN
                           x( 1, 1 ) = x( 1, 1 ) / xnorm
                           scale = scale / xnorm
                        END IF
                     END IF
*
*                    Scale if necessary
*
                     IF( scale.NE.one )
     $                  CALL dscal( ki, scale, work( 1+iv*n ), 1 )
                     work( j+iv*n ) = x( 1, 1 )
*
*                    Update right-hand side
*
                     CALL daxpy( j-1, -x( 1, 1 ), t( 1, j ), 1,
     $                           work( 1+iv*n ), 1 )
*
                  ELSE
*
*                    2-by-2 diagonal block
*
                     CALL dlaln2( .false., 2, 1, smin, one,
     $                            t( j-1, j-1 ), ldt, one, one,
     $                            work( j-1+iv*n ), n, wr, zero, x, 2,
     $                            scale, xnorm, ierr )
*
*                    Scale X(1,1) and X(2,1) to avoid overflow when
*                    updating the right-hand side.
*
                     IF( xnorm.GT.one ) THEN
                        beta = max( work( j-1 ), work( j ) )
                        IF( beta.GT.bignum / xnorm ) THEN
                           x( 1, 1 ) = x( 1, 1 ) / xnorm
                           x( 2, 1 ) = x( 2, 1 ) / xnorm
                           scale = scale / xnorm
                        END IF
                     END IF
*
*                    Scale if necessary
*
                     IF( scale.NE.one )
     $                  CALL dscal( ki, scale, work( 1+iv*n ), 1 )
                     work( j-1+iv*n ) = x( 1, 1 )
                     work( j  +iv*n ) = x( 2, 1 )
*
*                    Update right-hand side
*
                     CALL daxpy( j-2, -x( 1, 1 ), t( 1, j-1 ), 1,
     $                           work( 1+iv*n ), 1 )
                     CALL daxpy( j-2, -x( 2, 1 ), t( 1, j ), 1,
     $                           work( 1+iv*n ), 1 )
                  END IF
   60          CONTINUE
*
*              Copy the vector x or Q*x to VR and normalize.
*
               IF( .NOT.over ) THEN
*                 ------------------------------
*                 no back-transform: copy x to VR and normalize.
                  CALL dcopy( ki, work( 1 + iv*n ), 1, vr( 1, is ), 1 )
*
                  ii = idamax( ki, vr( 1, is ), 1 )
                  remax = one / abs( vr( ii, is ) )
                  CALL dscal( ki, remax, vr( 1, is ), 1 )
*
                  DO 70 k = ki + 1, n
                     vr( k, is ) = zero
   70             CONTINUE
*
               ELSE IF( nb.EQ.1 ) THEN
*                 ------------------------------
*                 version 1: back-transform each vector with GEMV, Q*x.
                  IF( ki.GT.1 )
     $               CALL dgemv( 'N', n, ki-1, one, vr, ldvr,
     $                           work( 1 + iv*n ), 1, work( ki + iv*n ),
     $                           vr( 1, ki ), 1 )
*
                  ii = idamax( n, vr( 1, ki ), 1 )
                  remax = one / abs( vr( ii, ki ) )
                  CALL dscal( n, remax, vr( 1, ki ), 1 )
*
               ELSE
*                 ------------------------------
*                 version 2: back-transform block of vectors with GEMM
*                 zero out below vector
                  DO k = ki + 1, n
                     work( k + iv*n ) = zero
                  END DO
                  iscomplex( iv ) = ip
*                 back-transform and normalization is done below
               END IF
            ELSE
*
*              --------------------------------------------------------
*              Complex right eigenvector.
*
*              Initial solve
*              [ ( T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I*WI) ]*X = 0.
*              [ ( T(KI,  KI-1) T(KI,  KI) )               ]
*
               IF( abs( t( ki-1, ki ) ).GE.abs( t( ki, ki-1 ) ) ) THEN
                  work( ki-1 + (iv-1)*n ) = one
                  work( ki   + (iv  )*n ) = wi / t( ki-1, ki )
               ELSE
                  work( ki-1 + (iv-1)*n ) = -wi / t( ki, ki-1 )
                  work( ki   + (iv  )*n ) = one
               END IF
               work( ki   + (iv-1)*n ) = zero
               work( ki-1 + (iv  )*n ) = zero
*
*              Form right-hand side.
*
               DO 80 k = 1, ki - 2
                  work( k+(iv-1)*n ) = -work( ki-1+(iv-1)*n )*t(k,ki-1)
                  work( k+(iv  )*n ) = -work( ki  +(iv  )*n )*t(k,ki  )
   80          CONTINUE
*
*              Solve upper quasi-triangular system:
*              [ T(1:KI-2,1:KI-2) - (WR+i*WI) ]*X = SCALE*(WORK+i*WORK2)
*
               jnxt = ki - 2
               DO 90 j = ki - 2, 1, -1
                  IF( j.GT.jnxt )
     $               GO TO 90
                  j1 = j
                  j2 = j
                  jnxt = j - 1
                  IF( j.GT.1 ) THEN
                     IF( t( j, j-1 ).NE.zero ) THEN
                        j1   = j - 1
                        jnxt = j - 2
                     END IF
                  END IF
*
                  IF( j1.EQ.j2 ) THEN
*
*                    1-by-1 diagonal block
*
                     CALL dlaln2( .false., 1, 2, smin, one, t( j, j ),
     $                            ldt, one, one, work( j+(iv-1)*n ), n,
     $                            wr, wi, x, 2, scale, xnorm, ierr )
*
*                    Scale X(1,1) and X(1,2) to avoid overflow when
*                    updating the right-hand side.
*
                     IF( xnorm.GT.one ) THEN
                        IF( work( j ).GT.bignum / xnorm ) THEN
                           x( 1, 1 ) = x( 1, 1 ) / xnorm
                           x( 1, 2 ) = x( 1, 2 ) / xnorm
                           scale = scale / xnorm
                        END IF
                     END IF
*
*                    Scale if necessary
*
                     IF( scale.NE.one ) THEN
                        CALL dscal( ki, scale, work( 1+(iv-1)*n ), 1 )
                        CALL dscal( ki, scale, work( 1+(iv  )*n ), 1 )
                     END IF
                     work( j+(iv-1)*n ) = x( 1, 1 )
                     work( j+(iv  )*n ) = x( 1, 2 )
*
*                    Update the right-hand side
*
                     CALL daxpy( j-1, -x( 1, 1 ), t( 1, j ), 1,
     $                           work( 1+(iv-1)*n ), 1 )
                     CALL daxpy( j-1, -x( 1, 2 ), t( 1, j ), 1,
     $                           work( 1+(iv  )*n ), 1 )
*
                  ELSE
*
*                    2-by-2 diagonal block
*
                     CALL dlaln2( .false., 2, 2, smin, one,
     $                            t( j-1, j-1 ), ldt, one, one,
     $                            work( j-1+(iv-1)*n ), n, wr, wi, x, 2,
     $                            scale, xnorm, ierr )
*
*                    Scale X to avoid overflow when updating
*                    the right-hand side.
*
                     IF( xnorm.GT.one ) THEN
                        beta = max( work( j-1 ), work( j ) )
                        IF( beta.GT.bignum / xnorm ) THEN
                           rec = one / xnorm
                           x( 1, 1 ) = x( 1, 1 )*rec
                           x( 1, 2 ) = x( 1, 2 )*rec
                           x( 2, 1 ) = x( 2, 1 )*rec
                           x( 2, 2 ) = x( 2, 2 )*rec
                           scale = scale*rec
                        END IF
                     END IF
*
*                    Scale if necessary
*
                     IF( scale.NE.one ) THEN
                        CALL dscal( ki, scale, work( 1+(iv-1)*n ), 1 )
                        CALL dscal( ki, scale, work( 1+(iv  )*n ), 1 )
                     END IF
                     work( j-1+(iv-1)*n ) = x( 1, 1 )
                     work( j  +(iv-1)*n ) = x( 2, 1 )
                     work( j-1+(iv  )*n ) = x( 1, 2 )
                     work( j  +(iv  )*n ) = x( 2, 2 )
*
*                    Update the right-hand side
*
                     CALL daxpy( j-2, -x( 1, 1 ), t( 1, j-1 ), 1,
     $                           work( 1+(iv-1)*n   ), 1 )
                     CALL daxpy( j-2, -x( 2, 1 ), t( 1, j ), 1,
     $                           work( 1+(iv-1)*n   ), 1 )
                     CALL daxpy( j-2, -x( 1, 2 ), t( 1, j-1 ), 1,
     $                           work( 1+(iv  )*n ), 1 )
                     CALL daxpy( j-2, -x( 2, 2 ), t( 1, j ), 1,
     $                           work( 1+(iv  )*n ), 1 )
                  END IF
   90          CONTINUE
*
*              Copy the vector x or Q*x to VR and normalize.
*
               IF( .NOT.over ) THEN
*                 ------------------------------
*                 no back-transform: copy x to VR and normalize.
                  CALL dcopy( ki, work( 1+(iv-1)*n ), 1, vr(1,is-1), 1 )
                  CALL dcopy( ki, work( 1+(iv  )*n ), 1, vr(1,is  ), 1 )
*
                  emax = zero
                  DO 100 k = 1, ki
                     emax = max( emax, abs( vr( k, is-1 ) )+
     $                                 abs( vr( k, is   ) ) )
  100             CONTINUE
                  remax = one / emax
                  CALL dscal( ki, remax, vr( 1, is-1 ), 1 )
                  CALL dscal( ki, remax, vr( 1, is   ), 1 )
*
                  DO 110 k = ki + 1, n
                     vr( k, is-1 ) = zero
                     vr( k, is   ) = zero
  110             CONTINUE
*
               ELSE IF( nb.EQ.1 ) THEN
*                 ------------------------------
*                 version 1: back-transform each vector with GEMV, Q*x.
                  IF( ki.GT.2 ) THEN
                     CALL dgemv( 'N', n, ki-2, one, vr, ldvr,
     $                           work( 1    + (iv-1)*n ), 1,
     $                           work( ki-1 + (iv-1)*n ), vr(1,ki-1), 1)
                     CALL dgemv( 'N', n, ki-2, one, vr, ldvr,
     $                           work( 1  + (iv)*n ), 1,
     $                           work( ki + (iv)*n ), vr( 1, ki ), 1 )
                  ELSE
                     CALL dscal( n, work(ki-1+(iv-1)*n), vr(1,ki-1), 1)
                     CALL dscal( n, work(ki  +(iv  )*n), vr(1,ki  ), 1)
                  END IF
*
                  emax = zero
                  DO 120 k = 1, n
                     emax = max( emax, abs( vr( k, ki-1 ) )+
     $                                 abs( vr( k, ki   ) ) )
  120             CONTINUE
                  remax = one / emax
                  CALL dscal( n, remax, vr( 1, ki-1 ), 1 )
                  CALL dscal( n, remax, vr( 1, ki   ), 1 )
*
               ELSE
*                 ------------------------------
*                 version 2: back-transform block of vectors with GEMM
*                 zero out below vector
                  DO k = ki + 1, n
                     work( k + (iv-1)*n ) = zero
                     work( k + (iv  )*n ) = zero
                  END DO
                  iscomplex( iv-1 ) = -ip
                  iscomplex( iv   ) =  ip
                  iv = iv - 1
*                 back-transform and normalization is done below
               END IF
            END IF
 
            IF( nb.GT.1 ) THEN
*              --------------------------------------------------------
*              Blocked version of back-transform
*              For complex case, KI2 includes both vectors (KI-1 and KI)
               IF( ip.EQ.0 ) THEN
                  ki2 = ki
               ELSE
                  ki2 = ki - 1
               END IF
 
*              Columns IV:NB of work are valid vectors.
*              When the number of vectors stored reaches NB-1 or NB,
*              or if this was last vector, do the GEMM
               IF( (iv.LE.2) .OR. (ki2.EQ.1) ) THEN
                  CALL dgemm( 'N', 'N', n, nb-iv+1, ki2+nb-iv, one,
     $                        vr, ldvr,
     $                        work( 1 + (iv)*n    ), n,
     $                        zero,
     $                        work( 1 + (nb+iv)*n ), n )
*                 normalize vectors
                  DO k = iv, nb
                     IF( iscomplex(k).EQ.0 ) THEN
*                       real eigenvector
                        ii = idamax( n, work( 1 + (nb+k)*n ), 1 )
                        remax = one / abs( work( ii + (nb+k)*n ) )
                     ELSE IF( iscomplex(k).EQ.1 ) THEN
*                       first eigenvector of conjugate pair
                        emax = zero
                        DO ii = 1, n
                           emax = max( emax,
     $                                 abs( work( ii + (nb+k  )*n ) )+
     $                                 abs( work( ii + (nb+k+1)*n ) ) )
                        END DO
                        remax = one / emax
*                    else if ISCOMPLEX(K).EQ.-1
*                       second eigenvector of conjugate pair
*                       reuse same REMAX as previous K
                     END IF
                     CALL dscal( n, remax, work( 1 + (nb+k)*n ), 1 )
                  END DO
                  CALL dlacpy( 'F', n, nb-iv+1,
     $                         work( 1 + (nb+iv)*n ), n,
     $                         vr( 1, ki2 ), ldvr )
                  iv = nb
               ELSE
                  iv = iv - 1
               END IF
            END IF ! blocked back-transform
*
            is = is - 1
            IF( ip.NE.0 )
     $         is = is - 1
  140    CONTINUE
      END IF
 
      IF( leftv ) THEN
*
*        ============================================================
*        Compute left eigenvectors.
*
*        IV is index of column in current block.
*        For complex left vector, uses IV for real part and IV+1 for complex part.
*        Non-blocked version always uses IV=1;
*        blocked     version starts with IV=1, goes up to NB-1 or NB.
*        (Note the "0-th" column is used for 1-norms computed above.)
         iv = 1
         ip = 0
         is = 1
         DO 260 ki = 1, n
            IF( ip.EQ.1 ) THEN
*              previous iteration (ki-1) was first of conjugate pair,
*              so this ki is second of conjugate pair; skip to end of loop
               ip = -1
               GO TO 260
            ELSE IF( ki.EQ.n ) THEN
*              last column, so this ki must be real eigenvalue
               ip = 0
            ELSE IF( t( ki+1, ki ).EQ.zero ) THEN
*              zero on sub-diagonal, so this ki is real eigenvalue
               ip = 0
            ELSE
*              non-zero on sub-diagonal, so this ki is first of conjugate pair
               ip = 1
            END IF
*
            IF( somev ) THEN
               IF( .NOT.SELECT( ki ) )
     $            GO TO 260
            END IF
*
*           Compute the KI-th eigenvalue (WR,WI).
*
            wr = t( ki, ki )
            wi = zero
            IF( ip.NE.0 )
     $         wi = sqrt( abs( t( ki, ki+1 ) ) )*
     $              sqrt( abs( t( ki+1, ki ) ) )
            smin = max( ulp*( abs( wr )+abs( wi ) ), smlnum )
*
            IF( ip.EQ.0 ) THEN
*
*              --------------------------------------------------------
*              Real left eigenvector
*
               work( ki + iv*n ) = one
*
*              Form right-hand side.
*
               DO 160 k = ki + 1, n
                  work( k + iv*n ) = -t( ki, k )
  160          CONTINUE
*
*              Solve transposed quasi-triangular system:
*              [ T(KI+1:N,KI+1:N) - WR ]**T * X = SCALE*WORK
*
               vmax = one
               vcrit = bignum
*
               jnxt = ki + 1
               DO 170 j = ki + 1, n
                  IF( j.LT.jnxt )
     $               GO TO 170
                  j1 = j
                  j2 = j
                  jnxt = j + 1
                  IF( j.LT.n ) THEN
                     IF( t( j+1, j ).NE.zero ) THEN
                        j2 = j + 1
                        jnxt = j + 2
                     END IF
                  END IF
*
                  IF( j1.EQ.j2 ) THEN
*
*                    1-by-1 diagonal block
*
*                    Scale if necessary to avoid overflow when forming
*                    the right-hand side.
*
                     IF( work( j ).GT.vcrit ) THEN
                        rec = one / vmax
                        CALL dscal( n-ki+1, rec, work( ki+iv*n ), 1 )
                        vmax = one
                        vcrit = bignum
                     END IF
*
                     work( j+iv*n ) = work( j+iv*n ) -
     $                                ddot( j-ki-1, t( ki+1, j ), 1,
     $                                      work( ki+1+iv*n ), 1 )
*
*                    Solve [ T(J,J) - WR ]**T * X = WORK
*
                     CALL dlaln2( .false., 1, 1, smin, one, t( j, j ),
     $                            ldt, one, one, work( j+iv*n ), n, wr,
     $                            zero, x, 2, scale, xnorm, ierr )
*
*                    Scale if necessary
*
                     IF( scale.NE.one )
     $                  CALL dscal( n-ki+1, scale, work( ki+iv*n ), 1 )
                     work( j+iv*n ) = x( 1, 1 )
                     vmax = max( abs( work( j+iv*n ) ), vmax )
                     vcrit = bignum / vmax
*
                  ELSE
*
*                    2-by-2 diagonal block
*
*                    Scale if necessary to avoid overflow when forming
*                    the right-hand side.
*
                     beta = max( work( j ), work( j+1 ) )
                     IF( beta.GT.vcrit ) THEN
                        rec = one / vmax
                        CALL dscal( n-ki+1, rec, work( ki+iv*n ), 1 )
                        vmax = one
                        vcrit = bignum
                     END IF
*
                     work( j+iv*n ) = work( j+iv*n ) -
     $                                ddot( j-ki-1, t( ki+1, j ), 1,
     $                                      work( ki+1+iv*n ), 1 )
*
                     work( j+1+iv*n ) = work( j+1+iv*n ) -
     $                                  ddot( j-ki-1, t( ki+1, j+1 ), 1,
     $                                        work( ki+1+iv*n ), 1 )
*
*                    Solve
*                    [ T(J,J)-WR   T(J,J+1)      ]**T * X = SCALE*( WORK1 )
*                    [ T(J+1,J)    T(J+1,J+1)-WR ]                ( WORK2 )
*
                     CALL dlaln2( .true., 2, 1, smin, one, t( j, j ),
     $                            ldt, one, one, work( j+iv*n ), n, wr,
     $                            zero, x, 2, scale, xnorm, ierr )
*
*                    Scale if necessary
*
                     IF( scale.NE.one )
     $                  CALL dscal( n-ki+1, scale, work( ki+iv*n ), 1 )
                     work( j  +iv*n ) = x( 1, 1 )
                     work( j+1+iv*n ) = x( 2, 1 )
*
                     vmax = max( abs( work( j  +iv*n ) ),
     $                           abs( work( j+1+iv*n ) ), vmax )
                     vcrit = bignum / vmax
*
                  END IF
  170          CONTINUE
*
*              Copy the vector x or Q*x to VL and normalize.
*
               IF( .NOT.over ) THEN
*                 ------------------------------
*                 no back-transform: copy x to VL and normalize.
                  CALL dcopy( n-ki+1, work( ki + iv*n ), 1,
     $                                vl( ki, is ), 1 )
*
                  ii = idamax( n-ki+1, vl( ki, is ), 1 ) + ki - 1
                  remax = one / abs( vl( ii, is ) )
                  CALL dscal( n-ki+1, remax, vl( ki, is ), 1 )
*
                  DO 180 k = 1, ki - 1
                     vl( k, is ) = zero
  180             CONTINUE
*
               ELSE IF( nb.EQ.1 ) THEN
*                 ------------------------------
*                 version 1: back-transform each vector with GEMV, Q*x.
                  IF( ki.LT.n )
     $               CALL dgemv( 'N', n, n-ki, one,
     $                           vl( 1, ki+1 ), ldvl,
     $                           work( ki+1 + iv*n ), 1,
     $                           work( ki   + iv*n ), vl( 1, ki ), 1 )
*
                  ii = idamax( n, vl( 1, ki ), 1 )
                  remax = one / abs( vl( ii, ki ) )
                  CALL dscal( n, remax, vl( 1, ki ), 1 )
*
               ELSE
*                 ------------------------------
*                 version 2: back-transform block of vectors with GEMM
*                 zero out above vector
*                 could go from KI-NV+1 to KI-1
                  DO k = 1, ki - 1
                     work( k + iv*n ) = zero
                  END DO
                  iscomplex( iv ) = ip
*                 back-transform and normalization is done below
               END IF
            ELSE
*
*              --------------------------------------------------------
*              Complex left eigenvector.
*
*              Initial solve:
*              [ ( T(KI,KI)    T(KI,KI+1)  )**T - (WR - I* WI) ]*X = 0.
*              [ ( T(KI+1,KI) T(KI+1,KI+1) )                   ]
*
               IF( abs( t( ki, ki+1 ) ).GE.abs( t( ki+1, ki ) ) ) THEN
                  work( ki   + (iv  )*n ) = wi / t( ki, ki+1 )
                  work( ki+1 + (iv+1)*n ) = one
               ELSE
                  work( ki   + (iv  )*n ) = one
                  work( ki+1 + (iv+1)*n ) = -wi / t( ki+1, ki )
               END IF
               work( ki+1 + (iv  )*n ) = zero
               work( ki   + (iv+1)*n ) = zero
*
*              Form right-hand side.
*
               DO 190 k = ki + 2, n
                  work( k+(iv  )*n ) = -work( ki  +(iv  )*n )*t(ki,  k)
                  work( k+(iv+1)*n ) = -work( ki+1+(iv+1)*n )*t(ki+1,k)
  190          CONTINUE
*
*              Solve transposed quasi-triangular system:
*              [ T(KI+2:N,KI+2:N)**T - (WR-i*WI) ]*X = WORK1+i*WORK2
*
               vmax = one
               vcrit = bignum
*
               jnxt = ki + 2
               DO 200 j = ki + 2, n
                  IF( j.LT.jnxt )
     $               GO TO 200
                  j1 = j
                  j2 = j
                  jnxt = j + 1
                  IF( j.LT.n ) THEN
                     IF( t( j+1, j ).NE.zero ) THEN
                        j2 = j + 1
                        jnxt = j + 2
                     END IF
                  END IF
*
                  IF( j1.EQ.j2 ) THEN
*
*                    1-by-1 diagonal block
*
*                    Scale if necessary to avoid overflow when
*                    forming the right-hand side elements.
*
                     IF( work( j ).GT.vcrit ) THEN
                        rec = one / vmax
                        CALL dscal( n-ki+1, rec, work(ki+(iv  )*n), 1 )
                        CALL dscal( n-ki+1, rec, work(ki+(iv+1)*n), 1 )
                        vmax = one
                        vcrit = bignum
                     END IF
*
                     work( j+(iv  )*n ) = work( j+(iv)*n ) -
     $                                  ddot( j-ki-2, t( ki+2, j ), 1,
     $                                        work( ki+2+(iv)*n ), 1 )
                     work( j+(iv+1)*n ) = work( j+(iv+1)*n ) -
     $                                  ddot( j-ki-2, t( ki+2, j ), 1,
     $                                        work( ki+2+(iv+1)*n ), 1 )
*
*                    Solve [ T(J,J)-(WR-i*WI) ]*(X11+i*X12)= WK+I*WK2
*
                     CALL dlaln2( .false., 1, 2, smin, one, t( j, j ),
     $                            ldt, one, one, work( j+iv*n ), n, wr,
     $                            -wi, x, 2, scale, xnorm, ierr )
*
*                    Scale if necessary
*
                     IF( scale.NE.one ) THEN
                        CALL dscal( n-ki+1, scale, work(ki+(iv  )*n), 1)
                        CALL dscal( n-ki+1, scale, work(ki+(iv+1)*n), 1)
                     END IF
                     work( j+(iv  )*n ) = x( 1, 1 )
                     work( j+(iv+1)*n ) = x( 1, 2 )
                     vmax = max( abs( work( j+(iv  )*n ) ),
     $                           abs( work( j+(iv+1)*n ) ), vmax )
                     vcrit = bignum / vmax
*
                  ELSE
*
*                    2-by-2 diagonal block
*
*                    Scale if necessary to avoid overflow when forming
*                    the right-hand side elements.
*
                     beta = max( work( j ), work( j+1 ) )
                     IF( beta.GT.vcrit ) THEN
                        rec = one / vmax
                        CALL dscal( n-ki+1, rec, work(ki+(iv  )*n), 1 )
                        CALL dscal( n-ki+1, rec, work(ki+(iv+1)*n), 1 )
                        vmax = one
                        vcrit = bignum
                     END IF
*
                     work( j  +(iv  )*n ) = work( j+(iv)*n ) -
     $                                ddot( j-ki-2, t( ki+2, j ), 1,
     $                                      work( ki+2+(iv)*n ), 1 )
*
                     work( j  +(iv+1)*n ) = work( j+(iv+1)*n ) -
     $                                ddot( j-ki-2, t( ki+2, j ), 1,
     $                                      work( ki+2+(iv+1)*n ), 1 )
*
                     work( j+1+(iv  )*n ) = work( j+1+(iv)*n ) -
     $                                ddot( j-ki-2, t( ki+2, j+1 ), 1,
     $                                      work( ki+2+(iv)*n ), 1 )
*
                     work( j+1+(iv+1)*n ) = work( j+1+(iv+1)*n ) -
     $                                ddot( j-ki-2, t( ki+2, j+1 ), 1,
     $                                      work( ki+2+(iv+1)*n ), 1 )
*
*                    Solve 2-by-2 complex linear equation
*                    [ (T(j,j)   T(j,j+1)  )**T - (wr-i*wi)*I ]*X = SCALE*B
*                    [ (T(j+1,j) T(j+1,j+1))                  ]
*
                     CALL dlaln2( .true., 2, 2, smin, one, t( j, j ),
     $                            ldt, one, one, work( j+iv*n ), n, wr,
     $                            -wi, x, 2, scale, xnorm, ierr )
*
*                    Scale if necessary
*
                     IF( scale.NE.one ) THEN
                        CALL dscal( n-ki+1, scale, work(ki+(iv  )*n), 1)
                        CALL dscal( n-ki+1, scale, work(ki+(iv+1)*n), 1)
                     END IF
                     work( j  +(iv  )*n ) = x( 1, 1 )
                     work( j  +(iv+1)*n ) = x( 1, 2 )
                     work( j+1+(iv  )*n ) = x( 2, 1 )
                     work( j+1+(iv+1)*n ) = x( 2, 2 )
                     vmax = max( abs( x( 1, 1 ) ), abs( x( 1, 2 ) ),
     $                           abs( x( 2, 1 ) ), abs( x( 2, 2 ) ),
     $                           vmax )
                     vcrit = bignum / vmax
*
                  END IF
  200          CONTINUE
*
*              Copy the vector x or Q*x to VL and normalize.
*
               IF( .NOT.over ) THEN
*                 ------------------------------
*                 no back-transform: copy x to VL and normalize.
                  CALL dcopy( n-ki+1, work( ki + (iv  )*n ), 1,
     $                        vl( ki, is   ), 1 )
                  CALL dcopy( n-ki+1, work( ki + (iv+1)*n ), 1,
     $                        vl( ki, is+1 ), 1 )
*
                  emax = zero
                  DO 220 k = ki, n
                     emax = max( emax, abs( vl( k, is   ) )+
     $                                 abs( vl( k, is+1 ) ) )
  220             CONTINUE
                  remax = one / emax
                  CALL dscal( n-ki+1, remax, vl( ki, is   ), 1 )
                  CALL dscal( n-ki+1, remax, vl( ki, is+1 ), 1 )
*
                  DO 230 k = 1, ki - 1
                     vl( k, is   ) = zero
                     vl( k, is+1 ) = zero
  230             CONTINUE
*
               ELSE IF( nb.EQ.1 ) THEN
*                 ------------------------------
*                 version 1: back-transform each vector with GEMV, Q*x.
                  IF( ki.LT.n-1 ) THEN
                     CALL dgemv( 'N', n, n-ki-1, one,
     $                           vl( 1, ki+2 ), ldvl,
     $                           work( ki+2 + (iv)*n ), 1,
     $                           work( ki   + (iv)*n ),
     $                           vl( 1, ki ), 1 )
                     CALL dgemv( 'N', n, n-ki-1, one,
     $                           vl( 1, ki+2 ), ldvl,
     $                           work( ki+2 + (iv+1)*n ), 1,
     $                           work( ki+1 + (iv+1)*n ),
     $                           vl( 1, ki+1 ), 1 )
                  ELSE
                     CALL dscal( n, work(ki+  (iv  )*n), vl(1, ki  ), 1)
                     CALL dscal( n, work(ki+1+(iv+1)*n), vl(1, ki+1), 1)
                  END IF
*
                  emax = zero
                  DO 240 k = 1, n
                     emax = max( emax, abs( vl( k, ki   ) )+
     $                                 abs( vl( k, ki+1 ) ) )
  240             CONTINUE
                  remax = one / emax
                  CALL dscal( n, remax, vl( 1, ki   ), 1 )
                  CALL dscal( n, remax, vl( 1, ki+1 ), 1 )
*
               ELSE
*                 ------------------------------
*                 version 2: back-transform block of vectors with GEMM
*                 zero out above vector
*                 could go from KI-NV+1 to KI-1
                  DO k = 1, ki - 1
                     work( k + (iv  )*n ) = zero
                     work( k + (iv+1)*n ) = zero
                  END DO
                  iscomplex( iv   ) =  ip
                  iscomplex( iv+1 ) = -ip
                  iv = iv + 1
*                 back-transform and normalization is done below
               END IF
            END IF
 
            IF( nb.GT.1 ) THEN
*              --------------------------------------------------------
*              Blocked version of back-transform
*              For complex case, KI2 includes both vectors (KI and KI+1)
               IF( ip.EQ.0 ) THEN
                  ki2 = ki
               ELSE
                  ki2 = ki + 1
               END IF
 
*              Columns 1:IV of work are valid vectors.
*              When the number of vectors stored reaches NB-1 or NB,
*              or if this was last vector, do the GEMM
               IF( (iv.GE.nb-1) .OR. (ki2.EQ.n) ) THEN
                  CALL dgemm( 'N', 'N', n, iv, n-ki2+iv, one,
     $                        vl( 1, ki2-iv+1 ), ldvl,
     $                        work( ki2-iv+1 + (1)*n ), n,
     $                        zero,
     $                        work( 1 + (nb+1)*n ), n )
*                 normalize vectors
                  DO k = 1, iv
                     IF( iscomplex(k).EQ.0) THEN
*                       real eigenvector
                        ii = idamax( n, work( 1 + (nb+k)*n ), 1 )
                        remax = one / abs( work( ii + (nb+k)*n ) )
                     ELSE IF( iscomplex(k).EQ.1) THEN
*                       first eigenvector of conjugate pair
                        emax = zero
                        DO ii = 1, n
                           emax = max( emax,
     $                                 abs( work( ii + (nb+k  )*n ) )+
     $                                 abs( work( ii + (nb+k+1)*n ) ) )
                        END DO
                        remax = one / emax
*                    else if ISCOMPLEX(K).EQ.-1
*                       second eigenvector of conjugate pair
*                       reuse same REMAX as previous K
                     END IF
                     CALL dscal( n, remax, work( 1 + (nb+k)*n ), 1 )
                  END DO
                  CALL dlacpy( 'F', n, iv,
     $                         work( 1 + (nb+1)*n ), n,
     $                         vl( 1, ki2-iv+1 ), ldvl )
                  iv = 1
               ELSE
                  iv = iv + 1
               END IF
            END IF ! blocked back-transform
*
            is = is + 1
            IF( ip.NE.0 )
     $         is = is + 1
  260    CONTINUE
      END IF
*
      RETURN
*
*     End of DTREVC3
*
      END
      SUBROUTINE zlatrs( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
     $                   CNORM, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   CNORM( * )
      COMPLEX*16         A( LDA, * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      parameter( zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0,
     $                   two = 2.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL,
     $                   xbnd, xj, xmax
      COMPLEX*16         CSUMJ, TJJS, USCAL, ZDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      COMPLEX*16         ZDOTC, ZDOTU, ZLADIV
      EXTERNAL           lsame, idamax, izamax, dlamch, dzasum, zdotc,
     $                   zdotu, zladiv
*     ..
*     .. External Subroutines ..
      EXTERNAL           dscal, xerbla, zaxpy, zdscal, ztrsv, dlabad
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dconjg, dimag, max, min
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1, CABS2
*     ..
*     .. Statement Function definitions ..
      cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )
      cabs2( zdum ) = abs( dble( zdum ) / 2.d0 ) +
     $                abs( dimag( zdum ) / 2.d0 )
*     ..
*     .. Executable Statements ..
*
      info = 0
      upper = lsame( uplo, 'U' )
      notran = lsame( trans, 'N' )
      nounit = lsame( diag, 'N' )
*
*     Test the input parameters.
*
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) .AND. .NOT.
     $         lsame( trans, 'C' ) ) THEN
         info = -2
      ELSE IF( .NOT.nounit .AND. .NOT.lsame( diag, 'U' ) ) THEN
         info = -3
      ELSE IF( .NOT.lsame( normin, 'Y' ) .AND. .NOT.
     $         lsame( normin, 'N' ) ) THEN
         info = -4
      ELSE IF( n.LT.0 ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -7
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZLATRS', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Determine machine dependent parameters to control overflow.
*
      smlnum = dlamch( 'Safe minimum' )
      bignum = one / smlnum
      CALL dlabad( smlnum, bignum )
      smlnum = smlnum / dlamch( 'Precision' )
      bignum = one / smlnum
      scale = one
*
      IF( lsame( normin, 'N' ) ) THEN
*
*        Compute the 1-norm of each column, not including the diagonal.
*
         IF( upper ) THEN
*
*           A is upper triangular.
*
            DO 10 j = 1, n
               cnorm( j ) = dzasum( j-1, a( 1, j ), 1 )
   10       CONTINUE
         ELSE
*
*           A is lower triangular.
*
            DO 20 j = 1, n - 1
               cnorm( j ) = dzasum( n-j, a( j+1, j ), 1 )
   20       CONTINUE
            cnorm( n ) = zero
         END IF
      END IF
*
*     Scale the column norms by TSCAL if the maximum element in CNORM is
*     greater than BIGNUM/2.
*
      imax = idamax( n, cnorm, 1 )
      tmax = cnorm( imax )
      IF( tmax.LE.bignum*half ) THEN
         tscal = one
      ELSE
         tscal = half / ( smlnum*tmax )
         CALL dscal( n, tscal, cnorm, 1 )
      END IF
*
*     Compute a bound on the computed solution vector to see if the
*     Level 2 BLAS routine ZTRSV can be used.
*
      xmax = zero
      DO 30 j = 1, n
         xmax = max( xmax, cabs2( x( j ) ) )
   30 CONTINUE
      xbnd = xmax
*
      IF( notran ) THEN
*
*        Compute the growth in A * x = b.
*
         IF( upper ) THEN
            jfirst = n
            jlast = 1
            jinc = -1
         ELSE
            jfirst = 1
            jlast = n
            jinc = 1
         END IF
*
         IF( tscal.NE.one ) THEN
            grow = zero
            GO TO 60
         END IF
*
         IF( nounit ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, G(0) = max{x(i), i=1,...,n}.
*
            grow = half / max( xbnd, smlnum )
            xbnd = grow
            DO 40 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 60
*
               tjjs = a( j, j )
               tjj = cabs1( tjjs )
*
               IF( tjj.GE.smlnum ) THEN
*
*                 M(j) = G(j-1) / abs(A(j,j))
*
                  xbnd = min( xbnd, min( one, tjj )*grow )
               ELSE
*
*                 M(j) could overflow, set XBND to 0.
*
                  xbnd = zero
               END IF
*
               IF( tjj+cnorm( j ).GE.smlnum ) THEN
*
*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
*
                  grow = grow*( tjj / ( tjj+cnorm( j ) ) )
               ELSE
*
*                 G(j) could overflow, set GROW to 0.
*
                  grow = zero
               END IF
   40       CONTINUE
            grow = xbnd
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            grow = min( one, half / max( xbnd, smlnum ) )
            DO 50 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 60
*
*              G(j) = G(j-1)*( 1 + CNORM(j) )
*
               grow = grow*( one / ( one+cnorm( j ) ) )
   50       CONTINUE
         END IF
   60    CONTINUE
*
      ELSE
*
*        Compute the growth in A**T * x = b  or  A**H * x = b.
*
         IF( upper ) THEN
            jfirst = 1
            jlast = n
            jinc = 1
         ELSE
            jfirst = n
            jlast = 1
            jinc = -1
         END IF
*
         IF( tscal.NE.one ) THEN
            grow = zero
            GO TO 90
         END IF
*
         IF( nounit ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, M(0) = max{x(i), i=1,...,n}.
*
            grow = half / max( xbnd, smlnum )
            xbnd = grow
            DO 70 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 90
*
*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
*
               xj = one + cnorm( j )
               grow = min( grow, xbnd / xj )
*
               tjjs = a( j, j )
               tjj = cabs1( tjjs )
*
               IF( tjj.GE.smlnum ) THEN
*
*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
*
                  IF( xj.GT.tjj )
     $               xbnd = xbnd*( tjj / xj )
               ELSE
*
*                 M(j) could overflow, set XBND to 0.
*
                  xbnd = zero
               END IF
   70       CONTINUE
            grow = min( grow, xbnd )
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            grow = min( one, half / max( xbnd, smlnum ) )
            DO 80 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 90
*
*              G(j) = ( 1 + CNORM(j) )*G(j-1)
*
               xj = one + cnorm( j )
               grow = grow / xj
   80       CONTINUE
         END IF
   90    CONTINUE
      END IF
*
      IF( ( grow*tscal ).GT.smlnum ) THEN
*
*        Use the Level 2 BLAS solve if the reciprocal of the bound on
*        elements of X is not too small.
*
         CALL ztrsv( uplo, trans, diag, n, a, lda, x, 1 )
      ELSE
*
*        Use a Level 1 BLAS solve, scaling intermediate results.
*
         IF( xmax.GT.bignum*half ) THEN
*
*           Scale X so that its components are less than or equal to
*           BIGNUM in absolute value.
*
            scale = ( bignum*half ) / xmax
            CALL zdscal( n, scale, x, 1 )
            xmax = bignum
         ELSE
            xmax = xmax*two
         END IF
*
         IF( notran ) THEN
*
*           Solve A * x = b
*
            DO 120 j = jfirst, jlast, jinc
*
*              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
*
               xj = cabs1( x( j ) )
               IF( nounit ) THEN
                  tjjs = a( j, j )*tscal
               ELSE
                  tjjs = tscal
                  IF( tscal.EQ.one )
     $               GO TO 110
               END IF
               tjj = cabs1( tjjs )
               IF( tjj.GT.smlnum ) THEN
*
*                    abs(A(j,j)) > SMLNUM:
*
                  IF( tjj.LT.one ) THEN
                     IF( xj.GT.tjj*bignum ) THEN
*
*                          Scale x by 1/b(j).
*
                        rec = one / xj
                        CALL zdscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     END IF
                  END IF
                  x( j ) = zladiv( x( j ), tjjs )
                  xj = cabs1( x( j ) )
               ELSE IF( tjj.GT.zero ) THEN
*
*                    0 < abs(A(j,j)) <= SMLNUM:
*
                  IF( xj.GT.tjj*bignum ) THEN
*
*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
*                       to avoid overflow when dividing by A(j,j).
*
                     rec = ( tjj*bignum ) / xj
                     IF( cnorm( j ).GT.one ) THEN
*
*                          Scale by 1/CNORM(j) to avoid overflow when
*                          multiplying x(j) times column j.
*
                        rec = rec / cnorm( j )
                     END IF
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  END IF
                  x( j ) = zladiv( x( j ), tjjs )
                  xj = cabs1( x( j ) )
               ELSE
*
*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to A*x = 0.
*
                  DO 100 i = 1, n
                     x( i ) = zero
  100             CONTINUE
                  x( j ) = one
                  xj = one
                  scale = zero
                  xmax = zero
               END IF
  110          CONTINUE
*
*              Scale x if necessary to avoid overflow when adding a
*              multiple of column j of A.
*
               IF( xj.GT.one ) THEN
                  rec = one / xj
                  IF( cnorm( j ).GT.( bignum-xmax )*rec ) THEN
*
*                    Scale x by 1/(2*abs(x(j))).
*
                     rec = rec*half
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                  END IF
               ELSE IF( xj*cnorm( j ).GT.( bignum-xmax ) ) THEN
*
*                 Scale x by 1/2.
*
                  CALL zdscal( n, half, x, 1 )
                  scale = scale*half
               END IF
*
               IF( upper ) THEN
                  IF( j.GT.1 ) THEN
*
*                    Compute the update
*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
*
                     CALL zaxpy( j-1, -x( j )*tscal, a( 1, j ), 1, x,
     $                           1 )
                     i = izamax( j-1, x, 1 )
                     xmax = cabs1( x( i ) )
                  END IF
               ELSE
                  IF( j.LT.n ) THEN
*
*                    Compute the update
*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
*
                     CALL zaxpy( n-j, -x( j )*tscal, a( j+1, j ), 1,
     $                           x( j+1 ), 1 )
                     i = j + izamax( n-j, x( j+1 ), 1 )
                     xmax = cabs1( x( i ) )
                  END IF
               END IF
  120       CONTINUE
*
         ELSE IF( lsame( trans, 'T' ) ) THEN
*
*           Solve A**T * x = b
*
            DO 170 j = jfirst, jlast, jinc
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               xj = cabs1( x( j ) )
               uscal = tscal
               rec = one / max( xmax, one )
               IF( cnorm( j ).GT.( bignum-xj )*rec ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  rec = rec*half
                  IF( nounit ) THEN
                     tjjs = a( j, j )*tscal
                  ELSE
                     tjjs = tscal
                  END IF
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.one ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     rec = min( one, rec*tjj )
                     uscal = zladiv( uscal, tjjs )
                  END IF
                  IF( rec.LT.one ) THEN
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  END IF
               END IF
*
               csumj = zero
               IF( uscal.EQ.dcmplx( one ) ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call ZDOTU to perform the dot product.
*
                  IF( upper ) THEN
                     csumj = zdotu( j-1, a( 1, j ), 1, x, 1 )
                  ELSE IF( j.LT.n ) THEN
                     csumj = zdotu( n-j, a( j+1, j ), 1, x( j+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( upper ) THEN
                     DO 130 i = 1, j - 1
                        csumj = csumj + ( a( i, j )*uscal )*x( i )
  130                CONTINUE
                  ELSE IF( j.LT.n ) THEN
                     DO 140 i = j + 1, n
                        csumj = csumj + ( a( i, j )*uscal )*x( i )
  140                CONTINUE
                  END IF
               END IF
*
               IF( uscal.EQ.dcmplx( tscal ) ) THEN
*
*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  x( j ) = x( j ) - csumj
                  xj = cabs1( x( j ) )
                  IF( nounit ) THEN
                     tjjs = a( j, j )*tscal
                  ELSE
                     tjjs = tscal
                     IF( tscal.EQ.one )
     $                  GO TO 160
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.smlnum ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( tjj.LT.one ) THEN
                        IF( xj.GT.tjj*bignum ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           rec = one / xj
                           CALL zdscal( n, rec, x, 1 )
                           scale = scale*rec
                           xmax = xmax*rec
                        END IF
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE IF( tjj.GT.zero ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( xj.GT.tjj*bignum ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        rec = ( tjj*bignum ) / xj
                        CALL zdscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0 and compute a solution to A**T *x = 0.
*
                     DO 150 i = 1, n
                        x( i ) = zero
  150                CONTINUE
                     x( j ) = one
                     scale = zero
                     xmax = zero
                  END IF
  160             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
*                 product has already been divided by 1/A(j,j).
*
                  x( j ) = zladiv( x( j ), tjjs ) - csumj
               END IF
               xmax = max( xmax, cabs1( x( j ) ) )
  170       CONTINUE
*
         ELSE
*
*           Solve A**H * x = b
*
            DO 220 j = jfirst, jlast, jinc
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               xj = cabs1( x( j ) )
               uscal = tscal
               rec = one / max( xmax, one )
               IF( cnorm( j ).GT.( bignum-xj )*rec ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  rec = rec*half
                  IF( nounit ) THEN
                     tjjs = dconjg( a( j, j ) )*tscal
                  ELSE
                     tjjs = tscal
                  END IF
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.one ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     rec = min( one, rec*tjj )
                     uscal = zladiv( uscal, tjjs )
                  END IF
                  IF( rec.LT.one ) THEN
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  END IF
               END IF
*
               csumj = zero
               IF( uscal.EQ.dcmplx( one ) ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call ZDOTC to perform the dot product.
*
                  IF( upper ) THEN
                     csumj = zdotc( j-1, a( 1, j ), 1, x, 1 )
                  ELSE IF( j.LT.n ) THEN
                     csumj = zdotc( n-j, a( j+1, j ), 1, x( j+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( upper ) THEN
                     DO 180 i = 1, j - 1
                        csumj = csumj + ( dconjg( a( i, j ) )*uscal )*
     $                          x( i )
  180                CONTINUE
                  ELSE IF( j.LT.n ) THEN
                     DO 190 i = j + 1, n
                        csumj = csumj + ( dconjg( a( i, j ) )*uscal )*
     $                          x( i )
  190                CONTINUE
                  END IF
               END IF
*
               IF( uscal.EQ.dcmplx( tscal ) ) THEN
*
*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  x( j ) = x( j ) - csumj
                  xj = cabs1( x( j ) )
                  IF( nounit ) THEN
                     tjjs = dconjg( a( j, j ) )*tscal
                  ELSE
                     tjjs = tscal
                     IF( tscal.EQ.one )
     $                  GO TO 210
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.smlnum ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( tjj.LT.one ) THEN
                        IF( xj.GT.tjj*bignum ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           rec = one / xj
                           CALL zdscal( n, rec, x, 1 )
                           scale = scale*rec
                           xmax = xmax*rec
                        END IF
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE IF( tjj.GT.zero ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( xj.GT.tjj*bignum ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        rec = ( tjj*bignum ) / xj
                        CALL zdscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0 and compute a solution to A**H *x = 0.
*
                     DO 200 i = 1, n
                        x( i ) = zero
  200                CONTINUE
                     x( j ) = one
                     scale = zero
                     xmax = zero
                  END IF
  210             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
*                 product has already been divided by 1/A(j,j).
*
                  x( j ) = zladiv( x( j ), tjjs ) - csumj
               END IF
               xmax = max( xmax, cabs1( x( j ) ) )
  220       CONTINUE
         END IF
         scale = scale / tscal
      END IF
*
*     Scale the column norms by 1/TSCAL for return.
*
      IF( tscal.NE.one ) THEN
         CALL dscal( n, one / tscal, cnorm, 1 )
      END IF
*
      RETURN
*
*     End of ZLATRS
*
      END
      SUBROUTINE dlaln2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B,
     $                   LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            LTRANS
      INTEGER            INFO, LDA, LDB, LDX, NA, NW
      DOUBLE PRECISION   CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
      DOUBLE PRECISION   TWO
      parameter( two = 2.0d0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICMAX, J
      DOUBLE PRECISION   BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21,
     $                   ci22, cmax, cnorm, cr21, cr22, csi, csr, li21,
     $                   lr21, smini, smlnum, temp, u22abs, ui11, ui11r,
     $                   ui12, ui12s, ui22, ur11, ur11r, ur12, ur12s,
     $                   ur22, xi1, xi2, xr1, xr2
*     ..
*     .. Local Arrays ..
      LOGICAL            RSWAP( 4 ), ZSWAP( 4 )
      INTEGER            IPIVOT( 4, 4 )
      DOUBLE PRECISION   CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           dladiv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max
*     ..
*     .. Equivalences ..
      equivalence( ci( 1, 1 ), civ( 1 ) ),
     $                   ( cr( 1, 1 ), crv( 1 ) )
*     ..
*     .. Data statements ..
      DATA               zswap / .false., .false., .true., .true. /
      DATA               rswap / .false., .true., .false., .true. /
      DATA               ipivot / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4,
     $                   3, 2, 1 /
*     ..
*     .. Executable Statements ..
*
*     Compute BIGNUM
*
      smlnum = two*dlamch( 'Safe minimum' )
      bignum = one / smlnum
      smini = max( smin, smlnum )
*
*     Don't check for input errors
*
      info = 0
*
*     Standard Initializations
*
      scale = one
*
      IF( na.EQ.1 ) THEN
*
*        1 x 1  (i.e., scalar) system   C X = B
*
         IF( nw.EQ.1 ) THEN
*
*           Real 1x1 system.
*
*           C = ca A - w D
*
            csr = ca*a( 1, 1 ) - wr*d1
            cnorm = abs( csr )
*
*           If | C | < SMINI, use C = SMINI
*
            IF( cnorm.LT.smini ) THEN
               csr = smini
               cnorm = smini
               info = 1
            END IF
*
*           Check scaling for  X = B / C
*
            bnorm = abs( b( 1, 1 ) )
            IF( cnorm.LT.one .AND. bnorm.GT.one ) THEN
               IF( bnorm.GT.bignum*cnorm )
     $            scale = one / bnorm
            END IF
*
*           Compute X
*
            x( 1, 1 ) = ( b( 1, 1 )*scale ) / csr
            xnorm = abs( x( 1, 1 ) )
         ELSE
*
*           Complex 1x1 system (w is complex)
*
*           C = ca A - w D
*
            csr = ca*a( 1, 1 ) - wr*d1
            csi = -wi*d1
            cnorm = abs( csr ) + abs( csi )
*
*           If | C | < SMINI, use C = SMINI
*
            IF( cnorm.LT.smini ) THEN
               csr = smini
               csi = zero
               cnorm = smini
               info = 1
            END IF
*
*           Check scaling for  X = B / C
*
            bnorm = abs( b( 1, 1 ) ) + abs( b( 1, 2 ) )
            IF( cnorm.LT.one .AND. bnorm.GT.one ) THEN
               IF( bnorm.GT.bignum*cnorm )
     $            scale = one / bnorm
            END IF
*
*           Compute X
*
            CALL dladiv( scale*b( 1, 1 ), scale*b( 1, 2 ), csr, csi,
     $                   x( 1, 1 ), x( 1, 2 ) )
            xnorm = abs( x( 1, 1 ) ) + abs( x( 1, 2 ) )
         END IF
*
      ELSE
*
*        2x2 System
*
*        Compute the real part of  C = ca A - w D  (or  ca A**T - w D )
*
         cr( 1, 1 ) = ca*a( 1, 1 ) - wr*d1
         cr( 2, 2 ) = ca*a( 2, 2 ) - wr*d2
         IF( ltrans ) THEN
            cr( 1, 2 ) = ca*a( 2, 1 )
            cr( 2, 1 ) = ca*a( 1, 2 )
         ELSE
            cr( 2, 1 ) = ca*a( 2, 1 )
            cr( 1, 2 ) = ca*a( 1, 2 )
         END IF
*
         IF( nw.EQ.1 ) THEN
*
*           Real 2x2 system  (w is real)
*
*           Find the largest element in C
*
            cmax = zero
            icmax = 0
*
            DO 10 j = 1, 4
               IF( abs( crv( j ) ).GT.cmax ) THEN
                  cmax = abs( crv( j ) )
                  icmax = j
               END IF
   10       CONTINUE
*
*           If norm(C) < SMINI, use SMINI*identity.
*
            IF( cmax.LT.smini ) THEN
               bnorm = max( abs( b( 1, 1 ) ), abs( b( 2, 1 ) ) )
               IF( smini.LT.one .AND. bnorm.GT.one ) THEN
                  IF( bnorm.GT.bignum*smini )
     $               scale = one / bnorm
               END IF
               temp = scale / smini
               x( 1, 1 ) = temp*b( 1, 1 )
               x( 2, 1 ) = temp*b( 2, 1 )
               xnorm = temp*bnorm
               info = 1
               RETURN
            END IF
*
*           Gaussian elimination with complete pivoting.
*
            ur11 = crv( icmax )
            cr21 = crv( ipivot( 2, icmax ) )
            ur12 = crv( ipivot( 3, icmax ) )
            cr22 = crv( ipivot( 4, icmax ) )
            ur11r = one / ur11
            lr21 = ur11r*cr21
            ur22 = cr22 - ur12*lr21
*
*           If smaller pivot < SMINI, use SMINI
*
            IF( abs( ur22 ).LT.smini ) THEN
               ur22 = smini
               info = 1
            END IF
            IF( rswap( icmax ) ) THEN
               br1 = b( 2, 1 )
               br2 = b( 1, 1 )
            ELSE
               br1 = b( 1, 1 )
               br2 = b( 2, 1 )
            END IF
            br2 = br2 - lr21*br1
            bbnd = max( abs( br1*( ur22*ur11r ) ), abs( br2 ) )
            IF( bbnd.GT.one .AND. abs( ur22 ).LT.one ) THEN
               IF( bbnd.GE.bignum*abs( ur22 ) )
     $            scale = one / bbnd
            END IF
*
            xr2 = ( br2*scale ) / ur22
            xr1 = ( scale*br1 )*ur11r - xr2*( ur11r*ur12 )
            IF( zswap( icmax ) ) THEN
               x( 1, 1 ) = xr2
               x( 2, 1 ) = xr1
            ELSE
               x( 1, 1 ) = xr1
               x( 2, 1 ) = xr2
            END IF
            xnorm = max( abs( xr1 ), abs( xr2 ) )
*
*           Further scaling if  norm(A) norm(X) > overflow
*
            IF( xnorm.GT.one .AND. cmax.GT.one ) THEN
               IF( xnorm.GT.bignum / cmax ) THEN
                  temp = cmax / bignum
                  x( 1, 1 ) = temp*x( 1, 1 )
                  x( 2, 1 ) = temp*x( 2, 1 )
                  xnorm = temp*xnorm
                  scale = temp*scale
               END IF
            END IF
         ELSE
*
*           Complex 2x2 system  (w is complex)
*
*           Find the largest element in C
*
            ci( 1, 1 ) = -wi*d1
            ci( 2, 1 ) = zero
            ci( 1, 2 ) = zero
            ci( 2, 2 ) = -wi*d2
            cmax = zero
            icmax = 0
*
            DO 20 j = 1, 4
               IF( abs( crv( j ) )+abs( civ( j ) ).GT.cmax ) THEN
                  cmax = abs( crv( j ) ) + abs( civ( j ) )
                  icmax = j
               END IF
   20       CONTINUE
*
*           If norm(C) < SMINI, use SMINI*identity.
*
            IF( cmax.LT.smini ) THEN
               bnorm = max( abs( b( 1, 1 ) )+abs( b( 1, 2 ) ),
     $                 abs( b( 2, 1 ) )+abs( b( 2, 2 ) ) )
               IF( smini.LT.one .AND. bnorm.GT.one ) THEN
                  IF( bnorm.GT.bignum*smini )
     $               scale = one / bnorm
               END IF
               temp = scale / smini
               x( 1, 1 ) = temp*b( 1, 1 )
               x( 2, 1 ) = temp*b( 2, 1 )
               x( 1, 2 ) = temp*b( 1, 2 )
               x( 2, 2 ) = temp*b( 2, 2 )
               xnorm = temp*bnorm
               info = 1
               RETURN
            END IF
*
*           Gaussian elimination with complete pivoting.
*
            ur11 = crv( icmax )
            ui11 = civ( icmax )
            cr21 = crv( ipivot( 2, icmax ) )
            ci21 = civ( ipivot( 2, icmax ) )
            ur12 = crv( ipivot( 3, icmax ) )
            ui12 = civ( ipivot( 3, icmax ) )
            cr22 = crv( ipivot( 4, icmax ) )
            ci22 = civ( ipivot( 4, icmax ) )
            IF( icmax.EQ.1 .OR. icmax.EQ.4 ) THEN
*
*              Code when off-diagonals of pivoted C are real
*
               IF( abs( ur11 ).GT.abs( ui11 ) ) THEN
                  temp = ui11 / ur11
                  ur11r = one / ( ur11*( one+temp**2 ) )
                  ui11r = -temp*ur11r
               ELSE
                  temp = ur11 / ui11
                  ui11r = -one / ( ui11*( one+temp**2 ) )
                  ur11r = -temp*ui11r
               END IF
               lr21 = cr21*ur11r
               li21 = cr21*ui11r
               ur12s = ur12*ur11r
               ui12s = ur12*ui11r
               ur22 = cr22 - ur12*lr21
               ui22 = ci22 - ur12*li21
            ELSE
*
*              Code when diagonals of pivoted C are real
*
               ur11r = one / ur11
               ui11r = zero
               lr21 = cr21*ur11r
               li21 = ci21*ur11r
               ur12s = ur12*ur11r
               ui12s = ui12*ur11r
               ur22 = cr22 - ur12*lr21 + ui12*li21
               ui22 = -ur12*li21 - ui12*lr21
            END IF
            u22abs = abs( ur22 ) + abs( ui22 )
*
*           If smaller pivot < SMINI, use SMINI
*
            IF( u22abs.LT.smini ) THEN
               ur22 = smini
               ui22 = zero
               info = 1
            END IF
            IF( rswap( icmax ) ) THEN
               br2 = b( 1, 1 )
               br1 = b( 2, 1 )
               bi2 = b( 1, 2 )
               bi1 = b( 2, 2 )
            ELSE
               br1 = b( 1, 1 )
               br2 = b( 2, 1 )
               bi1 = b( 1, 2 )
               bi2 = b( 2, 2 )
            END IF
            br2 = br2 - lr21*br1 + li21*bi1
            bi2 = bi2 - li21*br1 - lr21*bi1
            bbnd = max( ( abs( br1 )+abs( bi1 ) )*
     $             ( u22abs*( abs( ur11r )+abs( ui11r ) ) ),
     $             abs( br2 )+abs( bi2 ) )
            IF( bbnd.GT.one .AND. u22abs.LT.one ) THEN
               IF( bbnd.GE.bignum*u22abs ) THEN
                  scale = one / bbnd
                  br1 = scale*br1
                  bi1 = scale*bi1
                  br2 = scale*br2
                  bi2 = scale*bi2
               END IF
            END IF
*
            CALL dladiv( br2, bi2, ur22, ui22, xr2, xi2 )
            xr1 = ur11r*br1 - ui11r*bi1 - ur12s*xr2 + ui12s*xi2
            xi1 = ui11r*br1 + ur11r*bi1 - ui12s*xr2 - ur12s*xi2
            IF( zswap( icmax ) ) THEN
               x( 1, 1 ) = xr2
               x( 2, 1 ) = xr1
               x( 1, 2 ) = xi2
               x( 2, 2 ) = xi1
            ELSE
               x( 1, 1 ) = xr1
               x( 2, 1 ) = xr2
               x( 1, 2 ) = xi1
               x( 2, 2 ) = xi2
            END IF
            xnorm = max( abs( xr1 )+abs( xi1 ), abs( xr2 )+abs( xi2 ) )
*
*           Further scaling if  norm(A) norm(X) > overflow
*
            IF( xnorm.GT.one .AND. cmax.GT.one ) THEN
               IF( xnorm.GT.bignum / cmax ) THEN
                  temp = cmax / bignum
                  x( 1, 1 ) = temp*x( 1, 1 )
                  x( 2, 1 ) = temp*x( 2, 1 )
                  x( 1, 2 ) = temp*x( 1, 2 )
                  x( 2, 2 ) = temp*x( 2, 2 )
                  xnorm = temp*xnorm
                  scale = temp*scale
               END IF
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DLALN2
*
      END
      SUBROUTINE dlahqr( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, INFO )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
*     ..
*
*  =========================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      parameter( zero = 0.0d0, one = 1.0d0, two = 2.0d0 )
      DOUBLE PRECISION   DAT1, DAT2
      parameter( dat1 = 3.0d0 / 4.0d0, dat2 = -0.4375d0 )
      INTEGER            KEXSH
      parameter( kexsh = 10 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, AB, BA, BB, CS, DET, H11, H12, H21, H21S,
     $                   h22, rt1i, rt1r, rt2i, rt2r, rtdisc, s, safmax,
     $                   safmin, smlnum, sn, sum, t1, t2, t3, tr, tst,
     $                   ulp, v2, v3
      INTEGER            I, I1, I2, ITS, ITMAX, J, K, L, M, NH, NR, NZ,
     $                   kdefl 
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           dcopy, dlabad, dlanv2, dlarfg, drot
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, max, min, sqrt
*     ..
*     .. Executable Statements ..
*
      info = 0
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
      IF( ilo.EQ.ihi ) THEN
         wr( ilo ) = h( ilo, ilo )
         wi( ilo ) = zero
         RETURN
      END IF
*
*     ==== clear out the trash ====
      DO 10 j = ilo, ihi - 3
         h( j+2, j ) = zero
         h( j+3, j ) = zero
   10 CONTINUE
      IF( ilo.LE.ihi-2 )
     $   h( ihi, ihi-2 ) = zero
*
      nh = ihi - ilo + 1
      nz = ihiz - iloz + 1
*
*     Set machine-dependent constants for the stopping criterion.
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = one / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( nh ) / ulp )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( wantt ) THEN
         i1 = 1
         i2 = n
      END IF
*
*     ITMAX is the total number of QR iterations allowed.
*
      itmax = 30 * max( 10, nh )
*
*     KDEFL counts the number of iterations since a deflation
*
      kdefl = 0
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
*     H(L,L-1) is negligible so that the matrix splits.
*
      i = ihi
   20 CONTINUE
      l = ilo
      IF( i.LT.ilo )
     $   GO TO 160
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 or 2 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      DO 140 its = 0, itmax
*
*        Look for a single small subdiagonal element.
*
         DO 30 k = i, l + 1, -1
            IF( abs( h( k, k-1 ) ).LE.smlnum )
     $         GO TO 40
            tst = abs( h( k-1, k-1 ) ) + abs( h( k, k ) )
            IF( tst.EQ.zero ) THEN
               IF( k-2.GE.ilo )
     $            tst = tst + abs( h( k-1, k-2 ) )
               IF( k+1.LE.ihi )
     $            tst = tst + abs( h( k+1, k ) )
            END IF
*           ==== The following is a conservative small subdiagonal
*           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
*           .    1997). It has better mathematical foundation and
*           .    improves accuracy in some cases.  ====
            IF( abs( h( k, k-1 ) ).LE.ulp*tst ) THEN
               ab = max( abs( h( k, k-1 ) ), abs( h( k-1, k ) ) )
               ba = min( abs( h( k, k-1 ) ), abs( h( k-1, k ) ) )
               aa = max( abs( h( k, k ) ),
     $              abs( h( k-1, k-1 )-h( k, k ) ) )
               bb = min( abs( h( k, k ) ),
     $              abs( h( k-1, k-1 )-h( k, k ) ) )
               s = aa + ab
               IF( ba*( ab / s ).LE.max( smlnum,
     $             ulp*( bb*( aa / s ) ) ) )GO TO 40
            END IF
   30    CONTINUE
   40    CONTINUE
         l = k
         IF( l.GT.ilo ) THEN
*
*           H(L,L-1) is negligible
*
            h( l, l-1 ) = zero
         END IF
*
*        Exit from loop if a submatrix of order 1 or 2 has split off.
*
         IF( l.GE.i-1 )
     $      GO TO 150
         kdefl = kdefl + 1
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.wantt ) THEN
            i1 = l
            i2 = i
         END IF
*
         IF( mod(kdefl,2*kexsh).EQ.0 ) THEN
*
*           Exceptional shift.
*
            s = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
            h11 = dat1*s + h( i, i )
            h12 = dat2*s
            h21 = s
            h22 = h11
         ELSE IF( mod(kdefl,kexsh).EQ.0 ) THEN
*
*           Exceptional shift.
*
            s = abs( h( l+1, l ) ) + abs( h( l+2, l+1 ) )
            h11 = dat1*s + h( l, l )
            h12 = dat2*s
            h21 = s
            h22 = h11
         ELSE
*
*           Prepare to use Francis' double shift
*           (i.e. 2nd degree generalized Rayleigh quotient)
*
            h11 = h( i-1, i-1 )
            h21 = h( i, i-1 )
            h12 = h( i-1, i )
            h22 = h( i, i )
         END IF
         s = abs( h11 ) + abs( h12 ) + abs( h21 ) + abs( h22 )
         IF( s.EQ.zero ) THEN
            rt1r = zero
            rt1i = zero
            rt2r = zero
            rt2i = zero
         ELSE
            h11 = h11 / s
            h21 = h21 / s
            h12 = h12 / s
            h22 = h22 / s
            tr = ( h11+h22 ) / two
            det = ( h11-tr )*( h22-tr ) - h12*h21
            rtdisc = sqrt( abs( det ) )
            IF( det.GE.zero ) THEN
*
*              ==== complex conjugate shifts ====
*
               rt1r = tr*s
               rt2r = rt1r
               rt1i = rtdisc*s
               rt2i = -rt1i
            ELSE
*
*              ==== real shifts (use only one of them)  ====
*
               rt1r = tr + rtdisc
               rt2r = tr - rtdisc
               IF( abs( rt1r-h22 ).LE.abs( rt2r-h22 ) ) THEN
                  rt1r = rt1r*s
                  rt2r = rt1r
               ELSE
                  rt2r = rt2r*s
                  rt1r = rt2r
               END IF
               rt1i = zero
               rt2i = zero
            END IF
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 50 m = i - 2, l, -1
*           Determine the effect of starting the double-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.  (The following uses scaling to avoid
*           overflows and most underflows.)
*
            h21s = h( m+1, m )
            s = abs( h( m, m )-rt2r ) + abs( rt2i ) + abs( h21s )
            h21s = h( m+1, m ) / s
            v( 1 ) = h21s*h( m, m+1 ) + ( h( m, m )-rt1r )*
     $               ( ( h( m, m )-rt2r ) / s ) - rt1i*( rt2i / s )
            v( 2 ) = h21s*( h( m, m )+h( m+1, m+1 )-rt1r-rt2r )
            v( 3 ) = h21s*h( m+2, m+1 )
            s = abs( v( 1 ) ) + abs( v( 2 ) ) + abs( v( 3 ) )
            v( 1 ) = v( 1 ) / s
            v( 2 ) = v( 2 ) / s
            v( 3 ) = v( 3 ) / s
            IF( m.EQ.l )
     $         GO TO 60
            IF( abs( h( m, m-1 ) )*( abs( v( 2 ) )+abs( v( 3 ) ) ).LE.
     $          ulp*abs( v( 1 ) )*( abs( h( m-1, m-1 ) )+abs( h( m,
     $          m ) )+abs( h( m+1, m+1 ) ) ) )GO TO 60
   50    CONTINUE
   60    CONTINUE
*
*        Double-shift QR step
*
         DO 130 k = m, i - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix. NR is the order of G.
*
            nr = min( 3, i-k+1 )
            IF( k.GT.m )
     $         CALL dcopy( nr, h( k, k-1 ), 1, v, 1 )
            CALL dlarfg( nr, v( 1 ), v( 2 ), 1, t1 )
            IF( k.GT.m ) THEN
               h( k, k-1 ) = v( 1 )
               h( k+1, k-1 ) = zero
               IF( k.LT.i-1 )
     $            h( k+2, k-1 ) = zero
            ELSE IF( m.GT.l ) THEN
*               ==== Use the following instead of
*               .    H( K, K-1 ) = -H( K, K-1 ) to
*               .    avoid a bug when v(2) and v(3)
*               .    underflow. ====
               h( k, k-1 ) = h( k, k-1 )*( one-t1 )
            END IF
            v2 = v( 2 )
            t2 = t1*v2
            IF( nr.EQ.3 ) THEN
               v3 = v( 3 )
               t3 = t1*v3
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 70 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j ) + v3*h( k+2, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
                  h( k+2, j ) = h( k+2, j ) - sum*t3
   70          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 80 j = i1, min( k+3, i )
                  sum = h( j, k ) + v2*h( j, k+1 ) + v3*h( j, k+2 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
                  h( j, k+2 ) = h( j, k+2 ) - sum*t3
   80          CONTINUE
*
               IF( wantz ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 90 j = iloz, ihiz
                     sum = z( j, k ) + v2*z( j, k+1 ) + v3*z( j, k+2 )
                     z( j, k ) = z( j, k ) - sum*t1
                     z( j, k+1 ) = z( j, k+1 ) - sum*t2
                     z( j, k+2 ) = z( j, k+2 ) - sum*t3
   90             CONTINUE
               END IF
            ELSE IF( nr.EQ.2 ) THEN
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 100 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
  100          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 110 j = i1, i
                  sum = h( j, k ) + v2*h( j, k+1 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
  110          CONTINUE
*
               IF( wantz ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 120 j = iloz, ihiz
                     sum = z( j, k ) + v2*z( j, k+1 )
                     z( j, k ) = z( j, k ) - sum*t1
                     z( j, k+1 ) = z( j, k+1 ) - sum*t2
  120             CONTINUE
               END IF
            END IF
  130    CONTINUE
*
  140 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      info = i
      RETURN
*
  150 CONTINUE
*
      IF( l.EQ.i ) THEN
*
*        H(I,I-1) is negligible: one eigenvalue has converged.
*
         wr( i ) = h( i, i )
         wi( i ) = zero
      ELSE IF( l.EQ.i-1 ) THEN
*
*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
*
*        Transform the 2-by-2 submatrix to standard Schur form,
*        and compute and store the eigenvalues.
*
         CALL dlanv2( h( i-1, i-1 ), h( i-1, i ), h( i, i-1 ),
     $                h( i, i ), wr( i-1 ), wi( i-1 ), wr( i ), wi( i ),
     $                cs, sn )
*
         IF( wantt ) THEN
*
*           Apply the transformation to the rest of H.
*
            IF( i2.GT.i )
     $         CALL drot( i2-i, h( i-1, i+1 ), ldh, h( i, i+1 ), ldh,
     $                    cs, sn )
            CALL drot( i-i1-1, h( i1, i-1 ), 1, h( i1, i ), 1, cs, sn )
         END IF
         IF( wantz ) THEN
*
*           Apply the transformation to Z.
*
            CALL drot( nz, z( iloz, i-1 ), 1, z( iloz, i ), 1, cs, sn )
         END IF
      END IF
*     reset deflation counter
      kdefl = 0
*
*     return to start of the main loop with new value of I.
*
      i = l - 1
      GO TO 20
*
  160 CONTINUE
      RETURN
*
*     End of DLAHQR
*
      END
      SUBROUTINE dlaqr0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
     $                   z( ldz, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    DLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 15 )
*
*     ==== Exceptional deflation windows:  try to cure rare
*     .    slow convergence by varying the size of the
*     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      parameter( kexnw = 5 )
*
*     ==== Exceptional shifts: try to cure rare slow convergence
*     .    with ad-hoc exceptional shifts every KEXSH iterations.
*     .    ====
      INTEGER            KEXSH
      parameter( kexsh = 6 )
*
*     ==== The constants WILK1 and WILK2 are used to form the
*     .    exceptional shifts. ====
      DOUBLE PRECISION   WILK1, WILK2
      parameter( wilk1 = 0.75d0, wilk2 = -0.4375d0 )
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
     $                   kt, ktop, ku, kv, kwh, kwtop, kwv, ld, ls,
     $                   lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns,
     $                   nsmax, nsr, nve, nw, nwmax, nwr, nwupbd
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   ZDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlacpy, dlahqr, dlanv2, dlaqr3, dlaqr4, dlaqr5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, int, max, min, mod
*     ..
*     .. Executable Statements ..
      info = 0
*
*     ==== Quick return for N = 0: nothing to do. ====
*
      IF( n.EQ.0 ) THEN
         work( 1 ) = one
         RETURN
      END IF
*
      IF( n.LE.ntiny ) THEN
*
*        ==== Tiny matrices must use DLAHQR. ====
*
         lwkopt = 1
         IF( lwork.NE.-1 )
     $      CALL dlahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi,
     $                   iloz, ihiz, z, ldz, info )
      ELSE
*
*        ==== Use small bulge multi-shift QR with aggressive early
*        .    deflation on larger-than-tiny matrices. ====
*
*        ==== Hope for the best. ====
*
         info = 0
*
*        ==== Set up job flags for ILAENV. ====
*
         IF( wantt ) THEN
            jbcmpz( 1: 1 ) = 'S'
         ELSE
            jbcmpz( 1: 1 ) = 'E'
         END IF
         IF( wantz ) THEN
            jbcmpz( 2: 2 ) = 'V'
         ELSE
            jbcmpz( 2: 2 ) = 'N'
         END IF
*
*        ==== NWR = recommended deflation window size.  At this
*        .    point,  N .GT. NTINY = 15, so there is enough
*        .    subdiagonal workspace for NWR.GE.2 as required.
*        .    (In fact, there is enough subdiagonal space for
*        .    NWR.GE.4.) ====
*
         nwr = ilaenv( 13, 'DLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nwr = max( 2, nwr )
         nwr = min( ihi-ilo+1, ( n-1 ) / 3, nwr )
*
*        ==== NSR = recommended number of simultaneous shifts.
*        .    At this point N .GT. NTINY = 15, so there is at
*        .    enough subdiagonal workspace for NSR to be even
*        .    and greater than or equal to two as required. ====
*
         nsr = ilaenv( 15, 'DLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nsr = min( nsr, ( n-3 ) / 6, ihi-ilo )
         nsr = max( 2, nsr-mod( nsr, 2 ) )
*
*        ==== Estimate optimal workspace ====
*
*        ==== Workspace query call to DLAQR3 ====
*
         CALL dlaqr3( wantt, wantz, n, ilo, ihi, nwr+1, h, ldh, iloz,
     $                ihiz, z, ldz, ls, ld, wr, wi, h, ldh, n, h, ldh,
     $                n, h, ldh, work, -1 )
*
*        ==== Optimal workspace = MAX(DLAQR5, DLAQR3) ====
*
         lwkopt = max( 3*nsr / 2, int( work( 1 ) ) )
*
*        ==== Quick return in case of workspace query. ====
*
         IF( lwork.EQ.-1 ) THEN
            work( 1 ) = dble( lwkopt )
            RETURN
         END IF
*
*        ==== DLAHQR/DLAQR0 crossover point ====
*
         nmin = ilaenv( 12, 'DLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
*
*        ==== Nibble crossover point ====
*
         nibble = ilaenv( 14, 'DLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nibble = max( 0, nibble )
*
*        ==== Accumulate reflections during ttswp?  Use block
*        .    2-by-2 structure during matrix-matrix multiply? ====
*
         kacc22 = ilaenv( 16, 'DLAQR0', jbcmpz, n, ilo, ihi, lwork )
         kacc22 = max( 0, kacc22 )
         kacc22 = min( 2, kacc22 )
*
*        ==== NWMAX = the largest possible deflation window for
*        .    which there is sufficient workspace. ====
*
         nwmax = min( ( n-1 ) / 3, lwork / 2 )
         nw = nwmax
*
*        ==== NSMAX = the Largest number of simultaneous shifts
*        .    for which there is sufficient workspace. ====
*
         nsmax = min( ( n-3 ) / 6, 2*lwork / 3 )
         nsmax = nsmax - mod( nsmax, 2 )
*
*        ==== NDFL: an iteration count restarted at deflation. ====
*
         ndfl = 1
*
*        ==== ITMAX = iteration limit ====
*
         itmax = max( 30, 2*kexsh )*max( 10, ( ihi-ilo+1 ) )
*
*        ==== Last row and column in the active block ====
*
         kbot = ihi
*
*        ==== Main Loop ====
*
         DO 80 it = 1, itmax
*
*           ==== Done when KBOT falls below ILO ====
*
            IF( kbot.LT.ilo )
     $         GO TO 90
*
*           ==== Locate active block ====
*
            DO 10 k = kbot, ilo + 1, -1
               IF( h( k, k-1 ).EQ.zero )
     $            GO TO 20
   10       CONTINUE
            k = ilo
   20       CONTINUE
            ktop = k
*
*           ==== Select deflation window size:
*           .    Typical Case:
*           .      If possible and advisable, nibble the entire
*           .      active block.  If not, use size MIN(NWR,NWMAX)
*           .      or MIN(NWR+1,NWMAX) depending upon which has
*           .      the smaller corresponding subdiagonal entry
*           .      (a heuristic).
*           .
*           .    Exceptional Case:
*           .      If there have been no deflations in KEXNW or
*           .      more iterations, then vary the deflation window
*           .      size.   At first, because, larger windows are,
*           .      in general, more powerful than smaller ones,
*           .      rapidly increase the window to the maximum possible.
*           .      Then, gradually reduce the window size. ====
*
            nh = kbot - ktop + 1
            nwupbd = min( nh, nwmax )
            IF( ndfl.LT.kexnw ) THEN
               nw = min( nwupbd, nwr )
            ELSE
               nw = min( nwupbd, 2*nw )
            END IF
            IF( nw.LT.nwmax ) THEN
               IF( nw.GE.nh-1 ) THEN
                  nw = nh
               ELSE
                  kwtop = kbot - nw + 1
                  IF( abs( h( kwtop, kwtop-1 ) ).GT.
     $                abs( h( kwtop-1, kwtop-2 ) ) )nw = nw + 1
               END IF
            END IF
            IF( ndfl.LT.kexnw ) THEN
               ndec = -1
            ELSE IF( ndec.GE.0 .OR. nw.GE.nwupbd ) THEN
               ndec = ndec + 1
               IF( nw-ndec.LT.2 )
     $            ndec = 0
               nw = nw - ndec
            END IF
*
*           ==== Aggressive early deflation:
*           .    split workspace under the subdiagonal into
*           .      - an nw-by-nw work array V in the lower
*           .        left-hand-corner,
*           .      - an NW-by-at-least-NW-but-more-is-better
*           .        (NW-by-NHO) horizontal work array along
*           .        the bottom edge,
*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
*           .        vertical work array along the left-hand-edge.
*           .        ====
*
            kv = n - nw + 1
            kt = nw + 1
            nho = ( n-nw-1 ) - kt + 1
            kwv = nw + 2
            nve = ( n-nw ) - kwv + 1
*
*           ==== Aggressive early deflation ====
*
            CALL dlaqr3( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz,
     $                   ihiz, z, ldz, ls, ld, wr, wi, h( kv, 1 ), ldh,
     $                   nho, h( kv, kt ), ldh, nve, h( kwv, 1 ), ldh,
     $                   work, lwork )
*
*           ==== Adjust KBOT accounting for new deflations. ====
*
            kbot = kbot - ld
*
*           ==== KS points to the shifts. ====
*
            ks = kbot - ls + 1
*
*           ==== Skip an expensive QR sweep if there is a (partly
*           .    heuristic) reason to expect that many eigenvalues
*           .    will deflate without it.  Here, the QR sweep is
*           .    skipped if many eigenvalues have just been deflated
*           .    or if the remaining active block is small.
*
            IF( ( ld.EQ.0 ) .OR. ( ( 100*ld.LE.nw*nibble ) .AND. ( kbot-
     $          ktop+1.GT.min( nmin, nwmax ) ) ) ) THEN
*
*              ==== NS = nominal number of simultaneous shifts.
*              .    This may be lowered (slightly) if DLAQR3
*              .    did not provide that many shifts. ====
*
               ns = min( nsmax, nsr, max( 2, kbot-ktop ) )
               ns = ns - mod( ns, 2 )
*
*              ==== If there have been no deflations
*              .    in a multiple of KEXSH iterations,
*              .    then try exceptional shifts.
*              .    Otherwise use shifts provided by
*              .    DLAQR3 above or from the eigenvalues
*              .    of a trailing principal submatrix. ====
*
               IF( mod( ndfl, kexsh ).EQ.0 ) THEN
                  ks = kbot - ns + 1
                  DO 30 i = kbot, max( ks+1, ktop+2 ), -2
                     ss = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
                     aa = wilk1*ss + h( i, i )
                     bb = ss
                     cc = wilk2*ss
                     dd = aa
                     CALL dlanv2( aa, bb, cc, dd, wr( i-1 ), wi( i-1 ),
     $                            wr( i ), wi( i ), cs, sn )
   30             CONTINUE
                  IF( ks.EQ.ktop ) THEN
                     wr( ks+1 ) = h( ks+1, ks+1 )
                     wi( ks+1 ) = zero
                     wr( ks ) = wr( ks+1 )
                     wi( ks ) = wi( ks+1 )
                  END IF
               ELSE
*
*                 ==== Got NS/2 or fewer shifts? Use DLAQR4 or
*                 .    DLAHQR on a trailing principal submatrix to
*                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
*                 .    there is enough space below the subdiagonal
*                 .    to fit an NS-by-NS scratch array.) ====
*
                  IF( kbot-ks+1.LE.ns / 2 ) THEN
                     ks = kbot - ns + 1
                     kt = n - ns + 1
                     CALL dlacpy( 'A', ns, ns, h( ks, ks ), ldh,
     $                            h( kt, 1 ), ldh )
                     IF( ns.GT.nmin ) THEN
                        CALL dlaqr4( .false., .false., ns, 1, ns,
     $                               h( kt, 1 ), ldh, wr( ks ),
     $                               wi( ks ), 1, 1, zdum, 1, work,
     $                               lwork, inf )
                     ELSE
                        CALL dlahqr( .false., .false., ns, 1, ns,
     $                               h( kt, 1 ), ldh, wr( ks ),
     $                               wi( ks ), 1, 1, zdum, 1, inf )
                     END IF
                     ks = ks + inf
*
*                    ==== In case of a rare QR failure use
*                    .    eigenvalues of the trailing 2-by-2
*                    .    principal submatrix.  ====
*
                     IF( ks.GE.kbot ) THEN
                        aa = h( kbot-1, kbot-1 )
                        cc = h( kbot, kbot-1 )
                        bb = h( kbot-1, kbot )
                        dd = h( kbot, kbot )
                        CALL dlanv2( aa, bb, cc, dd, wr( kbot-1 ),
     $                               wi( kbot-1 ), wr( kbot ),
     $                               wi( kbot ), cs, sn )
                        ks = kbot - 1
                     END IF
                  END IF
*
                  IF( kbot-ks+1.GT.ns ) THEN
*
*                    ==== Sort the shifts (Helps a little)
*                    .    Bubble sort keeps complex conjugate
*                    .    pairs together. ====
*
                     sorted = .false.
                     DO 50 k = kbot, ks + 1, -1
                        IF( sorted )
     $                     GO TO 60
                        sorted = .true.
                        DO 40 i = ks, k - 1
                           IF( abs( wr( i ) )+abs( wi( i ) ).LT.
     $                         abs( wr( i+1 ) )+abs( wi( i+1 ) ) ) THEN
                              sorted = .false.
*
                              swap = wr( i )
                              wr( i ) = wr( i+1 )
                              wr( i+1 ) = swap
*
                              swap = wi( i )
                              wi( i ) = wi( i+1 )
                              wi( i+1 ) = swap
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
*
*                 ==== Shuffle shifts into pairs of real shifts
*                 .    and pairs of complex conjugate shifts
*                 .    assuming complex conjugate shifts are
*                 .    already adjacent to one another. (Yes,
*                 .    they are.)  ====
*
                  DO 70 i = kbot, ks + 2, -2
                     IF( wi( i ).NE.-wi( i-1 ) ) THEN
*
                        swap = wr( i )
                        wr( i ) = wr( i-1 )
                        wr( i-1 ) = wr( i-2 )
                        wr( i-2 ) = swap
*
                        swap = wi( i )
                        wi( i ) = wi( i-1 )
                        wi( i-1 ) = wi( i-2 )
                        wi( i-2 ) = swap
                     END IF
   70             CONTINUE
               END IF
*
*              ==== If there are only two shifts and both are
*              .    real, then use only one.  ====
*
               IF( kbot-ks+1.EQ.2 ) THEN
                  IF( wi( kbot ).EQ.zero ) THEN
                     IF( abs( wr( kbot )-h( kbot, kbot ) ).LT.
     $                   abs( wr( kbot-1 )-h( kbot, kbot ) ) ) THEN
                        wr( kbot-1 ) = wr( kbot )
                     ELSE
                        wr( kbot ) = wr( kbot-1 )
                     END IF
                  END IF
               END IF
*
*              ==== Use up to NS of the the smallest magnitude
*              .    shifts.  If there aren't NS shifts available,
*              .    then use them all, possibly dropping one to
*              .    make the number of shifts even. ====
*
               ns = min( ns, kbot-ks+1 )
               ns = ns - mod( ns, 2 )
               ks = kbot - ns + 1
*
*              ==== Small-bulge multi-shift QR sweep:
*              .    split workspace under the subdiagonal into
*              .    - a KDU-by-KDU work array U in the lower
*              .      left-hand-corner,
*              .    - a KDU-by-at-least-KDU-but-more-is-better
*              .      (KDU-by-NHo) horizontal work array WH along
*              .      the bottom edge,
*              .    - and an at-least-KDU-but-more-is-better-by-KDU
*              .      (NVE-by-KDU) vertical work WV arrow along
*              .      the left-hand-edge. ====
*
               kdu = 2*ns
               ku = n - kdu + 1
               kwh = kdu + 1
               nho = ( n-kdu+1-4 ) - ( kdu+1 ) + 1
               kwv = kdu + 4
               nve = n - kdu - kwv + 1
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               CALL dlaqr5( wantt, wantz, kacc22, n, ktop, kbot, ns,
     $                      wr( ks ), wi( ks ), h, ldh, iloz, ihiz, z,
     $                      ldz, work, 3, h( ku, 1 ), ldh, nve,
     $                      h( kwv, 1 ), ldh, nho, h( ku, kwh ), ldh )
            END IF
*
*           ==== Note progress (or the lack of it). ====
*
            IF( ld.GT.0 ) THEN
               ndfl = 1
            ELSE
               ndfl = ndfl + 1
            END IF
*
*           ==== End of main loop ====
   80    CONTINUE
*
*        ==== Iteration limit exceeded.  Set INFO to show where
*        .    the problem occurred and exit. ====
*
         info = kbot
   90    CONTINUE
      END IF
*
*     ==== Return the optimal value of LWORK. ====
*
      work( 1 ) = dble( lwkopt )
*
*     ==== End of DLAQR0 ====
*
      END
      SUBROUTINE dgehd2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlarf, dlarfg, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEHD2', -info )
         RETURN
      END IF
*
      DO 10 i = ilo, ihi - 1
*
*        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
*
         CALL dlarfg( ihi-i, a( i+1, i ), a( min( i+2, n ), i ), 1,
     $                tau( i ) )
         aii = a( i+1, i )
         a( i+1, i ) = one
*
*        Apply H(i) to A(1:ihi,i+1:ihi) from the right
*
         CALL dlarf( 'Right', ihi, ihi-i, a( i+1, i ), 1, tau( i ),
     $               a( 1, i+1 ), lda, work )
*
*        Apply H(i) to A(i+1:ihi,i+1:n) from the left
*
         CALL dlarf( 'Left', ihi-i, n-i, a( i+1, i ), 1, tau( i ),
     $               a( i+1, i+1 ), lda, work )
*
         a( i+1, i ) = aii
   10 CONTINUE
*
      RETURN
*
*     End of DGEHD2
*
      END
      SUBROUTINE dgebak( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * ), V( LDV, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      DOUBLE PRECISION   S
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dscal, dswap, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters
*
      rightv = lsame( side, 'R' )
      leftv = lsame( side, 'L' )
*
      info = 0
      IF( .NOT.lsame( job, 'N' ) .AND. .NOT.lsame( job, 'P' ) .AND.
     $    .NOT.lsame( job, 'S' ) .AND. .NOT.lsame( job, 'B' ) ) THEN
         info = -1
      ELSE IF( .NOT.rightv .AND. .NOT.leftv ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -5
      ELSE IF( m.LT.0 ) THEN
         info = -7
      ELSE IF( ldv.LT.max( 1, n ) ) THEN
         info = -9
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEBAK', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
      IF( m.EQ.0 )
     $   RETURN
      IF( lsame( job, 'N' ) )
     $   RETURN
*
      IF( ilo.EQ.ihi )
     $   GO TO 30
*
*     Backward balance
*
      IF( lsame( job, 'S' ) .OR. lsame( job, 'B' ) ) THEN
*
         IF( rightv ) THEN
            DO 10 i = ilo, ihi
               s = scale( i )
               CALL dscal( m, s, v( i, 1 ), ldv )
   10       CONTINUE
         END IF
*
         IF( leftv ) THEN
            DO 20 i = ilo, ihi
               s = one / scale( i )
               CALL dscal( m, s, v( i, 1 ), ldv )
   20       CONTINUE
         END IF
*
      END IF
*
*     Backward permutation
*
*     For  I = ILO-1 step -1 until 1,
*              IHI+1 step 1 until N do --
*
   30 CONTINUE
      IF( lsame( job, 'P' ) .OR. lsame( job, 'B' ) ) THEN
         IF( rightv ) THEN
            DO 40 ii = 1, n
               i = ii
               IF( i.GE.ilo .AND. i.LE.ihi )
     $            GO TO 40
               IF( i.LT.ilo )
     $            i = ilo - ii
               k = scale( i )
               IF( k.EQ.i )
     $            GO TO 40
               CALL dswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
   40       CONTINUE
         END IF
*
         IF( leftv ) THEN
            DO 50 ii = 1, n
               i = ii
               IF( i.GE.ilo .AND. i.LE.ihi )
     $            GO TO 50
               IF( i.LT.ilo )
     $            i = ilo - ii
               k = scale( i )
               IF( k.EQ.i )
     $            GO TO 50
               CALL dswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
   50       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEBAK
*
      END
      SUBROUTINE dgebal( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), SCALE( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      DOUBLE PRECISION   SCLFAC
      parameter( sclfac = 2.0d+0 )
      DOUBLE PRECISION   FACTOR
      parameter( factor = 0.95d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IEXC, IRA, J, K, L, M
      DOUBLE PRECISION   C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1,
     $                   SFMIN2
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN, LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DNRM2
      EXTERNAL           disnan, lsame, idamax, dlamch, dnrm2
*     ..
*     .. External Subroutines ..
      EXTERNAL           dscal, dswap, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min
*     ..
*     Test the input parameters
*
      info = 0
      IF( .NOT.lsame( job, 'N' ) .AND. .NOT.lsame( job, 'P' ) .AND.
     $    .NOT.lsame( job, 'S' ) .AND. .NOT.lsame( job, 'B' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEBAL', -info )
         RETURN
      END IF
*
      k = 1
      l = n
*
      IF( n.EQ.0 )
     $   GO TO 210
*
      IF( lsame( job, 'N' ) ) THEN
         DO 10 i = 1, n
            scale( i ) = one
   10    CONTINUE
         GO TO 210
      END IF
*
      IF( lsame( job, 'S' ) )
     $   GO TO 120
*
*     Permutation to isolate eigenvalues if possible
*
      GO TO 50
*
*     Row and column exchange.
*
   20 CONTINUE
      scale( m ) = j
      IF( j.EQ.m )
     $   GO TO 30
*
      CALL dswap( l, a( 1, j ), 1, a( 1, m ), 1 )
      CALL dswap( n-k+1, a( j, k ), lda, a( m, k ), lda )
*
   30 CONTINUE
      GO TO ( 40, 80 )iexc
*
*     Search for rows isolating an eigenvalue and push them down.
*
   40 CONTINUE
      IF( l.EQ.1 )
     $   GO TO 210
      l = l - 1
*
   50 CONTINUE
      DO 70 j = l, 1, -1
*
         DO 60 i = 1, l
            IF( i.EQ.j )
     $         GO TO 60
            IF( a( j, i ).NE.zero )
     $         GO TO 70
   60    CONTINUE
*
         m = l
         iexc = 1
         GO TO 20
   70 CONTINUE
*
      GO TO 90
*
*     Search for columns isolating an eigenvalue and push them left.
*
   80 CONTINUE
      k = k + 1
*
   90 CONTINUE
      DO 110 j = k, l
*
         DO 100 i = k, l
            IF( i.EQ.j )
     $         GO TO 100
            IF( a( i, j ).NE.zero )
     $         GO TO 110
  100    CONTINUE
*
         m = k
         iexc = 2
         GO TO 20
  110 CONTINUE
*
  120 CONTINUE
      DO 130 i = k, l
         scale( i ) = one
  130 CONTINUE
*
      IF( lsame( job, 'P' ) )
     $   GO TO 210
*
*     Balance the submatrix in rows K to L.
*
*     Iterative loop for norm reduction
*
      sfmin1 = dlamch( 'S' ) / dlamch( 'P' )
      sfmax1 = one / sfmin1
      sfmin2 = sfmin1*sclfac
      sfmax2 = one / sfmin2
*
  140 CONTINUE
      noconv = .false.
*
      DO 200 i = k, l
*
         c = dnrm2( l-k+1, a( k, i ), 1 )
         r = dnrm2( l-k+1, a( i, k ), lda )
         ica = idamax( l, a( 1, i ), 1 )
         ca = abs( a( ica, i ) )
         ira = idamax( n-k+1, a( i, k ), lda )
         ra = abs( a( i, ira+k-1 ) )
*
*        Guard against zero C or R due to underflow.
*
         IF( c.EQ.zero .OR. r.EQ.zero )
     $      GO TO 200
         g = r / sclfac
         f = one
         s = c + r
  160    CONTINUE
         IF( c.GE.g .OR. max( f, c, ca ).GE.sfmax2 .OR.
     $       min( r, g, ra ).LE.sfmin2 )GO TO 170
            IF( disnan( c+f+ca+r+g+ra ) ) THEN
*
*           Exit if NaN to avoid infinite loop
*
            info = -3
            CALL xerbla( 'DGEBAL', -info )
            RETURN
         END IF
         f = f*sclfac
         c = c*sclfac
         ca = ca*sclfac
         r = r / sclfac
         g = g / sclfac
         ra = ra / sclfac
         GO TO 160
*
  170    CONTINUE
         g = c / sclfac
  180    CONTINUE
         IF( g.LT.r .OR. max( r, ra ).GE.sfmax2 .OR.
     $       min( f, c, g, ca ).LE.sfmin2 )GO TO 190
         f = f / sclfac
         c = c / sclfac
         g = g / sclfac
         ca = ca / sclfac
         r = r*sclfac
         ra = ra*sclfac
         GO TO 180
*
*        Now balance.
*
  190    CONTINUE
         IF( ( c+r ).GE.factor*s )
     $      GO TO 200
         IF( f.LT.one .AND. scale( i ).LT.one ) THEN
            IF( f*scale( i ).LE.sfmin1 )
     $         GO TO 200
         END IF
         IF( f.GT.one .AND. scale( i ).GT.one ) THEN
            IF( scale( i ).GE.sfmax1 / f )
     $         GO TO 200
         END IF
         g = one / f
         scale( i ) = scale( i )*f
         noconv = .true.
*
         CALL dscal( n-k+1, g, a( i, k ), lda )
         CALL dscal( l, f, a( 1, i ), 1 )
*
  200 CONTINUE
*
      IF( noconv )
     $   GO TO 140
*
  210 CONTINUE
      ilo = k
      ihi = l
*
      RETURN
*
*     End of DGEBAL
*
      END
      DOUBLE PRECISION FUNCTION dzasum(N,ZX,INCX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER incx,n
*     ..
*     .. Array Arguments ..
      COMPLEX*16 zx(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION stemp
      INTEGER i,nincx
*     ..
*     .. External Functions ..
      DOUBLE PRECISION dcabs1
      EXTERNAL dcabs1
*     ..
      dzasum = 0.0d0
      stemp = 0.0d0
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO i = 1,n
            stemp = stemp + dcabs1(zx(i))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            stemp = stemp + dcabs1(zx(i))
         END DO
      END IF
      dzasum = stemp
      RETURN
*
*     End of DZASUM
*
      END
      LOGICAL FUNCTION disnan( DIN )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: din
*     ..
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL dlaisnan
      EXTERNAL dlaisnan
*  ..
*  .. Executable Statements ..
      disnan = dlaisnan(din,din)
      RETURN
      END
      SUBROUTINE dlanv2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      parameter( zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0,
     $                     two = 2.0d0 )
      DOUBLE PRECISION   MULTPL
      parameter( multpl = 4.0d+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BCMAX, BCMIS, CC, CS1, DD, EPS, P, SAB,
     $                   SAC, SCALE, SIGMA, SN1, TAU, TEMP, Z, SAFMIN, 
     $                   SAFMN2, SAFMX2
      INTEGER            COUNT
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           dlamch, dlapy2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min, sign, sqrt
*     ..
*     .. Executable Statements ..
*
      safmin = dlamch( 'S' )
      eps = dlamch( 'P' )
      safmn2 = dlamch( 'B' )**int( log( safmin / eps ) /
     $            log( dlamch( 'B' ) ) / two )
      safmx2 = one / safmn2
      IF( c.EQ.zero ) THEN
         cs = one
         sn = zero
*
      ELSE IF( b.EQ.zero ) THEN
*
*        Swap rows and columns
*
         cs = zero
         sn = one
         temp = d
         d = a
         a = temp
         b = -c
         c = zero
*
      ELSE IF( ( a-d ).EQ.zero .AND. sign( one, b ).NE.sign( one, c ) )
     $          THEN
         cs = one
         sn = zero
*
      ELSE
*
         temp = a - d
         p = half*temp
         bcmax = max( abs( b ), abs( c ) )
         bcmis = min( abs( b ), abs( c ) )*sign( one, b )*sign( one, c )
         scale = max( abs( p ), bcmax )
         z = ( p / scale )*p + ( bcmax / scale )*bcmis
*
*        If Z is of the order of the machine accuracy, postpone the
*        decision on the nature of eigenvalues
*
         IF( z.GE.multpl*eps ) THEN
*
*           Real eigenvalues. Compute A and D.
*
            z = p + sign( sqrt( scale )*sqrt( z ), p )
            a = d + z
            d = d - ( bcmax / z )*bcmis
*
*           Compute B and the rotation matrix
*
            tau = dlapy2( c, z )
            cs = z / tau
            sn = c / tau
            b = b - c
            c = zero
*
         ELSE
*
*           Complex eigenvalues, or real (almost) equal eigenvalues.
*           Make diagonal elements equal.
*
            count = 0
            sigma = b + c
   10       CONTINUE
            count = count + 1
            scale = max( abs(temp), abs(sigma) )
            IF( scale.GE.safmx2 ) THEN
               sigma = sigma * safmn2
               temp = temp * safmn2
               IF (count .LE. 20)
     $            GOTO 10
            END IF
            IF( scale.LE.safmn2 ) THEN
               sigma = sigma * safmx2
               temp = temp * safmx2
               IF (count .LE. 20)
     $            GOTO 10
            END IF
            p = half*temp
            tau = dlapy2( sigma, temp )
            cs = sqrt( half*( one+abs( sigma ) / tau ) )
            sn = -( p / ( tau*cs ) )*sign( one, sigma )
*
*           Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
*                   [ CC  DD ]   [ C  D ] [ SN  CS ]
*
            aa = a*cs + b*sn
            bb = -a*sn + b*cs
            cc = c*cs + d*sn
            dd = -c*sn + d*cs
*
*           Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
*                   [ C  D ]   [-SN  CS ] [ CC  DD ]
*
            a = aa*cs + cc*sn
            b = bb*cs + dd*sn
            c = -aa*sn + cc*cs
            d = -bb*sn + dd*cs
*
            temp = half*( a+d )
            a = temp
            d = temp
*
            IF( c.NE.zero ) THEN
               IF( b.NE.zero ) THEN
                  IF( sign( one, b ).EQ.sign( one, c ) ) THEN
*
*                    Real eigenvalues: reduce to upper triangular form
*
                     sab = sqrt( abs( b ) )
                     sac = sqrt( abs( c ) )
                     p = sign( sab*sac, c )
                     tau = one / sqrt( abs( b+c ) )
                     a = temp + p
                     d = temp - p
                     b = b - c
                     c = zero
                     cs1 = sab*tau
                     sn1 = sac*tau
                     temp = cs*cs1 - sn*sn1
                     sn = cs*sn1 + sn*cs1
                     cs = temp
                  END IF
               ELSE
                  b = -c
                  c = zero
                  temp = cs
                  cs = -sn
                  sn = temp
               END IF
            END IF
         END IF
*
      END IF
*
*     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
*
      rt1r = a
      rt2r = d
      IF( c.EQ.zero ) THEN
         rt1i = zero
         rt2i = zero
      ELSE
         rt1i = sqrt( abs( b ) )*sqrt( abs( c ) )
         rt2i = -rt1i
      END IF
      RETURN
*
*     End of DLANV2
*
      END
      SUBROUTINE dlaqr3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,
     $                   LDT, NV, WV, LDWV, WORK, LWORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
     $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),
     $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),
     $                   z( ldz, * )
*     ..
*
*  ================================================================
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S,
     $                   SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL,
     $                   KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3,
     $                   lwkopt, nmin
      LOGICAL            BULGE, SORTED
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            ILAENV
      EXTERNAL           dlamch, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           dcopy, dgehrd, dgemm, dlabad, dlacpy, dlahqr,
     $                   dlanv2, dlaqr4, dlarf, dlarfg, dlaset, dormhr,
     $                   dtrexc
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, int, max, min, sqrt
*     ..
*     .. Executable Statements ..
*
*     ==== Estimate optimal workspace. ====
*
      jw = min( nw, kbot-ktop+1 )
      IF( jw.LE.2 ) THEN
         lwkopt = 1
      ELSE
*
*        ==== Workspace query call to DGEHRD ====
*
         CALL dgehrd( jw, 1, jw-1, t, ldt, work, work, -1, info )
         lwk1 = int( work( 1 ) )
*
*        ==== Workspace query call to DORMHR ====
*
         CALL dormhr( 'R', 'N', jw, jw, 1, jw-1, t, ldt, work, v, ldv,
     $                work, -1, info )
         lwk2 = int( work( 1 ) )
*
*        ==== Workspace query call to DLAQR4 ====
*
         CALL dlaqr4( .true., .true., jw, 1, jw, t, ldt, sr, si, 1, jw,
     $                v, ldv, work, -1, infqr )
         lwk3 = int( work( 1 ) )
*
*        ==== Optimal workspace ====
*
         lwkopt = max( jw+max( lwk1, lwk2 ), lwk3 )
      END IF
*
*     ==== Quick return in case of workspace query. ====
*
      IF( lwork.EQ.-1 ) THEN
         work( 1 ) = dble( lwkopt )
         RETURN
      END IF
*
*     ==== Nothing to do ...
*     ... for an empty active block ... ====
      ns = 0
      nd = 0
      work( 1 ) = one
      IF( ktop.GT.kbot )
     $   RETURN
*     ... nor for an empty deflation window. ====
      IF( nw.LT.1 )
     $   RETURN
*
*     ==== Machine constants ====
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = one / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
*
*     ==== Setup deflation window ====
*
      jw = min( nw, kbot-ktop+1 )
      kwtop = kbot - jw + 1
      IF( kwtop.EQ.ktop ) THEN
         s = zero
      ELSE
         s = h( kwtop, kwtop-1 )
      END IF
*
      IF( kbot.EQ.kwtop ) THEN
*
*        ==== 1-by-1 deflation window: not much to do ====
*
         sr( kwtop ) = h( kwtop, kwtop )
         si( kwtop ) = zero
         ns = 1
         nd = 0
         IF( abs( s ).LE.max( smlnum, ulp*abs( h( kwtop, kwtop ) ) ) )
     $        THEN
            ns = 0
            nd = 1
            IF( kwtop.GT.ktop )
     $         h( kwtop, kwtop-1 ) = zero
         END IF
         work( 1 ) = one
         RETURN
      END IF
*
*     ==== Convert to spike-triangular form.  (In case of a
*     .    rare QR failure, this routine continues to do
*     .    aggressive early deflation using that part of
*     .    the deflation window that converged using INFQR
*     .    here and there to keep track.) ====
*
      CALL dlacpy( 'U', jw, jw, h( kwtop, kwtop ), ldh, t, ldt )
      CALL dcopy( jw-1, h( kwtop+1, kwtop ), ldh+1, t( 2, 1 ), ldt+1 )
*
      CALL dlaset( 'A', jw, jw, zero, one, v, ldv )
      nmin = ilaenv( 12, 'DLAQR3', 'SV', jw, 1, jw, lwork )
      IF( jw.GT.nmin ) THEN
         CALL dlaqr4( .true., .true., jw, 1, jw, t, ldt, sr( kwtop ),
     $                si( kwtop ), 1, jw, v, ldv, work, lwork, infqr )
      ELSE
         CALL dlahqr( .true., .true., jw, 1, jw, t, ldt, sr( kwtop ),
     $                si( kwtop ), 1, jw, v, ldv, infqr )
      END IF
*
*     ==== DTREXC needs a clean margin near the diagonal ====
*
      DO 10 j = 1, jw - 3
         t( j+2, j ) = zero
         t( j+3, j ) = zero
   10 CONTINUE
      IF( jw.GT.2 )
     $   t( jw, jw-2 ) = zero
*
*     ==== Deflation detection loop ====
*
      ns = jw
      ilst = infqr + 1
   20 CONTINUE
      IF( ilst.LE.ns ) THEN
         IF( ns.EQ.1 ) THEN
            bulge = .false.
         ELSE
            bulge = t( ns, ns-1 ).NE.zero
         END IF
*
*        ==== Small spike tip test for deflation ====
*
         IF( .NOT. bulge ) THEN
*
*           ==== Real eigenvalue ====
*
            foo = abs( t( ns, ns ) )
            IF( foo.EQ.zero )
     $         foo = abs( s )
            IF( abs( s*v( 1, ns ) ).LE.max( smlnum, ulp*foo ) ) THEN
*
*              ==== Deflatable ====
*
               ns = ns - 1
            ELSE
*
*              ==== Undeflatable.   Move it up out of the way.
*              .    (DTREXC can not fail in this case.) ====
*
               ifst = ns
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, work,
     $                      info )
               ilst = ilst + 1
            END IF
         ELSE
*
*           ==== Complex conjugate pair ====
*
            foo = abs( t( ns, ns ) ) + sqrt( abs( t( ns, ns-1 ) ) )*
     $            sqrt( abs( t( ns-1, ns ) ) )
            IF( foo.EQ.zero )
     $         foo = abs( s )
            IF( max( abs( s*v( 1, ns ) ), abs( s*v( 1, ns-1 ) ) ).LE.
     $          max( smlnum, ulp*foo ) ) THEN
*
*              ==== Deflatable ====
*
               ns = ns - 2
            ELSE
*
*              ==== Undeflatable. Move them up out of the way.
*              .    Fortunately, DTREXC does the right thing with
*              .    ILST in case of a rare exchange failure. ====
*
               ifst = ns
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, work,
     $                      info )
               ilst = ilst + 2
            END IF
         END IF
*
*        ==== End deflation detection loop ====
*
         GO TO 20
      END IF
*
*        ==== Return to Hessenberg form ====
*
      IF( ns.EQ.0 )
     $   s = zero
*
      IF( ns.LT.jw ) THEN
*
*        ==== sorting diagonal blocks of T improves accuracy for
*        .    graded matrices.  Bubble sort deals well with
*        .    exchange failures. ====
*
         sorted = .false.
         i = ns + 1
   30    CONTINUE
         IF( sorted )
     $      GO TO 50
         sorted = .true.
*
         kend = i - 1
         i = infqr + 1
         IF( i.EQ.ns ) THEN
            k = i + 1
         ELSE IF( t( i+1, i ).EQ.zero ) THEN
            k = i + 1
         ELSE
            k = i + 2
         END IF
   40    CONTINUE
         IF( k.LE.kend ) THEN
            IF( k.EQ.i+1 ) THEN
               evi = abs( t( i, i ) )
            ELSE
               evi = abs( t( i, i ) ) + sqrt( abs( t( i+1, i ) ) )*
     $               sqrt( abs( t( i, i+1 ) ) )
            END IF
*
            IF( k.EQ.kend ) THEN
               evk = abs( t( k, k ) )
            ELSE IF( t( k+1, k ).EQ.zero ) THEN
               evk = abs( t( k, k ) )
            ELSE
               evk = abs( t( k, k ) ) + sqrt( abs( t( k+1, k ) ) )*
     $               sqrt( abs( t( k, k+1 ) ) )
            END IF
*
            IF( evi.GE.evk ) THEN
               i = k
            ELSE
               sorted = .false.
               ifst = i
               ilst = k
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, work,
     $                      info )
               IF( info.EQ.0 ) THEN
                  i = ilst
               ELSE
                  i = k
               END IF
            END IF
            IF( i.EQ.kend ) THEN
               k = i + 1
            ELSE IF( t( i+1, i ).EQ.zero ) THEN
               k = i + 1
            ELSE
               k = i + 2
            END IF
            GO TO 40
         END IF
         GO TO 30
   50    CONTINUE
      END IF
*
*     ==== Restore shift/eigenvalue array from T ====
*
      i = jw
   60 CONTINUE
      IF( i.GE.infqr+1 ) THEN
         IF( i.EQ.infqr+1 ) THEN
            sr( kwtop+i-1 ) = t( i, i )
            si( kwtop+i-1 ) = zero
            i = i - 1
         ELSE IF( t( i, i-1 ).EQ.zero ) THEN
            sr( kwtop+i-1 ) = t( i, i )
            si( kwtop+i-1 ) = zero
            i = i - 1
         ELSE
            aa = t( i-1, i-1 )
            cc = t( i, i-1 )
            bb = t( i-1, i )
            dd = t( i, i )
            CALL dlanv2( aa, bb, cc, dd, sr( kwtop+i-2 ),
     $                   si( kwtop+i-2 ), sr( kwtop+i-1 ),
     $                   si( kwtop+i-1 ), cs, sn )
            i = i - 2
         END IF
         GO TO 60
      END IF
*
      IF( ns.LT.jw .OR. s.EQ.zero ) THEN
         IF( ns.GT.1 .AND. s.NE.zero ) THEN
*
*           ==== Reflect spike back into lower triangle ====
*
            CALL dcopy( ns, v, ldv, work, 1 )
            beta = work( 1 )
            CALL dlarfg( ns, beta, work( 2 ), 1, tau )
            work( 1 ) = one
*
            CALL dlaset( 'L', jw-2, jw-2, zero, zero, t( 3, 1 ), ldt )
*
            CALL dlarf( 'L', ns, jw, work, 1, tau, t, ldt,
     $                  work( jw+1 ) )
            CALL dlarf( 'R', ns, ns, work, 1, tau, t, ldt,
     $                  work( jw+1 ) )
            CALL dlarf( 'R', jw, ns, work, 1, tau, v, ldv,
     $                  work( jw+1 ) )
*
            CALL dgehrd( jw, 1, ns, t, ldt, work, work( jw+1 ),
     $                   lwork-jw, info )
         END IF
*
*        ==== Copy updated reduced window into place ====
*
         IF( kwtop.GT.1 )
     $      h( kwtop, kwtop-1 ) = s*v( 1, 1 )
         CALL dlacpy( 'U', jw, jw, t, ldt, h( kwtop, kwtop ), ldh )
         CALL dcopy( jw-1, t( 2, 1 ), ldt+1, h( kwtop+1, kwtop ),
     $               ldh+1 )
*
*        ==== Accumulate orthogonal matrix in order update
*        .    H and Z, if requested.  ====
*
         IF( ns.GT.1 .AND. s.NE.zero )
     $      CALL dormhr( 'R', 'N', jw, ns, 1, ns, t, ldt, work, v, ldv,
     $                   work( jw+1 ), lwork-jw, info )
*
*        ==== Update vertical slab in H ====
*
         IF( wantt ) THEN
            ltop = 1
         ELSE
            ltop = ktop
         END IF
         DO 70 krow = ltop, kwtop - 1, nv
            kln = min( nv, kwtop-krow )
            CALL dgemm( 'N', 'N', kln, jw, jw, one, h( krow, kwtop ),
     $                  ldh, v, ldv, zero, wv, ldwv )
            CALL dlacpy( 'A', kln, jw, wv, ldwv, h( krow, kwtop ), ldh )
   70    CONTINUE
*
*        ==== Update horizontal slab in H ====
*
         IF( wantt ) THEN
            DO 80 kcol = kbot + 1, n, nh
               kln = min( nh, n-kcol+1 )
               CALL dgemm( 'C', 'N', jw, kln, jw, one, v, ldv,
     $                     h( kwtop, kcol ), ldh, zero, t, ldt )
               CALL dlacpy( 'A', jw, kln, t, ldt, h( kwtop, kcol ),
     $                      ldh )
   80       CONTINUE
         END IF
*
*        ==== Update vertical slab in Z ====
*
         IF( wantz ) THEN
            DO 90 krow = iloz, ihiz, nv
               kln = min( nv, ihiz-krow+1 )
               CALL dgemm( 'N', 'N', kln, jw, jw, one, z( krow, kwtop ),
     $                     ldz, v, ldv, zero, wv, ldwv )
               CALL dlacpy( 'A', kln, jw, wv, ldwv, z( krow, kwtop ),
     $                      ldz )
   90       CONTINUE
         END IF
      END IF
*
*     ==== Return the number of deflations ... ====
*
      nd = jw - ns
*
*     ==== ... and the number of shifts. (Subtracting
*     .    INFQR from the spike length takes care
*     .    of the case of a rare QR failure while
*     .    calculating eigenvalues of the deflation
*     .    window.)  ====
*
      ns = ns - infqr
*
*      ==== Return optimal workspace. ====
*
      work( 1 ) = dble( lwkopt )
*
*     ==== End of DLAQR3 ====
*
      END
      SUBROUTINE dlaqr4( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
     $                   z( ldz, * )
*     ..
*
*  ================================================================
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    DLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 15 )
*
*     ==== Exceptional deflation windows:  try to cure rare
*     .    slow convergence by varying the size of the
*     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      parameter( kexnw = 5 )
*
*     ==== Exceptional shifts: try to cure rare slow convergence
*     .    with ad-hoc exceptional shifts every KEXSH iterations.
*     .    ====
      INTEGER            KEXSH
      parameter( kexsh = 6 )
*
*     ==== The constants WILK1 and WILK2 are used to form the
*     .    exceptional shifts. ====
      DOUBLE PRECISION   WILK1, WILK2
      parameter( wilk1 = 0.75d0, wilk2 = -0.4375d0 )
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
     $                   kt, ktop, ku, kv, kwh, kwtop, kwv, ld, ls,
     $                   lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns,
     $                   nsmax, nsr, nve, nw, nwmax, nwr, nwupbd
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   ZDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlacpy, dlahqr, dlanv2, dlaqr2, dlaqr5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, int, max, min, mod
*     ..
*     .. Executable Statements ..
      info = 0
*
*     ==== Quick return for N = 0: nothing to do. ====
*
      IF( n.EQ.0 ) THEN
         work( 1 ) = one
         RETURN
      END IF
*
      IF( n.LE.ntiny ) THEN
*
*        ==== Tiny matrices must use DLAHQR. ====
*
         lwkopt = 1
         IF( lwork.NE.-1 )
     $      CALL dlahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi,
     $                   iloz, ihiz, z, ldz, info )
      ELSE
*
*        ==== Use small bulge multi-shift QR with aggressive early
*        .    deflation on larger-than-tiny matrices. ====
*
*        ==== Hope for the best. ====
*
         info = 0
*
*        ==== Set up job flags for ILAENV. ====
*
         IF( wantt ) THEN
            jbcmpz( 1: 1 ) = 'S'
         ELSE
            jbcmpz( 1: 1 ) = 'E'
         END IF
         IF( wantz ) THEN
            jbcmpz( 2: 2 ) = 'V'
         ELSE
            jbcmpz( 2: 2 ) = 'N'
         END IF
*
*        ==== NWR = recommended deflation window size.  At this
*        .    point,  N .GT. NTINY = 15, so there is enough
*        .    subdiagonal workspace for NWR.GE.2 as required.
*        .    (In fact, there is enough subdiagonal space for
*        .    NWR.GE.4.) ====
*
         nwr = ilaenv( 13, 'DLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nwr = max( 2, nwr )
         nwr = min( ihi-ilo+1, ( n-1 ) / 3, nwr )
*
*        ==== NSR = recommended number of simultaneous shifts.
*        .    At this point N .GT. NTINY = 15, so there is at
*        .    enough subdiagonal workspace for NSR to be even
*        .    and greater than or equal to two as required. ====
*
         nsr = ilaenv( 15, 'DLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nsr = min( nsr, ( n-3 ) / 6, ihi-ilo )
         nsr = max( 2, nsr-mod( nsr, 2 ) )
*
*        ==== Estimate optimal workspace ====
*
*        ==== Workspace query call to DLAQR2 ====
*
         CALL dlaqr2( wantt, wantz, n, ilo, ihi, nwr+1, h, ldh, iloz,
     $                ihiz, z, ldz, ls, ld, wr, wi, h, ldh, n, h, ldh,
     $                n, h, ldh, work, -1 )
*
*        ==== Optimal workspace = MAX(DLAQR5, DLAQR2) ====
*
         lwkopt = max( 3*nsr / 2, int( work( 1 ) ) )
*
*        ==== Quick return in case of workspace query. ====
*
         IF( lwork.EQ.-1 ) THEN
            work( 1 ) = dble( lwkopt )
            RETURN
         END IF
*
*        ==== DLAHQR/DLAQR0 crossover point ====
*
         nmin = ilaenv( 12, 'DLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
*
*        ==== Nibble crossover point ====
*
         nibble = ilaenv( 14, 'DLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nibble = max( 0, nibble )
*
*        ==== Accumulate reflections during ttswp?  Use block
*        .    2-by-2 structure during matrix-matrix multiply? ====
*
         kacc22 = ilaenv( 16, 'DLAQR4', jbcmpz, n, ilo, ihi, lwork )
         kacc22 = max( 0, kacc22 )
         kacc22 = min( 2, kacc22 )
*
*        ==== NWMAX = the largest possible deflation window for
*        .    which there is sufficient workspace. ====
*
         nwmax = min( ( n-1 ) / 3, lwork / 2 )
         nw = nwmax
*
*        ==== NSMAX = the Largest number of simultaneous shifts
*        .    for which there is sufficient workspace. ====
*
         nsmax = min( ( n-3 ) / 6, 2*lwork / 3 )
         nsmax = nsmax - mod( nsmax, 2 )
*
*        ==== NDFL: an iteration count restarted at deflation. ====
*
         ndfl = 1
*
*        ==== ITMAX = iteration limit ====
*
         itmax = max( 30, 2*kexsh )*max( 10, ( ihi-ilo+1 ) )
*
*        ==== Last row and column in the active block ====
*
         kbot = ihi
*
*        ==== Main Loop ====
*
         DO 80 it = 1, itmax
*
*           ==== Done when KBOT falls below ILO ====
*
            IF( kbot.LT.ilo )
     $         GO TO 90
*
*           ==== Locate active block ====
*
            DO 10 k = kbot, ilo + 1, -1
               IF( h( k, k-1 ).EQ.zero )
     $            GO TO 20
   10       CONTINUE
            k = ilo
   20       CONTINUE
            ktop = k
*
*           ==== Select deflation window size:
*           .    Typical Case:
*           .      If possible and advisable, nibble the entire
*           .      active block.  If not, use size MIN(NWR,NWMAX)
*           .      or MIN(NWR+1,NWMAX) depending upon which has
*           .      the smaller corresponding subdiagonal entry
*           .      (a heuristic).
*           .
*           .    Exceptional Case:
*           .      If there have been no deflations in KEXNW or
*           .      more iterations, then vary the deflation window
*           .      size.   At first, because, larger windows are,
*           .      in general, more powerful than smaller ones,
*           .      rapidly increase the window to the maximum possible.
*           .      Then, gradually reduce the window size. ====
*
            nh = kbot - ktop + 1
            nwupbd = min( nh, nwmax )
            IF( ndfl.LT.kexnw ) THEN
               nw = min( nwupbd, nwr )
            ELSE
               nw = min( nwupbd, 2*nw )
            END IF
            IF( nw.LT.nwmax ) THEN
               IF( nw.GE.nh-1 ) THEN
                  nw = nh
               ELSE
                  kwtop = kbot - nw + 1
                  IF( abs( h( kwtop, kwtop-1 ) ).GT.
     $                abs( h( kwtop-1, kwtop-2 ) ) )nw = nw + 1
               END IF
            END IF
            IF( ndfl.LT.kexnw ) THEN
               ndec = -1
            ELSE IF( ndec.GE.0 .OR. nw.GE.nwupbd ) THEN
               ndec = ndec + 1
               IF( nw-ndec.LT.2 )
     $            ndec = 0
               nw = nw - ndec
            END IF
*
*           ==== Aggressive early deflation:
*           .    split workspace under the subdiagonal into
*           .      - an nw-by-nw work array V in the lower
*           .        left-hand-corner,
*           .      - an NW-by-at-least-NW-but-more-is-better
*           .        (NW-by-NHO) horizontal work array along
*           .        the bottom edge,
*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
*           .        vertical work array along the left-hand-edge.
*           .        ====
*
            kv = n - nw + 1
            kt = nw + 1
            nho = ( n-nw-1 ) - kt + 1
            kwv = nw + 2
            nve = ( n-nw ) - kwv + 1
*
*           ==== Aggressive early deflation ====
*
            CALL dlaqr2( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz,
     $                   ihiz, z, ldz, ls, ld, wr, wi, h( kv, 1 ), ldh,
     $                   nho, h( kv, kt ), ldh, nve, h( kwv, 1 ), ldh,
     $                   work, lwork )
*
*           ==== Adjust KBOT accounting for new deflations. ====
*
            kbot = kbot - ld
*
*           ==== KS points to the shifts. ====
*
            ks = kbot - ls + 1
*
*           ==== Skip an expensive QR sweep if there is a (partly
*           .    heuristic) reason to expect that many eigenvalues
*           .    will deflate without it.  Here, the QR sweep is
*           .    skipped if many eigenvalues have just been deflated
*           .    or if the remaining active block is small.
*
            IF( ( ld.EQ.0 ) .OR. ( ( 100*ld.LE.nw*nibble ) .AND. ( kbot-
     $          ktop+1.GT.min( nmin, nwmax ) ) ) ) THEN
*
*              ==== NS = nominal number of simultaneous shifts.
*              .    This may be lowered (slightly) if DLAQR2
*              .    did not provide that many shifts. ====
*
               ns = min( nsmax, nsr, max( 2, kbot-ktop ) )
               ns = ns - mod( ns, 2 )
*
*              ==== If there have been no deflations
*              .    in a multiple of KEXSH iterations,
*              .    then try exceptional shifts.
*              .    Otherwise use shifts provided by
*              .    DLAQR2 above or from the eigenvalues
*              .    of a trailing principal submatrix. ====
*
               IF( mod( ndfl, kexsh ).EQ.0 ) THEN
                  ks = kbot - ns + 1
                  DO 30 i = kbot, max( ks+1, ktop+2 ), -2
                     ss = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
                     aa = wilk1*ss + h( i, i )
                     bb = ss
                     cc = wilk2*ss
                     dd = aa
                     CALL dlanv2( aa, bb, cc, dd, wr( i-1 ), wi( i-1 ),
     $                            wr( i ), wi( i ), cs, sn )
   30             CONTINUE
                  IF( ks.EQ.ktop ) THEN
                     wr( ks+1 ) = h( ks+1, ks+1 )
                     wi( ks+1 ) = zero
                     wr( ks ) = wr( ks+1 )
                     wi( ks ) = wi( ks+1 )
                  END IF
               ELSE
*
*                 ==== Got NS/2 or fewer shifts? Use DLAHQR
*                 .    on a trailing principal submatrix to
*                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
*                 .    there is enough space below the subdiagonal
*                 .    to fit an NS-by-NS scratch array.) ====
*
                  IF( kbot-ks+1.LE.ns / 2 ) THEN
                     ks = kbot - ns + 1
                     kt = n - ns + 1
                     CALL dlacpy( 'A', ns, ns, h( ks, ks ), ldh,
     $                            h( kt, 1 ), ldh )
                     CALL dlahqr( .false., .false., ns, 1, ns,
     $                            h( kt, 1 ), ldh, wr( ks ), wi( ks ),
     $                            1, 1, zdum, 1, inf )
                     ks = ks + inf
*
*                    ==== In case of a rare QR failure use
*                    .    eigenvalues of the trailing 2-by-2
*                    .    principal submatrix.  ====
*
                     IF( ks.GE.kbot ) THEN
                        aa = h( kbot-1, kbot-1 )
                        cc = h( kbot, kbot-1 )
                        bb = h( kbot-1, kbot )
                        dd = h( kbot, kbot )
                        CALL dlanv2( aa, bb, cc, dd, wr( kbot-1 ),
     $                               wi( kbot-1 ), wr( kbot ),
     $                               wi( kbot ), cs, sn )
                        ks = kbot - 1
                     END IF
                  END IF
*
                  IF( kbot-ks+1.GT.ns ) THEN
*
*                    ==== Sort the shifts (Helps a little)
*                    .    Bubble sort keeps complex conjugate
*                    .    pairs together. ====
*
                     sorted = .false.
                     DO 50 k = kbot, ks + 1, -1
                        IF( sorted )
     $                     GO TO 60
                        sorted = .true.
                        DO 40 i = ks, k - 1
                           IF( abs( wr( i ) )+abs( wi( i ) ).LT.
     $                         abs( wr( i+1 ) )+abs( wi( i+1 ) ) ) THEN
                              sorted = .false.
*
                              swap = wr( i )
                              wr( i ) = wr( i+1 )
                              wr( i+1 ) = swap
*
                              swap = wi( i )
                              wi( i ) = wi( i+1 )
                              wi( i+1 ) = swap
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
*
*                 ==== Shuffle shifts into pairs of real shifts
*                 .    and pairs of complex conjugate shifts
*                 .    assuming complex conjugate shifts are
*                 .    already adjacent to one another. (Yes,
*                 .    they are.)  ====
*
                  DO 70 i = kbot, ks + 2, -2
                     IF( wi( i ).NE.-wi( i-1 ) ) THEN
*
                        swap = wr( i )
                        wr( i ) = wr( i-1 )
                        wr( i-1 ) = wr( i-2 )
                        wr( i-2 ) = swap
*
                        swap = wi( i )
                        wi( i ) = wi( i-1 )
                        wi( i-1 ) = wi( i-2 )
                        wi( i-2 ) = swap
                     END IF
   70             CONTINUE
               END IF
*
*              ==== If there are only two shifts and both are
*              .    real, then use only one.  ====
*
               IF( kbot-ks+1.EQ.2 ) THEN
                  IF( wi( kbot ).EQ.zero ) THEN
                     IF( abs( wr( kbot )-h( kbot, kbot ) ).LT.
     $                   abs( wr( kbot-1 )-h( kbot, kbot ) ) ) THEN
                        wr( kbot-1 ) = wr( kbot )
                     ELSE
                        wr( kbot ) = wr( kbot-1 )
                     END IF
                  END IF
               END IF
*
*              ==== Use up to NS of the the smallest magnitude
*              .    shifts.  If there aren't NS shifts available,
*              .    then use them all, possibly dropping one to
*              .    make the number of shifts even. ====
*
               ns = min( ns, kbot-ks+1 )
               ns = ns - mod( ns, 2 )
               ks = kbot - ns + 1
*
*              ==== Small-bulge multi-shift QR sweep:
*              .    split workspace under the subdiagonal into
*              .    - a KDU-by-KDU work array U in the lower
*              .      left-hand-corner,
*              .    - a KDU-by-at-least-KDU-but-more-is-better
*              .      (KDU-by-NHo) horizontal work array WH along
*              .      the bottom edge,
*              .    - and an at-least-KDU-but-more-is-better-by-KDU
*              .      (NVE-by-KDU) vertical work WV arrow along
*              .      the left-hand-edge. ====
*
               kdu = 2*ns
               ku = n - kdu + 1
               kwh = kdu + 1
               nho = ( n-kdu+1-4 ) - ( kdu+1 ) + 1
               kwv = kdu + 4
               nve = n - kdu - kwv + 1
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               CALL dlaqr5( wantt, wantz, kacc22, n, ktop, kbot, ns,
     $                      wr( ks ), wi( ks ), h, ldh, iloz, ihiz, z,
     $                      ldz, work, 3, h( ku, 1 ), ldh, nve,
     $                      h( kwv, 1 ), ldh, nho, h( ku, kwh ), ldh )
            END IF
*
*           ==== Note progress (or the lack of it). ====
*
            IF( ld.GT.0 ) THEN
               ndfl = 1
            ELSE
               ndfl = ndfl + 1
            END IF
*
*           ==== End of main loop ====
   80    CONTINUE
*
*        ==== Iteration limit exceeded.  Set INFO to show where
*        .    the problem occurred and exit. ====
*
         info = kbot
   90    CONTINUE
      END IF
*
*     ==== Return the optimal value of LWORK. ====
*
      work( 1 ) = dble( lwkopt )
*
*     ==== End of DLAQR4 ====
*
      END
      SUBROUTINE dlaqr5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS,
     $                   SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U,
     $                   LDU, NV, WV, LDWV, NH, WH, LDWH )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
     $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), U( LDU, * ),
     $                   V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ),
     $                   z( ldz, * )
*     ..
*
*  ================================================================
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ALPHA, BETA, H11, H12, H21, H22, REFSUM,
     $                   SAFMAX, SAFMIN, SCL, SMLNUM, SWAP, TST1, TST2,
     $                   ulp
      INTEGER            I, I2, I4, INCOL, J, JBOT, JCOL, JLEN,
     $                   JROW, JTOP, K, K1, KDU, KMS, KRCOL,
     $                   m, m22, mbot, mtop, nbmps, ndcol,
     $                   ns, nu
      LOGICAL            ACCUM, BMP22
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
*
      INTRINSIC          abs, dble, max, min, mod
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   VT( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgemm, dlabad, dlacpy, dlaqr1, dlarfg, dlaset,
     $                   dtrmm
*     ..
*     .. Executable Statements ..
*
*     ==== If there are no shifts, then there is nothing to do. ====
*
      IF( nshfts.LT.2 )
     $   RETURN
*
*     ==== If the active block is empty or 1-by-1, then there
*     .    is nothing to do. ====
*
      IF( ktop.GE.kbot )
     $   RETURN
*
*     ==== Shuffle shifts into pairs of real shifts and pairs
*     .    of complex conjugate shifts assuming complex
*     .    conjugate shifts are already adjacent to one
*     .    another. ====
*
      DO 10 i = 1, nshfts - 2, 2
         IF( si( i ).NE.-si( i+1 ) ) THEN
*
            swap = sr( i )
            sr( i ) = sr( i+1 )
            sr( i+1 ) = sr( i+2 )
            sr( i+2 ) = swap
*
            swap = si( i )
            si( i ) = si( i+1 )
            si( i+1 ) = si( i+2 )
            si( i+2 ) = swap
         END IF
   10 CONTINUE
*
*     ==== NSHFTS is supposed to be even, but if it is odd,
*     .    then simply reduce it by one.  The shuffle above
*     .    ensures that the dropped shift is real and that
*     .    the remaining shifts are paired. ====
*
      ns = nshfts - mod( nshfts, 2 )
*
*     ==== Machine constants for deflation ====
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = one / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
*
*     ==== Use accumulated reflections to update far-from-diagonal
*     .    entries ? ====
*
      accum = ( kacc22.EQ.1 ) .OR. ( kacc22.EQ.2 )
*
*     ==== clear trash ====
*
      IF( ktop+2.LE.kbot )
     $   h( ktop+2, ktop ) = zero
*
*     ==== NBMPS = number of 2-shift bulges in the chain ====
*
      nbmps = ns / 2
*
*     ==== KDU = width of slab ====
*
      kdu = 4*nbmps
*
*     ==== Create and chase chains of NBMPS bulges ====
*
      DO 180 incol = ktop - 2*nbmps + 1, kbot - 2, 2*nbmps
*
*        JTOP = Index from which updates from the right start.
*
         IF( accum ) THEN
            jtop = max( ktop, incol )
         ELSE IF( wantt ) THEN
            jtop = 1
         ELSE
            jtop = ktop
         END IF
*
         ndcol = incol + kdu
         IF( accum )
     $      CALL dlaset( 'ALL', kdu, kdu, zero, one, u, ldu )
*
*        ==== Near-the-diagonal bulge chase.  The following loop
*        .    performs the near-the-diagonal part of a small bulge
*        .    multi-shift QR sweep.  Each 4*NBMPS column diagonal
*        .    chunk extends from column INCOL to column NDCOL
*        .    (including both column INCOL and column NDCOL). The
*        .    following loop chases a 2*NBMPS+1 column long chain of
*        .    NBMPS bulges 2*NBMPS columns to the right.  (INCOL
*        .    may be less than KTOP and and NDCOL may be greater than
*        .    KBOT indicating phantom columns from which to chase
*        .    bulges before they are actually introduced or to which
*        .    to chase bulges beyond column KBOT.)  ====
*
         DO 145 krcol = incol, min( incol+2*nbmps-1, kbot-2 )
*
*           ==== Bulges number MTOP to MBOT are active double implicit
*           .    shift bulges.  There may or may not also be small
*           .    2-by-2 bulge, if there is room.  The inactive bulges
*           .    (if any) must wait until the active bulges have moved
*           .    down the diagonal to make room.  The phantom matrix
*           .    paradigm described above helps keep track.  ====
*
            mtop = max( 1, ( ktop-krcol ) / 2+1 )
            mbot = min( nbmps, ( kbot-krcol-1 ) / 2 )
            m22 = mbot + 1
            bmp22 = ( mbot.LT.nbmps ) .AND. ( krcol+2*( m22-1 ) ).EQ.
     $              ( kbot-2 )
*
*           ==== Generate reflections to chase the chain right
*           .    one column.  (The minimum value of K is KTOP-1.) ====
*
            IF ( bmp22 ) THEN
*
*              ==== Special case: 2-by-2 reflection at bottom treated
*              .    separately ====
*
               k = krcol + 2*( m22-1 )
               IF( k.EQ.ktop-1 ) THEN
                  CALL dlaqr1( 2, h( k+1, k+1 ), ldh, sr( 2*m22-1 ),
     $                         si( 2*m22-1 ), sr( 2*m22 ), si( 2*m22 ),
     $                         v( 1, m22 ) )
                  beta = v( 1, m22 )
                  CALL dlarfg( 2, beta, v( 2, m22 ), 1, v( 1, m22 ) )
               ELSE
                  beta = h( k+1, k )
                  v( 2, m22 ) = h( k+2, k )
                  CALL dlarfg( 2, beta, v( 2, m22 ), 1, v( 1, m22 ) )
                  h( k+1, k ) = beta
                  h( k+2, k ) = zero
               END IF
 
*
*              ==== Perform update from right within 
*              .    computational window. ====
*
               DO 30 j = jtop, min( kbot, k+3 )
                  refsum = v( 1, m22 )*( h( j, k+1 )+v( 2, m22 )*
     $                     h( j, k+2 ) )
                  h( j, k+1 ) = h( j, k+1 ) - refsum
                  h( j, k+2 ) = h( j, k+2 ) - refsum*v( 2, m22 )
   30          CONTINUE
*
*              ==== Perform update from left within 
*              .    computational window. ====
*
               IF( accum ) THEN
                  jbot = min( ndcol, kbot )
               ELSE IF( wantt ) THEN
                  jbot = n
               ELSE
                  jbot = kbot
               END IF
               DO 40 j = k+1, jbot
                  refsum = v( 1, m22 )*( h( k+1, j )+v( 2, m22 )*
     $                     h( k+2, j ) )
                  h( k+1, j ) = h( k+1, j ) - refsum
                  h( k+2, j ) = h( k+2, j ) - refsum*v( 2, m22 )
   40          CONTINUE
*
*              ==== The following convergence test requires that
*              .    the tradition small-compared-to-nearby-diagonals
*              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
*              .    criteria both be satisfied.  The latter improves
*              .    accuracy in some examples. Falling back on an
*              .    alternate convergence criterion when TST1 or TST2
*              .    is zero (as done here) is traditional but probably
*              .    unnecessary. ====
*
               IF( k.GE.ktop ) THEN
                  IF( h( k+1, k ).NE.zero ) THEN
                     tst1 = abs( h( k, k ) ) + abs( h( k+1, k+1 ) )
                     IF( tst1.EQ.zero ) THEN
                        IF( k.GE.ktop+1 )
     $                     tst1 = tst1 + abs( h( k, k-1 ) )
                        IF( k.GE.ktop+2 )
     $                     tst1 = tst1 + abs( h( k, k-2 ) )
                        IF( k.GE.ktop+3 )
     $                     tst1 = tst1 + abs( h( k, k-3 ) )
                        IF( k.LE.kbot-2 )
     $                     tst1 = tst1 + abs( h( k+2, k+1 ) )
                        IF( k.LE.kbot-3 )
     $                     tst1 = tst1 + abs( h( k+3, k+1 ) )
                        IF( k.LE.kbot-4 )
     $                     tst1 = tst1 + abs( h( k+4, k+1 ) )
                     END IF
                     IF( abs( h( k+1, k ) )
     $                   .LE.max( smlnum, ulp*tst1 ) ) THEN
                        h12 = max( abs( h( k+1, k ) ),
     $                             abs( h( k, k+1 ) ) )
                        h21 = min( abs( h( k+1, k ) ),
     $                             abs( h( k, k+1 ) ) )
                        h11 = max( abs( h( k+1, k+1 ) ),
     $                             abs( h( k, k )-h( k+1, k+1 ) ) )
                        h22 = min( abs( h( k+1, k+1 ) ),
     $                        abs( h( k, k )-h( k+1, k+1 ) ) )
                        scl = h11 + h12
                        tst2 = h22*( h11 / scl )
*
                        IF( tst2.EQ.zero .OR. h21*( h12 / scl ).LE.
     $                      max( smlnum, ulp*tst2 ) ) THEN
                           h( k+1, k ) = zero
                        END IF
                     END IF
                  END IF
               END IF
*
*              ==== Accumulate orthogonal transformations. ====
*
               IF( accum ) THEN
                  kms = k - incol
                  DO 50 j = max( 1, ktop-incol ), kdu
                     refsum = v( 1, m22 )*( u( j, kms+1 )+
     $                        v( 2, m22 )*u( j, kms+2 ) )
                     u( j, kms+1 ) = u( j, kms+1 ) - refsum
                     u( j, kms+2 ) = u( j, kms+2 ) - refsum*v( 2, m22 )
  50                 CONTINUE
               ELSE IF( wantz ) THEN
                  DO 60 j = iloz, ihiz
                     refsum = v( 1, m22 )*( z( j, k+1 )+v( 2, m22 )*
     $                        z( j, k+2 ) )
                     z( j, k+1 ) = z( j, k+1 ) - refsum
                     z( j, k+2 ) = z( j, k+2 ) - refsum*v( 2, m22 )
  60              CONTINUE
               END IF
            END IF
*
*           ==== Normal case: Chain of 3-by-3 reflections ====
*
            DO 80 m = mbot, mtop, -1
               k = krcol + 2*( m-1 )
               IF( k.EQ.ktop-1 ) THEN
                  CALL dlaqr1( 3, h( ktop, ktop ), ldh, sr( 2*m-1 ),
     $                         si( 2*m-1 ), sr( 2*m ), si( 2*m ),
     $                         v( 1, m ) )
                  alpha = v( 1, m )
                  CALL dlarfg( 3, alpha, v( 2, m ), 1, v( 1, m ) )
               ELSE
*
*                 ==== Perform delayed transformation of row below
*                 .    Mth bulge. Exploit fact that first two elements
*                 .    of row are actually zero. ====
*
                  refsum = v( 1, m )*v( 3, m )*h( k+3, k+2 )
                  h( k+3, k   ) = -refsum
                  h( k+3, k+1 ) = -refsum*v( 2, m )
                  h( k+3, k+2 ) = h( k+3, k+2 ) - refsum*v( 3, m )
*
*                 ==== Calculate reflection to move
*                 .    Mth bulge one step. ====
*
                  beta      = h( k+1, k )
                  v( 2, m ) = h( k+2, k )
                  v( 3, m ) = h( k+3, k )
                  CALL dlarfg( 3, beta, v( 2, m ), 1, v( 1, m ) )
*
*                 ==== A Bulge may collapse because of vigilant
*                 .    deflation or destructive underflow.  In the
*                 .    underflow case, try the two-small-subdiagonals
*                 .    trick to try to reinflate the bulge.  ====
*
                  IF( h( k+3, k ).NE.zero .OR. h( k+3, k+1 ).NE.
     $                zero .OR. h( k+3, k+2 ).EQ.zero ) THEN
*
*                    ==== Typical case: not collapsed (yet). ====
*
                     h( k+1, k ) = beta
                     h( k+2, k ) = zero
                     h( k+3, k ) = zero
                  ELSE
*
*                    ==== Atypical case: collapsed.  Attempt to
*                    .    reintroduce ignoring H(K+1,K) and H(K+2,K).
*                    .    If the fill resulting from the new
*                    .    reflector is too large, then abandon it.
*                    .    Otherwise, use the new one. ====
*
                     CALL dlaqr1( 3, h( k+1, k+1 ), ldh, sr( 2*m-1 ),
     $                            si( 2*m-1 ), sr( 2*m ), si( 2*m ),
     $                            vt )
                     alpha = vt( 1 )
                     CALL dlarfg( 3, alpha, vt( 2 ), 1, vt( 1 ) )
                     refsum = vt( 1 )*( h( k+1, k )+vt( 2 )*
     $                        h( k+2, k ) )
*
                     IF( abs( h( k+2, k )-refsum*vt( 2 ) )+
     $                   abs( refsum*vt( 3 ) ).GT.ulp*
     $                   ( abs( h( k, k ) )+abs( h( k+1,
     $                   k+1 ) )+abs( h( k+2, k+2 ) ) ) ) THEN
*
*                       ==== Starting a new bulge here would
*                       .    create non-negligible fill.  Use
*                       .    the old one with trepidation. ====
*
                        h( k+1, k ) = beta
                        h( k+2, k ) = zero
                        h( k+3, k ) = zero
                     ELSE
*
*                       ==== Starting a new bulge here would
*                       .    create only negligible fill.
*                       .    Replace the old reflector with
*                       .    the new one. ====
*
                        h( k+1, k ) = h( k+1, k ) - refsum
                        h( k+2, k ) = zero
                        h( k+3, k ) = zero
                        v( 1, m ) = vt( 1 )
                        v( 2, m ) = vt( 2 )
                        v( 3, m ) = vt( 3 )
                     END IF
                  END IF
               END IF
*
*              ====  Apply reflection from the right and
*              .     the first column of update from the left.
*              .     These updates are required for the vigilant
*              .     deflation check. We still delay most of the
*              .     updates from the left for efficiency. ====      
*
               DO 70 j = jtop, min( kbot, k+3 )
                  refsum = v( 1, m )*( h( j, k+1 )+v( 2, m )*
     $                     h( j, k+2 )+v( 3, m )*h( j, k+3 ) )
                  h( j, k+1 ) = h( j, k+1 ) - refsum
                  h( j, k+2 ) = h( j, k+2 ) - refsum*v( 2, m )
                  h( j, k+3 ) = h( j, k+3 ) - refsum*v( 3, m )
   70          CONTINUE
*
*              ==== Perform update from left for subsequent
*              .    column. ====
*
               refsum = v( 1, m )*( h( k+1, k+1 )+v( 2, m )*
     $                  h( k+2, k+1 )+v( 3, m )*h( k+3, k+1 ) )
               h( k+1, k+1 ) = h( k+1, k+1 ) - refsum
               h( k+2, k+1 ) = h( k+2, k+1 ) - refsum*v( 2, m )
               h( k+3, k+1 ) = h( k+3, k+1 ) - refsum*v( 3, m )
*
*              ==== The following convergence test requires that
*              .    the tradition small-compared-to-nearby-diagonals
*              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
*              .    criteria both be satisfied.  The latter improves
*              .    accuracy in some examples. Falling back on an
*              .    alternate convergence criterion when TST1 or TST2
*              .    is zero (as done here) is traditional but probably
*              .    unnecessary. ====
*
               IF( k.LT.ktop)
     $              cycle
               IF( h( k+1, k ).NE.zero ) THEN
                  tst1 = abs( h( k, k ) ) + abs( h( k+1, k+1 ) )
                  IF( tst1.EQ.zero ) THEN
                     IF( k.GE.ktop+1 )
     $                  tst1 = tst1 + abs( h( k, k-1 ) )
                     IF( k.GE.ktop+2 )
     $                  tst1 = tst1 + abs( h( k, k-2 ) )
                     IF( k.GE.ktop+3 )
     $                  tst1 = tst1 + abs( h( k, k-3 ) )
                     IF( k.LE.kbot-2 )
     $                  tst1 = tst1 + abs( h( k+2, k+1 ) )
                     IF( k.LE.kbot-3 )
     $                  tst1 = tst1 + abs( h( k+3, k+1 ) )
                     IF( k.LE.kbot-4 )
     $                  tst1 = tst1 + abs( h( k+4, k+1 ) )
                  END IF
                  IF( abs( h( k+1, k ) ).LE.max( smlnum, ulp*tst1 ) )
     $                 THEN
                     h12 = max( abs( h( k+1, k ) ), abs( h( k, k+1 ) ) )
                     h21 = min( abs( h( k+1, k ) ), abs( h( k, k+1 ) ) )
                     h11 = max( abs( h( k+1, k+1 ) ),
     $                     abs( h( k, k )-h( k+1, k+1 ) ) )
                     h22 = min( abs( h( k+1, k+1 ) ),
     $                     abs( h( k, k )-h( k+1, k+1 ) ) )
                     scl = h11 + h12
                     tst2 = h22*( h11 / scl )
*
                     IF( tst2.EQ.zero .OR. h21*( h12 / scl ).LE.
     $                   max( smlnum, ulp*tst2 ) ) THEN
                        h( k+1, k ) = zero
                     END IF
                  END IF
               END IF
   80       CONTINUE
*
*           ==== Multiply H by reflections from the left ====
*
            IF( accum ) THEN
               jbot = min( ndcol, kbot )
            ELSE IF( wantt ) THEN
               jbot = n
            ELSE
               jbot = kbot
            END IF
*
            DO 100 m = mbot, mtop, -1
               k = krcol + 2*( m-1 )
               DO 90 j = max( ktop, krcol + 2*m ), jbot
                  refsum = v( 1, m )*( h( k+1, j )+v( 2, m )*
     $                     h( k+2, j )+v( 3, m )*h( k+3, j ) )
                  h( k+1, j ) = h( k+1, j ) - refsum
                  h( k+2, j ) = h( k+2, j ) - refsum*v( 2, m )
                  h( k+3, j ) = h( k+3, j ) - refsum*v( 3, m )
   90          CONTINUE
  100       CONTINUE
*
*           ==== Accumulate orthogonal transformations. ====
*
            IF( accum ) THEN
*
*              ==== Accumulate U. (If needed, update Z later
*              .    with an efficient matrix-matrix
*              .    multiply.) ====
*
               DO 120 m = mbot, mtop, -1
                  k = krcol + 2*( m-1 )
                  kms = k - incol
                  i2 = max( 1, ktop-incol )
                  i2 = max( i2, kms-(krcol-incol)+1 )
                  i4 = min( kdu, krcol + 2*( mbot-1 ) - incol + 5 )
                  DO 110 j = i2, i4
                     refsum = v( 1, m )*( u( j, kms+1 )+v( 2, m )*
     $                        u( j, kms+2 )+v( 3, m )*u( j, kms+3 ) )
                     u( j, kms+1 ) = u( j, kms+1 ) - refsum
                     u( j, kms+2 ) = u( j, kms+2 ) - refsum*v( 2, m )
                     u( j, kms+3 ) = u( j, kms+3 ) - refsum*v( 3, m )
  110             CONTINUE
  120          CONTINUE
            ELSE IF( wantz ) THEN
*
*              ==== U is not accumulated, so update Z
*              .    now by multiplying by reflections
*              .    from the right. ====
*
               DO 140 m = mbot, mtop, -1
                  k = krcol + 2*( m-1 )
                  DO 130 j = iloz, ihiz
                     refsum = v( 1, m )*( z( j, k+1 )+v( 2, m )*
     $                        z( j, k+2 )+v( 3, m )*z( j, k+3 ) )
                     z( j, k+1 ) = z( j, k+1 ) - refsum
                     z( j, k+2 ) = z( j, k+2 ) - refsum*v( 2, m )
                     z( j, k+3 ) = z( j, k+3 ) - refsum*v( 3, m )
  130             CONTINUE
  140          CONTINUE
            END IF
*
*           ==== End of near-the-diagonal bulge chase. ====
*
  145    CONTINUE
*
*        ==== Use U (if accumulated) to update far-from-diagonal
*        .    entries in H.  If required, use U to update Z as
*        .    well. ====
*
         IF( accum ) THEN
            IF( wantt ) THEN
               jtop = 1
               jbot = n
            ELSE
               jtop = ktop
               jbot = kbot
            END IF
            k1 = max( 1, ktop-incol )
            nu = ( kdu-max( 0, ndcol-kbot ) ) - k1 + 1
*
*           ==== Horizontal Multiply ====
*
            DO 150 jcol = min( ndcol, kbot ) + 1, jbot, nh
               jlen = min( nh, jbot-jcol+1 )
               CALL dgemm( 'C', 'N', nu, jlen, nu, one, u( k1, k1 ),
     $                        ldu, h( incol+k1, jcol ), ldh, zero, wh,
     $                        ldwh )
               CALL dlacpy( 'ALL', nu, jlen, wh, ldwh,
     $                         h( incol+k1, jcol ), ldh )
  150       CONTINUE
*
*           ==== Vertical multiply ====
*
            DO 160 jrow = jtop, max( ktop, incol ) - 1, nv
               jlen = min( nv, max( ktop, incol )-jrow )
               CALL dgemm( 'N', 'N', jlen, nu, nu, one,
     $                     h( jrow, incol+k1 ), ldh, u( k1, k1 ),
     $                     ldu, zero, wv, ldwv )
               CALL dlacpy( 'ALL', jlen, nu, wv, ldwv,
     $                      h( jrow, incol+k1 ), ldh )
  160       CONTINUE
*
*           ==== Z multiply (also vertical) ====
*
            IF( wantz ) THEN
               DO 170 jrow = iloz, ihiz, nv
                  jlen = min( nv, ihiz-jrow+1 )
                  CALL dgemm( 'N', 'N', jlen, nu, nu, one,
     $                        z( jrow, incol+k1 ), ldz, u( k1, k1 ),
     $                        ldu, zero, wv, ldwv )
                  CALL dlacpy( 'ALL', jlen, nu, wv, ldwv,
     $                         z( jrow, incol+k1 ), ldz )
  170          CONTINUE
            END IF
         END IF
  180 CONTINUE
*
*     ==== End of DLAQR5 ====
*
      END
      LOGICAL FUNCTION dlaisnan( DIN1, DIN2 )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: din1, din2
*     ..
*
*  =====================================================================
*
*  .. Executable Statements ..
      dlaisnan = (din1.NE.din2)
      RETURN
      END
      SUBROUTINE dlaqr1( N, H, LDH, SR1, SI1, SR2, SI2, V )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   SI1, SI2, SR1, SR2
      INTEGER            LDH, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), V( * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   H21S, H31S, S
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( n.NE.2 .AND. n.NE.3 ) THEN
         RETURN
      END IF
*
      IF( n.EQ.2 ) THEN
         s = abs( h( 1, 1 )-sr2 ) + abs( si2 ) + abs( h( 2, 1 ) )
         IF( s.EQ.zero ) THEN
            v( 1 ) = zero
            v( 2 ) = zero
         ELSE
            h21s = h( 2, 1 ) / s
            v( 1 ) = h21s*h( 1, 2 ) + ( h( 1, 1 )-sr1 )*
     $               ( ( h( 1, 1 )-sr2 ) / s ) - si1*( si2 / s )
            v( 2 ) = h21s*( h( 1, 1 )+h( 2, 2 )-sr1-sr2 )
         END IF
      ELSE
         s = abs( h( 1, 1 )-sr2 ) + abs( si2 ) + abs( h( 2, 1 ) ) +
     $       abs( h( 3, 1 ) )
         IF( s.EQ.zero ) THEN
            v( 1 ) = zero
            v( 2 ) = zero
            v( 3 ) = zero
         ELSE
            h21s = h( 2, 1 ) / s
            h31s = h( 3, 1 ) / s
            v( 1 ) = ( h( 1, 1 )-sr1 )*( ( h( 1, 1 )-sr2 ) / s ) -
     $               si1*( si2 / s ) + h( 1, 2 )*h21s + h( 1, 3 )*h31s
            v( 2 ) = h21s*( h( 1, 1 )+h( 2, 2 )-sr1-sr2 ) +
     $               h( 2, 3 )*h31s
            v( 3 ) = h31s*( h( 1, 1 )+h( 3, 3 )-sr1-sr2 ) +
     $               h21s*h( 3, 2 )
         END IF
      END IF
      END
      SUBROUTINE dlaqr2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,
     $                   LDT, NV, WV, LDWV, WORK, LWORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
     $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),
     $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),
     $                   z( ldz, * )
*     ..
*
*  ================================================================
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S,
     $                   SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL,
     $                   KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2,
     $                   lwkopt
      LOGICAL            BULGE, SORTED
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           dcopy, dgehrd, dgemm, dlabad, dlacpy, dlahqr,
     $                   dlanv2, dlarf, dlarfg, dlaset, dormhr, dtrexc
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, int, max, min, sqrt
*     ..
*     .. Executable Statements ..
*
*     ==== Estimate optimal workspace. ====
*
      jw = min( nw, kbot-ktop+1 )
      IF( jw.LE.2 ) THEN
         lwkopt = 1
      ELSE
*
*        ==== Workspace query call to DGEHRD ====
*
         CALL dgehrd( jw, 1, jw-1, t, ldt, work, work, -1, info )
         lwk1 = int( work( 1 ) )
*
*        ==== Workspace query call to DORMHR ====
*
         CALL dormhr( 'R', 'N', jw, jw, 1, jw-1, t, ldt, work, v, ldv,
     $                work, -1, info )
         lwk2 = int( work( 1 ) )
*
*        ==== Optimal workspace ====
*
         lwkopt = jw + max( lwk1, lwk2 )
      END IF
*
*     ==== Quick return in case of workspace query. ====
*
      IF( lwork.EQ.-1 ) THEN
         work( 1 ) = dble( lwkopt )
         RETURN
      END IF
*
*     ==== Nothing to do ...
*     ... for an empty active block ... ====
      ns = 0
      nd = 0
      work( 1 ) = one
      IF( ktop.GT.kbot )
     $   RETURN
*     ... nor for an empty deflation window. ====
      IF( nw.LT.1 )
     $   RETURN
*
*     ==== Machine constants ====
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = one / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
*
*     ==== Setup deflation window ====
*
      jw = min( nw, kbot-ktop+1 )
      kwtop = kbot - jw + 1
      IF( kwtop.EQ.ktop ) THEN
         s = zero
      ELSE
         s = h( kwtop, kwtop-1 )
      END IF
*
      IF( kbot.EQ.kwtop ) THEN
*
*        ==== 1-by-1 deflation window: not much to do ====
*
         sr( kwtop ) = h( kwtop, kwtop )
         si( kwtop ) = zero
         ns = 1
         nd = 0
         IF( abs( s ).LE.max( smlnum, ulp*abs( h( kwtop, kwtop ) ) ) )
     $        THEN
            ns = 0
            nd = 1
            IF( kwtop.GT.ktop )
     $         h( kwtop, kwtop-1 ) = zero
         END IF
         work( 1 ) = one
         RETURN
      END IF
*
*     ==== Convert to spike-triangular form.  (In case of a
*     .    rare QR failure, this routine continues to do
*     .    aggressive early deflation using that part of
*     .    the deflation window that converged using INFQR
*     .    here and there to keep track.) ====
*
      CALL dlacpy( 'U', jw, jw, h( kwtop, kwtop ), ldh, t, ldt )
      CALL dcopy( jw-1, h( kwtop+1, kwtop ), ldh+1, t( 2, 1 ), ldt+1 )
*
      CALL dlaset( 'A', jw, jw, zero, one, v, ldv )
      CALL dlahqr( .true., .true., jw, 1, jw, t, ldt, sr( kwtop ),
     $             si( kwtop ), 1, jw, v, ldv, infqr )
*
*     ==== DTREXC needs a clean margin near the diagonal ====
*
      DO 10 j = 1, jw - 3
         t( j+2, j ) = zero
         t( j+3, j ) = zero
   10 CONTINUE
      IF( jw.GT.2 )
     $   t( jw, jw-2 ) = zero
*
*     ==== Deflation detection loop ====
*
      ns = jw
      ilst = infqr + 1
   20 CONTINUE
      IF( ilst.LE.ns ) THEN
         IF( ns.EQ.1 ) THEN
            bulge = .false.
         ELSE
            bulge = t( ns, ns-1 ).NE.zero
         END IF
*
*        ==== Small spike tip test for deflation ====
*
         IF( .NOT.bulge ) THEN
*
*           ==== Real eigenvalue ====
*
            foo = abs( t( ns, ns ) )
            IF( foo.EQ.zero )
     $         foo = abs( s )
            IF( abs( s*v( 1, ns ) ).LE.max( smlnum, ulp*foo ) ) THEN
*
*              ==== Deflatable ====
*
               ns = ns - 1
            ELSE
*
*              ==== Undeflatable.   Move it up out of the way.
*              .    (DTREXC can not fail in this case.) ====
*
               ifst = ns
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, work,
     $                      info )
               ilst = ilst + 1
            END IF
         ELSE
*
*           ==== Complex conjugate pair ====
*
            foo = abs( t( ns, ns ) ) + sqrt( abs( t( ns, ns-1 ) ) )*
     $            sqrt( abs( t( ns-1, ns ) ) )
            IF( foo.EQ.zero )
     $         foo = abs( s )
            IF( max( abs( s*v( 1, ns ) ), abs( s*v( 1, ns-1 ) ) ).LE.
     $          max( smlnum, ulp*foo ) ) THEN
*
*              ==== Deflatable ====
*
               ns = ns - 2
            ELSE
*
*              ==== Undeflatable. Move them up out of the way.
*              .    Fortunately, DTREXC does the right thing with
*              .    ILST in case of a rare exchange failure. ====
*
               ifst = ns
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, work,
     $                      info )
               ilst = ilst + 2
            END IF
         END IF
*
*        ==== End deflation detection loop ====
*
         GO TO 20
      END IF
*
*        ==== Return to Hessenberg form ====
*
      IF( ns.EQ.0 )
     $   s = zero
*
      IF( ns.LT.jw ) THEN
*
*        ==== sorting diagonal blocks of T improves accuracy for
*        .    graded matrices.  Bubble sort deals well with
*        .    exchange failures. ====
*
         sorted = .false.
         i = ns + 1
   30    CONTINUE
         IF( sorted )
     $      GO TO 50
         sorted = .true.
*
         kend = i - 1
         i = infqr + 1
         IF( i.EQ.ns ) THEN
            k = i + 1
         ELSE IF( t( i+1, i ).EQ.zero ) THEN
            k = i + 1
         ELSE
            k = i + 2
         END IF
   40    CONTINUE
         IF( k.LE.kend ) THEN
            IF( k.EQ.i+1 ) THEN
               evi = abs( t( i, i ) )
            ELSE
               evi = abs( t( i, i ) ) + sqrt( abs( t( i+1, i ) ) )*
     $               sqrt( abs( t( i, i+1 ) ) )
            END IF
*
            IF( k.EQ.kend ) THEN
               evk = abs( t( k, k ) )
            ELSE IF( t( k+1, k ).EQ.zero ) THEN
               evk = abs( t( k, k ) )
            ELSE
               evk = abs( t( k, k ) ) + sqrt( abs( t( k+1, k ) ) )*
     $               sqrt( abs( t( k, k+1 ) ) )
            END IF
*
            IF( evi.GE.evk ) THEN
               i = k
            ELSE
               sorted = .false.
               ifst = i
               ilst = k
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, work,
     $                      info )
               IF( info.EQ.0 ) THEN
                  i = ilst
               ELSE
                  i = k
               END IF
            END IF
            IF( i.EQ.kend ) THEN
               k = i + 1
            ELSE IF( t( i+1, i ).EQ.zero ) THEN
               k = i + 1
            ELSE
               k = i + 2
            END IF
            GO TO 40
         END IF
         GO TO 30
   50    CONTINUE
      END IF
*
*     ==== Restore shift/eigenvalue array from T ====
*
      i = jw
   60 CONTINUE
      IF( i.GE.infqr+1 ) THEN
         IF( i.EQ.infqr+1 ) THEN
            sr( kwtop+i-1 ) = t( i, i )
            si( kwtop+i-1 ) = zero
            i = i - 1
         ELSE IF( t( i, i-1 ).EQ.zero ) THEN
            sr( kwtop+i-1 ) = t( i, i )
            si( kwtop+i-1 ) = zero
            i = i - 1
         ELSE
            aa = t( i-1, i-1 )
            cc = t( i, i-1 )
            bb = t( i-1, i )
            dd = t( i, i )
            CALL dlanv2( aa, bb, cc, dd, sr( kwtop+i-2 ),
     $                   si( kwtop+i-2 ), sr( kwtop+i-1 ),
     $                   si( kwtop+i-1 ), cs, sn )
            i = i - 2
         END IF
         GO TO 60
      END IF
*
      IF( ns.LT.jw .OR. s.EQ.zero ) THEN
         IF( ns.GT.1 .AND. s.NE.zero ) THEN
*
*           ==== Reflect spike back into lower triangle ====
*
            CALL dcopy( ns, v, ldv, work, 1 )
            beta = work( 1 )
            CALL dlarfg( ns, beta, work( 2 ), 1, tau )
            work( 1 ) = one
*
            CALL dlaset( 'L', jw-2, jw-2, zero, zero, t( 3, 1 ), ldt )
*
            CALL dlarf( 'L', ns, jw, work, 1, tau, t, ldt,
     $                  work( jw+1 ) )
            CALL dlarf( 'R', ns, ns, work, 1, tau, t, ldt,
     $                  work( jw+1 ) )
            CALL dlarf( 'R', jw, ns, work, 1, tau, v, ldv,
     $                  work( jw+1 ) )
*
            CALL dgehrd( jw, 1, ns, t, ldt, work, work( jw+1 ),
     $                   lwork-jw, info )
         END IF
*
*        ==== Copy updated reduced window into place ====
*
         IF( kwtop.GT.1 )
     $      h( kwtop, kwtop-1 ) = s*v( 1, 1 )
         CALL dlacpy( 'U', jw, jw, t, ldt, h( kwtop, kwtop ), ldh )
         CALL dcopy( jw-1, t( 2, 1 ), ldt+1, h( kwtop+1, kwtop ),
     $               ldh+1 )
*
*        ==== Accumulate orthogonal matrix in order update
*        .    H and Z, if requested.  ====
*
         IF( ns.GT.1 .AND. s.NE.zero )
     $      CALL dormhr( 'R', 'N', jw, ns, 1, ns, t, ldt, work, v, ldv,
     $                   work( jw+1 ), lwork-jw, info )
*
*        ==== Update vertical slab in H ====
*
         IF( wantt ) THEN
            ltop = 1
         ELSE
            ltop = ktop
         END IF
         DO 70 krow = ltop, kwtop - 1, nv
            kln = min( nv, kwtop-krow )
            CALL dgemm( 'N', 'N', kln, jw, jw, one, h( krow, kwtop ),
     $                  ldh, v, ldv, zero, wv, ldwv )
            CALL dlacpy( 'A', kln, jw, wv, ldwv, h( krow, kwtop ), ldh )
   70    CONTINUE
*
*        ==== Update horizontal slab in H ====
*
         IF( wantt ) THEN
            DO 80 kcol = kbot + 1, n, nh
               kln = min( nh, n-kcol+1 )
               CALL dgemm( 'C', 'N', jw, kln, jw, one, v, ldv,
     $                     h( kwtop, kcol ), ldh, zero, t, ldt )
               CALL dlacpy( 'A', jw, kln, t, ldt, h( kwtop, kcol ),
     $                      ldh )
   80       CONTINUE
         END IF
*
*        ==== Update vertical slab in Z ====
*
         IF( wantz ) THEN
            DO 90 krow = iloz, ihiz, nv
               kln = min( nv, ihiz-krow+1 )
               CALL dgemm( 'N', 'N', kln, jw, jw, one, z( krow, kwtop ),
     $                     ldz, v, ldv, zero, wv, ldwv )
               CALL dlacpy( 'A', kln, jw, wv, ldwv, z( krow, kwtop ),
     $                      ldz )
   90       CONTINUE
         END IF
      END IF
*
*     ==== Return the number of deflations ... ====
*
      nd = jw - ns
*
*     ==== ... and the number of shifts. (Subtracting
*     .    INFQR from the spike length takes care
*     .    of the case of a rare QR failure while
*     .    calculating eigenvalues of the deflation
*     .    window.)  ====
*
      ns = ns - infqr
*
*      ==== Return optimal workspace. ====
*
      work( 1 ) = dble( lwkopt )
*
*     ==== End of DLAQR2 ====
*
      END
      SUBROUTINE dtrexc( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            HERE, NBF, NBL, NBNEXT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlaexc, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input arguments.
*
      info = 0
      wantq = lsame( compq, 'V' )
      IF( .NOT.wantq .AND. .NOT.lsame( compq, 'N' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( ldt.LT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ldq.LT.1 .OR. ( wantq .AND. ldq.LT.max( 1, n ) ) ) THEN
         info = -6
      ELSE IF(( ifst.LT.1 .OR. ifst.GT.n ).AND.( n.GT.0 )) THEN
         info = -7
      ELSE IF(( ilst.LT.1 .OR. ilst.GT.n ).AND.( n.GT.0 )) THEN
         info = -8
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DTREXC', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.LE.1 )
     $   RETURN
*
*     Determine the first row of specified block
*     and find out it is 1 by 1 or 2 by 2.
*
      IF( ifst.GT.1 ) THEN
         IF( t( ifst, ifst-1 ).NE.zero )
     $      ifst = ifst - 1
      END IF
      nbf = 1
      IF( ifst.LT.n ) THEN
         IF( t( ifst+1, ifst ).NE.zero )
     $      nbf = 2
      END IF
*
*     Determine the first row of the final block
*     and find out it is 1 by 1 or 2 by 2.
*
      IF( ilst.GT.1 ) THEN
         IF( t( ilst, ilst-1 ).NE.zero )
     $      ilst = ilst - 1
      END IF
      nbl = 1
      IF( ilst.LT.n ) THEN
         IF( t( ilst+1, ilst ).NE.zero )
     $      nbl = 2
      END IF
*
      IF( ifst.EQ.ilst )
     $   RETURN
*
      IF( ifst.LT.ilst ) THEN
*
*        Update ILST
*
         IF( nbf.EQ.2 .AND. nbl.EQ.1 )
     $      ilst = ilst - 1
         IF( nbf.EQ.1 .AND. nbl.EQ.2 )
     $      ilst = ilst + 1
*
         here = ifst
*
   10    CONTINUE
*
*        Swap block with next one below
*
         IF( nbf.EQ.1 .OR. nbf.EQ.2 ) THEN
*
*           Current block either 1 by 1 or 2 by 2
*
            nbnext = 1
            IF( here+nbf+1.LE.n ) THEN
               IF( t( here+nbf+1, here+nbf ).NE.zero )
     $            nbnext = 2
            END IF
            CALL dlaexc( wantq, n, t, ldt, q, ldq, here, nbf, nbnext,
     $                   work, info )
            IF( info.NE.0 ) THEN
               ilst = here
               RETURN
            END IF
            here = here + nbnext
*
*           Test if 2 by 2 block breaks into two 1 by 1 blocks
*
            IF( nbf.EQ.2 ) THEN
               IF( t( here+1, here ).EQ.zero )
     $            nbf = 3
            END IF
*
         ELSE
*
*           Current block consists of two 1 by 1 blocks each of which
*           must be swapped individually
*
            nbnext = 1
            IF( here+3.LE.n ) THEN
               IF( t( here+3, here+2 ).NE.zero )
     $            nbnext = 2
            END IF
            CALL dlaexc( wantq, n, t, ldt, q, ldq, here+1, 1, nbnext,
     $                   work, info )
            IF( info.NE.0 ) THEN
               ilst = here
               RETURN
            END IF
            IF( nbnext.EQ.1 ) THEN
*
*              Swap two 1 by 1 blocks, no problems possible
*
               CALL dlaexc( wantq, n, t, ldt, q, ldq, here, 1, nbnext,
     $                      work, info )
               here = here + 1
            ELSE
*
*              Recompute NBNEXT in case 2 by 2 split
*
               IF( t( here+2, here+1 ).EQ.zero )
     $            nbnext = 1
               IF( nbnext.EQ.2 ) THEN
*
*                 2 by 2 Block did not split
*
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here, 1,
     $                         nbnext, work, info )
                  IF( info.NE.0 ) THEN
                     ilst = here
                     RETURN
                  END IF
                  here = here + 2
               ELSE
*
*                 2 by 2 Block did split
*
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here, 1, 1,
     $                         work, info )
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here+1, 1, 1,
     $                         work, info )
                  here = here + 2
               END IF
            END IF
         END IF
         IF( here.LT.ilst )
     $      GO TO 10
*
      ELSE
*
         here = ifst
   20    CONTINUE
*
*        Swap block with next one above
*
         IF( nbf.EQ.1 .OR. nbf.EQ.2 ) THEN
*
*           Current block either 1 by 1 or 2 by 2
*
            nbnext = 1
            IF( here.GE.3 ) THEN
               IF( t( here-1, here-2 ).NE.zero )
     $            nbnext = 2
            END IF
            CALL dlaexc( wantq, n, t, ldt, q, ldq, here-nbnext, nbnext,
     $                   nbf, work, info )
            IF( info.NE.0 ) THEN
               ilst = here
               RETURN
            END IF
            here = here - nbnext
*
*           Test if 2 by 2 block breaks into two 1 by 1 blocks
*
            IF( nbf.EQ.2 ) THEN
               IF( t( here+1, here ).EQ.zero )
     $            nbf = 3
            END IF
*
         ELSE
*
*           Current block consists of two 1 by 1 blocks each of which
*           must be swapped individually
*
            nbnext = 1
            IF( here.GE.3 ) THEN
               IF( t( here-1, here-2 ).NE.zero )
     $            nbnext = 2
            END IF
            CALL dlaexc( wantq, n, t, ldt, q, ldq, here-nbnext, nbnext,
     $                   1, work, info )
            IF( info.NE.0 ) THEN
               ilst = here
               RETURN
            END IF
            IF( nbnext.EQ.1 ) THEN
*
*              Swap two 1 by 1 blocks, no problems possible
*
               CALL dlaexc( wantq, n, t, ldt, q, ldq, here, nbnext, 1,
     $                      work, info )
               here = here - 1
            ELSE
*
*              Recompute NBNEXT in case 2 by 2 split
*
               IF( t( here, here-1 ).EQ.zero )
     $            nbnext = 1
               IF( nbnext.EQ.2 ) THEN
*
*                 2 by 2 Block did not split
*
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here-1, 2, 1,
     $                         work, info )
                  IF( info.NE.0 ) THEN
                     ilst = here
                     RETURN
                  END IF
                  here = here - 2
               ELSE
*
*                 2 by 2 Block did split
*
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here, 1, 1,
     $                         work, info )
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here-1, 1, 1,
     $                         work, info )
                  here = here - 2
               END IF
            END IF
         END IF
         IF( here.GT.ilst )
     $      GO TO 20
      END IF
      ilst = here
*
      RETURN
*
*     End of DTREXC
*
      END
      SUBROUTINE dormhr( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
     $                   LDC, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           lsame, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           dormqr, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      nh = ihi - ilo
      left = lsame( side, 'L' )
      lquery = ( lwork.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( left ) THEN
         nq = m
         nw = max( 1, n )
      ELSE
         nq = n
         nw = max( 1, m )
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.lsame( trans, 'N' ) .AND. .NOT.lsame( trans, 'T' ) )
     $          THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, nq ) ) THEN
         info = -5
      ELSE IF( ihi.LT.min( ilo, nq ) .OR. ihi.GT.nq ) THEN
         info = -6
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -8
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -11
      ELSE IF( lwork.LT.nw .AND. .NOT.lquery ) THEN
         info = -13
      END IF
*
      IF( info.EQ.0 ) THEN
         IF( left ) THEN
            nb = ilaenv( 1, 'DORMQR', side // trans, nh, n, nh, -1 )
         ELSE
            nb = ilaenv( 1, 'DORMQR', side // trans, m, nh, nh, -1 )
         END IF
         lwkopt = nw*nb
         work( 1 ) = lwkopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORMHR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 .OR. nh.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
      IF( left ) THEN
         mi = nh
         ni = n
         i1 = ilo + 1
         i2 = 1
      ELSE
         mi = m
         ni = nh
         i1 = 1
         i2 = ilo + 1
      END IF
*
      CALL dormqr( side, trans, mi, ni, nh, a( ilo+1, ilo ), lda,
     $             tau( ilo ), c( i1, i2 ), ldc, work, lwork, iinfo )
*
      work( 1 ) = lwkopt
      RETURN
*
*     End of DORMHR
*
      END
      SUBROUTINE dlaexc( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,
     $                   INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ
      INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      DOUBLE PRECISION   TEN
      parameter( ten = 1.0d+1 )
      INTEGER            LDD, LDX
      parameter( ldd = 4, ldx = 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IERR, J2, J3, J4, K, ND
      DOUBLE PRECISION   CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22,
     $                   t33, tau, tau1, tau2, temp, thresh, wi1, wi2,
     $                   wr1, wr2, xnorm
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   D( LDD, 4 ), U( 3 ), U1( 3 ), U2( 3 ),
     $                   x( ldx, 2 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           dlamch, dlange
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlacpy, dlanv2, dlarfg, dlarfx, dlartg, dlasy2,
     $                   drot
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max
*     ..
*     .. Executable Statements ..
*
      info = 0
*
*     Quick return if possible
*
      IF( n.EQ.0 .OR. n1.EQ.0 .OR. n2.EQ.0 )
     $   RETURN
      IF( j1+n1.GT.n )
     $   RETURN
*
      j2 = j1 + 1
      j3 = j1 + 2
      j4 = j1 + 3
*
      IF( n1.EQ.1 .AND. n2.EQ.1 ) THEN
*
*        Swap two 1-by-1 blocks.
*
         t11 = t( j1, j1 )
         t22 = t( j2, j2 )
*
*        Determine the transformation to perform the interchange.
*
         CALL dlartg( t( j1, j2 ), t22-t11, cs, sn, temp )
*
*        Apply transformation to the matrix T.
*
         IF( j3.LE.n )
     $      CALL drot( n-j1-1, t( j1, j3 ), ldt, t( j2, j3 ), ldt, cs,
     $                 sn )
         CALL drot( j1-1, t( 1, j1 ), 1, t( 1, j2 ), 1, cs, sn )
*
         t( j1, j1 ) = t22
         t( j2, j2 ) = t11
*
         IF( wantq ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL drot( n, q( 1, j1 ), 1, q( 1, j2 ), 1, cs, sn )
         END IF
*
      ELSE
*
*        Swapping involves at least one 2-by-2 block.
*
*        Copy the diagonal block of order N1+N2 to the local array D
*        and compute its norm.
*
         nd = n1 + n2
         CALL dlacpy( 'Full', nd, nd, t( j1, j1 ), ldt, d, ldd )
         dnorm = dlange( 'Max', nd, nd, d, ldd, work )
*
*        Compute machine-dependent threshold for test for accepting
*        swap.
*
         eps = dlamch( 'P' )
         smlnum = dlamch( 'S' ) / eps
         thresh = max( ten*eps*dnorm, smlnum )
*
*        Solve T11*X - X*T22 = scale*T12 for X.
*
         CALL dlasy2( .false., .false., -1, n1, n2, d, ldd,
     $                d( n1+1, n1+1 ), ldd, d( 1, n1+1 ), ldd, scale, x,
     $                ldx, xnorm, ierr )
*
*        Swap the adjacent diagonal blocks.
*
         k = n1 + n1 + n2 - 3
         GO TO ( 10, 20, 30 )k
*
   10    CONTINUE
*
*        N1 = 1, N2 = 2: generate elementary reflector H so that:
*
*        ( scale, X11, X12 ) H = ( 0, 0, * )
*
         u( 1 ) = scale
         u( 2 ) = x( 1, 1 )
         u( 3 ) = x( 1, 2 )
         CALL dlarfg( 3, u( 3 ), u, 1, tau )
         u( 3 ) = one
         t11 = t( j1, j1 )
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL dlarfx( 'L', 3, 3, u, tau, d, ldd, work )
         CALL dlarfx( 'R', 3, 3, u, tau, d, ldd, work )
*
*        Test whether to reject swap.
*
         IF( max( abs( d( 3, 1 ) ), abs( d( 3, 2 ) ), abs( d( 3,
     $       3 )-t11 ) ).GT.thresh )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL dlarfx( 'L', 3, n-j1+1, u, tau, t( j1, j1 ), ldt, work )
         CALL dlarfx( 'R', j2, 3, u, tau, t( 1, j1 ), ldt, work )
*
         t( j3, j1 ) = zero
         t( j3, j2 ) = zero
         t( j3, j3 ) = t11
*
         IF( wantq ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL dlarfx( 'R', n, 3, u, tau, q( 1, j1 ), ldq, work )
         END IF
         GO TO 40
*
   20    CONTINUE
*
*        N1 = 2, N2 = 1: generate elementary reflector H so that:
*
*        H (  -X11 ) = ( * )
*          (  -X21 ) = ( 0 )
*          ( scale ) = ( 0 )
*
         u( 1 ) = -x( 1, 1 )
         u( 2 ) = -x( 2, 1 )
         u( 3 ) = scale
         CALL dlarfg( 3, u( 1 ), u( 2 ), 1, tau )
         u( 1 ) = one
         t33 = t( j3, j3 )
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL dlarfx( 'L', 3, 3, u, tau, d, ldd, work )
         CALL dlarfx( 'R', 3, 3, u, tau, d, ldd, work )
*
*        Test whether to reject swap.
*
         IF( max( abs( d( 2, 1 ) ), abs( d( 3, 1 ) ), abs( d( 1,
     $       1 )-t33 ) ).GT.thresh )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL dlarfx( 'R', j3, 3, u, tau, t( 1, j1 ), ldt, work )
         CALL dlarfx( 'L', 3, n-j1, u, tau, t( j1, j2 ), ldt, work )
*
         t( j1, j1 ) = t33
         t( j2, j1 ) = zero
         t( j3, j1 ) = zero
*
         IF( wantq ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL dlarfx( 'R', n, 3, u, tau, q( 1, j1 ), ldq, work )
         END IF
         GO TO 40
*
   30    CONTINUE
*
*        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
*        that:
*
*        H(2) H(1) (  -X11  -X12 ) = (  *  * )
*                  (  -X21  -X22 )   (  0  * )
*                  ( scale    0  )   (  0  0 )
*                  (    0  scale )   (  0  0 )
*
         u1( 1 ) = -x( 1, 1 )
         u1( 2 ) = -x( 2, 1 )
         u1( 3 ) = scale
         CALL dlarfg( 3, u1( 1 ), u1( 2 ), 1, tau1 )
         u1( 1 ) = one
*
         temp = -tau1*( x( 1, 2 )+u1( 2 )*x( 2, 2 ) )
         u2( 1 ) = -temp*u1( 2 ) - x( 2, 2 )
         u2( 2 ) = -temp*u1( 3 )
         u2( 3 ) = scale
         CALL dlarfg( 3, u2( 1 ), u2( 2 ), 1, tau2 )
         u2( 1 ) = one
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL dlarfx( 'L', 3, 4, u1, tau1, d, ldd, work )
         CALL dlarfx( 'R', 4, 3, u1, tau1, d, ldd, work )
         CALL dlarfx( 'L', 3, 4, u2, tau2, d( 2, 1 ), ldd, work )
         CALL dlarfx( 'R', 4, 3, u2, tau2, d( 1, 2 ), ldd, work )
*
*        Test whether to reject swap.
*
         IF( max( abs( d( 3, 1 ) ), abs( d( 3, 2 ) ), abs( d( 4, 1 ) ),
     $       abs( d( 4, 2 ) ) ).GT.thresh )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL dlarfx( 'L', 3, n-j1+1, u1, tau1, t( j1, j1 ), ldt, work )
         CALL dlarfx( 'R', j4, 3, u1, tau1, t( 1, j1 ), ldt, work )
         CALL dlarfx( 'L', 3, n-j1+1, u2, tau2, t( j2, j1 ), ldt, work )
         CALL dlarfx( 'R', j4, 3, u2, tau2, t( 1, j2 ), ldt, work )
*
         t( j3, j1 ) = zero
         t( j3, j2 ) = zero
         t( j4, j1 ) = zero
         t( j4, j2 ) = zero
*
         IF( wantq ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL dlarfx( 'R', n, 3, u1, tau1, q( 1, j1 ), ldq, work )
            CALL dlarfx( 'R', n, 3, u2, tau2, q( 1, j2 ), ldq, work )
         END IF
*
   40    CONTINUE
*
         IF( n2.EQ.2 ) THEN
*
*           Standardize new 2-by-2 block T11
*
            CALL dlanv2( t( j1, j1 ), t( j1, j2 ), t( j2, j1 ),
     $                   t( j2, j2 ), wr1, wi1, wr2, wi2, cs, sn )
            CALL drot( n-j1-1, t( j1, j1+2 ), ldt, t( j2, j1+2 ), ldt,
     $                 cs, sn )
            CALL drot( j1-1, t( 1, j1 ), 1, t( 1, j2 ), 1, cs, sn )
            IF( wantq )
     $         CALL drot( n, q( 1, j1 ), 1, q( 1, j2 ), 1, cs, sn )
         END IF
*
         IF( n1.EQ.2 ) THEN
*
*           Standardize new 2-by-2 block T22
*
            j3 = j1 + n2
            j4 = j3 + 1
            CALL dlanv2( t( j3, j3 ), t( j3, j4 ), t( j4, j3 ),
     $                   t( j4, j4 ), wr1, wi1, wr2, wi2, cs, sn )
            IF( j3+2.LE.n )
     $         CALL drot( n-j3-1, t( j3, j3+2 ), ldt, t( j4, j3+2 ),
     $                    ldt, cs, sn )
            CALL drot( j3-1, t( 1, j3 ), 1, t( 1, j4 ), 1, cs, sn )
            IF( wantq )
     $         CALL drot( n, q( 1, j3 ), 1, q( 1, j4 ), 1, cs, sn )
         END IF
*
      END IF
      RETURN
*
*     Exit with INFO = 1 if swap was rejected.
*
   50 CONTINUE
      info = 1
      RETURN
*
*     End of DLAEXC
*
      END
      SUBROUTINE dlahr2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  A( LDA, * ), T( LDT, NB ), TAU( NB ),
     $                   Y( LDY, NB )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      parameter( zero = 0.0d+0,
     $                     one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION  EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           daxpy, dcopy, dgemm, dgemv, dlacpy,
     $                   dlarfg, dscal, dtrmm, dtrmv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          min
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( n.LE.1 )
     $   RETURN
*
      DO 10 i = 1, nb
         IF( i.GT.1 ) THEN
*
*           Update A(K+1:N,I)
*
*           Update I-th column of A - Y * V**T
*
            CALL dgemv( 'NO TRANSPOSE', n-k, i-1, -one, y(k+1,1), ldy,
     $                  a( k+i-1, 1 ), lda, one, a( k+1, i ), 1 )
*
*           Apply I - V * T**T * V**T to this column (call it b) from the
*           left, using the last column of T as workspace
*
*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
*                    ( V2 )             ( b2 )
*
*           where V1 is unit lower triangular
*
*           w := V1**T * b1
*
            CALL dcopy( i-1, a( k+1, i ), 1, t( 1, nb ), 1 )
            CALL dtrmv( 'Lower', 'Transpose', 'UNIT',
     $                  i-1, a( k+1, 1 ),
     $                  lda, t( 1, nb ), 1 )
*
*           w := w + V2**T * b2
*
            CALL dgemv( 'Transpose', n-k-i+1, i-1,
     $                  one, a( k+i, 1 ),
     $                  lda, a( k+i, i ), 1, one, t( 1, nb ), 1 )
*
*           w := T**T * w
*
            CALL dtrmv( 'Upper', 'Transpose', 'NON-UNIT',
     $                  i-1, t, ldt,
     $                  t( 1, nb ), 1 )
*
*           b2 := b2 - V2*w
*
            CALL dgemv( 'NO TRANSPOSE', n-k-i+1, i-1, -one,
     $                  a( k+i, 1 ),
     $                  lda, t( 1, nb ), 1, one, a( k+i, i ), 1 )
*
*           b1 := b1 - V1*w
*
            CALL dtrmv( 'Lower', 'NO TRANSPOSE',
     $                  'UNIT', i-1,
     $                  a( k+1, 1 ), lda, t( 1, nb ), 1 )
            CALL daxpy( i-1, -one, t( 1, nb ), 1, a( k+1, i ), 1 )
*
            a( k+i-1, i-1 ) = ei
         END IF
*
*        Generate the elementary reflector H(I) to annihilate
*        A(K+I+1:N,I)
*
         CALL dlarfg( n-k-i+1, a( k+i, i ), a( min( k+i+1, n ), i ), 1,
     $                tau( i ) )
         ei = a( k+i, i )
         a( k+i, i ) = one
*
*        Compute  Y(K+1:N,I)
*
         CALL dgemv( 'NO TRANSPOSE', n-k, n-k-i+1,
     $               one, a( k+1, i+1 ),
     $               lda, a( k+i, i ), 1, zero, y( k+1, i ), 1 )
         CALL dgemv( 'Transpose', n-k-i+1, i-1,
     $               one, a( k+i, 1 ), lda,
     $               a( k+i, i ), 1, zero, t( 1, i ), 1 )
         CALL dgemv( 'NO TRANSPOSE', n-k, i-1, -one,
     $               y( k+1, 1 ), ldy,
     $               t( 1, i ), 1, one, y( k+1, i ), 1 )
         CALL dscal( n-k, tau( i ), y( k+1, i ), 1 )
*
*        Compute T(1:I,I)
*
         CALL dscal( i-1, -tau( i ), t( 1, i ), 1 )
         CALL dtrmv( 'Upper', 'No Transpose', 'NON-UNIT',
     $               i-1, t, ldt,
     $               t( 1, i ), 1 )
         t( i, i ) = tau( i )
*
   10 CONTINUE
      a( k+nb, nb ) = ei
*
*     Compute Y(1:K,1:NB)
*
      CALL dlacpy( 'ALL', k, nb, a( 1, 2 ), lda, y, ldy )
      CALL dtrmm( 'RIGHT', 'Lower', 'NO TRANSPOSE',
     $            'UNIT', k, nb,
     $            one, a( k+1, 1 ), lda, y, ldy )
      IF( n.GT.k+nb )
     $   CALL dgemm( 'NO TRANSPOSE', 'NO TRANSPOSE', k,
     $               nb, n-k-nb, one,
     $               a( 1, 2+nb ), lda, a( k+1+nb, 1 ), lda, one, y,
     $               ldy )
      CALL dtrmm( 'RIGHT', 'Upper', 'NO TRANSPOSE',
     $            'NON-UNIT', k, nb,
     $            one, t, ldt, y, ldy )
*
      RETURN
*
*     End of DLAHR2
*
      END
      SUBROUTINE dlarfx( SIDE, M, N, V, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            LDC, M, N
      DOUBLE PRECISION   TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9,
     $                   V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlarf
*     ..
*     .. Executable Statements ..
*
      IF( tau.EQ.zero )
     $   RETURN
      IF( lsame( side, 'L' ) ) THEN
*
*        Form  H * C, where H has order m.
*
         GO TO ( 10, 30, 50, 70, 90, 110, 130, 150,
     $           170, 190 )m
*
*        Code for general M
*
         CALL dlarf( side, m, n, v, 1, tau, c, ldc, work )
         GO TO 410
   10    CONTINUE
*
*        Special code for 1 x 1 Householder
*
         t1 = one - tau*v( 1 )*v( 1 )
         DO 20 j = 1, n
            c( 1, j ) = t1*c( 1, j )
   20    CONTINUE
         GO TO 410
   30    CONTINUE
*
*        Special code for 2 x 2 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         DO 40 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
   40    CONTINUE
         GO TO 410
   50    CONTINUE
*
*        Special code for 3 x 3 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         DO 60 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
   60    CONTINUE
         GO TO 410
   70    CONTINUE
*
*        Special code for 4 x 4 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         DO 80 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) +
     $            v4*c( 4, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
   80    CONTINUE
         GO TO 410
   90    CONTINUE
*
*        Special code for 5 x 5 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         DO 100 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) +
     $            v4*c( 4, j ) + v5*c( 5, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
  100    CONTINUE
         GO TO 410
  110    CONTINUE
*
*        Special code for 6 x 6 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         DO 120 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) +
     $            v4*c( 4, j ) + v5*c( 5, j ) + v6*c( 6, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
            c( 6, j ) = c( 6, j ) - sum*t6
  120    CONTINUE
         GO TO 410
  130    CONTINUE
*
*        Special code for 7 x 7 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         DO 140 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) +
     $            v4*c( 4, j ) + v5*c( 5, j ) + v6*c( 6, j ) +
     $            v7*c( 7, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
            c( 6, j ) = c( 6, j ) - sum*t6
            c( 7, j ) = c( 7, j ) - sum*t7
  140    CONTINUE
         GO TO 410
  150    CONTINUE
*
*        Special code for 8 x 8 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         DO 160 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) +
     $            v4*c( 4, j ) + v5*c( 5, j ) + v6*c( 6, j ) +
     $            v7*c( 7, j ) + v8*c( 8, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
            c( 6, j ) = c( 6, j ) - sum*t6
            c( 7, j ) = c( 7, j ) - sum*t7
            c( 8, j ) = c( 8, j ) - sum*t8
  160    CONTINUE
         GO TO 410
  170    CONTINUE
*
*        Special code for 9 x 9 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         v9 = v( 9 )
         t9 = tau*v9
         DO 180 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) +
     $            v4*c( 4, j ) + v5*c( 5, j ) + v6*c( 6, j ) +
     $            v7*c( 7, j ) + v8*c( 8, j ) + v9*c( 9, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
            c( 6, j ) = c( 6, j ) - sum*t6
            c( 7, j ) = c( 7, j ) - sum*t7
            c( 8, j ) = c( 8, j ) - sum*t8
            c( 9, j ) = c( 9, j ) - sum*t9
  180    CONTINUE
         GO TO 410
  190    CONTINUE
*
*        Special code for 10 x 10 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         v9 = v( 9 )
         t9 = tau*v9
         v10 = v( 10 )
         t10 = tau*v10
         DO 200 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) +
     $            v4*c( 4, j ) + v5*c( 5, j ) + v6*c( 6, j ) +
     $            v7*c( 7, j ) + v8*c( 8, j ) + v9*c( 9, j ) +
     $            v10*c( 10, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
            c( 6, j ) = c( 6, j ) - sum*t6
            c( 7, j ) = c( 7, j ) - sum*t7
            c( 8, j ) = c( 8, j ) - sum*t8
            c( 9, j ) = c( 9, j ) - sum*t9
            c( 10, j ) = c( 10, j ) - sum*t10
  200    CONTINUE
         GO TO 410
      ELSE
*
*        Form  C * H, where H has order n.
*
         GO TO ( 210, 230, 250, 270, 290, 310, 330, 350,
     $           370, 390 )n
*
*        Code for general N
*
         CALL dlarf( side, m, n, v, 1, tau, c, ldc, work )
         GO TO 410
  210    CONTINUE
*
*        Special code for 1 x 1 Householder
*
         t1 = one - tau*v( 1 )*v( 1 )
         DO 220 j = 1, m
            c( j, 1 ) = t1*c( j, 1 )
  220    CONTINUE
         GO TO 410
  230    CONTINUE
*
*        Special code for 2 x 2 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         DO 240 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
  240    CONTINUE
         GO TO 410
  250    CONTINUE
*
*        Special code for 3 x 3 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         DO 260 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
  260    CONTINUE
         GO TO 410
  270    CONTINUE
*
*        Special code for 4 x 4 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         DO 280 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) +
     $            v4*c( j, 4 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
  280    CONTINUE
         GO TO 410
  290    CONTINUE
*
*        Special code for 5 x 5 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         DO 300 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) +
     $            v4*c( j, 4 ) + v5*c( j, 5 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
  300    CONTINUE
         GO TO 410
  310    CONTINUE
*
*        Special code for 6 x 6 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         DO 320 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) +
     $            v4*c( j, 4 ) + v5*c( j, 5 ) + v6*c( j, 6 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
            c( j, 6 ) = c( j, 6 ) - sum*t6
  320    CONTINUE
         GO TO 410
  330    CONTINUE
*
*        Special code for 7 x 7 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         DO 340 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) +
     $            v4*c( j, 4 ) + v5*c( j, 5 ) + v6*c( j, 6 ) +
     $            v7*c( j, 7 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
            c( j, 6 ) = c( j, 6 ) - sum*t6
            c( j, 7 ) = c( j, 7 ) - sum*t7
  340    CONTINUE
         GO TO 410
  350    CONTINUE
*
*        Special code for 8 x 8 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         DO 360 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) +
     $            v4*c( j, 4 ) + v5*c( j, 5 ) + v6*c( j, 6 ) +
     $            v7*c( j, 7 ) + v8*c( j, 8 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
            c( j, 6 ) = c( j, 6 ) - sum*t6
            c( j, 7 ) = c( j, 7 ) - sum*t7
            c( j, 8 ) = c( j, 8 ) - sum*t8
  360    CONTINUE
         GO TO 410
  370    CONTINUE
*
*        Special code for 9 x 9 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         v9 = v( 9 )
         t9 = tau*v9
         DO 380 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) +
     $            v4*c( j, 4 ) + v5*c( j, 5 ) + v6*c( j, 6 ) +
     $            v7*c( j, 7 ) + v8*c( j, 8 ) + v9*c( j, 9 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
            c( j, 6 ) = c( j, 6 ) - sum*t6
            c( j, 7 ) = c( j, 7 ) - sum*t7
            c( j, 8 ) = c( j, 8 ) - sum*t8
            c( j, 9 ) = c( j, 9 ) - sum*t9
  380    CONTINUE
         GO TO 410
  390    CONTINUE
*
*        Special code for 10 x 10 Householder
*
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         v9 = v( 9 )
         t9 = tau*v9
         v10 = v( 10 )
         t10 = tau*v10
         DO 400 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) +
     $            v4*c( j, 4 ) + v5*c( j, 5 ) + v6*c( j, 6 ) +
     $            v7*c( j, 7 ) + v8*c( j, 8 ) + v9*c( j, 9 ) +
     $            v10*c( j, 10 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
            c( j, 6 ) = c( j, 6 ) - sum*t6
            c( j, 7 ) = c( j, 7 ) - sum*t7
            c( j, 8 ) = c( j, 8 ) - sum*t8
            c( j, 9 ) = c( j, 9 ) - sum*t9
            c( j, 10 ) = c( j, 10 ) - sum*t10
  400    CONTINUE
         GO TO 410
      END IF
  410 CONTINUE
      RETURN
*
*     End of DLARFX
*
      END
      SUBROUTINE dlasy2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,
     $                   LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            LTRANL, LTRANR
      INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
      DOUBLE PRECISION   SCALE, XNORM
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ),
     $                   x( ldx, * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      DOUBLE PRECISION   TWO, HALF, EIGHT
      parameter( two = 2.0d+0, half = 0.5d+0, eight = 8.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BSWAP, XSWAP
      INTEGER            I, IP, IPIV, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1,
     $                   temp, u11, u12, u22, xmax
*     ..
*     .. Local Arrays ..
      LOGICAL            BSWPIV( 4 ), XSWPIV( 4 )
      INTEGER            JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ),
     $                   locu22( 4 )
      DOUBLE PRECISION   BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 )
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           idamax, dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           dcopy, dswap
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max
*     ..
*     .. Data statements ..
      DATA               locu12 / 3, 4, 1, 2 / , locl21 / 2, 1, 4, 3 / ,
     $                   locu22 / 4, 3, 2, 1 /
      DATA               xswpiv / .false., .false., .true., .true. /
      DATA               bswpiv / .false., .true., .false., .true. /
*     ..
*     .. Executable Statements ..
*
*     Do not check the input parameters for errors
*
      info = 0
*
*     Quick return if possible
*
      IF( n1.EQ.0 .OR. n2.EQ.0 )
     $   RETURN
*
*     Set constants to control overflow
*
      eps = dlamch( 'P' )
      smlnum = dlamch( 'S' ) / eps
      sgn = isgn
*
      k = n1 + n1 + n2 - 2
      GO TO ( 10, 20, 30, 50 )k
*
*     1 by 1: TL11*X + SGN*X*TR11 = B11
*
   10 CONTINUE
      tau1 = tl( 1, 1 ) + sgn*tr( 1, 1 )
      bet = abs( tau1 )
      IF( bet.LE.smlnum ) THEN
         tau1 = smlnum
         bet = smlnum
         info = 1
      END IF
*
      scale = one
      gam = abs( b( 1, 1 ) )
      IF( smlnum*gam.GT.bet )
     $   scale = one / gam
*
      x( 1, 1 ) = ( b( 1, 1 )*scale ) / tau1
      xnorm = abs( x( 1, 1 ) )
      RETURN
*
*     1 by 2:
*     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
*                                       [TR21 TR22]
*
   20 CONTINUE
*
      smin = max( eps*max( abs( tl( 1, 1 ) ), abs( tr( 1, 1 ) ),
     $       abs( tr( 1, 2 ) ), abs( tr( 2, 1 ) ), abs( tr( 2, 2 ) ) ),
     $       smlnum )
      tmp( 1 ) = tl( 1, 1 ) + sgn*tr( 1, 1 )
      tmp( 4 ) = tl( 1, 1 ) + sgn*tr( 2, 2 )
      IF( ltranr ) THEN
         tmp( 2 ) = sgn*tr( 2, 1 )
         tmp( 3 ) = sgn*tr( 1, 2 )
      ELSE
         tmp( 2 ) = sgn*tr( 1, 2 )
         tmp( 3 ) = sgn*tr( 2, 1 )
      END IF
      btmp( 1 ) = b( 1, 1 )
      btmp( 2 ) = b( 1, 2 )
      GO TO 40
*
*     2 by 1:
*          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
*            [TL21 TL22] [X21]         [X21]         [B21]
*
   30 CONTINUE
      smin = max( eps*max( abs( tr( 1, 1 ) ), abs( tl( 1, 1 ) ),
     $       abs( tl( 1, 2 ) ), abs( tl( 2, 1 ) ), abs( tl( 2, 2 ) ) ),
     $       smlnum )
      tmp( 1 ) = tl( 1, 1 ) + sgn*tr( 1, 1 )
      tmp( 4 ) = tl( 2, 2 ) + sgn*tr( 1, 1 )
      IF( ltranl ) THEN
         tmp( 2 ) = tl( 1, 2 )
         tmp( 3 ) = tl( 2, 1 )
      ELSE
         tmp( 2 ) = tl( 2, 1 )
         tmp( 3 ) = tl( 1, 2 )
      END IF
      btmp( 1 ) = b( 1, 1 )
      btmp( 2 ) = b( 2, 1 )
   40 CONTINUE
*
*     Solve 2 by 2 system using complete pivoting.
*     Set pivots less than SMIN to SMIN.
*
      ipiv = idamax( 4, tmp, 1 )
      u11 = tmp( ipiv )
      IF( abs( u11 ).LE.smin ) THEN
         info = 1
         u11 = smin
      END IF
      u12 = tmp( locu12( ipiv ) )
      l21 = tmp( locl21( ipiv ) ) / u11
      u22 = tmp( locu22( ipiv ) ) - u12*l21
      xswap = xswpiv( ipiv )
      bswap = bswpiv( ipiv )
      IF( abs( u22 ).LE.smin ) THEN
         info = 1
         u22 = smin
      END IF
      IF( bswap ) THEN
         temp = btmp( 2 )
         btmp( 2 ) = btmp( 1 ) - l21*temp
         btmp( 1 ) = temp
      ELSE
         btmp( 2 ) = btmp( 2 ) - l21*btmp( 1 )
      END IF
      scale = one
      IF( ( two*smlnum )*abs( btmp( 2 ) ).GT.abs( u22 ) .OR.
     $    ( two*smlnum )*abs( btmp( 1 ) ).GT.abs( u11 ) ) THEN
         scale = half / max( abs( btmp( 1 ) ), abs( btmp( 2 ) ) )
         btmp( 1 ) = btmp( 1 )*scale
         btmp( 2 ) = btmp( 2 )*scale
      END IF
      x2( 2 ) = btmp( 2 ) / u22
      x2( 1 ) = btmp( 1 ) / u11 - ( u12 / u11 )*x2( 2 )
      IF( xswap ) THEN
         temp = x2( 2 )
         x2( 2 ) = x2( 1 )
         x2( 1 ) = temp
      END IF
      x( 1, 1 ) = x2( 1 )
      IF( n1.EQ.1 ) THEN
         x( 1, 2 ) = x2( 2 )
         xnorm = abs( x( 1, 1 ) ) + abs( x( 1, 2 ) )
      ELSE
         x( 2, 1 ) = x2( 2 )
         xnorm = max( abs( x( 1, 1 ) ), abs( x( 2, 1 ) ) )
      END IF
      RETURN
*
*     2 by 2:
*     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
*       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
*
*     Solve equivalent 4 by 4 system using complete pivoting.
*     Set pivots less than SMIN to SMIN.
*
   50 CONTINUE
      smin = max( abs( tr( 1, 1 ) ), abs( tr( 1, 2 ) ),
     $       abs( tr( 2, 1 ) ), abs( tr( 2, 2 ) ) )
      smin = max( smin, abs( tl( 1, 1 ) ), abs( tl( 1, 2 ) ),
     $       abs( tl( 2, 1 ) ), abs( tl( 2, 2 ) ) )
      smin = max( eps*smin, smlnum )
      btmp( 1 ) = zero
      CALL dcopy( 16, btmp, 0, t16, 1 )
      t16( 1, 1 ) = tl( 1, 1 ) + sgn*tr( 1, 1 )
      t16( 2, 2 ) = tl( 2, 2 ) + sgn*tr( 1, 1 )
      t16( 3, 3 ) = tl( 1, 1 ) + sgn*tr( 2, 2 )
      t16( 4, 4 ) = tl( 2, 2 ) + sgn*tr( 2, 2 )
      IF( ltranl ) THEN
         t16( 1, 2 ) = tl( 2, 1 )
         t16( 2, 1 ) = tl( 1, 2 )
         t16( 3, 4 ) = tl( 2, 1 )
         t16( 4, 3 ) = tl( 1, 2 )
      ELSE
         t16( 1, 2 ) = tl( 1, 2 )
         t16( 2, 1 ) = tl( 2, 1 )
         t16( 3, 4 ) = tl( 1, 2 )
         t16( 4, 3 ) = tl( 2, 1 )
      END IF
      IF( ltranr ) THEN
         t16( 1, 3 ) = sgn*tr( 1, 2 )
         t16( 2, 4 ) = sgn*tr( 1, 2 )
         t16( 3, 1 ) = sgn*tr( 2, 1 )
         t16( 4, 2 ) = sgn*tr( 2, 1 )
      ELSE
         t16( 1, 3 ) = sgn*tr( 2, 1 )
         t16( 2, 4 ) = sgn*tr( 2, 1 )
         t16( 3, 1 ) = sgn*tr( 1, 2 )
         t16( 4, 2 ) = sgn*tr( 1, 2 )
      END IF
      btmp( 1 ) = b( 1, 1 )
      btmp( 2 ) = b( 2, 1 )
      btmp( 3 ) = b( 1, 2 )
      btmp( 4 ) = b( 2, 2 )
*
*     Perform elimination
*
      DO 100 i = 1, 3
         xmax = zero
         DO 70 ip = i, 4
            DO 60 jp = i, 4
               IF( abs( t16( ip, jp ) ).GE.xmax ) THEN
                  xmax = abs( t16( ip, jp ) )
                  ipsv = ip
                  jpsv = jp
               END IF
   60       CONTINUE
   70    CONTINUE
         IF( ipsv.NE.i ) THEN
            CALL dswap( 4, t16( ipsv, 1 ), 4, t16( i, 1 ), 4 )
            temp = btmp( i )
            btmp( i ) = btmp( ipsv )
            btmp( ipsv ) = temp
         END IF
         IF( jpsv.NE.i )
     $      CALL dswap( 4, t16( 1, jpsv ), 1, t16( 1, i ), 1 )
         jpiv( i ) = jpsv
         IF( abs( t16( i, i ) ).LT.smin ) THEN
            info = 1
            t16( i, i ) = smin
         END IF
         DO 90 j = i + 1, 4
            t16( j, i ) = t16( j, i ) / t16( i, i )
            btmp( j ) = btmp( j ) - t16( j, i )*btmp( i )
            DO 80 k = i + 1, 4
               t16( j, k ) = t16( j, k ) - t16( j, i )*t16( i, k )
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      IF( abs( t16( 4, 4 ) ).LT.smin ) THEN
         info = 1
         t16( 4, 4 ) = smin
      END IF
      scale = one
      IF( ( eight*smlnum )*abs( btmp( 1 ) ).GT.abs( t16( 1, 1 ) ) .OR.
     $    ( eight*smlnum )*abs( btmp( 2 ) ).GT.abs( t16( 2, 2 ) ) .OR.
     $    ( eight*smlnum )*abs( btmp( 3 ) ).GT.abs( t16( 3, 3 ) ) .OR.
     $    ( eight*smlnum )*abs( btmp( 4 ) ).GT.abs( t16( 4, 4 ) ) ) THEN
         scale = ( one / eight ) / max( abs( btmp( 1 ) ),
     $           abs( btmp( 2 ) ), abs( btmp( 3 ) ), abs( btmp( 4 ) ) )
         btmp( 1 ) = btmp( 1 )*scale
         btmp( 2 ) = btmp( 2 )*scale
         btmp( 3 ) = btmp( 3 )*scale
         btmp( 4 ) = btmp( 4 )*scale
      END IF
      DO 120 i = 1, 4
         k = 5 - i
         temp = one / t16( k, k )
         tmp( k ) = btmp( k )*temp
         DO 110 j = k + 1, 4
            tmp( k ) = tmp( k ) - ( temp*t16( k, j ) )*tmp( j )
  110    CONTINUE
  120 CONTINUE
      DO 130 i = 1, 3
         IF( jpiv( 4-i ).NE.4-i ) THEN
            temp = tmp( 4-i )
            tmp( 4-i ) = tmp( jpiv( 4-i ) )
            tmp( jpiv( 4-i ) ) = temp
         END IF
  130 CONTINUE
      x( 1, 1 ) = tmp( 1 )
      x( 2, 1 ) = tmp( 2 )
      x( 1, 2 ) = tmp( 3 )
      x( 2, 2 ) = tmp( 4 )
      xnorm = max( abs( tmp( 1 ) )+abs( tmp( 3 ) ),
     $        abs( tmp( 2 ) )+abs( tmp( 4 ) ) )
      RETURN
*
*     End of DLASY2
*
      END
