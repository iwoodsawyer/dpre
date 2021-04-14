      SUBROUTINE ZPGORD( WANTQ, K, N, S, SELECT, A, LDA1, LDA2,
     $                   ALPHA, BETA, SCAL, Q, LDQ1, LDQ2, M,
     $                   ZWORK, LZWORK, INFO )
	IMPLICIT NONE
C
C     PURPOSE
C
C     ZPGORD reorders the periodic Schur decomposition of a complex
C     generalized matrix product
C
C                                 S(2)                 S(K)
C          A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K),
C
C     (in terms of unitary equivalence transformations), so that a
C     selected cluster of eigenvalues appears in the leading diagonal
C     blocks of the matrix product. The leading columns of the
C     orthogonal factors contained in Q form unitary bases of the
C     corresponding periodic eigenspaces (deflating subspaces).
C     A must be in periodic Schur form, that is, all factors of A must
C     be upper triangular.
C
C     If WANTQ = .TRUE., then the unitary factors are computed and
C     stored in the array Q so that for S(I) = 1,
C
C                         H
C             Q(:,:,I)(in)   A(:,:,I)(in)   Q(:,:,MOD(I,K)+1)(in)
C                          H                                        (1)
C         =   Q(:,:,I)(out)  A(:,:,I)(out)  Q(:,:,MOD(I,K)+1)(out),
C
C     and for S(I) = -1,
C
C                                  H
C             Q(:,:,MOD(I,K)+1)(in)   A(:,:,I)(in)   Q(:,:,I)(in)
C                                   H                               (2)
C         =   Q(:,:,MOD(I,K)+1)(out)  A(:,:,I)(out)  Q(:,:,I)(out).
C
C     ARGUMEMTS
C
C     Mode Parameters
C
C     WANTQ   (input) LOGICAL
C             = .FALSE.: do not modify Q;
C             = .TRUE. : modify the array Q by the unitary
C                        transformations that are applied to the
C                        matrices in A for reordering.
C
C     Input/Output Parameters
C
C     K       (input)  INTEGER
C             The number of factors.  K >= 1.
C
C     N       (input)  INTEGER
C             The order of each factor in A.  N >= 0.
C
C     S       (input)  INTEGER array, dimension (K)
C             The leading K elements of this array must contain the
C             signatures of the factors. Each entry in S must be either
C             1 or -1. By definition, S(1) must be set to 1.
C
C     SELECT  (input) LOGICAL array, dimension (N)
C             SELECT specifies the eigenvalues in the selected cluster.
C             To select the eigenvalue corresponding to the (j,j)
C             diagonal entries, SELECT(j) must be set to .TRUE..
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K)
C             On entry, the leading N-by-N-by-K part of this array
C             must contain the factors in periodic Schur form
C             form, that is, all factors are upper triangular.
C             On exit, if INFO = 0, the leading N-by-N-by-K part of
C             this array contains the factors of the reordered periodic
C             Schur form. Moreover, A(:,:,2), ..., A(:,:,K) are
C             normalized so that their diagonals contain nonnegative
C             real numbers.
C
C     LDA1    (input) INTEGER
C             The first leading dimension of A.  LDA1 >= MAX(1,N).
C
C     LDA2    (input) INTEGER
C             The second leading dimension of A.  LDA2 >= MAX(1,N).
C
C     Q       (input/output) COMPLEX*16 array, dimension (LDQ1,LDQ2,K)
C             On entry, if WANTQ = .TRUE., the leading N-by-N-by-K part
C             of this array must contain the initial unitary factors
C             as described in (1)-(2).
C             On exit, if WANTQ = .TRUE., the leading N-by-N-by-K part
C             of this array contains the modified orthogonal factors as
C             described in (1)-(2).
C
C     LDQ1    (input)  INTEGER
C             The first leading dimension of Q.
C             If WANTQ = .TRUE.,  LDQ1 >= MAX(1,N).
C
C     LDQ2    (input)  INTEGER
C             The second leading dimension of Q.
C             If WANTQ = .TRUE.,  LDQ2 >= MAX(1,N).
C
C     M       (output) INTEGER
C             The dimension of the specified periodic eigenspace.
C
C     Workspace
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C             On exit, if INFO = 0, ZWORK(1) returns the minimal value
C             of LZWORK.
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.  LZWORK >= MAX( K, N ).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0       : succesful exit;
C             < 0       : if INFO = -i, the i-th argument had an
C                         illegal value;
C             = 1       : the periodic QZ algorithm failed to converge.
C
C     METHOD
C
C     A complex version of the periodic QZ algorithm [1] with 
C     perfect shifts is used. For more details see [2]. This is not
C     a safe method. It is advisable to check whether the eigenvalues
C     are really reorderd.
C
C     REFERENCES
C
C     [1] Bojanczyk, A. and Golub, G. H. and Van Dooren, P.
C         The periodic Schur decomposition; algorithm and applications.
C         In Proc. SPIE Conference, pg. 31-42, vol. 1770, 1992.
C
C     [2] Kressner, D.
C         An efficient and reliable implementation of the periodic QZ
C         algorithm. In IFAC Workshop on Periodic Control Systems, 2001.
C
C     NUMERICAL ASPECTS
C
C     The implemented method is numerically backward stable.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, Dec. 2002.
C
C     ******************************************************************
C
C     .. Parameters ..
      COMPLEX*16        CONE, CZERO
      PARAMETER         ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                  CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Scalar Arguments ..
      LOGICAL           WANTQ
      INTEGER           K, INFO, LDA1, LDA2, LDQ1, LDQ2, LZWORK, M, N
C     .. Array Arguments ..
      LOGICAL           SELECT(*)
      INTEGER           S(*), SCAL(*)
      COMPLEX*16        A(LDA1, LDA2, *), ALPHA(*), BETA(*),
     $                  Q(LDQ1, LDQ2, *), ZWORK(*)
C     .. Local Scalars ..
      LOGICAL           SOK
      INTEGER           I, IERR, J, JS, L
	DOUBLE PRECISION  ABST, BASE, SAFMIN
	COMPLEX*16        TEMP
C     .. Local Arrays ..
	DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C     .. External Subroutines ..
      EXTERNAL          ZLAPR1, ZPGEX2, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         DCMPLX, DCONJG, MAX
C
C     .. Executable Statements ..
C
      INFO = 0
	SAFMIN = DLAMCH( 'SafeMinimum' )
      BASE = DLAMCH( 'Base' )
C
C     Check the scalar input parameters.
C
      IF ( K.LT.1 ) THEN
         INFO = -2
      ELSE IF ( N.LT.0 ) THEN
         INFO = -3
      ELSE
	   SOK = S(1).EQ.1
	   DO 10  L = 2, K
	      SOK = S(L).EQ.1 .OR. S(L).EQ.-1
   10    CONTINUE
         IF ( .NOT.SOK ) THEN
	      INFO = -4
         ELSE IF ( LDA1 .LT. MAX( 1, N ) ) THEN
            INFO = -7
         ELSE IF ( LDA2 .LT. MAX( 1, N ) ) THEN	
            INFO = -8
         ELSE IF ( WANTQ .AND. LDQ1.LT.MAX( 1, N ) ) THEN
            INFO = -13
         ELSE IF ( WANTQ .AND. LDQ2.LT.MAX( 1, N ) ) THEN
            INFO = -14
	   ELSE IF ( LZWORK.LT.MAX( K, N ) ) THEN
            INFO = -17
         END IF
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPGORD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         ZWORK(1) = CONE
         RETURN
      END IF
C
C     Set M to the dimension of the specified deflating subspace.
C
      M = 0
      DO 20 J = 1, N
         IF ( SELECT( J ) )
     $      M = M + 1
   20 CONTINUE
C
      JS = 0
      DO 40 J = 1, N
         IF ( SELECT( J ) ) THEN
            JS = JS + 1
C
C           Swap the J-th block to position JS.
C
            IF ( J.NE.JS ) THEN
	         DO 30  I = J-1, JS, -1
                  CALL ZPGEX2( .TRUE., K, N, I, S, A, LDA1, LDA2, Q,
     $                         LDQ1, LDQ2, ZWORK, IERR )
	            IF ( IERR.NE.0 ) THEN
	               INFO = 1
	               GO TO 120
	            END IF
   30          CONTINUE
	      END IF
         END IF
   40 CONTINUE
C
C     Rescale matrices and recompute eigenvalues.
C
      DO 100  L = K, 2, -1
	   IF ( S(L).EQ.1 )  THEN
	      DO 50 J = 1, N
               ABST = ABS( A(J,J,L) )
               IF ( ABST.GT.SAFMIN ) THEN
                  TEMP = DCONJG( A(J,J,L) / ABST )
                  A(J,J,L ) = ABST
                  CALL ZSCAL( N-J, TEMP, A(J,J+1,L), LDA1 )
	         ELSE
	            TEMP = CONE
	         END IF
               ZWORK(J) = TEMP
   50       CONTINUE	
	   ELSE
            DO 60  J = 1, N
               ABST = ABS( A(J,J,L) )
               IF ( ABST.GT.SAFMIN ) THEN
                  TEMP = DCONJG( A(J,J,L) / ABST )
                  A(J,J,L ) = ABST
                  CALL ZSCAL( J-1, TEMP, A(1,J,L), 1 )
	         ELSE
	            TEMP = CONE
	         END IF
               ZWORK(J) = DCONJG(TEMP)
   60       CONTINUE	
	   END IF
         IF ( WANTQ ) THEN
	      DO 70  J = 1, N
               CALL ZSCAL( N, DCONJG( ZWORK(J) ), Q(1,J,L), 1 )
   70       CONTINUE
         END IF
         IF ( S(L-1).EQ.1 )  THEN
	      DO 80  J = 1, N
	         CALL ZSCAL( J, DCONJG( ZWORK(J) ), A(1,J,L-1), 1 )
   80       CONTINUE
         ELSE
	      DO 90  J = 1, N
	         CALL ZSCAL( N-J+1, ZWORK(J), A(J,J,L-1), LDA1 )
   90       CONTINUE              
	   END IF
  100 CONTINUE
      DO 110  J = 1, N
         CALL ZLAPR1( BASE, K, S, A(J,J,1), LDA1*LDA2, ALPHA(J),
     $                BETA(J), SCAL(J) )
  110 CONTINUE
C
  120 CONTINUE
      ZWORK(1) = DCMPLX( MAX( K, N ), 0 )
      RETURN
C *** Last line of ZPGORD ***
      END
