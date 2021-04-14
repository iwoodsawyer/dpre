	SUBROUTINE ZPGEX2( WANTQ, K, N, J, S, A, LDA1, LDA2, Q, LDQ1,
     $                   LDQ2, ZWORK, INFO )
	IMPLICIT NONE
C
C     PURPOSE
C
C     ZPGEX2 swaps adjacent diagonal 1-by-1 blocks in a complex
C     generalized matrix product,
C
C                                 S(2)                 S(K)
C          A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K),
C
C     by unitary equivalence transformations. A must be in periodic
C     Schur form, that is, all factors of A must be upper triangular.
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
C             The order of each factor in A.  N >= 2.
C
C     J       (input) INTEGER
C             The index of the first block to be swapped.  1 <= J < N.
C
C     S       (input)  INTEGER array, dimension (K)
C             The leading K elements of this array must contain the
C             signatures of the factors. Each entry in S must be either
C             1 or -1. By definition, S(1) must be set to 1.
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K)
C             On entry, the leading N-by-N-by-K part of this array
C             must contain the factors in periodic Schur form
C             form, that is, all factors are upper triangular.
C             On exit, if INFO = 0, the leading N-by-N-by-K part of
C             this array contains the factors of the reordered periodic
C             Schur form.
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
C     Workspace
C
C     ZWORK   COMPLEX*16 array, dimension (K)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0       : succesful exit;
C             = 1       : the periodic QZ algorithm failed to converge.
C
C     METHOD
C
C     A complex version of the periodic QZ algorithm [1] with one
C     perfect shift is used. For more details see [2]. This is not
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
      DOUBLE PRECISION  ONE, ZERO
	PARAMETER         ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16        CONE, CZERO
      PARAMETER         ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                  CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Scalar Arguments ..
      LOGICAL           WANTQ
      INTEGER           J, K, INFO, LDA1, LDA2, LDQ1, LDQ2, N
C     .. Array Arguments ..
      INTEGER           S(*)
      COMPLEX*16        A(LDA1, LDA2, *), Q(LDQ1, LDQ2, *),
     $                  ZWORK(*)
C     .. Local Scalars ..
      LOGICAL           USEZQ
      INTEGER           JITER, L, LN
      DOUBLE PRECISION  CS, CSF, RHS, SAFMAX, SAFMIN, SMLNUM, ULP
	COMPLEX*16        SN, SNF, TEMP
C     .. Local Arrays ..
      INTEGER           ISEED(4)
	COMPLEX*16        RND(3)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
	EXTERNAL          DLAMCH
C     .. External Subroutines ..
      EXTERNAL          DLABAD, ZLARTG, ZROT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCONJG, MAX
C
C     .. Executable Statements ..
C
C     Save perfect shift.
C
      DO 10  L = 1, K
      ZWORK(L) = A(J,J,L)
   10 CONTINUE
      INFO = 0
      ISEED(1) = 1
      ISEED(2) = 0
      ISEED(3) = 0
      ISEED(4) = 0
      SAFMIN = DLAMCH( 'SafeMinimum' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'Precision' )
      CALL DLABAD( SAFMIN, SAFMAX )
      SMLNUM = SAFMIN*( 2 / ULP )
C
C     If there are infinite or zero eigenvalues in critical positions,
C     then unshifted ZQ iterations must be used.
C
      USEZQ = .FALSE.
      DO 20  L = 1, K
         USEZQ = ( USEZQ.OR.
     $             ( ( S(L).EQ.1 ).AND.( A(J+1,J+1,L).EQ.CZERO ) ).OR.
     $             ( ( S(L).EQ.-1 ).AND.( A(J,J,L).EQ.CZERO ) ) )
   20 CONTINUE
C
C     Destroy any triangular structure.
C
      CALL ZLARNV( 2, ISEED, 2, RND )
      CALL ZLARTG( RND(1), RND(2), CS, SN, TEMP )
      CSF = CS
      SNF = SN
      DO 30  L = 1, K
         IF ( WANTQ ) THEN
            CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, CS, DCONJG( SN ) )
         END IF
         IF ( S(L).EQ.1 ) THEN
            CALL ZROT( N-J+1, A(J,J,L), LDA1, A(J+1,J,L), LDA1, CS, SN )
         ELSE
            CALL ZROT( J+1, A(1,J,L), 1, A(1,J+1,L), 1, CS,
     $                 DCONJG( SN ) )	      
         END IF
         IF ( L.EQ.K ) THEN
            CS = CSF
            SN = SNF
         ELSE
            CALL ZLARNV( 2, ISEED, 2, RND )
            CALL ZLARTG( RND(1), RND(2), CS, SN, TEMP )
         END IF
         IF ( S(L).EQ.1 ) THEN
            CALL ZROT( J+1, A(1,J,L), 1, A(1,J+1,L), 1, CS,
     $                 DCONJG( SN ) )
         ELSE
            CALL ZROT( N-J+1, A(J,J,L), LDA1, A(J+1,J,L), LDA1, CS, SN )
         END IF
   30 CONTINUE
C
      IF ( USEZQ ) THEN
         DO 50  JITER = 1, 10
            DO 40  L = 1, K
               IF ( S(L).EQ.1 ) THEN
                  TEMP = A(J+1,J+1,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J+1,J+1,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( J, A(1,J+1,L), 1, A(1,J,L), 1, CS, SN )
                  SN = -SN
               ELSE
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J,J,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( N-J, A(J,J+1,L), LDA1, A(J+1,J+1,L), LDA1,
     $                       CS, SN )
               END IF
               LN = L + 1
               IF ( LN.GT.K )  LN = 1
               IF ( WANTQ ) THEN
                  CALL ZROT( N, Q(1,J,LN), 1, Q(1,J+1,LN), 1, CS,
     $                       DCONJG( SN ) )
               END IF
               IF ( S(LN).EQ.1 ) THEN
                  CALL ZROT( N-J+1, A(J,J,LN), LDA1, A(J+1,J,LN), LDA1,
     $                       CS, SN )
               ELSE
                  CALL ZROT( J+1, A(1,J,LN), 1, A(1,J+1,LN), 1, CS,
     $                       DCONJG( SN ) )
               END IF
   40       CONTINUE
            RHS = MAX( ABS( A(J,J,1) ), ABS( A(J+1,J+1,1) ) )
            IF ( RHS.EQ.ZERO )
     $         RHS = ABS( A(J,J+1,1) )
               RHS = MAX( ULP*RHS, SMLNUM )
            IF ( ABS( A(J+1,J,1) ).LE.RHS ) THEN
               A(J+1,J,1) = CZERO
               GO TO 90
            END IF
   50    CONTINUE
         INFO = 1
         GO TO 90
      END IF

      DO 80  JITER = 1, 42
C
C        Complex single shift.
C
         IF ( ( JITER / 10 )*10.EQ.JITER ) THEN
C
C           Random/exceptional shift.
C
            CALL ZLARNV( 2, ISEED, 2, RND )
            CALL ZLARTG( RND(1), RND(2), CS, SN, TEMP )
         ELSE
            CALL ZLARTG( CONE, CONE, CS, SN, TEMP )
            DO 60  L = K, 2, -1
               IF ( S(L).EQ.1 ) THEN
                  CALL ZLARTG( A(J,J,L)*CS, ZWORK(L)*DCONJG(SN),
     $                         CS, SN, TEMP )
               ELSE
                  CALL ZLARTG( ZWORK(L)*CS, -A(J,J,L)*DCONJG(SN),
     $                         CS, SN, TEMP )
                  SN = -SN
               END IF
   60       CONTINUE
            CALL ZLARTG( A(J,J,1)*CS - DCONJG(SN)*ZWORK(1),
     $                   A(J+1,J,1)*CS, CS, SN, TEMP )
         END IF
C
C        Do one QZ sweep.
C
         CALL ZROT( N-J+1, A(J,J,1), LDA1, A(J+1,J,1), LDA1,
     $              CS, SN )
         IF ( WANTQ ) THEN
            CALL ZROT( N, Q(1,J,1), 1, Q(1,J+1,1), 1, CS,
     $                 DCONJG( SN ) )
         END IF
C
C        Propagate rotation through AK, ..., A2 to A1.
C
         DO 70  L = K, 2, -1
            IF ( S(L).EQ.1 ) THEN
               CALL ZROT( J+1, A(1,J,L), 1, A(1,J+1,L), 1, CS,
     $                    DCONJG( SN ) )
               TEMP = A(J,J,L)
               CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J,J,L) )
               A(J+1,J,L) = CZERO
               CALL ZROT( N-J, A(J,J+1,L), LDA1, A(J+1,J+1,L), LDA1,
     $                    CS, SN )
            ELSE
               CALL ZROT( N-J+1, A(J,J,L), LDA1, A(J+1,J,L), LDA1, CS,
     $                    SN )
               TEMP = A(J+1,J+1,L)
               CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J+1,J+1,L) )
               A(J+1,J,L) = CZERO
               CALL ZROT( J, A(1,J+1,L), 1, A(1,J,L), 1, CS, SN )
               SN = -SN
            END IF
            IF ( WANTQ ) THEN
               CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, CS,
     $                    DCONJG( SN ) )
            END IF
   70    CONTINUE
         CALL ZROT( J+1, A(1,J,1), 1, A(1,J+1,1), 1, CS, DCONJG( SN ) )
C
C        Test for deflation.
C
         RHS = MAX( ABS( A(J,J,1) ), ABS( A(J+1,J+1,1) ) )
         IF ( RHS.EQ.ZERO )
     $      RHS = ABS( A(J,J+1,1) )
            RHS = MAX( ULP*RHS, SMLNUM )
         IF ( ABS( A(J+1,J,1) ).LE.RHS ) THEN
            A(J+1,J,1) = CZERO
            GO TO 90
         END IF
   80 CONTINUE
C     
C     Not converged.
C
      INFO = 1
C
   90 CONTINUE
C
C     Check for singular triangular factors.
C
      DO 100  L = 1, K
         RHS = MAX( ABS( A(J,J+1,L) ), ABS( A(J+1,J+1,L) ) )
         RHS = MAX( ULP*RHS, SMLNUM )
         IF ( ABS( A(J,J,L) ).LE.RHS ) THEN
            A(J,J,L) = CZERO
         END IF
         RHS = MAX( ABS( A(J,J+1,L) ), ABS( A(J,J,L) ) )
         RHS = MAX( ULP*RHS, SMLNUM )
         IF ( ABS( A(J+1,J+1,L) ).LE.RHS ) THEN
            A(J+1,J+1,L) = CZERO
         END IF
  100 CONTINUE
C
      RETURN
C
C *** Last line of ZPGEX2***
      END
