      SUBROUTINE ZPGHRD( COMPQ, K, N, ILO, IHI, S, A, LDA1, LDA2, Q,
     $                   LDQ1, LDQ2, DWORK, LDWORK, ZWORK, LZWORK,
     $                   INFO )
      IMPLICIT NONE
C
C     PURPOSE
C
C     To reduce the general complex product
C
C                               S(2)                 S(K)
C          A(:,:,1)  *  A(:,:,2)     * ... * A(:,:,K)
C
C     to upper Hessenberg-triangular form, where A is N-by-N-by-K and S
C     is the signature array with values 1 or -1. The matrix A(:,:,1)
C     is reduced to upper Hessenberg form while the other matrices are
C     triangularized.
C     Simple and unblocked version.
C
C     If COMPQ = 'V' or COMPQ = 'I', then the unitary factors are
C     computed and stored in the array Q so that for S(I) = 1,
C
C                         H
C             Q(:,:,I)(in)   A(:,:,I)(in)   Q(:,:,MOD(I,K)+1)(in)
C                          H
C         =   Q(:,:,I)(out)  A(:,:,I)(out)  Q(:,:,MOD(I,K)+1)(out),
C
C     and for S(I) = -1,
C
C                                  H
C             Q(:,:,MOD(I,K)+1)(in)   A(:,:,I)(in)   Q(:,:,I)(in)
C                                   H
C         =   Q(:,:,MOD(I,K)+1)(out)  A(:,:,I)(out)  Q(:,:,I)(out).
C
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COMPQ   (input) CHARACTER*1
C             = 'N': do not modify Q.
C             = 'V': modify the array Q by the unitary transformations
C                    that are applied to the matrices in A to reduce them
C                    to Hessenberg-triangular form.
C             = 'I': like COMPQ='V', except that each matrix in Q will
C                    be initialized to the identity first.
C
C     Input/Output Parameters
C
C     K       (input) INTEGER
C             The number of matrices in A.  K >= 1.
C
C     N       (input) INTEGER
C             Order of each factor in A.  N >= 0.
C
C     ILO     (input) INTEGER
C     IHI     (input) INTEGER
C             It is assumed that each factor in A is already upper
C             triangular in rows and columns 1:ILO-1 and IHI+1:N.
C             1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
C
C     S       (input) INTEGER array, dimension (K)
C             The leading K elements of this array must contain the
C             signatures of the factors. Each entry in S must be either
C             1 or -1. By definition, S(1) must be set to 1.
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K)
C             On entry, the leading N-by-N-by-K part of this array must
C             contain the factors of the general product to be reduced.
C             On exit, A(:,:,1) is overwritten by an upper Hessenberg
C             matrix and each A(:,:,I) for I not equal to 1 is
C             overwritten by an upper triangular matrix.
C
C     LDA1    (input) INTEGER
C             The first leading dimension of A. LDA1 >= max(1,N)
C
C     LDA2    (input) INTEGER
C             The second leading dimension of A. LDA2 >= max(1,N)
C
C     Q       (input/output) COMPLEX*16 array, dimension (LDQ1,LDQ2,K)
C             If COMPQ='N': Q is not referenced.
C             If COMPQ='I': On entry, Q need not to be set, and on exit
C                           it contains the unitary transformations.
C             If COMPQ='V': On entry, Q must contain unitary matrices,
C                           and on exit this is overwritten by the
C                           updated transformations.
C
C     LDQ1    (input) INTEGER
C             The first leading dimension of Q. LDQ1 >= max(1,N)
C
C     LDQ2    (input) INTEGER
C             The second leading dimension of Q. LDQ2 >= max(1,N)
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the minimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N). 
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C             On exit, if INFO = 0, ZWORK(1) returns the optimal value
C             of LZWORK.
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.  LZWORK >= MAX(1,2*N). 
C             For optimal performance this value should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C
C     METHOD
C
C     A slightly modified version of the periodic Hessenberg reduction
C     presented in [1] is used. For more details see [2].
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
C                                 3
C     The algorithm requires 0(K N ) floating point operations.
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
      CHARACTER*1       COMPQ
      INTEGER           K, ILO, IHI, INFO, LDA1, LDA2, LDQ1, LDQ2,
     $                  LDWORK, LZWORK, N
C     .. Array Arguments ..
      INTEGER           S(*)
      DOUBLE PRECISION  DWORK(*)
      COMPLEX*16        A(LDA1, LDA2, *), Q(LDQ1, LDQ2, *),
     $                  ZWORK(*)
C     .. Local Scalars ..
      LOGICAL           INITQ, SOK, WANTQ
      INTEGER           ICOLS, IERR, IROWS, JCOL, JROW, L, WRKOPT
      DOUBLE PRECISION  CS
      COMPLEX*16        SN, TEMP
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          ZGEQRF, ZGERQF, ZLACPY, ZLARTG, ZLASET, ZROT,
     $                  ZUNGQR, ZUNMQR, ZUNMRQ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, DCMPLX, DCONJG, INT, MAX
C
C     .. Executable Statements ..
C
      INFO   = 0
      WANTQ  = ( LSAME( COMPQ, 'V' ) ) .OR. ( LSAME( COMPQ, 'I' ) )
      INITQ  = LSAME( COMPQ, 'I' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.WANTQ .AND. .NOT.INITQ .AND. .NOT.LSAME(COMPQ,'N') ) THEN
         INFO = -1
      ELSE IF ( K .LT. 1 ) THEN
         INFO = -2
      ELSE IF ( N .LT. 0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
      ELSE
         SOK = S(1).EQ.1
         DO 10  L = 2, K
            SOK = S(L).EQ.1 .OR. S(L).EQ.-1
   10    CONTINUE
         IF ( .NOT.SOK ) THEN
            INFO = -6
         ELSE IF ( LDA1 .LT. MAX(1, N) ) THEN
            INFO = -8
         ELSE IF ( LDA2 .LT. MAX(1, N) ) THEN
            INFO = -9
         ELSE IF ( WANTQ .AND. LDQ1.LT.MAX(1, N) ) THEN
            INFO = -11
         ELSE IF ( WANTQ .AND. LDQ2.LT.MAX(1, N) ) THEN
            INFO = -12
         ELSE IF ( LDWORK.LT.MAX(1,N) ) THEN
            INFO = -14
         ELSE IF ( LZWORK.LT.MAX(1,2*N) ) THEN
            INFO = -16
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'ZPGHRD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = DBLE(1)
         ZWORK(1) = CONE
         RETURN
      END IF
      WRKOPT = 2*N
C
C     Initialize Q if desired.
C
      IF ( INITQ )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ1 )
C
C     Transform A(2,:,:),...,A(K,:,:) to upper triangular form.
C
      DO 30  L = K, 2, -1
         IF ( S(L).EQ.1 ) THEN
C
C           Compute a QR Decomposition of A(:,:,L).
C
            IROWS = IHI + 1 - ILO
            ICOLS = N + 1 - ILO
            CALL ZGEQRF( IROWS, ICOLS, A(ILO,ILO,L), LDA1, ZWORK,
     $                   ZWORK(N+1), LZWORK-N, IERR )
            WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
C
C           Apply transformation to A(:,:,L-1).
C
            IF ( S(L-1).EQ.1 ) THEN
               CALL ZUNMQR( 'Right', 'No transpose', IHI, IROWS, IROWS,
     $                      A(ILO,ILO,L), LDA1, ZWORK, A(1,ILO,L-1),
     $                      LDA1, ZWORK(N+1), LZWORK-N, IERR )
               WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            ELSE
               CALL ZUNMQR( 'Left', 'Complex transpose', IROWS, ICOLS,
     $                      IROWS, A(ILO,ILO,L), LDA1, ZWORK, 
     $                      A(ILO,ILO,L-1), LDA1, ZWORK(N+1), LZWORK-N,
     $                      IERR )
               WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            END IF
C
C           Update transformation matrix Q(:,:,L).
C
            IF ( INITQ ) THEN
C               CALL ZLASET( 'Full', N, N, CZERO, CONE, Q(1,1,L), LDQ1 )
               IF ( IROWS.GT.1 ) THEN
                  CALL ZLACPY( 'Lower', IROWS-1, IROWS-1,
     $                         A(ILO+1,ILO,L), LDA1, Q(ILO+1,ILO,L),
     $                         LDQ1 )
                  CALL ZUNGQR( IROWS, IROWS, IROWS, Q(ILO,ILO,L), LDQ1,
     $                         ZWORK, ZWORK(N+1), LZWORK-N, IERR)
                  WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
               END IF
            ELSE IF ( WANTQ ) THEN
               CALL ZUNMQR( 'Right', 'No transpose', N, IROWS, IROWS,
     $                       A(ILO,ILO,L), LDA1, ZWORK, Q(1,ILO,L),
     $                       LDQ1, ZWORK(N+1), LZWORK-N, IERR )
               WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            END IF
            IF ( IROWS.GT.1 )
     $         CALL ZLASET( 'Low', IROWS-1, IROWS-1, CZERO, CZERO,
     $                      A(ILO+1,ILO,L), LDA1)
         ELSE
C
C           Compute an RQ Decomposition of A(:,:,L).
C
            ICOLS = IHI + 1 - ILO
            CALL ZGERQF( IHI, ICOLS, A(1,ILO,L), LDA1, ZWORK,
     $                   ZWORK(N+1), LZWORK-N, IERR )
            WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
C
C           Apply transformation to A(:,:,L-1).
C
            IF ( S(L-1).EQ.1 ) THEN
               CALL ZUNMRQ( 'Right', 'Complex Transpose', IHI, ICOLS,
     $                      ICOLS, A(ILO,ILO,L), LDA1, ZWORK,
     $                      A(1,ILO,L-1), LDA1, ZWORK(N+1), LZWORK-N,
     $                      IERR )
               WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            ELSE
               CALL ZUNMRQ( 'Left', 'No transpose', ICOLS, N+1-ILO,
     $                      ICOLS, A(ILO,ILO,L), LDA1, ZWORK,
     $                      A(ILO,ILO,L-1), LDA1, ZWORK(N+1), LZWORK-N,
     $                      IERR )
               WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            END IF
C
C           Update transformation matrix Q(:,:,L).
C
            IF ( INITQ )
     $         CALL ZLASET( 'Full', N, N, CZERO, CONE, Q(1,1,L), LDQ1 )
            IF ( INITQ.OR.WANTQ )
     $         CALL ZUNMRQ( 'Right', 'Complex transpose', N, ICOLS,
     $                      ICOLS, A(1,1,L), LDA1, ZWORK, Q(1,ILO,L),
     $                      LDQ1, ZWORK(N+1), LZWORK-N, IERR )
            WRKOPT = MAX( WRKOPT, N + INT( ZWORK(N+1) ) )
            IF ( ICOLS.GT.1 )
     $         CALL ZLASET( 'Low', ICOLS-1, ICOLS-1, CZERO, CZERO,
     $                      A(ILO+1,ILO,L), LDA1 )
         END IF
   30 CONTINUE
C
C     Reduce A(:,:,1) to upper Hessenberg form.
C
      DO 110  JCOL = ILO, IHI - 2
C
C        Annihilate all elements below A(JCOL+1,JCOL,1).
C

         DO 40  JROW = IHI, JCOL + 2, -1
C
            TEMP = A(JROW-1,JCOL,1)
            CALL ZLARTG( TEMP, A(JROW,JCOL,1), CS, SN,
     $                   A(JROW-1,JCOL,1) )
            A(JROW,JCOL,1) = CZERO
            CALL ZROT( N-JCOL, A(JROW-1,JCOL+1,1), LDA1,
     $                 A(JROW,JCOL+1,1), LDA1, CS, SN )
            DWORK(JROW) = CS
            ZWORK(JROW) = SN
   40    CONTINUE
C
         IF ( WANTQ ) THEN
            DO 50  JROW = IHI, JCOL + 2, -1
               CALL ZROT( N, Q(1,JROW-1,1), 1, Q(1,JROW,1), 1,
     $                    DWORK(JROW), DCONJG( ZWORK(JROW) ) )
   50       CONTINUE
         END IF
C
C        Propagate transformations through A(:,:,K),...,A(:,:,2).
C
         DO 90  L = K, 2, -1
C
            IF ( S(L).EQ.1 ) THEN
               DO 60  JROW = IHI, JCOL + 2, -1
                  CALL ZROT( JROW, A(1,JROW-1,L), 1, A(1,JROW,L), 1,
     $                       DWORK(JROW), DCONJG( ZWORK(JROW) ) )
                  TEMP = A(JROW-1,JROW-1,L)
                  CALL ZLARTG( TEMP, A(JROW,JROW-1,L), CS, SN,
     $                         A(JROW-1,JROW-1,L) )
                  A(JROW,JROW-1,L) = CZERO
                  CALL ZROT( N-JROW+1, A(JROW-1,JROW,L), LDA1,
     $                       A(JROW,JROW,L), LDA1, CS, SN )
                  DWORK(JROW) = CS
                  ZWORK(JROW) = SN
   60          CONTINUE
            ELSE
               DO 70  JROW = IHI, JCOL + 2, -1
                  CALL ZROT( N+2-JROW, A(JROW-1,JROW-1,L), LDA1,
     $                       A(JROW,JROW-1,L), LDA1, DWORK(JROW),
     $                       ZWORK(JROW) )
                  TEMP = A(JROW,JROW,L)
                  CALL ZLARTG( TEMP, A(JROW,JROW-1,L), CS, SN,
     $                         A(JROW,JROW,L) )
                  A(JROW,JROW-1,L) = CZERO
                  CALL ZROT( JROW-1, A(1,JROW,L), 1, A(1,JROW-1,L), 1,
     $                       CS, SN )
                  DWORK(JROW) = CS
                  ZWORK(JROW) = -SN
   70          CONTINUE
            END IF
C
            IF ( WANTQ ) THEN
               DO 80  JROW = IHI, JCOL + 2, -1
                  CALL ZROT( N, Q(1,JROW-1,L), 1, Q(1,JROW,L), 1,
     $                       DWORK(JROW), DCONJG(ZWORK(JROW)) )
   80          CONTINUE
            END IF
   90    CONTINUE
C
C        Apply transformations to A(:,:,1).
C
         DO 100  JROW = IHI, JCOL + 2, -1
            CALL ZROT( IHI, A(1,JROW-1,L), 1, A(1,JROW,L), 1,
     $                 DWORK(JROW), DCONJG( ZWORK(JROW) ) )
  100    CONTINUE
  110 CONTINUE
      DWORK(1) = DBLE( N )
      ZWORK(1) = DCMPLX( WRKOPT, 0 )
      RETURN
C *** Last line of ZPGHRD ***
      END
