      SUBROUTINE ZPGEQZ( JOB, COMPQ, K, N, ILO, IHI, S, A, LDA1, LDA2, 
     $                   ALPHA, BETA, SCAL, Q, LDQ1, LDQ2, DWORK,
     $                   LDWORK, ZWORK, LZWORK, INFO )
	IMPLICIT NONE
C
C     PURPOSE
C
C     ZPGEQZ implements a single-shift version of the periodic QZ
C     method for finding the eigenvalues of the complex generalized
C     matrix product
C
C                                 S(2)                 S(K)
C          A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K).
C
C     In addition, A may be reduced to periodic Schur form by unitary
C     transformations: all factors A(:,:,i) become upper triangular.
C
C     If COMPQ = 'V' or COMPZ = 'I', then the unitary factors are
C     computed and stored in the array Q so that for S(I) = 1,
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
C     JOB     (input) CHARACTER*1
C             = 'E': compute only the eigenvalues; A will not
C                    necessarily be put into periodic Schur form.
C             = 'S': put A into periodic Schur form, as well
C                    as computing the eigenvalues contained in ALPHAR,
C                    ALPHAI, BETA and SCAL.
C
C     COMPQ   (input) CHARACTER*1
C             = 'N': do not modify Q.
C             = 'V': modify the array Q by the unitary transformations
C                    that are applied to the matrices in A to reduce them
C                    to periodic Schur form.
C             = 'I': like COMPQ='V', except that each matrix in Q will
C                    be initialized to the identity first.
C
C     Input/Output Parameters
C
C     K       (input)  INTEGER
C             The number of factors.  K >= 1.
C
C     N       (input)  INTEGER
C             The order of each factor in A.  N >= 0.
C
C     ILO     (input)  INTEGER
C     IHI     (input)  INTEGER
C             It is assumed that each factor in A is already upper
C             triangular in rows and columns 1:ILO-1 and IHI+1:N.
C             1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
C
C     S       (input)  INTEGER array, dimension (K)
C             The leading K elements of this array must contain the
C             signatures of the factors. Each entry in S must be either
C             1 or -1. By definition, S(1) must be set to 1.
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K)
C             On entry, the leading N-by-N-by-K part of this array
C             must contain the factors in upper Hessenberg-triangular
C             form, that is, A(:,:,1) is upper Hessenberg and the other
C             factors are upper triangular.
C             On exit, if JOB = 'S' and INFO = 0, the leading
C             N-by-N-by-K part of this array contains the factors of
C             A in periodic Schur form. All factors are reduched to
C             upper triangular form and, moreover, A(:,:,2),...,
C             A(:,:,K) are normalized so that their diagonals contain
C             nonnegative real numbers.
C             On exit, if JOB = 'E', then the leading N-by-N-by-K part
C             of this array contains meaningless elements.
C
C     LDA1    (input) INTEGER
C             The first leading dimension of A.  LDA1 >= MAX(1,N).
C
C     LDA2    (input) INTEGER
C             The second leading dimension of A.  LDA2 >= MAX(1,N).
C
C     ALPHA   (output) COMPLEX*16 array, dimension (N)
C             On exit, if INFO = 0, the leading N elements of this
C             array contain the scaled eigenvalues of A. The i-th
C             eigenvalue of A is given by
C
C             ALPHA(I) / BETA(I) * BASE**(SCAL(I)),
C
C             where 1.0 <= ABS(ALPHA(I)) < BASE and BASE is the machine
C             base (normally 2.0).
C
C     BETA    (output) COMPLEX*16 array, dimension (N)
C             On exit, if INFO = 0, the leading N elements of this
C             array contain indicators for infinite eigenvalues. That
C             is, if BETA(I) = 0.0, then the i-th eigenvalue is
C             infinite. Otherwise BETA(I) is set to 1.0.
C
C     SCAL    (output) INTEGER array, dimension (N)
C             On exit, if INFO = 0, the leading N elements of this
C             array contain the scaling parameters for the eigenvalues
C             of A.
C
C     Q       (input/output) COMPLEX*16 array, dimension (LDQ1,LDQ2,K)
C             On entry, if COMPQ = 'V', the leading N-by-N-by-K part
C             of this array must contain the initial unitary factors
C             as described in (1)-(2).
C             On exit, if COMPQ = 'V' or COMPQ = 'I', the leading
C             N-by-N-by-K part of this array contains the modified
C             orthogonal factors as described in (1)-(2).
C
C     LDQ1    (input)  INTEGER
C             The first leading dimension of Q.  LDQ1 >= MAX(1,N).
C
C     LDQ2    (input)  INTEGER
C             The second leading dimension of Q.  LDQ2 >= MAX(1,N).
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
C             On exit, if INFO = 0, ZWORK(1) returns the minimal value
C             of LZWORK.
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.  LZWORK >= MAX(1,N). 
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0       : succesful exit;
C             < 0       : if INFO = -i, the i-th argument had an
C                         illegal value;
C             = 1,..,N  : the periodic QZ iteration did not converge.
C                         A is not in periodic Schur form, but
C                         ALPHA(I), BETA(I) and SCAL(I), for
C                         I = INFO+1,...,N should be correct.
C
C     METHOD
C
C     A slightly modified version of the periodic QZ algorithm is
C     used. For more details see [2].
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
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16        CONE, CZERO
      PARAMETER         ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                  CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Scalar Arguments ..
      CHARACTER*1       COMPQ, JOB
      INTEGER           K, ILO, IHI, INFO, LDA1, LDA2, LDQ1, LDQ2,
     $                  LDWORK, LZWORK, N
C     .. Array Arguments ..
      INTEGER           S(*), SCAL(*)
	DOUBLE PRECISION  DWORK(*)
      COMPLEX*16        A(LDA1, LDA2, *), ALPHA(*), BETA(*),
     $                  Q(LDQ1, LDQ2, *), ZWORK(*)
C     .. Local Scalars ..
      LOGICAL           INITQ, LSCHR, WANTQ, SOK
      INTEGER           IFIRST, IFRSTM, IITER, ILAST, ILASTM, IN, J, J1,
     $                  JDEF, JITER, JLO, L, LDEF, LN, MAXIT, NTRA,
     $                  ZITER
      DOUBLE PRECISION  ABST, BASE, CS, SAFMIN, SAFMAX, SMLNUM, ULP,
     $                  TOL
	COMPLEX*16        SN, TEMP
C     .. Local Arrays ..
      INTEGER           ISEED(4)
	COMPLEX*16        RND(4)
C     .. External Functions ..
      LOGICAL           LSAME
	DOUBLE PRECISION  DLAMCH, ZLANHS
      EXTERNAL          DLAMCH, LSAME, ZLANHS
C     .. External Subroutines ..
      EXTERNAL          DLABAD, ZLAPR1, ZLARNV, ZLARTG, ZROT, ZSCAL
      INTRINSIC         ABS, DBLE, DCMPLX, DCONJG, INT, LOG, MAX
C
C     .. Executable Statements ..
C
      INFO = 0
      LSCHR = LSAME( JOB,'S' )
      WANTQ = LSAME( COMPQ,'V' ).OR.LSAME( COMPQ,'I' )
      INITQ = LSAME( COMPQ,'I' )
C
C     Check the scalar input parameters.
C
      IF ( .NOT. ( LSCHR .OR. LSAME( JOB,'E' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.WANTQ .AND. .NOT.INITQ .AND.
     $         .NOT.LSAME(COMPQ,'N') ) THEN
         INFO = -2
      ELSE IF ( K.LT.1 ) THEN
         INFO = -3
      ELSE IF ( N.LT.0 ) THEN
         INFO = -4
      ELSE IF ( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF ( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN   
         INFO = -6
      ELSE
	   SOK = S(1).EQ.1
	   DO 10  L = 2, K
	      SOK = S(L).EQ.1 .OR. S(L).EQ.-1
   10    CONTINUE
         IF ( .NOT.SOK ) THEN
	      INFO = -7
         ELSE IF ( LDA1 .LT. MAX(1, N) ) THEN
            INFO = -9
         ELSE IF ( LDA2 .LT. MAX(1, N) ) THEN	
            INFO = -10
         ELSE IF ( WANTQ .AND. LDQ1.LT.MAX(1, N) ) THEN
            INFO = -15
         ELSE IF ( WANTQ .AND. LDQ2.LT.MAX(1, N) ) THEN
            INFO = -16
         ELSE IF ( LDWORK.LT.MAX(1,N) ) THEN
            INFO = -18
	   ELSE IF ( LZWORK.LT.MAX(1,N) ) THEN
            INFO = -20
         END IF
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPGEQZ', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ONE
         ZWORK(1) = CONE
         RETURN
      END IF
C
C     Initialize Q.
C
      IF ( INITQ ) THEN
         DO 20  L = 1, K
            CALL ZLASET( 'Full', N, N, CZERO, CONE, Q(1,1,L), LDQ1 )
   20    CONTINUE
      END IF
C
C     Machine Constants
C
      IN = IHI + 1 - ILO
      SAFMIN = DLAMCH( 'SafeMinimum' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'Precision' )
      CALL DLABAD( SAFMIN, SAFMAX )
      SMLNUM = SAFMIN*( IN / ULP )
      BASE = DLAMCH( 'Base' )
      IF ( K.GE.INT( LOG( DLAMCH('Underflow') ) / LOG( ULP ) ) ) THEN
C
C        Start Iteration with a controlled zero shift.
C
         ZITER = -1
      ELSE
         ZITER = 0
      END IF
C
C     Set Eigenvalues IHI+1:N
C
      DO 30  J = IHI + 1, N
	   CALL ZLAPR1( BASE, K, S, A(J,J,1), LDA1*LDA2, ALPHA(J),
     $                BETA(J), SCAL(J) )
   30 CONTINUE
C
C     If IHI < ILO, skip QZ steps
C
      If ( IHI.LT.ILO )  GOTO 470
C
C     MAIN PERIODIC QZ ITERATION LOOP
C     
C     Initialize dynamic indices
C
C     Eigenvalues ILAST+1:N have been found.
C        Column operations modify rows IFRSTM:whatever.
C        Row operations modify columns whatever:ILASTM.
C
C     If only eigenvalues are being computed, then
C        IFRSTM is the row of the last splitting row above row ILAST;
C        this is always at least ILO.
C     IITER counts iterations since the last eigenvalue was found,
C        to tell when to use an observed zero or random shift.
C     MAXIT is the maximum number of QZ sweeps allowed.
C
      ILAST = IHI
      IF ( LSCHR ) THEN
         IFRSTM = 1
         ILASTM = N
      ELSE
         IFRSTM = ILO
         ILASTM = IHI
      END IF
      IITER = 0
      ISEED(1) = 1
      ISEED(2) = 0
      ISEED(3) = 0
      ISEED(4) = 0
      MAXIT = 30 * IN
C
      DO  460 JITER = 1, MAXIT
C
C        Special Case: ILAST = ILO
C
         IF ( ILAST.EQ.ILO )  GOTO 400
C
C        **************************************************************
C        *                     CHECK FOR DEFLATION                    *
C        **************************************************************
C
C        Test 1:  Deflation in the Hessenberg matrix.
C
         JLO = ILO
         DO 40  J = ILAST, ILO+1, -1
            TOL = ABS( A(J-1,J-1,1) ) + ABS( A(J,J,1) )
            IF ( TOL.EQ.ZERO )
     $         TOL = ZLANHS( '1',J-ILO+1, A(ILO,ILO,1), LDA1,
     $                       DWORK(1) )
            TOL = MAX( ULP*TOL, SMLNUM )
            IF ( ABS( A(J,J-1,1) ).LE.TOL ) THEN
               A(J,J-1,1) = CZERO
               JLO = J
               IF ( J.EQ.ILAST )  GOTO 400
               GOTO 50
            END IF
   40    CONTINUE
C
   50    CONTINUE
C
C        Test 2:  Deflation in the triangular matrices with index 1.
C
         DO 70  LDEF = 2,K
            IF ( S(LDEF).EQ.1 ) THEN
               DO 60  J = ILAST, JLO, -1
                  IF ( J.EQ.ILAST ) THEN
                     TOL = ABS( A(J-1,J,LDEF) )
                  ELSE IF ( J.EQ.JLO ) THEN
                     TOL = ABS( A(J,J+1,LDEF) )
                  ELSE
                     TOL = ABS( A(J-1,J,LDEF) ) + ABS( A(J,J+1,LDEF) )
                  END IF
                  IF ( TOL.EQ.ZERO )
     $               TOL = ZLANHS( '1', J-JLO+1, A(JLO,JLO,LDEF), LDA1,
     $                             DWORK(1) )
                  TOL = MAX( ULP*TOL, SMLNUM )
                  IF ( ABS( A(J,J,LDEF) ).LE.TOL ) THEN
                     A(J,J,LDEF) = CZERO
                     GOTO 180
                  END IF          
   60          CONTINUE
            END IF
   70    CONTINUE
C
C        Test 3:  Deflation in the triangular matrices with index -1.
C
         DO 90  LDEF = 2,K
            IF ( S(LDEF).EQ.-1 ) THEN
               DO 80  J = ILAST, JLO, -1
                  IF ( J.EQ.ILAST ) THEN
                     TOL = ABS( A(J-1,J,LDEF) )
                  ELSE IF ( J.EQ.JLO ) THEN
                     TOL = ABS( A(J,J+1,LDEF) )
                  ELSE
                     TOL = ABS( A(J-1,J,LDEF) ) + ABS( A(J,J+1,LDEF) )
                  END IF
                  IF ( TOL.EQ.ZERO )
     $                TOL = ZLANHS( '1', J-JLO+1, A(JLO,JLO,LDEF), LDA1,
     $                              DWORK(1) )
                  TOL = MAX( ULP*TOL, SMLNUM )
                  IF ( ABS( A(J,J,LDEF) ).LE.TOL ) THEN
                     A(J,J,LDEF) = CZERO
                     GOTO 330
                  END IF          
   80          CONTINUE
            END IF
   90    CONTINUE
C
C        Test 4:  Controlled zero shift.
C
         IF ( ZITER.GE.7 .OR. ZITER.LT.0 ) THEN
C
C           Make Hessenberg matrix upper triangular.
C
            DO 100 J = JLO, ILAST-1
               TEMP = A(J,J,1)
               CALL ZLARTG( TEMP, A(J+1,J,1), CS, SN, A(J,J,1) )
               A(J+1,J,1) = CZERO
               CALL ZROT( ILASTM-J, A(J,J+1,1), LDA1,
     $                    A(J+1,J+1,1), LDA1, CS, SN )
               DWORK(J) = CS
               ZWORK(J) = SN
  100       CONTINUE
            IF ( WANTQ ) THEN
               DO 110  J = JLO, ILAST-1
                  CALL ZROT( N, Q(1,J,1), 1, Q(1,J+1,1), 1, DWORK(J),
     $                      DCONJG( ZWORK(J) ) )
  110          CONTINUE
            END IF
C
C           Propagate Transformations back to A_1.
C
            DO 150  L = K, 2, -1
               IF ( S(L).EQ.1 ) THEN
                  DO 120  J = JLO, ILAST-1
                     CS = DWORK(J)
                     SN = ZWORK(J)
                     IF ( SN.NE.CZERO ) THEN
                        CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,L), 1,
     $                             A(IFRSTM,J+1,L), 1, CS,
     $                             DCONJG( SN ) )
C
C                       Check for deflation
C
                        TOL = ABS( A(J,J,L) ) + ABS( A(J+1,J+1,L) )
                        IF ( TOL.EQ.ZERO )
     $                     TOL = ZLANHS( '1',J-JLO+2, A(JLO,JLO,L),
     $                                   LDA1, DWORK(1) )
                        TOL = MAX( ULP*TOL, SMLNUM )
                        IF ( ABS( A(J+1,J,L) ).LE.TOL ) THEN
                           CS = ONE
                           SN = CZERO
                           A(J+1,J,L) = CZERO
                        END IF
C
                        TEMP = A(J,J,L)
                        CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN,
     $                               A(J,J,L) )
                        A(J+1,J,L) = CZERO
                        CALL ZROT( ILASTM-J, A(J,J+1,L), LDA1,
     $                             A(J+1,J+1,L), LDA1, CS, SN )
                     END IF
                     DWORK(J) = CS
                     ZWORK(J) = SN
  120             CONTINUE
               ELSE
                  DO 130  J = JLO, ILAST-1
                     CS = DWORK(J)
                     SN = ZWORK(J)
                     IF ( SN.NE.CZERO ) THEN
                        CALL ZROT( ILASTM-J+1, A(J,J,L), LDA1,
     $                             A(J+1,J,L), LDA1, CS, SN )
C
C                       Check for deflation
C
                        TOL = ABS( A(J,J,L) ) + ABS( A(J+1,J+1,L) )
                        IF ( TOL.EQ.ZERO )
     $                     TOL = ZLANHS( '1',J-JLO+2, A(JLO,JLO,L),
     $                                   LDA1, DWORK(1) )
                        TOL = MAX( ULP*TOL, SMLNUM )
                        IF ( ABS( A(J+1,J,L) ).LE.TOL ) THEN
                           CS = ONE
                           SN = CZERO
                           A(J+1,J,L) = CZERO
                        END IF
C
                        TEMP = A(J+1,J+1,L)
                        CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN,
     $                               A(J+1,J+1,L) )
                        A(J+1,J,L) = CZERO
                        CALL ZROT( J+1-IFRSTM, A(IFRSTM,J+1,L), 1,
     $                             A(IFRSTM,J,L), 1, CS, SN )
                        DWORK(J) = CS
                        ZWORK(J) = -SN
                     END IF
  130             CONTINUE
               END IF
C
               IF ( WANTQ ) THEN
                  DO 140  J = JLO, ILAST-1
                     CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, DWORK(J),
     $                          DCONJG( ZWORK(J) ) )
  140             CONTINUE
               END IF
  150       CONTINUE

C
C           Apply the transformations to the right hand side of the
C           Hessenberg factor.
C
            ZITER = 0
            DO 160  J = JLO, ILAST-1
               CS = DWORK(J)
               SN = ZWORK(J)
               CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,1), 1,
     $                    A(IFRSTM,J+1,1), 1, CS, DCONJG( SN ) )
	         IF ( SN.EQ.CZERO )
     $            ZITER = 1
  160       CONTINUE
C
C           No QZ iteration.
C
            GOTO 450
         END IF
C
C        **************************************************************
C        *                     HANDLE DEFLATIONS                      *
C        **************************************************************
C
C        Case I: Deflation occurs in the Hessenberg matrix. The QZ
C                iteration is only applied to the JLO:ILAST part.
C
  170    CONTINUE
         IFIRST = JLO
C
C        Go to the periodic QZ steps
C
         GOTO 410
C
C        Case II: Deflation occurs in a triangular matrix with index 1.
C
C        Do an unshifted periodic QZ step.
C
  180    JDEF = J
         DO 190  J = JLO, JDEF-1
            TEMP = A(J,J,1)
            CALL ZLARTG( TEMP, A(J+1,J,1), CS, SN, A(J,J,1) )
            A(J+1,J,1) = CZERO
            CALL ZROT( ILASTM-J, A(J,J+1,1), LDA1, A(J+1,J+1,1), LDA1,
     $                 CS, SN )
            DWORK(J) = CS
            ZWORK(J) = SN
  190    CONTINUE
         IF ( WANTQ ) THEN
            DO 200  J = JLO, JDEF-1
               CALL ZROT( N, Q(1,J,1), 1, Q(1,J+1,1), 1, DWORK(J),
     $                    DCONJG( ZWORK(J) ) )
  200       CONTINUE
         END IF
C
C        Propagate the transformations through the triangular matrices.
C        Due to the zero element on the diagonal of the LDEF-th
C        factor the number of transformations drops by one.
C
         DO 240  L = K, 2, -1
            IF ( L.LT.LDEF ) THEN
               NTRA = JDEF-2
            ELSE
               NTRA = JDEF-1
            END IF
            IF ( S(L).EQ.1 ) THEN
               DO 210  J = JLO, NTRA
                  CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,L), 1,
     $                       A(IFRSTM,J+1,L), 1, DWORK(J),
     $                       DCONJG( ZWORK (J) ) )
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J,J,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( ILASTM-J, A(J,J+1,L), LDA1,
     $                       A(J+1,J+1,L), LDA1, CS, SN )
                  DWORK(J) = CS
                  ZWORK(J) = SN
  210          CONTINUE
            ELSE
               DO 220  J = JLO, NTRA
                  CALL ZROT( ILASTM-J+1, A(J,J,L), LDA1, A(J+1,J,L),
     $                       LDA1, DWORK(J), ZWORK(J) )
                  TEMP = A(J+1,J+1,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J+1,J+1,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( J+1-IFRSTM, A(IFRSTM,J+1,L), 1,
     $                       A(IFRSTM,J,L), 1, CS, SN )
                  DWORK(J) = CS
                  ZWORK(J) = -SN
  220          CONTINUE
            END IF
            IF ( WANTQ ) THEN
               DO 230  J = JLO, NTRA
                  CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, DWORK(J),
     $                       DCONJG( ZWORK(J) ) )
  230          CONTINUE
            END IF
  240    CONTINUE
C
C        Apply the transformations to the right hand side of the
C        Hessenberg factor.
C
         DO 250  J = JLO, JDEF-2
            CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,1), 1, A(IFRSTM,J+1,1),
     $                 1, DWORK(J), DCONJG( ZWORK(J) ) )
  250    CONTINUE
C
C        Do an unshifted periodic ZQ step.
C
         DO 260  J = ILAST, JDEF+1, -1
            TEMP = A(J,J,1)
            CALL ZLARTG( TEMP, A(J,J-1,1), CS, SN, A(J,J,1) )
            A(J,J-1,1) = CZERO
            CALL ZROT( J-IFRSTM, A(IFRSTM,J,1), 1,
     $                 A(IFRSTM,J-1,1), 1, CS, SN )
            DWORK(J) = CS
            ZWORK(J) = -SN
  260    CONTINUE
         IF ( WANTQ ) THEN
            DO 270  J = ILAST, JDEF+1, -1
               CALL ZROT( N, Q(1,J-1,2), 1, Q(1,J,2),
     $                    1, DWORK(J), DCONJG( ZWORK(J) ) )
  270       CONTINUE
         END IF

C
C        Propagate the transformations through the triangular matrices.
C
         DO 310  L = 2, K
            IF ( L.GT.LDEF ) THEN
               NTRA = JDEF+2
            ELSE
               NTRA = JDEF+1
            END IF
            IF ( S(L).EQ.-1 ) THEN
               DO 280  J = ILAST, NTRA, -1
                  CS = DWORK(J)
                  SN = ZWORK(J)
                  CALL ZROT( J+1-IFRSTM, A(IFRSTM,J-1,L), 1,
     $                       A(IFRSTM,J,L), 1, CS, DCONJG( SN ) )
                  TEMP = A(J-1,J-1,L)
                  CALL ZLARTG( TEMP, A(J,J-1,L), CS, SN, A(J-1,J-1,L) )
                  A(J,J-1,L) = CZERO
                  CALL ZROT( ILASTM-J+1, A(J-1,J,L), LDA1, A(J,J,L),
     $                       LDA1, CS, SN )
                  DWORK(J) = CS
                  ZWORK(J) = SN
  280          CONTINUE
            ELSE
               DO 290  J = ILAST, NTRA, -1
                  CALL ZROT( ILASTM-J+2, A(J-1,J-1,L), LDA1,
     $                       A(J,J-1,L), LDA1, DWORK(J), ZWORK(J) )
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J,J-1,L), CS, SN, A(J,J,L) )
                  A(J,J-1,L) = CZERO
                  CALL ZROT( J-IFRSTM, A(IFRSTM,J,L), 1,
     $                       A(IFRSTM,J-1,L), 1, CS, SN )
                  DWORK(J) = CS
                  ZWORK(J) = -SN
  290          CONTINUE
            END IF
            IF ( WANTQ ) THEN
               LN = L+1
               IF ( L.EQ.K )  LN = 1
               DO 300  J = ILAST, NTRA, -1
                  CALL ZROT( N, Q(1,J-1,LN), 1, Q(1,J,LN), 1, DWORK(J),
     $                       DCONJG( ZWORK(J) ) )
  300          CONTINUE
            END IF
  310    CONTINUE
C
C        Apply the transformations to the left hand side of the
C        Hessenberg factor.
C
         DO 320  J = ILAST, JDEF+2, -1
            CALL ZROT( ILASTM-J+2, A(J-1,J-1,1), LDA1, A(J,J-1,1),
     $                 LDA1, DWORK(J), ZWORK(J) )
  320    CONTINUE
C
C        No QZ iteration.
C
         GOTO 450
C
C        Case III: Deflation occurs in a triangular matrix with
C                  index -1.
C
  330    CONTINUE
         JDEF = J
         IF ( JDEF.GT.( (ILAST-JLO+1)/2 ) ) THEN
C
C           Chase the zero downwards to the last position
C
            DO 350  J1 = JDEF, ILAST-1
               J = J1
               TEMP = A(J,J+1,LDEF)
               CALL ZLARTG( TEMP, A(J+1,J+1,LDEF), CS, SN,
     $                      A(J,J+1,LDEF) )
               A(J+1,J+1,LDEF) = CZERO
               CALL ZROT( ILASTM-J-1, A(J,J+2,LDEF), LDA1,
     $                    A(J+1,J+2,LDEF), LDA1, CS, SN )
               LN = LDEF+1
               IF ( LDEF.EQ.K )  LN = 1
               IF ( WANTQ ) THEN
                  CALL ZROT( N, Q(1,J,LN), 1, Q(1,J+1,LN), 1, CS,
     $                       DCONJG( SN ) )
	         END IF
               DO 340  L = 1, K-1
                  IF ( LN.EQ.1 ) THEN
                     CALL ZROT( ILASTM-J+2, A(J,J-1,LN), LDA1,
     $                          A(J+1,J-1,LN), LDA1, CS, SN )
                     TEMP = A(J+1,J,LN)
                     CALL ZLARTG( TEMP, A(J+1,J-1,LN), CS, SN,
     $                            A(J+1,J,LN) )
                     A(J+1,J-1,LN) = CZERO
                     CALL ZROT( J-IFRSTM+1, A(IFRSTM,J,LN), 1,
     $                          A(IFRSTM,J-1,LN), 1, CS, SN )
	               SN = -SN
                     J = J - 1
                  ELSE IF ( S(LN).EQ.1 ) THEN
                     CALL ZROT( ILASTM-J+1, A(J,J,LN), LDA1,
     $                          A(J+1,J,LN), LDA1, CS, SN )
                     TEMP = A(J+1,J+1,LN)
                     CALL ZLARTG( TEMP, A(J+1,J,LN), CS, SN,
     $                            A(J+1,J+1,LN) )
                     A(J+1,J,LN) = CZERO
                     CALL ZROT( J-IFRSTM+1, A(IFRSTM,J+1,LN), 1,
     $                          A(IFRSTM,J,LN), 1, CS, SN )
                     SN = -SN
                  ELSE
                     CALL ZROT( J-IFRSTM+2, A(IFRSTM,J,LN), 1,
     $                          A(IFRSTM,J+1,LN), 1, CS, DCONJG( SN ) )
                     TEMP = A(J,J,LN)
                     CALL ZLARTG( TEMP, A(J+1,J,LN), CS, SN, A(J,J,LN) )
                     A(J+1,J,LN) = CZERO
                     CALL ZROT( ILASTM-J, A(J,J+1,LN), LDA1,
     $                          A(J+1,J+1,LN), LDA1, CS, SN )
                  END IF
                  LN = LN+1
                  IF ( LN.GT.K )  LN = 1
                  IF ( WANTQ ) THEN
                     CALL ZROT( N, Q(1,J,LN), 1, Q(1,J+1,LN), 1, CS,
     $                          DCONJG( SN ) )
	            END IF
  340          CONTINUE
               CALL ZROT( J-IFRSTM+1, A(IFRSTM,J,LDEF), 1,
     $                    A(IFRSTM,J+1,LDEF), 1, CS, DCONJG( SN ) )
  350       CONTINUE
C
C           Deflate the last element in the Hessenberg matrix.
C
            J = ILAST
            TEMP = A(J,J,1)
            CALL ZLARTG( TEMP, A(J,J-1,1), CS, SN, A(J,J,1) )
            A(J,J-1,1) = CZERO
            CALL ZROT( J-IFRSTM, A(IFRSTM,J,1), 1,
     $                 A(IFRSTM,J-1,1), 1, CS, SN )
            SN = -SN
   !         IF ( WANTQ.NE.0 ) THEN
            IF ( WANTQ ) THEN
               CALL ZROT( N, Q(1,J-1,2), 1, Q(1,J,2), 1, CS,
     $                    DCONJG( SN ) )
            END IF
            DO 360  L = 2, LDEF-1
               IF ( S(L).EQ.-1 ) THEN
                  CALL ZROT( J+1-IFRSTM, A(IFRSTM,J-1,L), 1,
     $                       A(IFRSTM,J,L), 1, CS, DCONJG( SN ) )
                  TEMP = A(J-1,J-1,L)
                  CALL ZLARTG( TEMP, A(J,J-1,L), CS, SN,
     $                         A(J-1,J-1,L) )
                  A(J,J-1,L) = CZERO
                  CALL ZROT( ILASTM-J+1, A(J-1,J,L), LDA1,
     $                       A(J,J,L), LDA1, CS, SN )
               ELSE
                  CALL ZROT( ILASTM-J+2, A(J-1,J-1,L), LDA1,
     $                       A(J,J-1,L), LDA1, CS, SN )
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J,J-1,L), CS, SN,
     $                         A(J,J,L) )
                  A(J,J-1,L) = CZERO
                  CALL ZROT( J-IFRSTM, A(IFRSTM,J,L), 1,
     $                       A(IFRSTM,J-1,L), 1, CS, SN )
                  SN = -SN
               END IF
               IF ( WANTQ ) THEN
                  LN = L+1
                  IF ( L.EQ.K )  LN = 1
                  CALL ZROT( N, Q(1,J-1,LN), 1, Q(1,J,LN), 1, CS,
     $                       DCONJG( SN ) )
               END IF
 360        CONTINUE
            CALL ZROT( J+1-IFRSTM, A(IFRSTM,J-1,LDEF), 1,
     $                 A(IFRSTM,J,LDEF), 1, CS, DCONJG( SN ) )
         ELSE
C
C           Chase the zero upwards to the first position.
C
            DO 380  J1 = JDEF, JLO+1,-1
               J = J1
               TEMP = A(J-1,J,LDEF)
               CALL ZLARTG( TEMP, A(J-1,J-1,LDEF), CS, SN,
     $                      A(J-1,J,LDEF) )
               A(J-1,J-1,LDEF) = CZERO
               CALL ZROT( J-IFRSTM-1, A(IFRSTM,J,LDEF), 1,
     $                    A(IFRSTM,J-1,LDEF), 1, CS, SN )
               SN = -SN
               IF ( WANTQ ) THEN
                  CALL ZROT( N, Q(1,J-1,LDEF), 1, Q(1,J,LDEF), 1, CS,
     $                       DCONJG( SN ) )
	         END IF
               LN = LDEF - 1
               DO 370  L = 1, K-1
                  IF ( LN.EQ.1 ) THEN
                     CALL ZROT( J-IFRSTM+2, A(IFRSTM,J-1,LN), 1,
     $                          A(IFRSTM,J,LN), 1, CS, DCONJG( SN ) )
                     TEMP = A(J,J-1,LN)
                     CALL ZLARTG( TEMP, A(J+1,J-1,LN), CS, SN,
     $                            A(J,J-1,LN) )
                     A(J+1,J-1,LN) = CZERO
                     CALL ZROT( ILASTM-J+1, A(J,J,LN), LDA1,
     $                          A(J+1,J,LN), LDA1, CS, SN )
                     J = J + 1
                  ELSE IF ( S(LN).EQ.-1 ) THEN
                     CALL ZROT( ILASTM-J+2, A(J-1,J-1,LN), LDA1,
     $                          A(J,J-1,LN), LDA1, CS, SN )
                     TEMP = A(J,J,LN)
                     CALL ZLARTG( TEMP, A(J,J-1,LN), CS, SN,
     $                            A(J,J,LN) )
                     A(J,J-1,LN) = CZERO
                     CALL ZROT( J-IFRSTM, A(IFRSTM,J,LN), 1,
     $                          A(IFRSTM,J-1,LN), 1, CS, SN )
                     SN = -SN
                  ELSE
                     CALL ZROT( J-IFRSTM+1, A(IFRSTM,J-1,LN), 1,
     $                          A(IFRSTM,J,LN), 1, CS, DCONJG( SN ) )
                     TEMP = A(J-1,J-1,LN)
                     CALL ZLARTG( TEMP, A(J,J-1,LN), CS, SN,
     $                            A(J-1,J-1,LN) )
                     A(J,J-1,LN) = CZERO
                     CALL ZROT( ILASTM-J+1, A(J-1,J,LN), LDA1,
     $                          A(J,J,LN), LDA1, CS, SN )
                  END IF
	            IF ( WANTQ ) THEN
                     CALL ZROT( N, Q(1,J-1,LN), 1, Q(1,J,LN), 1, CS,
     $                          DCONJG( SN ) )
	            END IF
                  LN = LN - 1
                  IF ( LN.LE.0 )  LN = K
  370          CONTINUE
               CALL ZROT( ILASTM-J+1, A(J-1,J,LDEF), LDA1, A(J,J,LDEF),
     $                    LDA1, CS, SN )
  380       CONTINUE
C
C           Deflate the first element in the Hessenberg matrix.
C
            J = JLO
            TEMP = A(J,J,1)
            CALL ZLARTG( TEMP, A(J+1,J,1), CS, SN, A(J,J,1) )
            A(J+1,J,1) = CZERO
            CALL ZROT( ILASTM-J, A(J,J+1,1), LDA1, A(J+1,J+1,1),
     $                 LDA1, CS, SN )
            IF ( WANTQ ) THEN
               CALL ZROT( N, Q(1,J,1), 1, Q(1,J+1,1), 1, CS,
     $                    DCONJG( SN ) )
            END IF
            DO 390  L = K, LDEF+1, -1
               IF ( S(L).EQ.1 ) THEN
                  CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,L), 1,
     $                       A(IFRSTM,J+1,L), 1, CS, DCONJG( SN ) )
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J,J,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( ILASTM-J, A(J,J+1,L), LDA1,
     $                       A(J+1,J+1,L), LDA1, CS, SN )
               ELSE
                  CALL ZROT( ILASTM-J+1, A(J,J,L), LDA1,
     $                       A(J+1,J,L), LDA1, CS, SN )
                  TEMP = A(J+1,J+1,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN,
     $                         A(J+1,J+1,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( J+1-IFRSTM, A(IFRSTM,J+1,L), 1,
     $                       A(IFRSTM,J,L), 1, CS, SN )
                  SN = -SN
               END IF
               IF ( WANTQ ) THEN
                  CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, CS,
     $                       DCONJG( SN ) )
               END IF
  390       CONTINUE
            CALL ZROT( ILASTM-J, A(J,J+1,LDEF), LDA1, A(J+1,J+1,LDEF),
     $                 LDA1, CS, SN )
         END IF
C
C        No QZ iteration.
C
         GOTO 450
C
C        Special case: A 1x1 block splits off at the bottom
C
  400    CONTINUE
         CALL ZLAPR1( BASE, K, S, A(ILAST,ILAST,1), LDA1*LDA2,
     $                ALPHA(ILAST), BETA(ILAST), SCAL(ILAST) )
C
C        Go to next block - exit if finished.
C
         ILAST = ILAST - 1
         IF ( ILAST.LT.ILO )  GOTO 470
C
C        Reset iteration counters.
C
         IITER = 0
         IF ( ZITER.NE.-1 )  ZITER = 0
         IF ( .NOT.LSCHR ) THEN
            ILASTM = ILAST
            IF ( IFRSTM.GT.ILAST )  IFRSTM = ILO
         END IF
C
C        No QZ iteration.
C
         GOTO 450        
C  
C        **************************************************************
C        *                      PERIODIC QZ STEP                      *
C        **************************************************************
C
C        It is assumed that IFIRST < ILAST.
C
  410    CONTINUE
C
         IITER = IITER + 1
         ZITER = ZITER + 1
         IF( .NOT.LSCHR ) THEN
            IFRSTM = IFIRST
         END IF
C
C        Complex single shift.
C
         IF ( ( IITER / 10 )*10.EQ.IITER ) THEN
C
C           Exceptional shift.
C
            CALL ZLARNV( 2, ISEED, 2, RND )
            CALL ZLARTG( RND(1), RND(2), CS, SN, TEMP )
         ELSE
	      CALL ZLARTG( CONE, CONE, CS, SN, TEMP )
            DO 420  L = K, 2, -1
               IF ( S(L).EQ.1 ) THEN
	            CALL ZLARTG( A(IFIRST,IFIRST,L)*CS,
     $                         A(ILAST,ILAST,L)*DCONJG(SN),
     $                         CS, SN, TEMP )
               ELSE
	            CALL ZLARTG( A(ILAST,ILAST,L)*CS,
     $                         -A(IFIRST,IFIRST,L)*DCONJG(SN),
     $                         CS, SN, TEMP )
	            SN = -SN
               END IF
  420       CONTINUE
            CALL ZLARTG( A(IFIRST,IFIRST,1)*CS
     $                   -DCONJG(SN)*A(ILAST,ILAST,1),
     $                   A(IFIRST+1,IFIRST,1)*CS, CS, SN, TEMP )
         END IF
C
C        Do the sweeps.
C
         DO 440  J1 = IFIRST-1, ILAST-2
            J = J1 + 1
C
C           Create bulge if J1 = IFIRST - 1, otherwise chase bulge.
C
            IF ( J1.LT.IFIRST ) THEN
               CALL ZROT( ILASTM-J+1, A(J,J,1), LDA1, A(J+1,J,1), LDA1,
     $                    CS, SN )
            ELSE
               TEMP = A(J,J-1,1)
               CALL ZLARTG( TEMP, A(J+1,J-1,1), CS, SN, A(J,J-1,1) )
               A(J+1,J-1,1) = CZERO
               CALL ZROT( ILASTM-J+1, A(J,J,1), LDA1, A(J+1,J,1),
     $                    LDA1, CS, SN )
            END IF
            IF ( WANTQ ) THEN
               CALL ZROT( N, Q(1,J,1), 1, Q(1,J+1,1), 1, CS,
     $                    DCONJG( SN ) )
            END IF
C
C           Propagate rotation through AK, ..., A2 to A1.
C
            DO 430  L = K, 2, -1
               IF ( S(L).EQ.1 ) THEN
                  CALL ZROT( J+2-IFRSTM, A(IFRSTM,J,L), 1,
     $                       A(IFRSTM,J+1,L), 1, CS, DCONJG( SN ) )
                  TEMP = A(J,J,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J,J,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( ILASTM-J, A(J,J+1,L), LDA1,
     $                       A(J+1,J+1,L), LDA1, CS, SN )
               ELSE
                  CALL ZROT( ILASTM-J+1, A(J,J,L), LDA1, A(J+1,J,L),
     $                       LDA1, CS, SN )
                  TEMP = A(J+1,J+1,L)
                  CALL ZLARTG( TEMP, A(J+1,J,L), CS, SN, A(J+1,J+1,L) )
                  A(J+1,J,L) = CZERO
                  CALL ZROT( J+1-IFRSTM, A(IFRSTM,J+1,L), 1,
     $                       A(IFRSTM,J,L), 1, CS, SN )
                  SN = -SN
               END IF
               IF ( WANTQ ) THEN
                  CALL ZROT( N, Q(1,J,L), 1, Q(1,J+1,L), 1, CS,
     $                       DCONJG( SN ) )
               END IF
  430       CONTINUE
            CALL ZROT( MIN(J+2,ILASTM)-IFRSTM+1, A(IFRSTM,J,1), 1,
     $                 A(IFRSTM,J+1,1), 1, CS, DCONJG( SN ) )
  440    CONTINUE
C
C        End of iteration loop.
C
  450    CONTINUE
  460 CONTINUE
C
C     Drop through = non-convergence
C
      INFO = ILAST
      GO TO 550
C
C     Successful completion of all QZ steps
C
  470 CONTINUE
C
C     Set eigenvalues 1:ILO-1
C
      DO 480  J = 1,ILO-1
         CALL ZLAPR1( BASE, K, S, A(J,J,1), LDA1*LDA2, ALPHA(J),
     $                BETA(J), SCAL(J) )
  480 CONTINUE
      IF ( LSCHR ) THEN
C
C        Scale A(2,:,:) .. A(K,:,:).
C
         DO 540  L = K, 2, -1
	      IF ( S(L).EQ.1 )  THEN
	         DO 490 J = 1, N
                  ABST = ABS( A(J,J,L) )
                  IF ( ABST.GT.SAFMIN ) THEN
                     TEMP = DCONJG( A(J,J,L) / ABST )
                     A(J,J,L ) = ABST
                     IF ( J.LT.N )
     $                  CALL ZSCAL( N-J, TEMP, A(J,J+1,L), LDA1 )
	            ELSE
	               TEMP = CONE
	            END IF
                  ZWORK(J) = TEMP
  490          CONTINUE	
	      ELSE
               DO 500  J = 1, N
                  ABST = ABS( A(J,J,L) )
                  IF ( ABST.GT.SAFMIN ) THEN
                     TEMP = DCONJG( A(J,J,L) / ABST )
                     A(J,J,L ) = ABST
                     CALL ZSCAL( J-1, TEMP, A(1,J,L), 1 )
	            ELSE
	               TEMP = CONE
	            END IF
                  ZWORK(J) = DCONJG(TEMP)
  500          CONTINUE	
	      END IF
            IF ( WANTQ ) THEN
	         DO 510  J = 1, N
                  CALL ZSCAL( N, DCONJG( ZWORK(J) ), Q(1,J,L), 1 )
  510          CONTINUE
            END IF
         	IF ( S(L-1).EQ.1 )  THEN
	         DO 520  J = 1, N
	            CALL ZSCAL( J, DCONJG( ZWORK(J) ), A(1,J,L-1), 1 )
  520          CONTINUE
            ELSE
	         DO 530  J = 1, N
	            CALL ZSCAL( N-J+1, ZWORK(J), A(J,J,L-1), LDA1 )
  530          CONTINUE              
	      END IF
  540    CONTINUE
	END IF
      INFO = 0
C
  550 CONTINUE
C
      DWORK(1) = DBLE( N )
      ZWORK(1) = DCMPLX( N, 0 )
      RETURN     
C *** Last line of ZPGEQZ ***
      END
