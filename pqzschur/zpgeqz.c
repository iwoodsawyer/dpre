/* zpgeqz.f -- translated by f2c (version 20041007).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int zpgeqz(char *job, char *compq, integer *k, integer *n, 
	integer *ilo, integer *ihi, integer *s, doublecomplex *a, integer *
	lda1, integer *lda2, doublecomplex *alpha, doublecomplex *beta, 
	integer *scal, doublecomplex *q, integer *ldq1, integer *ldq2, 
	doublereal *dwork, integer *ldwork, doublecomplex *zwork, integer *
	lzwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_dim2, a_offset, q_dim1, q_dim2, q_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    double log(doublereal), z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, l, j1;
    static doublereal cs;
    static integer in, ln;
    static doublecomplex sn;
    static integer jlo;
    static doublecomplex rnd[4];
    static logical sok;
    static doublereal tol, ulp;
    static integer jdef, ldef;
    static doublereal base, abst;
    static integer ntra;
    static doublecomplex temp;
    extern /* Subroutine */ int zrot(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static integer iseed[4];
    extern logical lsame(char *, char *);
    static logical lschr;
    static integer iiter, ilast, jiter;
    extern /* Subroutine */ int zscal(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer maxit;
    static logical initq, wantq;
    static integer ziter;
    extern /* Subroutine */ int dlabad(doublereal *, doublereal *), zlapr1(
	    doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *);
    extern doublereal dlamch(char *);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla(char *, integer *);
    static doublereal safmax;
    static integer ilastm;
    extern doublereal zlanhs(char *, integer *, doublecomplex *, integer *, 
	    doublereal *);
    static integer ifirst;
    extern /* Subroutine */ int zlartg(doublecomplex *, doublecomplex *, 
	    doublereal *, doublecomplex *, doublecomplex *);
    static integer ifrstm;
    extern /* Subroutine */ int zlaset(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int zlarnv(integer *, integer *, integer *, 
	    doublecomplex *);


/*     PURPOSE */

/*     ZPGEQZ implements a single-shift version of the periodic QZ */
/*     method for finding the eigenvalues of the complex generalized */
/*     matrix product */

/*                                 S(2)                 S(K) */
/*          A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K). */

/*     In addition, A may be reduced to periodic Schur form by unitary */
/*     transformations: all factors A(:,:,i) become upper triangular. */

/*     If COMPQ = 'V' or COMPZ = 'I', then the unitary factors are */
/*     computed and stored in the array Q so that for S(I) = 1, */

/*                         H */
/*             Q(:,:,I)(in)   A(:,:,I)(in)   Q(:,:,MOD(I,K)+1)(in) */
/*                          H                                        (1) */
/*         =   Q(:,:,I)(out)  A(:,:,I)(out)  Q(:,:,MOD(I,K)+1)(out), */

/*     and for S(I) = -1, */

/*                                  H */
/*             Q(:,:,MOD(I,K)+1)(in)   A(:,:,I)(in)   Q(:,:,I)(in) */
/*                                   H                               (2) */
/*         =   Q(:,:,MOD(I,K)+1)(out)  A(:,:,I)(out)  Q(:,:,I)(out). */

/*     ARGUMEMTS */

/*     Mode Parameters */

/*     JOB     (input) CHARACTER*1 */
/*             = 'E': compute only the eigenvalues; A will not */
/*                    necessarily be put into periodic Schur form. */
/*             = 'S': put A into periodic Schur form, as well */
/*                    as computing the eigenvalues contained in ALPHAR, */
/*                    ALPHAI, BETA and SCAL. */

/*     COMPQ   (input) CHARACTER*1 */
/*             = 'N': do not modify Q. */
/*             = 'V': modify the array Q by the unitary transformations */
/*                    that are applied to the matrices in A to reduce them */
/*                    to periodic Schur form. */
/*             = 'I': like COMPQ='V', except that each matrix in Q will */
/*                    be initialized to the identity first. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of factors.  K >= 1. */

/*     N       (input)  INTEGER */
/*             The order of each factor in A.  N >= 0. */

/*     ILO     (input)  INTEGER */
/*     IHI     (input)  INTEGER */
/*             It is assumed that each factor in A is already upper */
/*             triangular in rows and columns 1:ILO-1 and IHI+1:N. */
/*             1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */

/*     S       (input)  INTEGER array, dimension (K) */
/*             The leading K elements of this array must contain the */
/*             signatures of the factors. Each entry in S must be either */
/*             1 or -1. By definition, S(1) must be set to 1. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K) */
/*             On entry, the leading N-by-N-by-K part of this array */
/*             must contain the factors in upper Hessenberg-triangular */
/*             form, that is, A(:,:,1) is upper Hessenberg and the other */
/*             factors are upper triangular. */
/*             On exit, if JOB = 'S' and INFO = 0, the leading */
/*             N-by-N-by-K part of this array contains the factors of */
/*             A in periodic Schur form. All factors are reduched to */
/*             upper triangular form and, moreover, A(:,:,2),..., */
/*             A(:,:,K) are normalized so that their diagonals contain */
/*             nonnegative real numbers. */
/*             On exit, if JOB = 'E', then the leading N-by-N-by-K part */
/*             of this array contains meaningless elements. */

/*     LDA1    (input) INTEGER */
/*             The first leading dimension of A.  LDA1 >= MAX(1,N). */

/*     LDA2    (input) INTEGER */
/*             The second leading dimension of A.  LDA2 >= MAX(1,N). */

/*     ALPHA   (output) COMPLEX*16 array, dimension (N) */
/*             On exit, if INFO = 0, the leading N elements of this */
/*             array contain the scaled eigenvalues of A. The i-th */
/*             eigenvalue of A is given by */

/*             ALPHA(I) / BETA(I) * BASE**(SCAL(I)), */

/*             where 1.0 <= ABS(ALPHA(I)) < BASE and BASE is the machine */
/*             base (normally 2.0). */

/*     BETA    (output) COMPLEX*16 array, dimension (N) */
/*             On exit, if INFO = 0, the leading N elements of this */
/*             array contain indicators for infinite eigenvalues. That */
/*             is, if BETA(I) = 0.0, then the i-th eigenvalue is */
/*             infinite. Otherwise BETA(I) is set to 1.0. */

/*     SCAL    (output) INTEGER array, dimension (N) */
/*             On exit, if INFO = 0, the leading N elements of this */
/*             array contain the scaling parameters for the eigenvalues */
/*             of A. */

/*     Q       (input/output) COMPLEX*16 array, dimension (LDQ1,LDQ2,K) */
/*             On entry, if COMPQ = 'V', the leading N-by-N-by-K part */
/*             of this array must contain the initial unitary factors */
/*             as described in (1)-(2). */
/*             On exit, if COMPQ = 'V' or COMPQ = 'I', the leading */
/*             N-by-N-by-K part of this array contains the modified */
/*             orthogonal factors as described in (1)-(2). */

/*     LDQ1    (input)  INTEGER */
/*             The first leading dimension of Q.  LDQ1 >= MAX(1,N). */

/*     LDQ2    (input)  INTEGER */
/*             The second leading dimension of Q.  LDQ2 >= MAX(1,N). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the minimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) returns the minimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK.  LZWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0       : succesful exit; */
/*             < 0       : if INFO = -i, the i-th argument had an */
/*                         illegal value; */
/*             = 1,..,N  : the periodic QZ iteration did not converge. */
/*                         A is not in periodic Schur form, but */
/*                         ALPHA(I), BETA(I) and SCAL(I), for */
/*                         I = INFO+1,...,N should be correct. */

/*     METHOD */

/*     A slightly modified version of the periodic QZ algorithm is */
/*     used. For more details see [2]. */

/*     REFERENCES */

/*     [1] Bojanczyk, A. and Golub, G. H. and Van Dooren, P. */
/*         The periodic Schur decomposition; algorithm and applications. */
/*         In Proc. SPIE Conference, pg. 31-42, vol. 1770, 1992. */

/*     [2] Kressner, D. */
/*         An efficient and reliable implementation of the periodic QZ */
/*         algorithm. In IFAC Workshop on Periodic Control Systems, 2001. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically backward stable. */
/*                                 3 */
/*     The algorithm requires 0(K N ) floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, Dec. 2002. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --s;
    a_dim1 = *lda1;
    a_dim2 = *lda2;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    --alpha;
    --beta;
    --scal;
    q_dim1 = *ldq1;
    q_dim2 = *ldq2;
    q_offset = 1 + q_dim1 * (1 + q_dim2);
    q -= q_offset;
    --dwork;
    --zwork;

    /* Function Body */
    *info = 0;
    lschr = lsame(job, "S");
    wantq = lsame(compq, "V") || lsame(compq, "I");
    initq = lsame(compq, "I",);

/*     Check the scalar input parameters. */

    if (! (lschr || lsame(job, "E"))) {
	*info = -1;
    } else if (! wantq && ! initq && ! lsame(compq, "N") {
	*info = -2;
    } else if (*k < 1) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ilo < 1) {
	*info = -5;
    } else if (*ihi > *n || *ihi < *ilo - 1) {
	*info = -6;
    } else {
	sok = s[1] == 1;
	i__1 = *k;
	for (l = 2; l <= i__1; ++l) {
	    sok = s[l] == 1 || s[l] == -1;
/* L10: */
	}
	if (! sok) {
	    *info = -7;
	} else if (*lda1 < max(1,*n)) {
	    *info = -9;
	} else if (*lda2 < max(1,*n)) {
	    *info = -10;
	} else if (wantq && *ldq1 < max(1,*n)) {
	    *info = -15;
	} else if (wantq && *ldq2 < max(1,*n)) {
	    *info = -16;
	} else if (*ldwork < max(1,*n)) {
	    *info = -18;
	} else if (*lzwork < max(1,*n)) {
	    *info = -20;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla("ZPGEQZ", &i__1);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	zwork[1].r = 1., zwork[1].i = 0.;
	return 0;
    }

/*     Initialize Q. */

    if (initq) {
	i__1 = *k;
	for (l = 1; l <= i__1; ++l) {
	    zlaset("Full", n, n, &c_b2, &c_b1, &q[(l * q_dim2 + 1) * q_dim1 
		    + 1], ldq1);
/* L20: */
	}
    }

/*     Machine Constants */

    in = *ihi + 1 - *ilo;
    safmin = dlamch("SafeMinimum");
    safmax = 1. / safmin;
    ulp = dlamch("Precision");
    dlabad(&safmin, &safmax);
    smlnum = safmin * (in / ulp);
    base = dlamch("Base");
    if (*k >= (integer) (log(dlamch("Underflow")) / log(ulp))) {

/*        Start Iteration with a controlled zero shift. */

	ziter = -1;
    } else {
	ziter = 0;
    }

/*     Set Eigenvalues IHI+1:N */

    i__1 = *n;
    for (j = *ihi + 1; j <= i__1; ++j) {
	i__2 = *lda1 * *lda2;
	zlapr1(&base, k, &s[1], &a[j + (j + a_dim2) * a_dim1], &i__2, &alpha[
		j], &beta[j], &scal[j]);
/* L30: */
    }

/*     If IHI < ILO, skip QZ steps */

    if (*ihi < *ilo) {
	goto L470;
    }

/*     MAIN PERIODIC QZ ITERATION LOOP */

/*     Initialize dynamic indices */

/*     Eigenvalues ILAST+1:N have been found. */
/*        Column operations modify rows IFRSTM:whatever. */
/*        Row operations modify columns whatever:ILASTM. */

/*     If only eigenvalues are being computed, then */
/*        IFRSTM is the row of the last splitting row above row ILAST; */
/*        this is always at least ILO. */
/*     IITER counts iterations since the last eigenvalue was found, */
/*        to tell when to use an observed zero or random shift. */
/*     MAXIT is the maximum number of QZ sweeps allowed. */

    ilast = *ihi;
    if (lschr) {
	ifrstm = 1;
	ilastm = *n;
    } else {
	ifrstm = *ilo;
	ilastm = *ihi;
    }
    iiter = 0;
    iseed[0] = 1;
    iseed[1] = 0;
    iseed[2] = 0;
    iseed[3] = 0;
    maxit = in * 30;

    i__1 = maxit;
    for (jiter = 1; jiter <= i__1; ++jiter) {

/*        Special Case: ILAST = ILO */

	if (ilast == *ilo) {
	    goto L400;
	}

/*        ************************************************************** */
/*        *                     CHECK FOR DEFLATION                    * */
/*        ************************************************************** */

/*        Test 1:  Deflation in the Hessenberg matrix. */

	jlo = *ilo;
	i__2 = *ilo + 1;
	for (j = ilast; j >= i__2; --j) {
	    tol = z_abs(&a[j - 1 + (j - 1 + a_dim2) * a_dim1]) + z_abs(&a[j + 
		    (j + a_dim2) * a_dim1]);
	    if (tol == 0.) {
		i__3 = j - *ilo + 1;
		tol = zlanhs("1", &i__3, &a[*ilo + (*ilo + a_dim2) * a_dim1],
			 lda1, &dwork[1]);
	    }
/* Computing MAX */
	    d__1 = ulp * tol;
	    tol = max(d__1,smlnum);
	    if (z_abs(&a[j + (j - 1 + a_dim2) * a_dim1]) <= tol) {
		i__3 = j + (j - 1 + a_dim2) * a_dim1;
		a[i__3].r = 0., a[i__3].i = 0.;
		jlo = j;
		if (j == ilast) {
		    goto L400;
		}
		goto L50;
	    }
/* L40: */
	}

L50:

/*        Test 2:  Deflation in the triangular matrices with index 1. */

	i__2 = *k;
	for (ldef = 2; ldef <= i__2; ++ldef) {
	    if (s[ldef] == 1) {
		i__3 = jlo;
		for (j = ilast; j >= i__3; --j) {
		    if (j == ilast) {
			tol = z_abs(&a[j - 1 + (j + ldef * a_dim2) * a_dim1]);
		    } else if (j == jlo) {
			tol = z_abs(&a[j + (j + 1 + ldef * a_dim2) * a_dim1]);
		    } else {
			tol = z_abs(&a[j - 1 + (j + ldef * a_dim2) * a_dim1]) 
				+ z_abs(&a[j + (j + 1 + ldef * a_dim2) * 
				a_dim1]);
		    }
		    if (tol == 0.) {
			i__4 = j - jlo + 1;
			tol = zlanhs("1", &i__4, &a[jlo + (jlo + ldef * 
				a_dim2) * a_dim1], lda1, &dwork[1])
				;
		    }
/* Computing MAX */
		    d__1 = ulp * tol;
		    tol = max(d__1,smlnum);
		    if (z_abs(&a[j + (j + ldef * a_dim2) * a_dim1]) <= tol) {
			i__4 = j + (j + ldef * a_dim2) * a_dim1;
			a[i__4].r = 0., a[i__4].i = 0.;
			goto L180;
		    }
/* L60: */
		}
	    }
/* L70: */
	}

/*        Test 3:  Deflation in the triangular matrices with index -1. */

	i__2 = *k;
	for (ldef = 2; ldef <= i__2; ++ldef) {
	    if (s[ldef] == -1) {
		i__3 = jlo;
		for (j = ilast; j >= i__3; --j) {
		    if (j == ilast) {
			tol = z_abs(&a[j - 1 + (j + ldef * a_dim2) * a_dim1]);
		    } else if (j == jlo) {
			tol = z_abs(&a[j + (j + 1 + ldef * a_dim2) * a_dim1]);
		    } else {
			tol = z_abs(&a[j - 1 + (j + ldef * a_dim2) * a_dim1]) 
				+ z_abs(&a[j + (j + 1 + ldef * a_dim2) * 
				a_dim1]);
		    }
		    if (tol == 0.) {
			i__4 = j - jlo + 1;
			tol = zlanhs("1", &i__4, &a[jlo + (jlo + ldef * 
				a_dim2) * a_dim1], lda1, &dwork[1])
				;
		    }
/* Computing MAX */
		    d__1 = ulp * tol;
		    tol = max(d__1,smlnum);
		    if (z_abs(&a[j + (j + ldef * a_dim2) * a_dim1]) <= tol) {
			i__4 = j + (j + ldef * a_dim2) * a_dim1;
			a[i__4].r = 0., a[i__4].i = 0.;
			goto L330;
		    }
/* L80: */
		}
	    }
/* L90: */
	}

/*        Test 4:  Controlled zero shift. */

	if (ziter >= 7 || ziter < 0) {

/*           Make Hessenberg matrix upper triangular. */

	    i__2 = ilast - 1;
	    for (j = jlo; j <= i__2; ++j) {
		i__3 = j + (j + a_dim2) * a_dim1;
		temp.r = a[i__3].r, temp.i = a[i__3].i;
		zlartg(&temp, &a[j + 1 + (j + a_dim2) * a_dim1], &cs, &sn, &
			a[j + (j + a_dim2) * a_dim1]);
		i__3 = j + 1 + (j + a_dim2) * a_dim1;
		a[i__3].r = 0., a[i__3].i = 0.;
		i__3 = ilastm - j;
		zrot(&i__3, &a[j + (j + 1 + a_dim2) * a_dim1], lda1, &a[j + 
			1 + (j + 1 + a_dim2) * a_dim1], lda1, &cs, &sn);
		dwork[j] = cs;
		i__3 = j;
		zwork[i__3].r = sn.r, zwork[i__3].i = sn.i;
/* L100: */
	    }
	    if (wantq) {
		i__2 = ilast - 1;
		for (j = jlo; j <= i__2; ++j) {
		    d_cnjg(&z__1, &zwork[j]);
		    zrot(n, &q[(j + q_dim2) * q_dim1 + 1], &c__1, &q[(j + 1 
			    + q_dim2) * q_dim1 + 1], &c__1, &dwork[j], &z__1);
/* L110: */
		}
	    }

/*           Propagate Transformations back to A_1. */

	    for (l = *k; l >= 2; --l) {
		if (s[l] == 1) {
		    i__2 = ilast - 1;
		    for (j = jlo; j <= i__2; ++j) {
			cs = dwork[j];
			i__3 = j;
			sn.r = zwork[i__3].r, sn.i = zwork[i__3].i;
			if (sn.r != 0. || sn.i != 0.) {
			    i__3 = j + 2 - ifrstm;
			    d_cnjg(&z__1, &sn);
			    zrot(&i__3, &a[ifrstm + (j + l * a_dim2) * 
				    a_dim1], &c__1, &a[ifrstm + (j + 1 + l * 
				    a_dim2) * a_dim1], &c__1, &cs, &z__1);

/*                       Check for deflation */

			    tol = z_abs(&a[j + (j + l * a_dim2) * a_dim1]) + 
				    z_abs(&a[j + 1 + (j + 1 + l * a_dim2) * 
				    a_dim1]);
			    if (tol == 0.) {
				i__3 = j - jlo + 2;
				tol = zlanhs("1", &i__3, &a[jlo + (jlo + l * 
					a_dim2) * a_dim1], lda1, &dwork[1], (
					ftnlen)1);
			    }
/* Computing MAX */
			    d__1 = ulp * tol;
			    tol = max(d__1,smlnum);
			    if (z_abs(&a[j + 1 + (j + l * a_dim2) * a_dim1]) 
				    <= tol) {
				cs = 1.;
				sn.r = 0., sn.i = 0.;
				i__3 = j + 1 + (j + l * a_dim2) * a_dim1;
				a[i__3].r = 0., a[i__3].i = 0.;
			    }

			    i__3 = j + (j + l * a_dim2) * a_dim1;
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
			    zlartg(&temp, &a[j + 1 + (j + l * a_dim2) * 
				    a_dim1], &cs, &sn, &a[j + (j + l * a_dim2)
				     * a_dim1]);
			    i__3 = j + 1 + (j + l * a_dim2) * a_dim1;
			    a[i__3].r = 0., a[i__3].i = 0.;
			    i__3 = ilastm - j;
			    zrot(&i__3, &a[j + (j + 1 + l * a_dim2) * a_dim1]
				    , lda1, &a[j + 1 + (j + 1 + l * a_dim2) * 
				    a_dim1], lda1, &cs, &sn);
			}
			dwork[j] = cs;
			i__3 = j;
			zwork[i__3].r = sn.r, zwork[i__3].i = sn.i;
/* L120: */
		    }
		} else {
		    i__2 = ilast - 1;
		    for (j = jlo; j <= i__2; ++j) {
			cs = dwork[j];
			i__3 = j;
			sn.r = zwork[i__3].r, sn.i = zwork[i__3].i;
			if (sn.r != 0. || sn.i != 0.) {
			    i__3 = ilastm - j + 1;
			    zrot(&i__3, &a[j + (j + l * a_dim2) * a_dim1], 
				    lda1, &a[j + 1 + (j + l * a_dim2) * 
				    a_dim1], lda1, &cs, &sn);

/*                       Check for deflation */

			    tol = z_abs(&a[j + (j + l * a_dim2) * a_dim1]) + 
				    z_abs(&a[j + 1 + (j + 1 + l * a_dim2) * 
				    a_dim1]);
			    if (tol == 0.) {
				i__3 = j - jlo + 2;
				tol = zlanhs("1", &i__3, &a[jlo + (jlo + l * 
					a_dim2) * a_dim1], lda1, &dwork[1], (
					ftnlen)1);
			    }
/* Computing MAX */
			    d__1 = ulp * tol;
			    tol = max(d__1,smlnum);
			    if (z_abs(&a[j + 1 + (j + l * a_dim2) * a_dim1]) 
				    <= tol) {
				cs = 1.;
				sn.r = 0., sn.i = 0.;
				i__3 = j + 1 + (j + l * a_dim2) * a_dim1;
				a[i__3].r = 0., a[i__3].i = 0.;
			    }

			    i__3 = j + 1 + (j + 1 + l * a_dim2) * a_dim1;
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
			    zlartg(&temp, &a[j + 1 + (j + l * a_dim2) * 
				    a_dim1], &cs, &sn, &a[j + 1 + (j + 1 + l *
				     a_dim2) * a_dim1]);
			    i__3 = j + 1 + (j + l * a_dim2) * a_dim1;
			    a[i__3].r = 0., a[i__3].i = 0.;
			    i__3 = j + 1 - ifrstm;
			    zrot(&i__3, &a[ifrstm + (j + 1 + l * a_dim2) * 
				    a_dim1], &c__1, &a[ifrstm + (j + l * 
				    a_dim2) * a_dim1], &c__1, &cs, &sn);
			    dwork[j] = cs;
			    i__3 = j;
			    z__1.r = -sn.r, z__1.i = -sn.i;
			    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
			}
/* L130: */
		    }
		}

		if (wantq) {
		    i__2 = ilast - 1;
		    for (j = jlo; j <= i__2; ++j) {
			d_cnjg(&z__1, &zwork[j]);
			zrot(n, &q[(j + l * q_dim2) * q_dim1 + 1], &c__1, &q[
				(j + 1 + l * q_dim2) * q_dim1 + 1], &c__1, &
				dwork[j], &z__1);
/* L140: */
		    }
		}
/* L150: */
	    }

/*           Apply the transformations to the right hand side of the */
/*           Hessenberg factor. */

	    ziter = 0;
	    i__2 = ilast - 1;
	    for (j = jlo; j <= i__2; ++j) {
		cs = dwork[j];
		i__3 = j;
		sn.r = zwork[i__3].r, sn.i = zwork[i__3].i;
		i__3 = j + 2 - ifrstm;
		d_cnjg(&z__1, &sn);
		zrot(&i__3, &a[ifrstm + (j + a_dim2) * a_dim1], &c__1, &a[
			ifrstm + (j + 1 + a_dim2) * a_dim1], &c__1, &cs, &
			z__1);
		if (sn.r == 0. && sn.i == 0.) {
		    ziter = 1;
		}
/* L160: */
	    }

/*           No QZ iteration. */

	    goto L450;
	}

/*        ************************************************************** */
/*        *                     HANDLE DEFLATIONS                      * */
/*        ************************************************************** */

/*        Case I: Deflation occurs in the Hessenberg matrix. The QZ */
/*                iteration is only applied to the JLO:ILAST part. */

/* L170: */
	ifirst = jlo;

/*        Go to the periodic QZ steps */

	goto L410;

/*        Case II: Deflation occurs in a triangular matrix with index 1. */

/*        Do an unshifted periodic QZ step. */

L180:
	jdef = j;
	i__2 = jdef - 1;
	for (j = jlo; j <= i__2; ++j) {
	    i__3 = j + (j + a_dim2) * a_dim1;
	    temp.r = a[i__3].r, temp.i = a[i__3].i;
	    zlartg(&temp, &a[j + 1 + (j + a_dim2) * a_dim1], &cs, &sn, &a[j 
		    + (j + a_dim2) * a_dim1]);
	    i__3 = j + 1 + (j + a_dim2) * a_dim1;
	    a[i__3].r = 0., a[i__3].i = 0.;
	    i__3 = ilastm - j;
	    zrot(&i__3, &a[j + (j + 1 + a_dim2) * a_dim1], lda1, &a[j + 1 + (
		    j + 1 + a_dim2) * a_dim1], lda1, &cs, &sn);
	    dwork[j] = cs;
	    i__3 = j;
	    zwork[i__3].r = sn.r, zwork[i__3].i = sn.i;
/* L190: */
	}
	if (wantq) {
	    i__2 = jdef - 1;
	    for (j = jlo; j <= i__2; ++j) {
		d_cnjg(&z__1, &zwork[j]);
		zrot(n, &q[(j + q_dim2) * q_dim1 + 1], &c__1, &q[(j + 1 + 
			q_dim2) * q_dim1 + 1], &c__1, &dwork[j], &z__1);
/* L200: */
	    }
	}

/*        Propagate the transformations through the triangular matrices. */
/*        Due to the zero element on the diagonal of the LDEF-th */
/*        factor the number of transformations drops by one. */

	for (l = *k; l >= 2; --l) {
	    if (l < ldef) {
		ntra = jdef - 2;
	    } else {
		ntra = jdef - 1;
	    }
	    if (s[l] == 1) {
		i__2 = ntra;
		for (j = jlo; j <= i__2; ++j) {
		    i__3 = j + 2 - ifrstm;
		    d_cnjg(&z__1, &zwork[j]);
		    zrot(&i__3, &a[ifrstm + (j + l * a_dim2) * a_dim1], &
			    c__1, &a[ifrstm + (j + 1 + l * a_dim2) * a_dim1], 
			    &c__1, &dwork[j], &z__1);
		    i__3 = j + (j + l * a_dim2) * a_dim1;
		    temp.r = a[i__3].r, temp.i = a[i__3].i;
		    zlartg(&temp, &a[j + 1 + (j + l * a_dim2) * a_dim1], &cs,
			     &sn, &a[j + (j + l * a_dim2) * a_dim1]);
		    i__3 = j + 1 + (j + l * a_dim2) * a_dim1;
		    a[i__3].r = 0., a[i__3].i = 0.;
		    i__3 = ilastm - j;
		    zrot(&i__3, &a[j + (j + 1 + l * a_dim2) * a_dim1], lda1, 
			    &a[j + 1 + (j + 1 + l * a_dim2) * a_dim1], lda1, &
			    cs, &sn);
		    dwork[j] = cs;
		    i__3 = j;
		    zwork[i__3].r = sn.r, zwork[i__3].i = sn.i;
/* L210: */
		}
	    } else {
		i__2 = ntra;
		for (j = jlo; j <= i__2; ++j) {
		    i__3 = ilastm - j + 1;
		    zrot(&i__3, &a[j + (j + l * a_dim2) * a_dim1], lda1, &a[
			    j + 1 + (j + l * a_dim2) * a_dim1], lda1, &dwork[
			    j], &zwork[j]);
		    i__3 = j + 1 + (j + 1 + l * a_dim2) * a_dim1;
		    temp.r = a[i__3].r, temp.i = a[i__3].i;
		    zlartg(&temp, &a[j + 1 + (j + l * a_dim2) * a_dim1], &cs,
			     &sn, &a[j + 1 + (j + 1 + l * a_dim2) * a_dim1]);
		    i__3 = j + 1 + (j + l * a_dim2) * a_dim1;
		    a[i__3].r = 0., a[i__3].i = 0.;
		    i__3 = j + 1 - ifrstm;
		    zrot(&i__3, &a[ifrstm + (j + 1 + l * a_dim2) * a_dim1], &
			    c__1, &a[ifrstm + (j + l * a_dim2) * a_dim1], &
			    c__1, &cs, &sn);
		    dwork[j] = cs;
		    i__3 = j;
		    z__1.r = -sn.r, z__1.i = -sn.i;
		    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L220: */
		}
	    }
	    if (wantq) {
		i__2 = ntra;
		for (j = jlo; j <= i__2; ++j) {
		    d_cnjg(&z__1, &zwork[j]);
		    zrot(n, &q[(j + l * q_dim2) * q_dim1 + 1], &c__1, &q[(j 
			    + 1 + l * q_dim2) * q_dim1 + 1], &c__1, &dwork[j],
			     &z__1);
/* L230: */
		}
	    }
/* L240: */
	}

/*        Apply the transformations to the right hand side of the */
/*        Hessenberg factor. */

	i__2 = jdef - 2;
	for (j = jlo; j <= i__2; ++j) {
	    i__3 = j + 2 - ifrstm;
	    d_cnjg(&z__1, &zwork[j]);
	    zrot(&i__3, &a[ifrstm + (j + a_dim2) * a_dim1], &c__1, &a[ifrstm 
		    + (j + 1 + a_dim2) * a_dim1], &c__1, &dwork[j], &z__1);
/* L250: */
	}

/*        Do an unshifted periodic ZQ step. */

	i__2 = jdef + 1;
	for (j = ilast; j >= i__2; --j) {
	    i__3 = j + (j + a_dim2) * a_dim1;
	    temp.r = a[i__3].r, temp.i = a[i__3].i;
	    zlartg(&temp, &a[j + (j - 1 + a_dim2) * a_dim1], &cs, &sn, &a[j 
		    + (j + a_dim2) * a_dim1]);
	    i__3 = j + (j - 1 + a_dim2) * a_dim1;
	    a[i__3].r = 0., a[i__3].i = 0.;
	    i__3 = j - ifrstm;
	    zrot(&i__3, &a[ifrstm + (j + a_dim2) * a_dim1], &c__1, &a[ifrstm 
		    + (j - 1 + a_dim2) * a_dim1], &c__1, &cs, &sn);
	    dwork[j] = cs;
	    i__3 = j;
	    z__1.r = -sn.r, z__1.i = -sn.i;
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L260: */
	}
	if (wantq) {
	    i__2 = jdef + 1;
	    for (j = ilast; j >= i__2; --j) {
		d_cnjg(&z__1, &zwork[j]);
		zrot(n, &q[(j - 1 + (q_dim2 << 1)) * q_dim1 + 1], &c__1, &q[(
			j + (q_dim2 << 1)) * q_dim1 + 1], &c__1, &dwork[j], &
			z__1);
/* L270: */
	    }
	}

/*        Propagate the transformations through the triangular matrices. */

	i__2 = *k;
	for (l = 2; l <= i__2; ++l) {
	    if (l > ldef) {
		ntra = jdef + 2;
	    } else {
		ntra = jdef + 1;
	    }
	    if (s[l] == -1) {
		i__3 = ntra;
		for (j = ilast; j >= i__3; --j) {
		    cs = dwork[j];
		    i__4 = j;
		    sn.r = zwork[i__4].r, sn.i = zwork[i__4].i;
		    i__4 = j + 1 - ifrstm;
		    d_cnjg(&z__1, &sn);
		    zrot(&i__4, &a[ifrstm + (j - 1 + l * a_dim2) * a_dim1], &
			    c__1, &a[ifrstm + (j + l * a_dim2) * a_dim1], &
			    c__1, &cs, &z__1);
		    i__4 = j - 1 + (j - 1 + l * a_dim2) * a_dim1;
		    temp.r = a[i__4].r, temp.i = a[i__4].i;
		    zlartg(&temp, &a[j + (j - 1 + l * a_dim2) * a_dim1], &cs,
			     &sn, &a[j - 1 + (j - 1 + l * a_dim2) * a_dim1]);
		    i__4 = j + (j - 1 + l * a_dim2) * a_dim1;
		    a[i__4].r = 0., a[i__4].i = 0.;
		    i__4 = ilastm - j + 1;
		    zrot(&i__4, &a[j - 1 + (j + l * a_dim2) * a_dim1], lda1, 
			    &a[j + (j + l * a_dim2) * a_dim1], lda1, &cs, &sn)
			    ;
		    dwork[j] = cs;
		    i__4 = j;
		    zwork[i__4].r = sn.r, zwork[i__4].i = sn.i;
/* L280: */
		}
	    } else {
		i__3 = ntra;
		for (j = ilast; j >= i__3; --j) {
		    i__4 = ilastm - j + 2;
		    zrot(&i__4, &a[j - 1 + (j - 1 + l * a_dim2) * a_dim1], 
			    lda1, &a[j + (j - 1 + l * a_dim2) * a_dim1], lda1,
			     &dwork[j], &zwork[j]);
		    i__4 = j + (j + l * a_dim2) * a_dim1;
		    temp.r = a[i__4].r, temp.i = a[i__4].i;
		    zlartg(&temp, &a[j + (j - 1 + l * a_dim2) * a_dim1], &cs,
			     &sn, &a[j + (j + l * a_dim2) * a_dim1]);
		    i__4 = j + (j - 1 + l * a_dim2) * a_dim1;
		    a[i__4].r = 0., a[i__4].i = 0.;
		    i__4 = j - ifrstm;
		    zrot(&i__4, &a[ifrstm + (j + l * a_dim2) * a_dim1], &
			    c__1, &a[ifrstm + (j - 1 + l * a_dim2) * a_dim1], 
			    &c__1, &cs, &sn);
		    dwork[j] = cs;
		    i__4 = j;
		    z__1.r = -sn.r, z__1.i = -sn.i;
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
/* L290: */
		}
	    }
	    if (wantq) {
		ln = l + 1;
		if (l == *k) {
		    ln = 1;
		}
		i__3 = ntra;
		for (j = ilast; j >= i__3; --j) {
		    d_cnjg(&z__1, &zwork[j]);
		    zrot(n, &q[(j - 1 + ln * q_dim2) * q_dim1 + 1], &c__1, &
			    q[(j + ln * q_dim2) * q_dim1 + 1], &c__1, &dwork[
			    j], &z__1);
/* L300: */
		}
	    }
/* L310: */
	}

/*        Apply the transformations to the left hand side of the */
/*        Hessenberg factor. */

	i__2 = jdef + 2;
	for (j = ilast; j >= i__2; --j) {
	    i__3 = ilastm - j + 2;
	    zrot(&i__3, &a[j - 1 + (j - 1 + a_dim2) * a_dim1], lda1, &a[j + (
		    j - 1 + a_dim2) * a_dim1], lda1, &dwork[j], &zwork[j]);
/* L320: */
	}

/*        No QZ iteration. */

	goto L450;

/*        Case III: Deflation occurs in a triangular matrix with */
/*                  index -1. */

L330:
	jdef = j;
	if (jdef > (ilast - jlo + 1) / 2) {

/*           Chase the zero downwards to the last position */

	    i__2 = ilast - 1;
	    for (j1 = jdef; j1 <= i__2; ++j1) {
		j = j1;
		i__3 = j + (j + 1 + ldef * a_dim2) * a_dim1;
		temp.r = a[i__3].r, temp.i = a[i__3].i;
		zlartg(&temp, &a[j + 1 + (j + 1 + ldef * a_dim2) * a_dim1], &
			cs, &sn, &a[j + (j + 1 + ldef * a_dim2) * a_dim1]);
		i__3 = j + 1 + (j + 1 + ldef * a_dim2) * a_dim1;
		a[i__3].r = 0., a[i__3].i = 0.;
		i__3 = ilastm - j - 1;
		zrot(&i__3, &a[j + (j + 2 + ldef * a_dim2) * a_dim1], lda1, &
			a[j + 1 + (j + 2 + ldef * a_dim2) * a_dim1], lda1, &
			cs, &sn);
		ln = ldef + 1;
		if (ldef == *k) {
		    ln = 1;
		}
		if (wantq) {
		    d_cnjg(&z__1, &sn);
		    zrot(n, &q[(j + ln * q_dim2) * q_dim1 + 1], &c__1, &q[(j 
			    + 1 + ln * q_dim2) * q_dim1 + 1], &c__1, &cs, &
			    z__1);
		}
		i__3 = *k - 1;
		for (l = 1; l <= i__3; ++l) {
		    if (ln == 1) {
			i__4 = ilastm - j + 2;
			zrot(&i__4, &a[j + (j - 1 + ln * a_dim2) * a_dim1], 
				lda1, &a[j + 1 + (j - 1 + ln * a_dim2) * 
				a_dim1], lda1, &cs, &sn);
			i__4 = j + 1 + (j + ln * a_dim2) * a_dim1;
			temp.r = a[i__4].r, temp.i = a[i__4].i;
			zlartg(&temp, &a[j + 1 + (j - 1 + ln * a_dim2) * 
				a_dim1], &cs, &sn, &a[j + 1 + (j + ln * 
				a_dim2) * a_dim1]);
			i__4 = j + 1 + (j - 1 + ln * a_dim2) * a_dim1;
			a[i__4].r = 0., a[i__4].i = 0.;
			i__4 = j - ifrstm + 1;
			zrot(&i__4, &a[ifrstm + (j + ln * a_dim2) * a_dim1], 
				&c__1, &a[ifrstm + (j - 1 + ln * a_dim2) * 
				a_dim1], &c__1, &cs, &sn);
			z__1.r = -sn.r, z__1.i = -sn.i;
			sn.r = z__1.r, sn.i = z__1.i;
			--j;
		    } else if (s[ln] == 1) {
			i__4 = ilastm - j + 1;
			zrot(&i__4, &a[j + (j + ln * a_dim2) * a_dim1], lda1,
				 &a[j + 1 + (j + ln * a_dim2) * a_dim1], lda1,
				 &cs, &sn);
			i__4 = j + 1 + (j + 1 + ln * a_dim2) * a_dim1;
			temp.r = a[i__4].r, temp.i = a[i__4].i;
			zlartg(&temp, &a[j + 1 + (j + ln * a_dim2) * a_dim1],
				 &cs, &sn, &a[j + 1 + (j + 1 + ln * a_dim2) * 
				a_dim1]);
			i__4 = j + 1 + (j + ln * a_dim2) * a_dim1;
			a[i__4].r = 0., a[i__4].i = 0.;
			i__4 = j - ifrstm + 1;
			zrot(&i__4, &a[ifrstm + (j + 1 + ln * a_dim2) * 
				a_dim1], &c__1, &a[ifrstm + (j + ln * a_dim2) 
				* a_dim1], &c__1, &cs, &sn);
			z__1.r = -sn.r, z__1.i = -sn.i;
			sn.r = z__1.r, sn.i = z__1.i;
		    } else {
			i__4 = j - ifrstm + 2;
			d_cnjg(&z__1, &sn);
			zrot(&i__4, &a[ifrstm + (j + ln * a_dim2) * a_dim1], 
				&c__1, &a[ifrstm + (j + 1 + ln * a_dim2) * 
				a_dim1], &c__1, &cs, &z__1);
			i__4 = j + (j + ln * a_dim2) * a_dim1;
			temp.r = a[i__4].r, temp.i = a[i__4].i;
			zlartg(&temp, &a[j + 1 + (j + ln * a_dim2) * a_dim1],
				 &cs, &sn, &a[j + (j + ln * a_dim2) * a_dim1])
				;
			i__4 = j + 1 + (j + ln * a_dim2) * a_dim1;
			a[i__4].r = 0., a[i__4].i = 0.;
			i__4 = ilastm - j;
			zrot(&i__4, &a[j + (j + 1 + ln * a_dim2) * a_dim1], 
				lda1, &a[j + 1 + (j + 1 + ln * a_dim2) * 
				a_dim1], lda1, &cs, &sn);
		    }
		    ++ln;
		    if (ln > *k) {
			ln = 1;
		    }
		    if (wantq) {
			d_cnjg(&z__1, &sn);
			zrot(n, &q[(j + ln * q_dim2) * q_dim1 + 1], &c__1, &
				q[(j + 1 + ln * q_dim2) * q_dim1 + 1], &c__1, 
				&cs, &z__1);
		    }
/* L340: */
		}
		i__3 = j - ifrstm + 1;
		d_cnjg(&z__1, &sn);
		zrot(&i__3, &a[ifrstm + (j + ldef * a_dim2) * a_dim1], &c__1,
			 &a[ifrstm + (j + 1 + ldef * a_dim2) * a_dim1], &c__1,
			 &cs, &z__1);
/* L350: */
	    }

/*           Deflate the last element in the Hessenberg matrix. */

	    j = ilast;
	    i__2 = j + (j + a_dim2) * a_dim1;
	    temp.r = a[i__2].r, temp.i = a[i__2].i;
	    zlartg(&temp, &a[j + (j - 1 + a_dim2) * a_dim1], &cs, &sn, &a[j 
		    + (j + a_dim2) * a_dim1]);
	    i__2 = j + (j - 1 + a_dim2) * a_dim1;
	    a[i__2].r = 0., a[i__2].i = 0.;
	    i__2 = j - ifrstm;
	    zrot(&i__2, &a[ifrstm + (j + a_dim2) * a_dim1], &c__1, &a[ifrstm 
		    + (j - 1 + a_dim2) * a_dim1], &c__1, &cs, &sn);
	    z__1.r = -sn.r, z__1.i = -sn.i;
	    sn.r = z__1.r, sn.i = z__1.i;
/*         IF ( WANTQ.NE.0 ) THEN */
	    if (wantq) {
		d_cnjg(&z__1, &sn);
		zrot(n, &q[(j - 1 + (q_dim2 << 1)) * q_dim1 + 1], &c__1, &q[(
			j + (q_dim2 << 1)) * q_dim1 + 1], &c__1, &cs, &z__1);
	    }
	    i__2 = ldef - 1;
	    for (l = 2; l <= i__2; ++l) {
		if (s[l] == -1) {
		    i__3 = j + 1 - ifrstm;
		    d_cnjg(&z__1, &sn);
		    zrot(&i__3, &a[ifrstm + (j - 1 + l * a_dim2) * a_dim1], &
			    c__1, &a[ifrstm + (j + l * a_dim2) * a_dim1], &
			    c__1, &cs, &z__1);
		    i__3 = j - 1 + (j - 1 + l * a_dim2) * a_dim1;
		    temp.r = a[i__3].r, temp.i = a[i__3].i;
		    zlartg(&temp, &a[j + (j - 1 + l * a_dim2) * a_dim1], &cs,
			     &sn, &a[j - 1 + (j - 1 + l * a_dim2) * a_dim1]);
		    i__3 = j + (j - 1 + l * a_dim2) * a_dim1;
		    a[i__3].r = 0., a[i__3].i = 0.;
		    i__3 = ilastm - j + 1;
		    zrot(&i__3, &a[j - 1 + (j + l * a_dim2) * a_dim1], lda1, 
			    &a[j + (j + l * a_dim2) * a_dim1], lda1, &cs, &sn)
			    ;
		} else {
		    i__3 = ilastm - j + 2;
		    zrot(&i__3, &a[j - 1 + (j - 1 + l * a_dim2) * a_dim1], 
			    lda1, &a[j + (j - 1 + l * a_dim2) * a_dim1], lda1,
			     &cs, &sn);
		    i__3 = j + (j + l * a_dim2) * a_dim1;
		    temp.r = a[i__3].r, temp.i = a[i__3].i;
		    zlartg(&temp, &a[j + (j - 1 + l * a_dim2) * a_dim1], &cs,
			     &sn, &a[j + (j + l * a_dim2) * a_dim1]);
		    i__3 = j + (j - 1 + l * a_dim2) * a_dim1;
		    a[i__3].r = 0., a[i__3].i = 0.;
		    i__3 = j - ifrstm;
		    zrot(&i__3, &a[ifrstm + (j + l * a_dim2) * a_dim1], &
			    c__1, &a[ifrstm + (j - 1 + l * a_dim2) * a_dim1], 
			    &c__1, &cs, &sn);
		    z__1.r = -sn.r, z__1.i = -sn.i;
		    sn.r = z__1.r, sn.i = z__1.i;
		}
		if (wantq) {
		    ln = l + 1;
		    if (l == *k) {
			ln = 1;
		    }
		    d_cnjg(&z__1, &sn);
		    zrot(n, &q[(j - 1 + ln * q_dim2) * q_dim1 + 1], &c__1, &
			    q[(j + ln * q_dim2) * q_dim1 + 1], &c__1, &cs, &
			    z__1);
		}
/* L360: */
	    }
	    i__2 = j + 1 - ifrstm;
	    d_cnjg(&z__1, &sn);
	    zrot(&i__2, &a[ifrstm + (j - 1 + ldef * a_dim2) * a_dim1], &c__1,
		     &a[ifrstm + (j + ldef * a_dim2) * a_dim1], &c__1, &cs, &
		    z__1);
	} else {

/*           Chase the zero upwards to the first position. */

	    i__2 = jlo + 1;
	    for (j1 = jdef; j1 >= i__2; --j1) {
		j = j1;
		i__3 = j - 1 + (j + ldef * a_dim2) * a_dim1;
		temp.r = a[i__3].r, temp.i = a[i__3].i;
		zlartg(&temp, &a[j - 1 + (j - 1 + ldef * a_dim2) * a_dim1], &
			cs, &sn, &a[j - 1 + (j + ldef * a_dim2) * a_dim1]);
		i__3 = j - 1 + (j - 1 + ldef * a_dim2) * a_dim1;
		a[i__3].r = 0., a[i__3].i = 0.;
		i__3 = j - ifrstm - 1;
		zrot(&i__3, &a[ifrstm + (j + ldef * a_dim2) * a_dim1], &c__1,
			 &a[ifrstm + (j - 1 + ldef * a_dim2) * a_dim1], &c__1,
			 &cs, &sn);
		z__1.r = -sn.r, z__1.i = -sn.i;
		sn.r = z__1.r, sn.i = z__1.i;
		if (wantq) {
		    d_cnjg(&z__1, &sn);
		    zrot(n, &q[(j - 1 + ldef * q_dim2) * q_dim1 + 1], &c__1, 
			    &q[(j + ldef * q_dim2) * q_dim1 + 1], &c__1, &cs, 
			    &z__1);
		}
		ln = ldef - 1;
		i__3 = *k - 1;
		for (l = 1; l <= i__3; ++l) {
		    if (ln == 1) {
			i__4 = j - ifrstm + 2;
			d_cnjg(&z__1, &sn);
			zrot(&i__4, &a[ifrstm + (j - 1 + ln * a_dim2) * 
				a_dim1], &c__1, &a[ifrstm + (j + ln * a_dim2) 
				* a_dim1], &c__1, &cs, &z__1);
			i__4 = j + (j - 1 + ln * a_dim2) * a_dim1;
			temp.r = a[i__4].r, temp.i = a[i__4].i;
			zlartg(&temp, &a[j + 1 + (j - 1 + ln * a_dim2) * 
				a_dim1], &cs, &sn, &a[j + (j - 1 + ln * 
				a_dim2) * a_dim1]);
			i__4 = j + 1 + (j - 1 + ln * a_dim2) * a_dim1;
			a[i__4].r = 0., a[i__4].i = 0.;
			i__4 = ilastm - j + 1;
			zrot(&i__4, &a[j + (j + ln * a_dim2) * a_dim1], lda1,
				 &a[j + 1 + (j + ln * a_dim2) * a_dim1], lda1,
				 &cs, &sn);
			++j;
		    } else if (s[ln] == -1) {
			i__4 = ilastm - j + 2;
			zrot(&i__4, &a[j - 1 + (j - 1 + ln * a_dim2) * 
				a_dim1], lda1, &a[j + (j - 1 + ln * a_dim2) * 
				a_dim1], lda1, &cs, &sn);
			i__4 = j + (j + ln * a_dim2) * a_dim1;
			temp.r = a[i__4].r, temp.i = a[i__4].i;
			zlartg(&temp, &a[j + (j - 1 + ln * a_dim2) * a_dim1],
				 &cs, &sn, &a[j + (j + ln * a_dim2) * a_dim1])
				;
			i__4 = j + (j - 1 + ln * a_dim2) * a_dim1;
			a[i__4].r = 0., a[i__4].i = 0.;
			i__4 = j - ifrstm;
			zrot(&i__4, &a[ifrstm + (j + ln * a_dim2) * a_dim1], 
				&c__1, &a[ifrstm + (j - 1 + ln * a_dim2) * 
				a_dim1], &c__1, &cs, &sn);
			z__1.r = -sn.r, z__1.i = -sn.i;
			sn.r = z__1.r, sn.i = z__1.i;
		    } else {
			i__4 = j - ifrstm + 1;
			d_cnjg(&z__1, &sn);
			zrot(&i__4, &a[ifrstm + (j - 1 + ln * a_dim2) * 
				a_dim1], &c__1, &a[ifrstm + (j + ln * a_dim2) 
				* a_dim1], &c__1, &cs, &z__1);
			i__4 = j - 1 + (j - 1 + ln * a_dim2) * a_dim1;
			temp.r = a[i__4].r, temp.i = a[i__4].i;
			zlartg(&temp, &a[j + (j - 1 + ln * a_dim2) * a_dim1],
				 &cs, &sn, &a[j - 1 + (j - 1 + ln * a_dim2) * 
				a_dim1]);
			i__4 = j + (j - 1 + ln * a_dim2) * a_dim1;
			a[i__4].r = 0., a[i__4].i = 0.;
			i__4 = ilastm - j + 1;
			zrot(&i__4, &a[j - 1 + (j + ln * a_dim2) * a_dim1], 
				lda1, &a[j + (j + ln * a_dim2) * a_dim1], 
				lda1, &cs, &sn);
		    }
		    if (wantq) {
			d_cnjg(&z__1, &sn);
			zrot(n, &q[(j - 1 + ln * q_dim2) * q_dim1 + 1], &
				c__1, &q[(j + ln * q_dim2) * q_dim1 + 1], &
				c__1, &cs, &z__1);
		    }
		    --ln;
		    if (ln <= 0) {
			ln = *k;
		    }
/* L370: */
		}
		i__3 = ilastm - j + 1;
		zrot(&i__3, &a[j - 1 + (j + ldef * a_dim2) * a_dim1], lda1, &
			a[j + (j + ldef * a_dim2) * a_dim1], lda1, &cs, &sn);
/* L380: */
	    }

/*           Deflate the first element in the Hessenberg matrix. */

	    j = jlo;
	    i__2 = j + (j + a_dim2) * a_dim1;
	    temp.r = a[i__2].r, temp.i = a[i__2].i;
	    zlartg(&temp, &a[j + 1 + (j + a_dim2) * a_dim1], &cs, &sn, &a[j 
		    + (j + a_dim2) * a_dim1]);
	    i__2 = j + 1 + (j + a_dim2) * a_dim1;
	    a[i__2].r = 0., a[i__2].i = 0.;
	    i__2 = ilastm - j;
	    zrot(&i__2, &a[j + (j + 1 + a_dim2) * a_dim1], lda1, &a[j + 1 + (
		    j + 1 + a_dim2) * a_dim1], lda1, &cs, &sn);
	    if (wantq) {
		d_cnjg(&z__1, &sn);
		zrot(n, &q[(j + q_dim2) * q_dim1 + 1], &c__1, &q[(j + 1 + 
			q_dim2) * q_dim1 + 1], &c__1, &cs, &z__1);
	    }
	    i__2 = ldef + 1;
	    for (l = *k; l >= i__2; --l) {
		if (s[l] == 1) {
		    i__3 = j + 2 - ifrstm;
		    d_cnjg(&z__1, &sn);
		    zrot(&i__3, &a[ifrstm + (j + l * a_dim2) * a_dim1], &
			    c__1, &a[ifrstm + (j + 1 + l * a_dim2) * a_dim1], 
			    &c__1, &cs, &z__1);
		    i__3 = j + (j + l * a_dim2) * a_dim1;
		    temp.r = a[i__3].r, temp.i = a[i__3].i;
		    zlartg(&temp, &a[j + 1 + (j + l * a_dim2) * a_dim1], &cs,
			     &sn, &a[j + (j + l * a_dim2) * a_dim1]);
		    i__3 = j + 1 + (j + l * a_dim2) * a_dim1;
		    a[i__3].r = 0., a[i__3].i = 0.;
		    i__3 = ilastm - j;
		    zrot(&i__3, &a[j + (j + 1 + l * a_dim2) * a_dim1], lda1, 
			    &a[j + 1 + (j + 1 + l * a_dim2) * a_dim1], lda1, &
			    cs, &sn);
		} else {
		    i__3 = ilastm - j + 1;
		    zrot(&i__3, &a[j + (j + l * a_dim2) * a_dim1], lda1, &a[
			    j + 1 + (j + l * a_dim2) * a_dim1], lda1, &cs, &
			    sn);
		    i__3 = j + 1 + (j + 1 + l * a_dim2) * a_dim1;
		    temp.r = a[i__3].r, temp.i = a[i__3].i;
		    zlartg(&temp, &a[j + 1 + (j + l * a_dim2) * a_dim1], &cs,
			     &sn, &a[j + 1 + (j + 1 + l * a_dim2) * a_dim1]);
		    i__3 = j + 1 + (j + l * a_dim2) * a_dim1;
		    a[i__3].r = 0., a[i__3].i = 0.;
		    i__3 = j + 1 - ifrstm;
		    zrot(&i__3, &a[ifrstm + (j + 1 + l * a_dim2) * a_dim1], &
			    c__1, &a[ifrstm + (j + l * a_dim2) * a_dim1], &
			    c__1, &cs, &sn);
		    z__1.r = -sn.r, z__1.i = -sn.i;
		    sn.r = z__1.r, sn.i = z__1.i;
		}
		if (wantq) {
		    d_cnjg(&z__1, &sn);
		    zrot(n, &q[(j + l * q_dim2) * q_dim1 + 1], &c__1, &q[(j 
			    + 1 + l * q_dim2) * q_dim1 + 1], &c__1, &cs, &
			    z__1);
		}
/* L390: */
	    }
	    i__2 = ilastm - j;
	    zrot(&i__2, &a[j + (j + 1 + ldef * a_dim2) * a_dim1], lda1, &a[j 
		    + 1 + (j + 1 + ldef * a_dim2) * a_dim1], lda1, &cs, &sn);
	}

/*        No QZ iteration. */

	goto L450;

/*        Special case: A 1x1 block splits off at the bottom */

L400:
	i__2 = *lda1 * *lda2;
	zlapr1(&base, k, &s[1], &a[ilast + (ilast + a_dim2) * a_dim1], &i__2,
		 &alpha[ilast], &beta[ilast], &scal[ilast]);

/*        Go to next block - exit if finished. */

	--ilast;
	if (ilast < *ilo) {
	    goto L470;
	}

/*        Reset iteration counters. */

	iiter = 0;
	if (ziter != -1) {
	    ziter = 0;
	}
	if (! lschr) {
	    ilastm = ilast;
	    if (ifrstm > ilast) {
		ifrstm = *ilo;
	    }
	}

/*        No QZ iteration. */

	goto L450;

/*        ************************************************************** */
/*        *                      PERIODIC QZ STEP                      * */
/*        ************************************************************** */

/*        It is assumed that IFIRST < ILAST. */

L410:

	++iiter;
	++ziter;
	if (! lschr) {
	    ifrstm = ifirst;
	}

/*        Complex single shift. */

	if (iiter / 10 * 10 == iiter) {

/*           Exceptional shift. */

	    zlarnv(&c__2, iseed, &c__2, rnd);
	    zlartg(rnd, &rnd[1], &cs, &sn, &temp);
	} else {
	    zlartg(&c_b1, &c_b1, &cs, &sn, &temp);
	    for (l = *k; l >= 2; --l) {
		if (s[l] == 1) {
		    i__2 = ifirst + (ifirst + l * a_dim2) * a_dim1;
		    z__1.r = cs * a[i__2].r, z__1.i = cs * a[i__2].i;
		    i__3 = ilast + (ilast + l * a_dim2) * a_dim1;
		    d_cnjg(&z__3, &sn);
		    z__2.r = a[i__3].r * z__3.r - a[i__3].i * z__3.i, z__2.i =
			     a[i__3].r * z__3.i + a[i__3].i * z__3.r;
		    zlartg(&z__1, &z__2, &cs, &sn, &temp);
		} else {
		    i__2 = ilast + (ilast + l * a_dim2) * a_dim1;
		    z__1.r = cs * a[i__2].r, z__1.i = cs * a[i__2].i;
		    i__3 = ifirst + (ifirst + l * a_dim2) * a_dim1;
		    z__3.r = -a[i__3].r, z__3.i = -a[i__3].i;
		    d_cnjg(&z__4, &sn);
		    z__2.r = z__3.r * z__4.r - z__3.i * z__4.i, z__2.i = 
			    z__3.r * z__4.i + z__3.i * z__4.r;
		    zlartg(&z__1, &z__2, &cs, &sn, &temp);
		    z__1.r = -sn.r, z__1.i = -sn.i;
		    sn.r = z__1.r, sn.i = z__1.i;
		}
/* L420: */
	    }
	    i__2 = ifirst + (ifirst + a_dim2) * a_dim1;
	    z__2.r = cs * a[i__2].r, z__2.i = cs * a[i__2].i;
	    d_cnjg(&z__4, &sn);
	    i__3 = ilast + (ilast + a_dim2) * a_dim1;
	    z__3.r = z__4.r * a[i__3].r - z__4.i * a[i__3].i, z__3.i = z__4.r 
		    * a[i__3].i + z__4.i * a[i__3].r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    i__4 = ifirst + 1 + (ifirst + a_dim2) * a_dim1;
	    z__5.r = cs * a[i__4].r, z__5.i = cs * a[i__4].i;
	    zlartg(&z__1, &z__5, &cs, &sn, &temp);
	}

/*        Do the sweeps. */

	i__2 = ilast - 2;
	for (j1 = ifirst - 1; j1 <= i__2; ++j1) {
	    j = j1 + 1;

/*           Create bulge if J1 = IFIRST - 1, otherwise chase bulge. */

	    if (j1 < ifirst) {
		i__3 = ilastm - j + 1;
		zrot(&i__3, &a[j + (j + a_dim2) * a_dim1], lda1, &a[j + 1 + (
			j + a_dim2) * a_dim1], lda1, &cs, &sn);
	    } else {
		i__3 = j + (j - 1 + a_dim2) * a_dim1;
		temp.r = a[i__3].r, temp.i = a[i__3].i;
		zlartg(&temp, &a[j + 1 + (j - 1 + a_dim2) * a_dim1], &cs, &
			sn, &a[j + (j - 1 + a_dim2) * a_dim1]);
		i__3 = j + 1 + (j - 1 + a_dim2) * a_dim1;
		a[i__3].r = 0., a[i__3].i = 0.;
		i__3 = ilastm - j + 1;
		zrot(&i__3, &a[j + (j + a_dim2) * a_dim1], lda1, &a[j + 1 + (
			j + a_dim2) * a_dim1], lda1, &cs, &sn);
	    }
	    if (wantq) {
		d_cnjg(&z__1, &sn);
		zrot(n, &q[(j + q_dim2) * q_dim1 + 1], &c__1, &q[(j + 1 + 
			q_dim2) * q_dim1 + 1], &c__1, &cs, &z__1);
	    }

/*           Propagate rotation through AK, ..., A2 to A1. */

	    for (l = *k; l >= 2; --l) {
		if (s[l] == 1) {
		    i__3 = j + 2 - ifrstm;
		    d_cnjg(&z__1, &sn);
		    zrot(&i__3, &a[ifrstm + (j + l * a_dim2) * a_dim1], &
			    c__1, &a[ifrstm + (j + 1 + l * a_dim2) * a_dim1], 
			    &c__1, &cs, &z__1);
		    i__3 = j + (j + l * a_dim2) * a_dim1;
		    temp.r = a[i__3].r, temp.i = a[i__3].i;
		    zlartg(&temp, &a[j + 1 + (j + l * a_dim2) * a_dim1], &cs,
			     &sn, &a[j + (j + l * a_dim2) * a_dim1]);
		    i__3 = j + 1 + (j + l * a_dim2) * a_dim1;
		    a[i__3].r = 0., a[i__3].i = 0.;
		    i__3 = ilastm - j;
		    zrot(&i__3, &a[j + (j + 1 + l * a_dim2) * a_dim1], lda1, 
			    &a[j + 1 + (j + 1 + l * a_dim2) * a_dim1], lda1, &
			    cs, &sn);
		} else {
		    i__3 = ilastm - j + 1;
		    zrot(&i__3, &a[j + (j + l * a_dim2) * a_dim1], lda1, &a[
			    j + 1 + (j + l * a_dim2) * a_dim1], lda1, &cs, &
			    sn);
		    i__3 = j + 1 + (j + 1 + l * a_dim2) * a_dim1;
		    temp.r = a[i__3].r, temp.i = a[i__3].i;
		    zlartg(&temp, &a[j + 1 + (j + l * a_dim2) * a_dim1], &cs,
			     &sn, &a[j + 1 + (j + 1 + l * a_dim2) * a_dim1]);
		    i__3 = j + 1 + (j + l * a_dim2) * a_dim1;
		    a[i__3].r = 0., a[i__3].i = 0.;
		    i__3 = j + 1 - ifrstm;
		    zrot(&i__3, &a[ifrstm + (j + 1 + l * a_dim2) * a_dim1], &
			    c__1, &a[ifrstm + (j + l * a_dim2) * a_dim1], &
			    c__1, &cs, &sn);
		    z__1.r = -sn.r, z__1.i = -sn.i;
		    sn.r = z__1.r, sn.i = z__1.i;
		}
		if (wantq) {
		    d_cnjg(&z__1, &sn);
		    zrot(n, &q[(j + l * q_dim2) * q_dim1 + 1], &c__1, &q[(j 
			    + 1 + l * q_dim2) * q_dim1 + 1], &c__1, &cs, &
			    z__1);
		}
/* L430: */
	    }
/* Computing MIN */
	    i__4 = j + 2;
	    i__3 = min(i__4,ilastm) - ifrstm + 1;
	    d_cnjg(&z__1, &sn);
	    zrot(&i__3, &a[ifrstm + (j + a_dim2) * a_dim1], &c__1, &a[ifrstm 
		    + (j + 1 + a_dim2) * a_dim1], &c__1, &cs, &z__1);
/* L440: */
	}

/*        End of iteration loop. */

L450:
/* L460: */
	;
    }

/*     Drop through = non-convergence */

    *info = ilast;
    goto L550;

/*     Successful completion of all QZ steps */

L470:

/*     Set eigenvalues 1:ILO-1 */

    i__1 = *ilo - 1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *lda1 * *lda2;
	zlapr1(&base, k, &s[1], &a[j + (j + a_dim2) * a_dim1], &i__2, &alpha[
		j], &beta[j], &scal[j]);
/* L480: */
    }
    if (lschr) {

/*        Scale A(2,:,:) .. A(K,:,:). */

	for (l = *k; l >= 2; --l) {
	    if (s[l] == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    abst = z_abs(&a[j + (j + l * a_dim2) * a_dim1]);
		    if (abst > safmin) {
			i__2 = j + (j + l * a_dim2) * a_dim1;
			z__2.r = a[i__2].r / abst, z__2.i = a[i__2].i / abst;
			d_cnjg(&z__1, &z__2);
			temp.r = z__1.r, temp.i = z__1.i;
			i__2 = j + (j + l * a_dim2) * a_dim1;
			a[i__2].r = abst, a[i__2].i = 0.;
			if (j < *n) {
			    i__2 = *n - j;
			    zscal(&i__2, &temp, &a[j + (j + 1 + l * a_dim2) *
				     a_dim1], lda1);
			}
		    } else {
			temp.r = 1., temp.i = 0.;
		    }
		    i__2 = j;
		    zwork[i__2].r = temp.r, zwork[i__2].i = temp.i;
/* L490: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    abst = z_abs(&a[j + (j + l * a_dim2) * a_dim1]);
		    if (abst > safmin) {
			i__2 = j + (j + l * a_dim2) * a_dim1;
			z__2.r = a[i__2].r / abst, z__2.i = a[i__2].i / abst;
			d_cnjg(&z__1, &z__2);
			temp.r = z__1.r, temp.i = z__1.i;
			i__2 = j + (j + l * a_dim2) * a_dim1;
			a[i__2].r = abst, a[i__2].i = 0.;
			i__2 = j - 1;
			zscal(&i__2, &temp, &a[(j + l * a_dim2) * a_dim1 + 1]
				, &c__1);
		    } else {
			temp.r = 1., temp.i = 0.;
		    }
		    i__2 = j;
		    d_cnjg(&z__1, &temp);
		    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L500: */
		}
	    }
	    if (wantq) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    d_cnjg(&z__1, &zwork[j]);
		    zscal(n, &z__1, &q[(j + l * q_dim2) * q_dim1 + 1], &c__1)
			    ;
/* L510: */
		}
	    }
	    if (s[l - 1] == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    d_cnjg(&z__1, &zwork[j]);
		    zscal(&j, &z__1, &a[(j + (l - 1) * a_dim2) * a_dim1 + 1],
			     &c__1);
/* L520: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n - j + 1;
		    zscal(&i__2, &zwork[j], &a[j + (j + (l - 1) * a_dim2) * 
			    a_dim1], lda1);
/* L530: */
		}
	    }
/* L540: */
	}
    }
    *info = 0;

L550:

    dwork[1] = (doublereal) (*n);
    d__1 = (doublereal) (*n);
    z__1.r = d__1, z__1.i = 0.;
    zwork[1].r = z__1.r, zwork[1].i = z__1.i;
    return 0;
/* *** Last line of ZPGEQZ *** */
} /* zpgeqz_ */

