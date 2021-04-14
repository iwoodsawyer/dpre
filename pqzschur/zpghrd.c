/* zpghrd.f -- translated by f2c (version 20041007).
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

/* Subroutine */ int zpghrd(char *compq, integer *k, integer *n, integer *
	ilo, integer *ihi, integer *s, doublecomplex *a, integer *lda1, 
	integer *lda2, doublecomplex *q, integer *ldq1, integer *ldq2, 
	doublereal *dwork, integer *ldwork, doublecomplex *zwork, integer *
	lzwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_dim2, a_offset, q_dim1, q_dim2, q_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer l;
    static doublereal cs;
    static doublecomplex sn;
    static logical sok;
    static integer jcol, ierr;
    static doublecomplex temp;
    static integer jrow;
    extern /* Subroutine */ int zrot(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    extern logical lsame(char *, char *);
    static integer icols;
    static logical initq, wantq;
    static integer irows;
    extern /* Subroutine */ int xerbla(char *, integer *), zgeqrf(
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	    doublecomplex *, integer *, integer *), zgerqf(integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, integer *), zlacpy(char *, integer *,
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *
	    ), zlartg(doublecomplex *, doublecomplex *, doublereal *,
	     doublecomplex *, doublecomplex *), zlaset(char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *);
    static integer wrkopt;
    extern /* Subroutine */ int zungqr(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zunmqr(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	    ), zunmrq(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	    );


/*     PURPOSE */

/*     To reduce the general complex product */

/*                               S(2)                 S(K) */
/*          A(:,:,1)  *  A(:,:,2)     * ... * A(:,:,K) */

/*     to upper Hessenberg-triangular form, where A is N-by-N-by-K and S */
/*     is the signature array with values 1 or -1. The matrix A(:,:,1) */
/*     is reduced to upper Hessenberg form while the other matrices are */
/*     triangularized. */
/*     Simple and unblocked version. */

/*     If COMPQ = 'V' or COMPQ = 'I', then the unitary factors are */
/*     computed and stored in the array Q so that for S(I) = 1, */

/*                         H */
/*             Q(:,:,I)(in)   A(:,:,I)(in)   Q(:,:,MOD(I,K)+1)(in) */
/*                          H */
/*         =   Q(:,:,I)(out)  A(:,:,I)(out)  Q(:,:,MOD(I,K)+1)(out), */

/*     and for S(I) = -1, */

/*                                  H */
/*             Q(:,:,MOD(I,K)+1)(in)   A(:,:,I)(in)   Q(:,:,I)(in) */
/*                                   H */
/*         =   Q(:,:,MOD(I,K)+1)(out)  A(:,:,I)(out)  Q(:,:,I)(out). */


/*     ARGUMENTS */

/*     Mode Parameters */

/*     COMPQ   (input) CHARACTER*1 */
/*             = 'N': do not modify Q. */
/*             = 'V': modify the array Q by the unitary transformations */
/*                    that are applied to the matrices in A to reduce them */
/*                    to Hessenberg-triangular form. */
/*             = 'I': like COMPQ='V', except that each matrix in Q will */
/*                    be initialized to the identity first. */

/*     Input/Output Parameters */

/*     K       (input) INTEGER */
/*             The number of matrices in A.  K >= 1. */

/*     N       (input) INTEGER */
/*             Order of each factor in A.  N >= 0. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             It is assumed that each factor in A is already upper */
/*             triangular in rows and columns 1:ILO-1 and IHI+1:N. */
/*             1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */

/*     S       (input) INTEGER array, dimension (K) */
/*             The leading K elements of this array must contain the */
/*             signatures of the factors. Each entry in S must be either */
/*             1 or -1. By definition, S(1) must be set to 1. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K) */
/*             On entry, the leading N-by-N-by-K part of this array must */
/*             contain the factors of the general product to be reduced. */
/*             On exit, A(:,:,1) is overwritten by an upper Hessenberg */
/*             matrix and each A(:,:,I) for I not equal to 1 is */
/*             overwritten by an upper triangular matrix. */

/*     LDA1    (input) INTEGER */
/*             The first leading dimension of A. LDA1 >= max(1,N) */

/*     LDA2    (input) INTEGER */
/*             The second leading dimension of A. LDA2 >= max(1,N) */

/*     Q       (input/output) COMPLEX*16 array, dimension (LDQ1,LDQ2,K) */
/*             If COMPQ='N': Q is not referenced. */
/*             If COMPQ='I': On entry, Q need not to be set, and on exit */
/*                           it contains the unitary transformations. */
/*             If COMPQ='V': On entry, Q must contain unitary matrices, */
/*                           and on exit this is overwritten by the */
/*                           updated transformations. */

/*     LDQ1    (input) INTEGER */
/*             The first leading dimension of Q. LDQ1 >= max(1,N) */

/*     LDQ2    (input) INTEGER */
/*             The second leading dimension of Q. LDQ2 >= max(1,N) */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the minimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) returns the optimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK.  LZWORK >= MAX(1,2*N). */
/*             For optimal performance this value should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */

/*     METHOD */

/*     A slightly modified version of the periodic Hessenberg reduction */
/*     presented in [1] is used. For more details see [2]. */

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
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --s;
    a_dim1 = *lda1;
    a_dim2 = *lda2;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    q_dim1 = *ldq1;
    q_dim2 = *ldq2;
    q_offset = 1 + q_dim1 * (1 + q_dim2);
    q -= q_offset;
    --dwork;
    --zwork;

    /* Function Body */
    *info = 0;
    wantq = lsame(compq, "V") || lsame(compq, "I");
    initq = lsame(compq, "I");

/*     Test the input scalar arguments. */

    if (! wantq && ! initq && ! lsame(compq, "N")) {
	*info = -1;
    } else if (*k < 1) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1) {
	*info = -4;
    } else if (*ihi > *n || *ihi < *ilo - 1) {
	*info = -5;
    } else {
	sok = s[1] == 1;
	i__1 = *k;
	for (l = 2; l <= i__1; ++l) {
	    sok = s[l] == 1 || s[l] == -1;
/* L10: */
	}
	if (! sok) {
	    *info = -6;
	} else if (*lda1 < max(1,*n)) {
	    *info = -8;
	} else if (*lda2 < max(1,*n)) {
	    *info = -9;
	} else if (wantq && *ldq1 < max(1,*n)) {
	    *info = -11;
	} else if (wantq && *ldq2 < max(1,*n)) {
	    *info = -12;
	} else if (*ldwork < max(1,*n)) {
	    *info = -14;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = 1, i__2 = *n << 1;
	    if (*lzwork < max(i__1,i__2)) {
		*info = -16;
	    }
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla("ZPGHRD", &i__1);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	zwork[1].r = 1., zwork[1].i = 0.;
	return 0;
    }
    wrkopt = *n << 1;

/*     Initialize Q if desired. */

    if (initq) {
	zlaset("Full", n, n, &c_b2, &c_b1, &q[q_offset], ldq1);
    }

/*     Transform A(2,:,:),...,A(K,:,:) to upper triangular form. */

    for (l = *k; l >= 2; --l) {
	if (s[l] == 1) {

/*           Compute a QR Decomposition of A(:,:,L). */

	    irows = *ihi + 1 - *ilo;
	    icols = *n + 1 - *ilo;
	    i__1 = *lzwork - *n;
	    zgeqrf(&irows, &icols, &a[*ilo + (*ilo + l * a_dim2) * a_dim1], 
		    lda1, &zwork[1], &zwork[*n + 1], &i__1, &ierr);
/* Computing MAX */
	    i__3 = *n + 1;
	    i__1 = wrkopt, i__2 = *n + (integer) zwork[i__3].r;
	    wrkopt = max(i__1,i__2);

/*           Apply transformation to A(:,:,L-1). */

	    if (s[l - 1] == 1) {
		i__1 = *lzwork - *n;
		zunmqr("Right", "No transpose", ihi, &irows, &irows, &a[*ilo 
			+ (*ilo + l * a_dim2) * a_dim1], lda1, &zwork[1], &a[(
			*ilo + (l - 1) * a_dim2) * a_dim1 + 1], lda1, &zwork[*
			n + 1], &i__1, &ierr);
/* Computing MAX */
		i__3 = *n + 1;
		i__1 = wrkopt, i__2 = *n + (integer) zwork[i__3].r;
		wrkopt = max(i__1,i__2);
	    } else {
		i__1 = *lzwork - *n;
		zunmqr("Left", "Complex transpose", &irows, &icols, &irows, &
			a[*ilo + (*ilo + l * a_dim2) * a_dim1], lda1, &zwork[
			1], &a[*ilo + (*ilo + (l - 1) * a_dim2) * a_dim1], 
			lda1, &zwork[*n + 1], &i__1, &ierr);
/* Computing MAX */
		i__3 = *n + 1;
		i__1 = wrkopt, i__2 = *n + (integer) zwork[i__3].r;
		wrkopt = max(i__1,i__2);
	    }

/*           Update transformation matrix Q(:,:,L). */

	    if (initq) {
/*               CALL ZLASET( 'Full', N, N, CZERO, CONE, Q(1,1,L), LDQ1 ) */
		if (irows > 1) {
		    i__1 = irows - 1;
		    i__2 = irows - 1;
		    zlacpy("Lower", &i__1, &i__2, &a[*ilo + 1 + (*ilo + l * 
			    a_dim2) * a_dim1], lda1, &q[*ilo + 1 + (*ilo + l *
			     q_dim2) * q_dim1], ldq1);
		    i__1 = *lzwork - *n;
		    zungqr(&irows, &irows, &irows, &q[*ilo + (*ilo + l * 
			    q_dim2) * q_dim1], ldq1, &zwork[1], &zwork[*n + 1]
			    , &i__1, &ierr);
/* Computing MAX */
		    i__3 = *n + 1;
		    i__1 = wrkopt, i__2 = *n + (integer) zwork[i__3].r;
		    wrkopt = max(i__1,i__2);
		}
	    } else if (wantq) {
		i__1 = *lzwork - *n;
		zunmqr("Right", "No transpose", n, &irows, &irows, &a[*ilo + 
			(*ilo + l * a_dim2) * a_dim1], lda1, &zwork[1], &q[(*
			ilo + l * q_dim2) * q_dim1 + 1], ldq1, &zwork[*n + 1],
			 &i__1, &ierr);
/* Computing MAX */
		i__3 = *n + 1;
		i__1 = wrkopt, i__2 = *n + (integer) zwork[i__3].r;
		wrkopt = max(i__1,i__2);
	    }
	    if (irows > 1) {
		i__1 = irows - 1;
		i__2 = irows - 1;
		zlaset("Low", &i__1, &i__2, &c_b2, &c_b2, &a[*ilo + 1 + (*
			ilo + l * a_dim2) * a_dim1], lda1);
	    }
	} else {

/*           Compute an RQ Decomposition of A(:,:,L). */

	    icols = *ihi + 1 - *ilo;
	    i__1 = *lzwork - *n;
	    zgerqf(ihi, &icols, &a[(*ilo + l * a_dim2) * a_dim1 + 1], lda1, &
		    zwork[1], &zwork[*n + 1], &i__1, &ierr);
/* Computing MAX */
	    i__3 = *n + 1;
	    i__1 = wrkopt, i__2 = *n + (integer) zwork[i__3].r;
	    wrkopt = max(i__1,i__2);

/*           Apply transformation to A(:,:,L-1). */

	    if (s[l - 1] == 1) {
		i__1 = *lzwork - *n;
		zunmrq("Right", "Complex Transpose", ihi, &icols, &icols, &a[
			*ilo + (*ilo + l * a_dim2) * a_dim1], lda1, &zwork[1],
			 &a[(*ilo + (l - 1) * a_dim2) * a_dim1 + 1], lda1, &
			zwork[*n + 1], &i__1, &ierr);
/* Computing MAX */
		i__3 = *n + 1;
		i__1 = wrkopt, i__2 = *n + (integer) zwork[i__3].r;
		wrkopt = max(i__1,i__2);
	    } else {
		i__1 = *n + 1 - *ilo;
		i__2 = *lzwork - *n;
		zunmrq("Left", "No transpose", &icols, &i__1, &icols, &a[*
			ilo + (*ilo + l * a_dim2) * a_dim1], lda1, &zwork[1], 
			&a[*ilo + (*ilo + (l - 1) * a_dim2) * a_dim1], lda1, &
			zwork[*n + 1], &i__2, &ierr);
/* Computing MAX */
		i__3 = *n + 1;
		i__1 = wrkopt, i__2 = *n + (integer) zwork[i__3].r;
		wrkopt = max(i__1,i__2);
	    }

/*           Update transformation matrix Q(:,:,L). */

	    if (initq) {
		zlaset("Full", n, n, &c_b2, &c_b1, &q[(l * q_dim2 + 1) * 
			q_dim1 + 1], ldq1);
	    }
	    if (initq || wantq) {
		i__1 = *lzwork - *n;
		zunmrq("Right", "Complex transpose", n, &icols, &icols, &a[(
			l * a_dim2 + 1) * a_dim1 + 1], lda1, &zwork[1], &q[(*
			ilo + l * q_dim2) * q_dim1 + 1], ldq1, &zwork[*n + 1],
			 &i__1, &ierr);
	    }
/* Computing MAX */
	    i__3 = *n + 1;
	    i__1 = wrkopt, i__2 = *n + (integer) zwork[i__3].r;
	    wrkopt = max(i__1,i__2);
	    if (icols > 1) {
		i__1 = icols - 1;
		i__2 = icols - 1;
		zlaset("Low", &i__1, &i__2, &c_b2, &c_b2, &a[*ilo + 1 + (*
			ilo + l * a_dim2) * a_dim1], lda1);
	    }
	}
/* L30: */
    }

/*     Reduce A(:,:,1) to upper Hessenberg form. */

    i__1 = *ihi - 2;
    for (jcol = *ilo; jcol <= i__1; ++jcol) {

/*        Annihilate all elements below A(JCOL+1,JCOL,1). */

	i__2 = jcol + 2;
	for (jrow = *ihi; jrow >= i__2; --jrow) {

	    i__3 = jrow - 1 + (jcol + a_dim2) * a_dim1;
	    temp.r = a[i__3].r, temp.i = a[i__3].i;
	    zlartg(&temp, &a[jrow + (jcol + a_dim2) * a_dim1], &cs, &sn, &a[
		    jrow - 1 + (jcol + a_dim2) * a_dim1]);
	    i__3 = jrow + (jcol + a_dim2) * a_dim1;
	    a[i__3].r = 0., a[i__3].i = 0.;
	    i__3 = *n - jcol;
	    zrot(&i__3, &a[jrow - 1 + (jcol + 1 + a_dim2) * a_dim1], lda1, &
		    a[jrow + (jcol + 1 + a_dim2) * a_dim1], lda1, &cs, &sn);
	    dwork[jrow] = cs;
	    i__3 = jrow;
	    zwork[i__3].r = sn.r, zwork[i__3].i = sn.i;
/* L40: */
	}

	if (wantq) {
	    i__2 = jcol + 2;
	    for (jrow = *ihi; jrow >= i__2; --jrow) {
		d_cnjg(&z__1, &zwork[jrow]);
		zrot(n, &q[(jrow - 1 + q_dim2) * q_dim1 + 1], &c__1, &q[(
			jrow + q_dim2) * q_dim1 + 1], &c__1, &dwork[jrow], &
			z__1);
/* L50: */
	    }
	}

/*        Propagate transformations through A(:,:,K),...,A(:,:,2). */

	for (l = *k; l >= 2; --l) {

	    if (s[l] == 1) {
		i__2 = jcol + 2;
		for (jrow = *ihi; jrow >= i__2; --jrow) {
		    d_cnjg(&z__1, &zwork[jrow]);
		    zrot(&jrow, &a[(jrow - 1 + l * a_dim2) * a_dim1 + 1], &
			    c__1, &a[(jrow + l * a_dim2) * a_dim1 + 1], &c__1,
			     &dwork[jrow], &z__1);
		    i__3 = jrow - 1 + (jrow - 1 + l * a_dim2) * a_dim1;
		    temp.r = a[i__3].r, temp.i = a[i__3].i;
		    zlartg(&temp, &a[jrow + (jrow - 1 + l * a_dim2) * a_dim1]
			    , &cs, &sn, &a[jrow - 1 + (jrow - 1 + l * a_dim2) 
			    * a_dim1]);
		    i__3 = jrow + (jrow - 1 + l * a_dim2) * a_dim1;
		    a[i__3].r = 0., a[i__3].i = 0.;
		    i__3 = *n - jrow + 1;
		    zrot(&i__3, &a[jrow - 1 + (jrow + l * a_dim2) * a_dim1], 
			    lda1, &a[jrow + (jrow + l * a_dim2) * a_dim1], 
			    lda1, &cs, &sn);
		    dwork[jrow] = cs;
		    i__3 = jrow;
		    zwork[i__3].r = sn.r, zwork[i__3].i = sn.i;
/* L60: */
		}
	    } else {
		i__2 = jcol + 2;
		for (jrow = *ihi; jrow >= i__2; --jrow) {
		    i__3 = *n + 2 - jrow;
		    zrot(&i__3, &a[jrow - 1 + (jrow - 1 + l * a_dim2) * 
			    a_dim1], lda1, &a[jrow + (jrow - 1 + l * a_dim2) *
			     a_dim1], lda1, &dwork[jrow], &zwork[jrow]);
		    i__3 = jrow + (jrow + l * a_dim2) * a_dim1;
		    temp.r = a[i__3].r, temp.i = a[i__3].i;
		    zlartg(&temp, &a[jrow + (jrow - 1 + l * a_dim2) * a_dim1]
			    , &cs, &sn, &a[jrow + (jrow + l * a_dim2) * 
			    a_dim1]);
		    i__3 = jrow + (jrow - 1 + l * a_dim2) * a_dim1;
		    a[i__3].r = 0., a[i__3].i = 0.;
		    i__3 = jrow - 1;
		    zrot(&i__3, &a[(jrow + l * a_dim2) * a_dim1 + 1], &c__1, 
			    &a[(jrow - 1 + l * a_dim2) * a_dim1 + 1], &c__1, &
			    cs, &sn);
		    dwork[jrow] = cs;
		    i__3 = jrow;
		    z__1.r = -sn.r, z__1.i = -sn.i;
		    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L70: */
		}
	    }

	    if (wantq) {
		i__2 = jcol + 2;
		for (jrow = *ihi; jrow >= i__2; --jrow) {
		    d_cnjg(&z__1, &zwork[jrow]);
		    zrot(n, &q[(jrow - 1 + l * q_dim2) * q_dim1 + 1], &c__1, 
			    &q[(jrow + l * q_dim2) * q_dim1 + 1], &c__1, &
			    dwork[jrow], &z__1);
/* L80: */
		}
	    }
/* L90: */
	}

/*        Apply transformations to A(:,:,1). */

	i__2 = jcol + 2;
	for (jrow = *ihi; jrow >= i__2; --jrow) {
	    d_cnjg(&z__1, &zwork[jrow]);
	    zrot(ihi, &a[(jrow - 1 + l * a_dim2) * a_dim1 + 1], &c__1, &a[(
		    jrow + l * a_dim2) * a_dim1 + 1], &c__1, &dwork[jrow], &
		    z__1);
/* L100: */
	}
/* L110: */
    }
    dwork[1] = (doublereal) (*n);
    d__1 = (doublereal) wrkopt;
    z__1.r = d__1, z__1.i = 0.;
    zwork[1].r = z__1.r, zwork[1].i = z__1.i;
    return 0;
/* *** Last line of ZPGHRD *** */
} /* zpghrd_ */

