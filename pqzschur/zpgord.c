/* zpgord.f -- translated by f2c (version 20041007).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "dprex.h"

/* Table of constant values */

static logical c_true = TRUE_;
static integer c__1 = 1;

/* Subroutine */ int zpgord(logical *wantq, integer *k, integer *n, integer *
	s, logical *select, doublecomplex *a, integer *lda1, integer *lda2, 
	doublecomplex *alpha, doublecomplex *beta, integer *scal, 
	doublecomplex *q, integer *ldq1, integer *ldq2, integer *m, 
	doublecomplex *zwork, integer *lzwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_dim2, a_offset, q_dim1, q_dim2, q_offset, i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, l, js;
    static logical sok;
    static doublereal base, abst;
    static integer ierr;
    static doublecomplex temp;
    static doublereal safmin;


/*     PURPOSE */

/*     ZPGORD reorders the periodic Schur decomposition of a complex */
/*     generalized matrix product */

/*                                 S(2)                 S(K) */
/*          A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K), */

/*     (in terms of unitary equivalence transformations), so that a */
/*     selected cluster of eigenvalues appears in the leading diagonal */
/*     blocks of the matrix product. The leading columns of the */
/*     orthogonal factors contained in Q form unitary bases of the */
/*     corresponding periodic eigenspaces (deflating subspaces). */
/*     A must be in periodic Schur form, that is, all factors of A must */
/*     be upper triangular. */

/*     If WANTQ = .TRUE., then the unitary factors are computed and */
/*     stored in the array Q so that for S(I) = 1, */

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

/*     WANTQ   (input) LOGICAL */
/*             = .FALSE.: do not modify Q; */
/*             = .TRUE. : modify the array Q by the unitary */
/*                        transformations that are applied to the */
/*                        matrices in A for reordering. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of factors.  K >= 1. */

/*     N       (input)  INTEGER */
/*             The order of each factor in A.  N >= 0. */

/*     S       (input)  INTEGER array, dimension (K) */
/*             The leading K elements of this array must contain the */
/*             signatures of the factors. Each entry in S must be either */
/*             1 or -1. By definition, S(1) must be set to 1. */

/*     SELECT  (input) LOGICAL array, dimension (N) */
/*             SELECT specifies the eigenvalues in the selected cluster. */
/*             To select the eigenvalue corresponding to the (j,j) */
/*             diagonal entries, SELECT(j) must be set to .TRUE.. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K) */
/*             On entry, the leading N-by-N-by-K part of this array */
/*             must contain the factors in periodic Schur form */
/*             form, that is, all factors are upper triangular. */
/*             On exit, if INFO = 0, the leading N-by-N-by-K part of */
/*             this array contains the factors of the reordered periodic */
/*             Schur form. Moreover, A(:,:,2), ..., A(:,:,K) are */
/*             normalized so that their diagonals contain nonnegative */
/*             real numbers. */

/*     LDA1    (input) INTEGER */
/*             The first leading dimension of A.  LDA1 >= MAX(1,N). */

/*     LDA2    (input) INTEGER */
/*             The second leading dimension of A.  LDA2 >= MAX(1,N). */

/*     Q       (input/output) COMPLEX*16 array, dimension (LDQ1,LDQ2,K) */
/*             On entry, if WANTQ = .TRUE., the leading N-by-N-by-K part */
/*             of this array must contain the initial unitary factors */
/*             as described in (1)-(2). */
/*             On exit, if WANTQ = .TRUE., the leading N-by-N-by-K part */
/*             of this array contains the modified orthogonal factors as */
/*             described in (1)-(2). */

/*     LDQ1    (input)  INTEGER */
/*             The first leading dimension of Q. */
/*             If WANTQ = .TRUE.,  LDQ1 >= MAX(1,N). */

/*     LDQ2    (input)  INTEGER */
/*             The second leading dimension of Q. */
/*             If WANTQ = .TRUE.,  LDQ2 >= MAX(1,N). */

/*     M       (output) INTEGER */
/*             The dimension of the specified periodic eigenspace. */

/*     Workspace */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) returns the minimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK.  LZWORK >= MAX( K, N ). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0       : succesful exit; */
/*             < 0       : if INFO = -i, the i-th argument had an */
/*                         illegal value; */
/*             = 1       : the periodic QZ algorithm failed to converge. */

/*     METHOD */

/*     A complex version of the periodic QZ algorithm [1] with */
/*     perfect shifts is used. For more details see [2]. This is not */
/*     a safe method. It is advisable to check whether the eigenvalues */
/*     are really reorderd. */

/*     REFERENCES */

/*     [1] Bojanczyk, A. and Golub, G. H. and Van Dooren, P. */
/*         The periodic Schur decomposition; algorithm and applications. */
/*         In Proc. SPIE Conference, pg. 31-42, vol. 1770, 1992. */

/*     [2] Kressner, D. */
/*         An efficient and reliable implementation of the periodic QZ */
/*         algorithm. In IFAC Workshop on Periodic Control Systems, 2001. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically backward stable. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, Dec. 2002. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --s;
    --select;
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
    --zwork;

    /* Function Body */
    *info = 0;
    safmin = dlamch("SafeMinimum");
    base = dlamch("Base");

/*     Check the scalar input parameters. */

    if (*k < 1) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else {
	sok = s[1] == 1;
	i__1 = *k;
	for (l = 2; l <= i__1; ++l) {
	    sok = s[l] == 1 || s[l] == -1;
/* L10: */
	}
	if (! sok) {
	    *info = -4;
	} else if (*lda1 < max(1,*n)) {
	    *info = -7;
	} else if (*lda2 < max(1,*n)) {
	    *info = -8;
	} else if (*wantq && *ldq1 < max(1,*n)) {
	    *info = -13;
	} else if (*wantq && *ldq2 < max(1,*n)) {
	    *info = -14;
	} else if (*lzwork < max(*k,*n)) {
	    *info = -17;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla("ZPGORD", &i__1);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	zwork[1].r = 1., zwork[1].i = 0.;
	return 0;
    }

/*     Set M to the dimension of the specified deflating subspace. */

    *m = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (select[j]) {
	    ++(*m);
	}
/* L20: */
    }

    js = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (select[j]) {
	    ++js;

/*           Swap the J-th block to position JS. */

	    if (j != js) {
		i__2 = js;
		for (i__ = j - 1; i__ >= i__2; --i__) {
		    zpgex2(&c_true, k, n, &i__, &s[1], &a[a_offset], lda1, 
			    lda2, &q[q_offset], ldq1, ldq2, &zwork[1], &ierr);
		    if (ierr != 0) {
			*info = 1;
			goto L120;
		    }
/* L30: */
		}
	    }
	}
/* L40: */
    }

/*     Rescale matrices and recompute eigenvalues. */

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
		    i__2 = *n - j;
		    zscal(&i__2, &temp, &a[j + (j + 1 + l * a_dim2) * a_dim1]
			    , lda1);
		} else {
		    temp.r = 1., temp.i = 0.;
		}
		i__2 = j;
		zwork[i__2].r = temp.r, zwork[i__2].i = temp.i;
/* L50: */
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
		    zscal(&i__2, &temp, &a[(j + l * a_dim2) * a_dim1 + 1], &
			    c__1);
		} else {
		    temp.r = 1., temp.i = 0.;
		}
		i__2 = j;
		d_cnjg(&z__1, &temp);
		zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L60: */
	    }
	}
	if (*wantq) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		d_cnjg(&z__1, &zwork[j]);
		zscal(n, &z__1, &q[(j + l * q_dim2) * q_dim1 + 1], &c__1);
/* L70: */
	    }
	}
	if (s[l - 1] == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		d_cnjg(&z__1, &zwork[j]);
		zscal(&j, &z__1, &a[(j + (l - 1) * a_dim2) * a_dim1 + 1], &
			c__1);
/* L80: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j + 1;
		zscal(&i__2, &zwork[j], &a[j + (j + (l - 1) * a_dim2) * 
			a_dim1], lda1);
/* L90: */
	    }
	}
/* L100: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *lda1 * *lda2;
	zlapr1(&base, k, &s[1], &a[j + (j + a_dim2) * a_dim1], &i__2, &alpha[
		j], &beta[j], &scal[j]);
/* L110: */
    }

L120:
    d__1 = (doublereal) max(*k,*n);
    z__1.r = d__1, z__1.i = 0.;
    zwork[1].r = z__1.r, zwork[1].i = z__1.i;
    return 0;
/* *** Last line of ZPGORD *** */
} /* zpgord_ */

