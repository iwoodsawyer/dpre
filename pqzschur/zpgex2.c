/* zpgex2.f -- translated by f2c (version 20041007).
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

static doublecomplex c_b1 = {1.,0.};
static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int zpgex2(logical *wantq, integer *k, integer *n, integer *
	j, integer *s, doublecomplex *a, integer *lda1, integer *lda2, 
	doublecomplex *q, integer *ldq1, integer *ldq2, doublecomplex *zwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_dim2, a_offset, q_dim1, q_dim2, q_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer l;
    static doublereal cs;
    static integer ln;
    static doublecomplex sn;
    static doublereal csf;
    static doublecomplex rnd[3], snf;
    static doublereal rhs, ulp;
    static doublecomplex temp;
    static integer iseed[4], jiter;
    static logical usezq;
    static doublereal safmin, safmax;
    static doublereal smlnum;


/*     PURPOSE */

/*     ZPGEX2 swaps adjacent diagonal 1-by-1 blocks in a complex */
/*     generalized matrix product, */

/*                                 S(2)                 S(K) */
/*          A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K), */

/*     by unitary equivalence transformations. A must be in periodic */
/*     Schur form, that is, all factors of A must be upper triangular. */

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
/*             The order of each factor in A.  N >= 2. */

/*     J       (input) INTEGER */
/*             The index of the first block to be swapped.  1 <= J < N. */

/*     S       (input)  INTEGER array, dimension (K) */
/*             The leading K elements of this array must contain the */
/*             signatures of the factors. Each entry in S must be either */
/*             1 or -1. By definition, S(1) must be set to 1. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA1,LDA2,K) */
/*             On entry, the leading N-by-N-by-K part of this array */
/*             must contain the factors in periodic Schur form */
/*             form, that is, all factors are upper triangular. */
/*             On exit, if INFO = 0, the leading N-by-N-by-K part of */
/*             this array contains the factors of the reordered periodic */
/*             Schur form. */

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

/*     Workspace */

/*     ZWORK   COMPLEX*16 array, dimension (K) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0       : succesful exit; */
/*             = 1       : the periodic QZ algorithm failed to converge. */

/*     METHOD */

/*     A complex version of the periodic QZ algorithm [1] with one */
/*     perfect shift is used. For more details see [2]. This is not */
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
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Save perfect shift. */

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
    --zwork;

    /* Function Body */
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	i__2 = l;
	i__3 = *j + (*j + l * a_dim2) * a_dim1;
	zwork[i__2].r = a[i__3].r, zwork[i__2].i = a[i__3].i;
/* L10: */
    }
    *info = 0;
    iseed[0] = 1;
    iseed[1] = 0;
    iseed[2] = 0;
    iseed[3] = 0;
    safmin = dlamch("SafeMinimum");
    safmax = 1. / safmin;
    ulp = dlamch("Precision");
    dlabad(&safmin, &safmax);
    smlnum = safmin * (2 / ulp);

/*     If there are infinite or zero eigenvalues in critical positions, */
/*     then unshifted ZQ iterations must be used. */

    usezq = FALSE_;
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *j + 1 + (*j + 1 + l * a_dim2) * a_dim1;
	i__3 = *j + (*j + l * a_dim2) * a_dim1;
	usezq = usezq || s[l] == 1 && (a[i__2].r == 0. && a[i__2].i == 0.) || 
		s[l] == -1 && (a[i__3].r == 0. && a[i__3].i == 0.);
/* L20: */
    }

/*     Destroy any triangular structure. */

    zlarnv(&c__2, iseed, &c__2, rnd);
    zlartg(rnd, &rnd[1], &cs, &sn, &temp);
    csf = cs;
    snf.r = sn.r, snf.i = sn.i;
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	if (*wantq) {
	    d_cnjg(&z__1, &sn);
	    zrot(n, &q[(*j + l * q_dim2) * q_dim1 + 1], &c__1, &q[(*j + 1 + 
		    l * q_dim2) * q_dim1 + 1], &c__1, &cs, &z__1);
	}
	if (s[l] == 1) {
	    i__2 = *n - *j + 1;
	    zrot(&i__2, &a[*j + (*j + l * a_dim2) * a_dim1], lda1, &a[*j + 1 
		    + (*j + l * a_dim2) * a_dim1], lda1, &cs, &sn);
	} else {
	    i__2 = *j + 1;
	    d_cnjg(&z__1, &sn);
	    zrot(&i__2, &a[(*j + l * a_dim2) * a_dim1 + 1], &c__1, &a[(*j + 
		    1 + l * a_dim2) * a_dim1 + 1], &c__1, &cs, &z__1);
	}
	if (l == *k) {
	    cs = csf;
	    sn.r = snf.r, sn.i = snf.i;
	} else {
	    zlarnv(&c__2, iseed, &c__2, rnd);
	    zlartg(rnd, &rnd[1], &cs, &sn, &temp);
	}
	if (s[l] == 1) {
	    i__2 = *j + 1;
	    d_cnjg(&z__1, &sn);
	    zrot(&i__2, &a[(*j + l * a_dim2) * a_dim1 + 1], &c__1, &a[(*j + 
		    1 + l * a_dim2) * a_dim1 + 1], &c__1, &cs, &z__1);
	} else {
	    i__2 = *n - *j + 1;
	    zrot(&i__2, &a[*j + (*j + l * a_dim2) * a_dim1], lda1, &a[*j + 1 
		    + (*j + l * a_dim2) * a_dim1], lda1, &cs, &sn);
	}
/* L30: */
    }

    if (usezq) {
	for (jiter = 1; jiter <= 10; ++jiter) {
	    i__1 = *k;
	    for (l = 1; l <= i__1; ++l) {
		if (s[l] == 1) {
		    i__2 = *j + 1 + (*j + 1 + l * a_dim2) * a_dim1;
		    temp.r = a[i__2].r, temp.i = a[i__2].i;
		    zlartg(&temp, &a[*j + 1 + (*j + l * a_dim2) * a_dim1], &
			    cs, &sn, &a[*j + 1 + (*j + 1 + l * a_dim2) * 
			    a_dim1]);
		    i__2 = *j + 1 + (*j + l * a_dim2) * a_dim1;
		    a[i__2].r = 0., a[i__2].i = 0.;
		    zrot(j, &a[(*j + 1 + l * a_dim2) * a_dim1 + 1], &c__1, &
			    a[(*j + l * a_dim2) * a_dim1 + 1], &c__1, &cs, &
			    sn);
		    z__1.r = -sn.r, z__1.i = -sn.i;
		    sn.r = z__1.r, sn.i = z__1.i;
		} else {
		    i__2 = *j + (*j + l * a_dim2) * a_dim1;
		    temp.r = a[i__2].r, temp.i = a[i__2].i;
		    zlartg(&temp, &a[*j + 1 + (*j + l * a_dim2) * a_dim1], &
			    cs, &sn, &a[*j + (*j + l * a_dim2) * a_dim1]);
		    i__2 = *j + 1 + (*j + l * a_dim2) * a_dim1;
		    a[i__2].r = 0., a[i__2].i = 0.;
		    i__2 = *n - *j;
		    zrot(&i__2, &a[*j + (*j + 1 + l * a_dim2) * a_dim1], 
			    lda1, &a[*j + 1 + (*j + 1 + l * a_dim2) * a_dim1],
			     lda1, &cs, &sn);
		}
		ln = l + 1;
		if (ln > *k) {
		    ln = 1;
		}
		if (*wantq) {
		    d_cnjg(&z__1, &sn);
		    zrot(n, &q[(*j + ln * q_dim2) * q_dim1 + 1], &c__1, &q[(*
			    j + 1 + ln * q_dim2) * q_dim1 + 1], &c__1, &cs, &
			    z__1);
		}
		if (s[ln] == 1) {
		    i__2 = *n - *j + 1;
		    zrot(&i__2, &a[*j + (*j + ln * a_dim2) * a_dim1], lda1, &
			    a[*j + 1 + (*j + ln * a_dim2) * a_dim1], lda1, &
			    cs, &sn);
		} else {
		    i__2 = *j + 1;
		    d_cnjg(&z__1, &sn);
		    zrot(&i__2, &a[(*j + ln * a_dim2) * a_dim1 + 1], &c__1, &
			    a[(*j + 1 + ln * a_dim2) * a_dim1 + 1], &c__1, &
			    cs, &z__1);
		}
/* L40: */
	    }
/* Computing MAX */
	    d__1 = z_abs(&a[*j + (*j + a_dim2) * a_dim1]), d__2 = z_abs(&a[*j 
		    + 1 + (*j + 1 + a_dim2) * a_dim1]);
	    rhs = max(d__1,d__2);
	    if (rhs == 0.) {
		rhs = z_abs(&a[*j + (*j + 1 + a_dim2) * a_dim1]);
	    }
/* Computing MAX */
	    d__1 = ulp * rhs;
	    rhs = max(d__1,smlnum);
	    if (z_abs(&a[*j + 1 + (*j + a_dim2) * a_dim1]) <= rhs) {
		i__1 = *j + 1 + (*j + a_dim2) * a_dim1;
		a[i__1].r = 0., a[i__1].i = 0.;
		goto L90;
	    }
/* L50: */
	}
	*info = 1;
	goto L90;
    }
    for (jiter = 1; jiter <= 42; ++jiter) {

/*        Complex single shift. */

	if (jiter / 10 * 10 == jiter) {

/*           Random/exceptional shift. */

	    zlarnv(&c__2, iseed, &c__2, rnd);
	    zlartg(rnd, &rnd[1], &cs, &sn, &temp);
	} else {
	    zlartg(&c_b1, &c_b1, &cs, &sn, &temp);
	    for (l = *k; l >= 2; --l) {
		if (s[l] == 1) {
		    i__1 = *j + (*j + l * a_dim2) * a_dim1;
		    z__1.r = cs * a[i__1].r, z__1.i = cs * a[i__1].i;
		    i__2 = l;
		    d_cnjg(&z__3, &sn);
		    z__2.r = zwork[i__2].r * z__3.r - zwork[i__2].i * z__3.i, 
			    z__2.i = zwork[i__2].r * z__3.i + zwork[i__2].i * 
			    z__3.r;
		    zlartg(&z__1, &z__2, &cs, &sn, &temp);
		} else {
		    i__1 = l;
		    z__1.r = cs * zwork[i__1].r, z__1.i = cs * zwork[i__1].i;
		    i__2 = *j + (*j + l * a_dim2) * a_dim1;
		    z__3.r = -a[i__2].r, z__3.i = -a[i__2].i;
		    d_cnjg(&z__4, &sn);
		    z__2.r = z__3.r * z__4.r - z__3.i * z__4.i, z__2.i = 
			    z__3.r * z__4.i + z__3.i * z__4.r;
		    zlartg(&z__1, &z__2, &cs, &sn, &temp);
		    z__1.r = -sn.r, z__1.i = -sn.i;
		    sn.r = z__1.r, sn.i = z__1.i;
		}
/* L60: */
	    }
	    i__1 = *j + (*j + a_dim2) * a_dim1;
	    z__2.r = cs * a[i__1].r, z__2.i = cs * a[i__1].i;
	    d_cnjg(&z__4, &sn);
	    z__3.r = z__4.r * zwork[1].r - z__4.i * zwork[1].i, z__3.i = 
		    z__4.r * zwork[1].i + z__4.i * zwork[1].r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    i__2 = *j + 1 + (*j + a_dim2) * a_dim1;
	    z__5.r = cs * a[i__2].r, z__5.i = cs * a[i__2].i;
	    zlartg(&z__1, &z__5, &cs, &sn, &temp);
	}

/*        Do one QZ sweep. */

	i__1 = *n - *j + 1;
	zrot(&i__1, &a[*j + (*j + a_dim2) * a_dim1], lda1, &a[*j + 1 + (*j + 
		a_dim2) * a_dim1], lda1, &cs, &sn);
	if (*wantq) {
	    d_cnjg(&z__1, &sn);
	    zrot(n, &q[(*j + q_dim2) * q_dim1 + 1], &c__1, &q[(*j + 1 + 
		    q_dim2) * q_dim1 + 1], &c__1, &cs, &z__1);
	}

/*        Propagate rotation through AK, ..., A2 to A1. */

	for (l = *k; l >= 2; --l) {
	    if (s[l] == 1) {
		i__1 = *j + 1;
		d_cnjg(&z__1, &sn);
		zrot(&i__1, &a[(*j + l * a_dim2) * a_dim1 + 1], &c__1, &a[(*
			j + 1 + l * a_dim2) * a_dim1 + 1], &c__1, &cs, &z__1);
		i__1 = *j + (*j + l * a_dim2) * a_dim1;
		temp.r = a[i__1].r, temp.i = a[i__1].i;
		zlartg(&temp, &a[*j + 1 + (*j + l * a_dim2) * a_dim1], &cs, &
			sn, &a[*j + (*j + l * a_dim2) * a_dim1]);
		i__1 = *j + 1 + (*j + l * a_dim2) * a_dim1;
		a[i__1].r = 0., a[i__1].i = 0.;
		i__1 = *n - *j;
		zrot(&i__1, &a[*j + (*j + 1 + l * a_dim2) * a_dim1], lda1, &
			a[*j + 1 + (*j + 1 + l * a_dim2) * a_dim1], lda1, &cs,
			 &sn);
	    } else {
		i__1 = *n - *j + 1;
		zrot(&i__1, &a[*j + (*j + l * a_dim2) * a_dim1], lda1, &a[*j 
			+ 1 + (*j + l * a_dim2) * a_dim1], lda1, &cs, &sn);
		i__1 = *j + 1 + (*j + 1 + l * a_dim2) * a_dim1;
		temp.r = a[i__1].r, temp.i = a[i__1].i;
		zlartg(&temp, &a[*j + 1 + (*j + l * a_dim2) * a_dim1], &cs, &
			sn, &a[*j + 1 + (*j + 1 + l * a_dim2) * a_dim1]);
		i__1 = *j + 1 + (*j + l * a_dim2) * a_dim1;
		a[i__1].r = 0., a[i__1].i = 0.;
		zrot(j, &a[(*j + 1 + l * a_dim2) * a_dim1 + 1], &c__1, &a[(*
			j + l * a_dim2) * a_dim1 + 1], &c__1, &cs, &sn);
		z__1.r = -sn.r, z__1.i = -sn.i;
		sn.r = z__1.r, sn.i = z__1.i;
	    }
	    if (*wantq) {
		d_cnjg(&z__1, &sn);
		zrot(n, &q[(*j + l * q_dim2) * q_dim1 + 1], &c__1, &q[(*j + 
			1 + l * q_dim2) * q_dim1 + 1], &c__1, &cs, &z__1);
	    }
/* L70: */
	}
	i__1 = *j + 1;
	d_cnjg(&z__1, &sn);
	zrot(&i__1, &a[(*j + a_dim2) * a_dim1 + 1], &c__1, &a[(*j + 1 + 
		a_dim2) * a_dim1 + 1], &c__1, &cs, &z__1);

/*        Test for deflation. */

/* Computing MAX */
	d__1 = z_abs(&a[*j + (*j + a_dim2) * a_dim1]), d__2 = z_abs(&a[*j + 1 
		+ (*j + 1 + a_dim2) * a_dim1]);
	rhs = max(d__1,d__2);
	if (rhs == 0.) {
	    rhs = z_abs(&a[*j + (*j + 1 + a_dim2) * a_dim1]);
	}
/* Computing MAX */
	d__1 = ulp * rhs;
	rhs = max(d__1,smlnum);
	if (z_abs(&a[*j + 1 + (*j + a_dim2) * a_dim1]) <= rhs) {
	    i__1 = *j + 1 + (*j + a_dim2) * a_dim1;
	    a[i__1].r = 0., a[i__1].i = 0.;
	    goto L90;
	}
/* L80: */
    }

/*     Not converged. */

    *info = 1;

L90:

/*     Check for singular triangular factors. */

    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
/* Computing MAX */
	d__1 = z_abs(&a[*j + (*j + 1 + l * a_dim2) * a_dim1]), d__2 = z_abs(&
		a[*j + 1 + (*j + 1 + l * a_dim2) * a_dim1]);
	rhs = max(d__1,d__2);
/* Computing MAX */
	d__1 = ulp * rhs;
	rhs = max(d__1,smlnum);
	if (z_abs(&a[*j + (*j + l * a_dim2) * a_dim1]) <= rhs) {
	    i__2 = *j + (*j + l * a_dim2) * a_dim1;
	    a[i__2].r = 0., a[i__2].i = 0.;
	}
/* Computing MAX */
	d__1 = z_abs(&a[*j + (*j + 1 + l * a_dim2) * a_dim1]), d__2 = z_abs(&
		a[*j + (*j + l * a_dim2) * a_dim1]);
	rhs = max(d__1,d__2);
/* Computing MAX */
	d__1 = ulp * rhs;
	rhs = max(d__1,smlnum);
	if (z_abs(&a[*j + 1 + (*j + 1 + l * a_dim2) * a_dim1]) <= rhs) {
	    i__2 = *j + 1 + (*j + 1 + l * a_dim2) * a_dim1;
	    a[i__2].r = 0., a[i__2].i = 0.;
	}
/* L100: */
    }

    return 0;

/* *** Last line of ZPGEX2*** */
} /* zpgex2_ */

