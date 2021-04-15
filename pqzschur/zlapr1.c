/* zlapr1.f -- translated by f2c (version 20041007).
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

/* Subroutine */ int zlapr1(doublereal *base, integer *k, integer *s, 
	doublecomplex *a, integer *inca, doublecomplex *alpha, doublecomplex *
	beta, integer *scal)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, inda;


/*     PURPOSE */

/*     Computes the general product of K complex scalars trying to avoid */
/*     over- and underflow. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     BASE    (input)  DOUBLE PRECISION */
/*             Machine base. */

/*     K       (input)  INTEGER */
/*             The number of scalars.  K >= 0. */

/*     S       (input)  INTEGER array, dimension (K) */
/*             The signature array. Each entry of S must be 1 or -1. */

/*     A       (input)  COMPLEX*16 array, dimension (K) */
/*             Vector of real scalars. */

/*     INCA    (input)  INTEGER */
/*             Increment for the array A. incA <> 0. */

/*     ALPHA   (output)  COMPLEX*16 */
/*             ALPHA is a real scalar with 1.0 <= ABS(ALPHA) < BASE such */
/*             that */

/*                  ALPHA / BETA * BASE**(SCAL) */

/*             is the general product of A. */

/*     BETA    (output)  COMPLEX*16 */
/*             BETA is either 0.0 or 1.0. */
/*             See also the description of ALPHA. */

/*     SCAL    (output)  INTEGER */
/*             Scaling factor, see ALPHA. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, Dec. 2002. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Array Arguments .. */
/*     .. Scalar Arguments .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --a;
    --s;

    /* Function Body */
    alpha->r = 1., alpha->i = 0.;
    beta->r = 1., beta->i = 0.;
    *scal = 0;
    inda = 1;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s[i__] == 1) {
	    i__2 = inda;
	    z__1.r = alpha->r * a[i__2].r - alpha->i * a[i__2].i, z__1.i = 
		    alpha->r * a[i__2].i + alpha->i * a[i__2].r;
	    alpha->r = z__1.r, alpha->i = z__1.i;
	} else {
	    i__2 = inda;
	    if (a[i__2].r == 0. && a[i__2].i == 0.) {
		beta->r = 0., beta->i = 0.;
	    } else {
		z_div(&z__1, alpha, &a[inda]);
		alpha->r = z__1.r, alpha->i = z__1.i;
	    }
	}
	if (z_abs(alpha) == 0.) {
	    *scal = 0;
	    alpha->r = 0., alpha->i = 0.;
	} else {
L10:
	    if (z_abs(alpha) >= 1.) {
		goto L20;
	    }
	    z__2.r = *base, z__2.i = 0.;
	    z__1.r = z__2.r * alpha->r - z__2.i * alpha->i, z__1.i = z__2.r * 
		    alpha->i + z__2.i * alpha->r;
	    alpha->r = z__1.r, alpha->i = z__1.i;
	    --(*scal);
	    goto L10;
L20:
L30:
	    if (z_abs(alpha) < *base) {
		goto L40;
	    }
	    z__2.r = *base, z__2.i = 0.;
	    z_div(&z__1, alpha, &z__2);
	    alpha->r = z__1.r, alpha->i = z__1.i;
	    ++(*scal);
	    goto L30;
L40:
	    ;
	}
	inda += *inca;
/* L50: */
    }

    return 0;

/* *** Last line of ZLAPR1 *** */
} /* zlapr1_ */

