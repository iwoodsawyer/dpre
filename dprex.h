/* Starting from version 7.8, MATLAB BLAS expects ptrdiff_t arguments for integers */
#if MATLAB_VERSION >= 0x0708
#include <stddef.h>
#include <stdlib.h>
#endif
#include <string.h>

/* Define MX_HAS_INTERLEAVED_COMPLEX for version <9.4 */
#ifndef MX_HAS_INTERLEAVED_COMPLEX
#define MX_HAS_INTERLEAVED_COMPLEX 0
#endif

/* Starting from version 7.6, MATLAB BLAS is seperated */
#if MATLAB_VERSION >= 0x0705
#include <blas.h>
#endif
#include <lapack.h>
#include "f2c.h"

// Conversion of optimal problems with coupling weighting terms to standard problems
#define sb02mt FORTRAN_WRAPPER(sb02mt)
extern void sb02mt(
        const char *jobg,
        const char *jobl,
        const char *fact,
        const char *uplo,
        const ptrdiff_t *n,
        const ptrdiff_t *m,
        double *a,
        const ptrdiff_t *lda,
        double *b,
        const ptrdiff_t *ldb,
        double *q,
        const ptrdiff_t *ldq,
        double *r,
        const ptrdiff_t *ldr,
        double *l,
        const ptrdiff_t *ldl,
        ptrdiff_t *ipiv,
        ptrdiff_t *oufact,
        double *g,
        const ptrdiff_t *ldg,
        ptrdiff_t *iwork,
        double *dwork,
        const ptrdiff_t *ldwork,
        ptrdiff_t *info
        );

// Conversion of optimal problems with coupling weighting terms to standard problems (more flexibility)
#define sb02mx FORTRAN_WRAPPER(sb02mx)
extern void sb02mx(
        const char *jobg,
        const char *jobl,
        const char *fact,
        const char *uplo,
        const char *trans,
        const char *flag,
        const char *def,
        const ptrdiff_t *n,
        const ptrdiff_t *m,
        double *a,
        const ptrdiff_t *lda,
        double *b,
        const ptrdiff_t *ldb,
        double *q,
        const ptrdiff_t *ldq,
        double *r,
        const ptrdiff_t *ldr,
        double *l,
        const ptrdiff_t *ldl,
        ptrdiff_t *ipiv,
        ptrdiff_t *oufact,
        double *g,
        const ptrdiff_t *ldg,
        ptrdiff_t *iwork,
        double *dwork,
        const ptrdiff_t *ldwork,
        ptrdiff_t *info
        );

// Constructing the 2n-by-2n Hamiltonian or symplectic matrix for linear-quadratic optimization problems
#define sb02mu FORTRAN_WRAPPER(sb02mu)
extern void sb02mu(
        const char *dico,
        const char *hinv,
        const char *uplo,
        const ptrdiff_t *n,
        double *a,
        const ptrdiff_t *lda,
        const double *g,
        const ptrdiff_t *ldg,
        const double *q,
        const ptrdiff_t *ldq,
        double *s,
        const ptrdiff_t *lds,
        ptrdiff_t *iwork,
        double *dwork,
        const ptrdiff_t *ldwork,
        ptrdiff_t *info
        );

// Constructing the 2n-by-2n Hamiltonian or symplectic matrix for linear-quadratic optimization problems (improved)
#define sb02ru FORTRAN_WRAPPER(sb02ru)
extern void sb02ru(
        const char *dico,
        const char *hinv,
        const char *trana,
        const char *uplo,
        const ptrdiff_t *n,
        const double *a,
        const ptrdiff_t *lda,
        double *g,
        const ptrdiff_t *ldg,
        double *q,
        const ptrdiff_t *ldq,
        double *s,
        const ptrdiff_t *lds,
        ptrdiff_t *iwork,
        double *dwork,
        const ptrdiff_t *ldwork,
        ptrdiff_t *info
        );

// Constructing the extended Hamiltonian or symplectic matrix pairs for linear-quadratic optimization problems, and compressing them to 2N-by-2N matrices
#define sb02oy FORTRAN_WRAPPER(sb02oy)
extern void sb02oy(
        const char *type,
        const char *dico,
        const char *jobb,
        const char *fact,
        const char *uplo,
        const char *jobl,
        const char *jobe,
        const ptrdiff_t *n,
        const ptrdiff_t *m,        
        const ptrdiff_t *p,        
        const double *a,
        const ptrdiff_t *lda,
        const double *b,
        const ptrdiff_t *ldb,        
        const double *q,
        const ptrdiff_t *ldq,
        const double *r,
        const ptrdiff_t *ldr,
        const double *l,
        const ptrdiff_t *ldl,
        const double *e,
        const ptrdiff_t *lde,
        double *af,
        const ptrdiff_t *ldaf,
        double *bf,
        const ptrdiff_t *ldbf,
        double *tol,
        ptrdiff_t *iwork,
        double *dwork,
        const ptrdiff_t *ldwork,
        ptrdiff_t *info
        );

// Optimal state feedback matrix for an optimal control problem
#define sb02nd FORTRAN_WRAPPER(sb02nd)
extern void sb02nd(
        const char *dico,
        const char *fact,
        const char *uplo,
        const char *jobl,
        const ptrdiff_t *n,
        const ptrdiff_t *m,
        const ptrdiff_t *p,
        const double *a,
        const ptrdiff_t *lda,
        double *b,
        const ptrdiff_t *ldb,
        double *r,
        const ptrdiff_t *ldr,
        ptrdiff_t *ipiv,
        const double *l,
        const ptrdiff_t *ldl,
        double *x,
        const ptrdiff_t *ldx,
        const double *rnorm,
        double *f,
        const ptrdiff_t *ldf,
        ptrdiff_t *oufact,
        ptrdiff_t *iwork,
        double *dwork,
        const ptrdiff_t *ldwork,
        ptrdiff_t *info
        );

// Solution of continuous- or discrete-time algebraic Riccati equations for descriptor systems
#define sg02nd FORTRAN_WRAPPER(sg02nd)
extern void sg02nd(
        const char *dico,
        const char *jobe,
        const char *job,
        const char *jobx,
        const char *fact,
        const char *uplo,
        const char *jobl,
        const char *trans,
        const ptrdiff_t *n,
        const ptrdiff_t *m,        
        const ptrdiff_t *p,
        const double *a,
        const ptrdiff_t *lda,
        const double *e,
        const ptrdiff_t *lde,
        double *b,
        const ptrdiff_t *ldb,
        double *r,
        const ptrdiff_t *ldr,
        ptrdiff_t *ipiv,
        const double *l,
        const ptrdiff_t *ldl,
        double *x,
        const ptrdiff_t *ldx,
        const double *rnorm,
        double *k,
        const ptrdiff_t *ldk,
        double *h,
        const ptrdiff_t *ldh,
        double *xe,
        const ptrdiff_t *ldxe,
        ptrdiff_t *oufact,
        ptrdiff_t *iwork,
        double *dwork,
        const ptrdiff_t *ldwork,
        ptrdiff_t *info
        );

// Solution of linear equations X op(A) = B
#define mb02vd FORTRAN_WRAPPER(mb02vd)
extern void mb02vd(
        const char *trans,
        const ptrdiff_t *m,
        const ptrdiff_t *n,
        double *a,
        const ptrdiff_t *lda,
        ptrdiff_t *ipiv,
        double *b,
        const ptrdiff_t *ldb,
        ptrdiff_t *info
        );

// Periodic Hessenberg form of a product of p matrices using orthogonal similarity transformations
#define mb03vd FORTRAN_WRAPPER(mb03vd)
extern void mb03vd(
        const ptrdiff_t *n,
        const ptrdiff_t *p,
        const ptrdiff_t *ilo,
        const ptrdiff_t *ihi,
        double *a,
        const ptrdiff_t *lda1,
        const ptrdiff_t *lda2,
        double *tau,
        const ptrdiff_t *ldtau,
        double *dwork,
        ptrdiff_t *info
        );

// Orthogonal matrices for reduction to periodic Hessenberg form of a product of matrices
#define mb03vy FORTRAN_WRAPPER(mb03vy)
extern void mb03vy(
        const ptrdiff_t *n,
        const ptrdiff_t *p,
        const ptrdiff_t *ilo,
        const ptrdiff_t *ihi,
        double *a,
        const ptrdiff_t *lda1,
        const ptrdiff_t *lda2,
        double *tau,
        const ptrdiff_t *ldtau,
        double *dwork,
        const ptrdiff_t *ldwork,
        ptrdiff_t *info
        );

// Schur decomposition and eigenvalues of a product of matrices in periodic Hessenberg form
#define mb03wd FORTRAN_WRAPPER(mb03wd)
extern void mb03wd(
        const char *job,
        const char *compz,
        const ptrdiff_t *n,
        const ptrdiff_t *p,
        const ptrdiff_t *ilo,
        const ptrdiff_t *ihi,
        const ptrdiff_t *iloz,
        const ptrdiff_t *ihiz,
        double *h,
        const ptrdiff_t *ldh1,
        const ptrdiff_t *ldh2,
        double *z,
        const ptrdiff_t *ldz1,
        const ptrdiff_t *ldz2,
        double *wr,
        double *wi,
        double *dwork,
        const ptrdiff_t *ldwork,
        ptrdiff_t *info
        );


// Computes the general product of K complex scalars trying to avoid over- and underflow. 
#define zlapr1 FORTRAN_WRAPPER(zlapr1)
extern int zlapr1(
        doublereal *base, 
        integer *k, 
        integer *s, 
        doublecomplex *a, 
        integer *inca, 
        doublecomplex *alpha, 
        doublecomplex *beta, 
        integer *scal
        );

// Finding the eigenvalues of the complex generalized matrix product. 
#define zpgeqz FORTRAN_WRAPPER(zpgeqz)
extern int zpgeqz(
        char *job, 
        char *compq, 
        integer *k, 
        integer *n,	
        integer *ilo, 
        integer *ihi, 
        integer *s, 
        doublecomplex *a, 
        integer *lda1, 
        integer *lda2, 
        doublecomplex *alpha, 
        doublecomplex *beta,
        integer *scal,
        doublecomplex *q,
        integer *ldq1,
        integer *ldq2,
        doublereal *dwork,
        integer *ldwork,
        doublecomplex *zwork,
        integer *lzwork, 
        integer *info
        );

// Swaps adjacent diagonal 1-by-1 blocks in a complex generalized matrix product. 
#define zpgex2 FORTRAN_WRAPPER(zpgex2)
extern int zpgex2(
        logical *wantq, 
        integer *k, 
        integer *n, 
        integer *j, 
        integer *s, 
        doublecomplex *a, 
        integer *lda1, 
        integer *lda2, 
        doublecomplex *q, 
        integer *ldq1, 
        integer *ldq2, 
        doublecomplex *zwork,
        integer *info
        );

// Swaps adjacent diagonal 1-by-1 blocks in a complex generalized matrix product. 
#define zpghrd FORTRAN_WRAPPER(zpghrd)
extern int zpghrd(
        char *compq, 
        integer *k, 
        integer *n, 
        integer *ilo, 
        integer *ihi, 
        integer *s, 
        doublecomplex *a, 
        integer *lda1, 
        integer *lda2, 
        doublecomplex *q, 
        integer *ldq1, 
        integer *ldq2, 
        doublereal *dwork, 
        integer *ldwork, 
        doublecomplex *zwork, 
        integer *lzwork, 
        integer *info
        );

// Reorders the periodic Schur decomposition of a complex generalized matrix product. 
#define zpgord FORTRAN_WRAPPER(zpgord)
extern int zpgord(
        logical *wantq, 
        integer *k, 
        integer *n, 
        integer *s, 
        logical *select, 
        doublecomplex *a, 
        integer *lda1, 
        integer *lda2, 
        doublecomplex *alpha,
        doublecomplex *beta,
        integer *scal, 
        doublecomplex *q, 
        integer *ldq1, 
        integer *ldq2, 
        integer *m,
        doublecomplex *zwork,
        integer *lzwork,
        integer *info
        );

