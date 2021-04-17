/*
 * Periodic Discrete-time Algebraic Riccati Equation (DPRE)
 *
 * [X,K,L] = dprex(A,B,Q,R,S,E,method)
 *
 * compile command:
 * mex -O dpre.c libmwblas.lib libmwlapack.lib libmwslicot.lib (>= R2012A)
 *
 * calls the MB03VD/MB03VY and MB03WD named SLICOT functions
 * or the PQZSCHUR functions (default)
 *
 * Ivo Houtzager
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "dprex.h"
#include <stdio.h>
#include <string.h>

void dprex_periodicqr(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize ndim, dims[3] = {1, 1, 1};
    mwSignedIndex m, n, p, nn, ilo, ihi, nre, iscl =0, info = 0; 
    mwSignedIndex ldwork, ldtau, ldmin = -1, ldzero = 0, ldone = 1;
    mwSignedIndex i, k;
    double one = 1.0, half = 0.5, qnorm = 0.0, gnorm = 0.0; 
    double *Apr, *Bpr, *Qpr, *Rpr, *Spr, *Xpr, *Kpr;
    double *Alpr, *Blpr, *Qlpr, *Rlpr, *Slpr, *Glpr, *Hlpr, *Zlpr, *Xlpr, *Tlpr;
    double *dwork, *tau, *wr, *wi, *rnorm, twork[2] = {0,0};
    mwSignedIndex *ipivr, *ipiv, *iwork, *oufact;
    char job, jobg, jobs, fact, uplo, dico, hinv, compz, trana, trans;
    char flag, def, full = 'F', full2 = 'G', lower = 'L', norm = '1';
    
    
    /* check first input argument */
    if (!mxIsNumeric(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt( "Input A must be a full double-precision real matrix." );
    }
    ndim = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *sa = mxGetDimensions(prhs[0]);
    if (ndim == 3) {
        p = sa[2];
    }
    else if (ndim == 2) {
        p = 1;
    }
    else {
        mexErrMsgTxt( "Input array A must have 2 or 3 dimensions." );
    }
    n = sa[0];
    nn = 2*n;
    if (n != sa[1]) {
        mexErrMsgTxt( "The number of rows and colums of array A must be the same size." );
    }
    Apr = mxGetData(prhs[0]);
    Alpr = mxMalloc(n*n*p*sizeof(double));
    memcpy(Alpr,Apr,n*n*p*sizeof(double));
    
    
    /* check second input argument */
    if (!mxIsNumeric(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgTxt( "Input B must be a full double-precision real matrix." );
    }
    mwSize db = mxGetNumberOfDimensions(prhs[1]);
    if (db != ndim) {
        mexErrMsgTxt( "Input array A and B must have the same number of dimensions." );
    }
    const mwSize *sb = mxGetDimensions(prhs[1]);
    m = sb[1];
    if (n != sb[0]) {
        mexErrMsgTxt( "The number of rows of arrays A and B must be the same size." );
    }
    if (ndim == 3) {
        if (p != sb[2]) {
            mexErrMsgTxt( "The number of periods must be the same for arrays A and B." );
        }
    }
    Bpr = mxGetData(prhs[1]);
    Blpr = mxMalloc(m*n*p*sizeof(double));
    memcpy(Blpr,Bpr,m*n*p*sizeof(double));
    
    
    /* check third input argument */
    if (!mxIsNumeric(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
        mexErrMsgTxt( "Input Q must be a full double-precision real matrix." );
    }
    mwSize dq = mxGetNumberOfDimensions(prhs[2]);
    if (dq != ndim) {
        mexErrMsgTxt( "Input array A and Q must have the same number of dimensions." );
    }
    const mwSize *sq = mxGetDimensions(prhs[2]);
    if (n != sq[0]) {
        mexErrMsgTxt( "The number of rows of arrays A and Q must be the same size." );
    }
    if (n != sq[1]) {
        mexErrMsgTxt( "The number of columns of arrays A and Q must be the same size." );
    }
    if (ndim == 3) {
        if (p != sq[2]) {
            mexErrMsgTxt( "The number of periods must be the same for arrays A and B." );
        }
    }
    Qpr = mxGetData(prhs[2]);   
    Qlpr = mxMalloc(n*n*p*sizeof(double));
    memcpy(Qlpr,Qpr,n*n*p*sizeof(double));
    
    
    /* check fourth input argument */
    if (nrhs > 3) {
        if (!mxIsNumeric(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
            mexErrMsgTxt( "Input R must be a full double-precision real matrix." );
        }
        mwSize dr = mxGetNumberOfDimensions(prhs[3]);
        if (dr != ndim) {
            mexErrMsgTxt( "Input array A and R must have the same number of dimensions." );
        }
        const mwSize *sr = mxGetDimensions(prhs[3]);
        if (m != sr[0]) {
            mexErrMsgTxt( "The number of rows of array R must be the same size as the number of columns of array B." );
        }
        if (m != sr[1]) {
            mexErrMsgTxt( "The number of columns of arrays R and B must be the same size." );
        }
        if (ndim == 3) {
            if (p != sr[2]) {
                mexErrMsgTxt( "The number of periods must be the same for arrays A and R." );
            }
        }
        Rpr = mxGetData(prhs[3]);        
        Rlpr = mxMalloc(m*m*p*sizeof(double));
        memcpy(Rlpr,Rpr,m*m*p*sizeof(double));
    }
    else {
        /* allocate matrix R=I if not provided */
        Rpr = mxCalloc(m*m*p,sizeof(double));
        for (k=0; k<p; k++) {
            for (i=0; i<m; i++) {
                Rpr[k*m*m+i*m+i] = 1;
            }
        }
        Rlpr = mxMalloc(m*m*p*sizeof(double));
        memcpy(Rlpr,Rpr,m*m*p*sizeof(double));
    }
    
    
    /* check fifth input argument */
    if (nrhs > 4) {
        if (!mxIsNumeric(prhs[4]) || mxIsSparse(prhs[4]) || !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])) {
            mexErrMsgTxt( "Input S must be a full double-precision real matrix." );
        }
        mwSize ds = mxGetNumberOfDimensions(prhs[4]);
        if (ds != ndim) {
            mexErrMsgTxt( "Input array A and S must have the same number of dimensions." );
        }
        const mwSize *ss = mxGetDimensions(prhs[4]);
        if (n != ss[0]) {
            mexErrMsgTxt( "The number of rows of arrays S and B must be the same size." );
        }
        if (m != ss[1]) {
            mexErrMsgTxt( "The number of columns of arrays S and B must be the same size." );
        }
        if (ndim == 3) {
            if (p != ss[2]) {
                mexErrMsgTxt( "The number of periods must be the same for arrays A and S." );
            }
        }
        Spr = mxGetData(prhs[4]);  
        Slpr = mxMalloc(n*m*p*sizeof(double));
        memcpy(Slpr,Spr,n*m*p*sizeof(double));
        jobs = 'N';
    }
    else {
        /* allocate matrix S=0 if not provided */
        Spr = mxCalloc(n*m*p,sizeof(double));
        Slpr = mxCalloc(n*m*p,sizeof(double));
        jobs = 'Z';
    }
    
    /* check sixth input argument */
    if (nrhs > 5) {
        const mwSize *se = mxGetDimensions(prhs[5]);
        if (se[0] != 0) {
            mexErrMsgTxt( "The periodic QR method does not support generalized state space matrix E." );
        }
    }
    
    /* allocate output arrays */
    dims[0] = n;
    dims[1] = n;
    dims[2] = p;
    plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    Xpr = mxGetData(plhs[0]);
    if (nlhs > 1) {
        dims[0] = m;
        dims[1] = n;
        dims[2] = p;
        plhs[1] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
        Kpr = mxGetData(plhs[1]);
    }
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(nn,1,mxCOMPLEX);
    }

    // return empty arrays if any size is zero
    if ((n == 0) || (m == 0)) {
        return;
    }
    
    
    // Perform SB02MT
    trans = 'N';
    jobg = 'G';
    fact = 'N';
    uplo = 'U';
    flag = 'M';
    def = 'D';
    Glpr = mxCalloc(n*n*p,sizeof(double));
    oufact = mxCalloc(p*2,sizeof(mwSignedIndex));
    ipivr = mxCalloc(p*m,sizeof(mwSignedIndex));
    rnorm = mxCalloc(p,sizeof(double));
    iwork = mxCalloc(m,sizeof(mwSignedIndex));
    ldwork = max(3*m,m*n);
    dwork = mxMalloc(ldwork*sizeof(double));
    for (k=0; k<p; k++) {
        sb02mt( &jobg, &jobs, &fact, &uplo, &n, &m, 
                &Alpr[k*n*n], &n, &Blpr[k*m*n], &n, &Qlpr[k*n*n], &n, &Rlpr[k*m*m], &m, 
                &Slpr[k*m*n], &n, &ipivr[k*m], &oufact[2*k], &Glpr[k*n*n], &n,
                iwork, twork, &ldmin, &info);
        // sb02mx( &jobg, &jobs, &fact, &uplo, &trans, &flag, &def, &n, &m, 
        //         &Alpr[k*n*n], &n, &Blpr[k*m*n], &n, &Qlpr[k*n*n], &n, &Rlpr[k*m*m], &m, 
        //         &Slpr[k*m*n], &n, &ipivr[k*m], &oufact[2*k], &Glpr[k*n*n], &n,
        //         iwork, twork, &ldwork, &info);
        
        if (info == 0) {
            ldwork = max(ldwork,max((mwSignedIndex)twork[0],max(3*m,m*n)));
            dwork = mxRealloc(dwork,ldwork*sizeof(double));
            sb02mt( &jobg, &jobs, &fact, &uplo, &n, &m,
                    &Alpr[k*n*n], &n, &Blpr[k*m*n], &n, &Qlpr[k*n*n], &n, &Rlpr[k*m*m], &m,
                    &Slpr[k*m*n], &n, &ipivr[k*m], &oufact[2*k], &Glpr[k*n*n], &n,
                    iwork, dwork, &ldwork, &info);
            // sb02mx( &jobg, &jobs, &fact, &uplo, &trans, &flag, &def, &n, &m,
            //         &Alpr[k*n*n], &n, &Blpr[k*m*n], &n, &Qlpr[k*n*n], &n, &Rlpr[k*m*m], &m,
            //         &Slpr[k*m*n], &n, &ipivr[k*m], &oufact[2*k], &Glpr[k*n*n], &n,
            //         iwork, dwork, &ldwork, &info);
        }
        if (info != 0) {
            mxFree(Alpr);
            mxFree(Blpr);
            mxFree(Qlpr);
            mxFree(Rlpr);
            mxFree(Slpr);
            mxFree(Glpr);
            mxFree(oufact);
            mxFree(ipivr);
            mxFree(rnorm);
            mxFree(iwork);
            mxFree(dwork);
            if (info == m+1) {
                mexErrMsgTxt("The matrix R is numerically singular.\n");
            }
            else {
                mexPrintf("SB02MT returned INFO=%d.\n",info);
                mexErrMsgTxt("SB02MT not successful.");
            }
        }
        rnorm[k] = dwork[1];
        
        // Compute the norms of the matrices Q and G.
        qnorm = qnorm + dlansy( &norm, &uplo, &n, &Qlpr[k*n*n], &n, dwork )/((double)p);
        gnorm = gnorm + dlansy( &norm, &uplo, &n, &Glpr[k*n*n], &n, dwork )/((double)p);
    }
    mxFree(dwork);
    mxFree(iwork);

    if ((qnorm > gnorm) && (gnorm > 0.0)) {
        iscl = 1; // Do scaling
    }

    
    // Perform SB02RU
    dico = 'D';
    hinv = 'D';
    trana = 'N';
    Hlpr = mxCalloc(nn*nn*p,sizeof(double));
    iwork = mxCalloc(2*n,sizeof(mwSignedIndex));
    ldwork = 6*n;
    dwork = mxMalloc(ldwork*sizeof(double));
    for (k=0; k<p; k++) {
        if (iscl) {
            // apply scaling
            dlascl( &full2, &ldzero, &ldzero, &qnorm, &gnorm, &n, &n, &Qlpr[k*n*n], &n, &info );
            dlascl( &full2, &ldzero, &ldzero, &gnorm, &qnorm, &n, &n, &Glpr[k*n*n], &n, &info );
        }
        
        sb02ru( &dico, &hinv, &trana, &uplo, &n, &Alpr[k*n*n], &n, 
                &Glpr[k*n*n], &n, &Qlpr[k*n*n], &n, &Hlpr[k*nn*nn], &nn, 
                iwork, dwork, &ldwork, &info);
        //sb02mu( &dico, &hinv, &uplo, &n, &Alpr[k*n*n], &n, 
        //        &Glpr[k*n*n], &n, &Qlpr[k*n*n], &n, &Hlpr[k*nn*nn], &nn, 
        //        iwork, dwork, &ldwork, &info);

        if (info != 0) {
            mxFree(Alpr);
            mxFree(Blpr);
            mxFree(Qlpr);
            mxFree(Rlpr);
            mxFree(Slpr);
            mxFree(Glpr);
            mxFree(Hlpr);
            mxFree(oufact);
            mxFree(ipivr);
            mxFree(rnorm);
            mxFree(iwork);
            mxFree(dwork);
            mexPrintf("SB02RU returned INFO=%d.\n",info);
            mexErrMsgTxt("SB02RU not successful.");
        }
    }
    mxFree(iwork);
    mxFree(dwork);
    mxFree(Qlpr);
    mxFree(Glpr);
    

    // Perform MB03VD
    ilo = 1;
    ihi = nn;
    ldtau = nn-1;
    tau = mxCalloc(p*ldtau,sizeof(double));
    dwork = mxMalloc(p*n*sizeof(double));
    mb03vd( &nn, &p, &ilo, &ihi, Hlpr, &nn, &nn, tau, &ldtau, dwork, &info );
    if (info != 0) {
        mxFree(Alpr);
        mxFree(Blpr);
        mxFree(Rlpr);
        mxFree(Slpr);
        mxFree(Hlpr);
        mxFree(oufact);
        mxFree(ipivr);
        mxFree(rnorm);
        mxFree(dwork);
        mxFree(tau);
        mexPrintf("MB03VD returned INFO=%d.\n",info);
        mexErrMsgTxt("MB03VD not successful.");
    }
    mxFree(dwork);
    
    // copy lower part H to Z
    Zlpr = mxCalloc(nn*nn*p,sizeof(double));
    for (k=0; k<p; k++) {
        dlacpy( &lower, &nn, &nn, &Hlpr[k*nn*nn], &nn, &Zlpr[k*nn*nn], &nn );
    }
    
    
    // Perform MB03VY
    mb03vy( &nn, &p, &ilo, &ihi, Zlpr, &nn, &nn, tau, &ldtau, twork, &ldmin, &info );
    if (info == 0) {
        ldwork = max(ldwork,max((mwSignedIndex)twork[0],n));
        dwork = mxMalloc(ldwork*sizeof(double));
        mb03vy( &nn, &p, &ilo, &ihi, Zlpr, &nn, &nn, tau, &ldtau, dwork, &ldwork, &info );
    }
    if (info != 0) {
        mxFree(Alpr);
        mxFree(Blpr);
        mxFree(Rlpr);
        mxFree(Slpr);
        mxFree(Hlpr);
        mxFree(Zlpr);
        mxFree(oufact);
        mxFree(ipivr);
        mxFree(rnorm);
        mxFree(dwork);
        mxFree(tau);
        mexPrintf("MB03VY returned INFO=%d.\n",info);
        mexErrMsgTxt("MB03VY not successful.");
    }
    mxFree(tau);
    mxFree(dwork);
    
    
    // Perform MB03WD
    job = 'S';
    compz = 'V';
    wr = mxCalloc(nn,sizeof(double));
    wi = mxCalloc(nn,sizeof(double));
    ldwork = ihi-ilo+p-1;
    dwork = mxMalloc(ldwork*sizeof(double));
    mb03wd( &job, &compz, &nn, &p, &ilo, &ihi, &ilo, &ihi, Hlpr, &nn, &nn, Zlpr, &nn, &nn, wr, wi, dwork, &ldwork, &info );
    if (info != 0) {
        mxFree(Alpr);
        mxFree(Blpr);
        mxFree(Rlpr);
        mxFree(Slpr);
        mxFree(Hlpr);
        mxFree(Zlpr);
        mxFree(oufact);
        mxFree(ipivr);
        mxFree(rnorm);
        mxFree(dwork);
        if (info > 0) {
            mexPrintf("MB03WD failed to converge INFO=%d.\n",info);
        }
        else {
            mexPrintf("MB03WD returned INFO=%d.\n",info);
        }
        mexErrMsgTxt("MB03WD not successful.");
    }
    mxFree(dwork);
    mxFree(Hlpr);
    
    if (nlhs > 2) {
        #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble *pc = mxGetComplexDoubles(plhs[2]);
        for (i=0; i<nn; i++) {
            pc[i].real = wr[i];
            pc[i].imag = wi[i];
        }
        #else
        double *pr = mxGetData(nlhs[2]);
        memcpy(pr,wr,nn*sizeof(double));
        double *pi = mxGetImagData(nlhs[2]);
        memcpy(pi,wi,nn*sizeof(double));
        #endif
    }
    mxFree(wr);
    mxFree(wi);
        
    // copy part Z to T and Z
    Tlpr = mxMalloc(p*n*n*sizeof(double));
    for (k=0; k<p; k++) {
        dlacpy( &full, &n, &n, &Zlpr[k*nn*nn], &nn, &Tlpr[k*n*n], &n );
        dlacpy( &full, &n, &n, &Zlpr[k*nn*nn+n], &nn, &Xpr[k*n*n], &n );
    }
    mxFree(Zlpr);

    
    // Perform MB02VD
    ipiv = mxCalloc(n,sizeof(mwSignedIndex));
    for (k=0; k<p; k++) {
        mb02vd( &trans, &n, &n, &Tlpr[k*n*n], &n, ipiv, &Xpr[k*n*n], &n, &info );
        if (info != 0) {
            mxFree(Alpr);
            mxFree(Blpr);
            mxFree(Rlpr);
            mxFree(Slpr);
            mxFree(Tlpr);
            mxFree(oufact);
            mxFree(ipivr);
            mxFree(rnorm);
            mxFree(ipiv);
            mexPrintf("MB02VD returned INFO=%d.\n",info);
            mexErrMsgTxt("MB02VD not successful.");
        }
    }
    mxFree(ipiv);
    mxFree(Tlpr);  
    
    // Make sure the solution matrix X is symmetric.
    for (k=0; k<p; k++) {
        for (i=0; i<n; i++) {
            nre = n-i;
            daxpy( &nre, &one,  &Xpr[k*n*n], &n, &Xpr[k*n*n], &ldone);
            dscal( &nre, &half, &Xpr[k*n*n], &ldone);
            dcopy( &nre, &Xpr[k*n*n], &ldone, &Xpr[k*n*n], &n);
        }
        if (iscl) {
            // undo scaling
            dlascl( &full2, &ldzero, &ldzero, &gnorm, &qnorm, &n, &n, &Xpr[k*n*n], &n, &info );
        }
    }

    
    // perform SB02ND
    if (nlhs > 1) {
        Xlpr = mxCalloc(n*n,sizeof(double));
        iwork = mxCalloc(m,sizeof(mwSignedIndex));
        ldwork = max(m*n,max(n+3*m+2,4*n+1));
        dwork = mxMalloc(ldwork*sizeof(double));
        for (k=0; k<p; k++) {
            switch (oufact[2*k]) {
                case 0:
                    fact = 'N';
                    break;
                case 1:
                    fact = 'C';
                    break;
                case 2:
                    fact = 'U';
                    mexErrMsgTxt("SG02ND fact=U is not supported.\n");
                    break;
            }
            
            if (fact == 'C')
            {
                dlacpy( &full, &n, &m, &Bpr[k*n*m], &n, &Blpr[k*n*m], &n );
                if (k<(p-1)) {
                    dlacpy( &uplo, &n, &n, &Xpr[(k+1)*n*n], &n, Xlpr, &n );
                }
                else {
                    dlacpy( &uplo, &n, &n, &Xpr[0], &n, Xlpr, &n );
                }
                
                sb02nd( &dico, &fact, &uplo, &jobs, &n, &m, &m,
                        &Apr[k*n*n], &n, &Blpr[k*m*n], &n,
                        &Rlpr[k*m*m], &m, &ipivr[k*m], &Spr[k*m*n], &n,
                        Xlpr, &n, &rnorm[k], &Kpr[k*m*n], &m,
                        &oufact[2*k], iwork, twork, &ldmin, &info);
                
                if (info == 0) {
                    ldwork = max(ldwork,max((mwSignedIndex)twork[0],m*n));
                    dwork = mxRealloc(dwork,ldwork*sizeof(double));
                    sb02nd( &dico, &fact, &uplo, &jobs, &n, &m, &m,
                            &Apr[k*n*n], &n, &Blpr[k*m*n], &n,
                            &Rlpr[k*m*m], &m, &ipivr[k*m], &Spr[k*m*n], &n,
                            Xlpr, &n, &rnorm[k], &Kpr[k*m*n], &m,
                            &oufact[2*k], iwork, dwork, &ldwork, &info);
                }
            }
            if ((fact == 'N') || (info > m+1)) {
                // call the routine again with the unfactored matrix R.
                if (fact != 'N') {
                    mexPrintf("SB02ND again with the unfactored matrix R (p=%d).\n",k);
                    fact = 'N';
                }
                dlacpy( &full, &n, &m, &Bpr[k*n*m], &n, &Blpr[k*n*m], &n );
                dlacpy( &full, &m, &m, &Rpr[k*m*m], &m, &Rlpr[k*m*m], &m );
                if (k<(p-1)) {
                    dlacpy( &full, &n, &n, &Xpr[(k+1)*n*n], &n, Xlpr, &n );
                }
                else {
                    dlacpy( &full, &n, &n, &Xpr[0], &n, Xlpr, &n );
                }
                
                sb02nd( &dico, &fact, &uplo, &jobs, &n, &m, &m,
                        &Apr[k*n*n], &n, &Blpr[k*m*n], &n,
                        &Rlpr[k*m*m], &m, &ipivr[k*m], &Spr[k*m*n], &n,
                        Xlpr, &n, &rnorm[k], &Kpr[k*m*n], &m,
                        &oufact[2*k], iwork, twork, &ldmin, &info);
                
                if (info == 0) {
                    ldwork = max(ldwork,max((mwSignedIndex)twork[0],max(n+3*m+2,4*n+1)));
                    dwork = mxRealloc(dwork,ldwork*sizeof(double));
                    sb02nd( &dico, &fact, &uplo, &jobs, &n, &m, &m,
                            &Apr[k*n*n], &n, &Blpr[k*m*n], &n,
                            &Rlpr[k*m*m], &m, &ipivr[k*m], &Spr[k*m*n], &n,
                            Xlpr, &n, &rnorm[k], &Kpr[k*m*n], &m,
                            &oufact[2*k], iwork, dwork, &ldwork, &info);
                }
            }
            if (info != 0) {
                mxFree(Alpr);
                mxFree(Blpr);
                mxFree(Rlpr);
                mxFree(Slpr);
                mxFree(Xlpr);
                mxFree(oufact);
                mxFree(ipivr);
                mxFree(rnorm);
                mxFree(iwork);
                mxFree(dwork);
                mexPrintf("SB02ND returned INFO=%d.\n",info);
                mexErrMsgTxt("SB02ND not successful.");
            }
        }
        mxFree(Xlpr);
        mxFree(dwork);
        mxFree(iwork);
    }
    
    mxFree(Alpr);
    mxFree(Blpr);
    mxFree(Rlpr);
    mxFree(Slpr);
    mxFree(oufact);
    mxFree(ipivr);
    mxFree(rnorm);
}


void dprex_complexqz(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize ndim, dims[3] = {1, 1, 1};
    mwSignedIndex m, n, p, nn, ilo, ihi, nre, info = 0; 
    mwSignedIndex ldwork, lzwork, ldmin = -1, ldzero = 0, ldone = 1;
    mwSignedIndex i, k;
    double one = 1.0, half = 0.5, tol = 0.0; 
    double *Apr, *Bpr, *Qpr, *Rpr, *Spr, *Epr, *Xpr, *Kpr, *alpha, *beta;
    double *Blpr, *Rlpr, *AFlpr, *BFlpr, *Hlpr, *Zlpr, *Xlpr, *Tlpr;
    double *dwork, *zwork, *rnorm, twork[2] = {0,0};
    mwSignedIndex *Jlpr, *ipiv, *iwork, *oufact, *scal;
    char type, compq, job, jobb, jobs, jobe, jobx, fact, uplo, dico, trans;
    char full = 'F';
    
    
    /* check first input argument */
    if (!mxIsNumeric(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt( "Input A must be a full double-precision real matrix." );
    }
    ndim = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *sa = mxGetDimensions(prhs[0]);
    if (ndim == 3) {
        p = sa[2];
    }
    else if (ndim == 2) {
        p = 1;
    }
    else {
        mexErrMsgTxt( "Input array A must have 2 or 3 dimensions." );
    }
    n = sa[0];
    if (n != sa[1]) {
        mexErrMsgTxt( "The number of rows and colums of array A must be the same size." );
    }
    Apr = mxGetData(prhs[0]);
    
    
    /* check second input argument */
    if (!mxIsNumeric(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgTxt( "Input B must be a full double-precision real matrix." );
    }
    mwSize db = mxGetNumberOfDimensions(prhs[1]);
    if (db != ndim) {
        mexErrMsgTxt( "Input array A and B must have the same number of dimensions." );
    }
    const mwSize *sb = mxGetDimensions(prhs[1]);
    m = sb[1];
    nn = 2*n+m;
    if (n != sb[0]) {
        mexErrMsgTxt( "The number of rows of arrays A and B must be the same size." );
    }
    if (ndim == 3) {
        if (p != sb[2]) {
            mexErrMsgTxt( "The number of periods must be the same for arrays A and B." );
        }
    }
    Bpr = mxGetData(prhs[1]);
    Blpr = mxMalloc(m*n*p*sizeof(double));
    memcpy(Blpr,Bpr,m*n*p*sizeof(double));
    
    
    /* check third input argument */
    if (!mxIsNumeric(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
        mexErrMsgTxt( "Input Q must be a full double-precision real matrix." );
    }
    mwSize dq = mxGetNumberOfDimensions(prhs[2]);
    if (dq != ndim) {
        mexErrMsgTxt( "Input array A and Q must have the same number of dimensions." );
    }
    const mwSize *sq = mxGetDimensions(prhs[2]);
    if (n != sq[0]) {
        mexErrMsgTxt( "The number of rows of arrays A and Q must be the same size." );
    }
    if (n != sq[1]) {
        mexErrMsgTxt( "The number of columns of arrays A and Q must be the same size." );
    }
    if (ndim == 3) {
        if (p != sq[2]) {
            mexErrMsgTxt( "The number of periods must be the same for arrays A and B." );
        }
    }
    Qpr = mxGetData(prhs[2]);
    
    
    /* check fourth input argument */
    if (nrhs > 3) {
        if (!mxIsNumeric(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
            mexErrMsgTxt( "Input R must be a full double-precision real matrix." );
        }
        mwSize dr = mxGetNumberOfDimensions(prhs[3]);
        if (dr != ndim) {
            mexErrMsgTxt( "Input array A and R must have the same number of dimensions." );
        }
        const mwSize *sr = mxGetDimensions(prhs[3]);
        if (m != sr[0]) {
            mexErrMsgTxt( "The number of rows of array R must be the same size as the number of columns of array B." );
        }
        if (m != sr[1]) {
            mexErrMsgTxt( "The number of columns of arrays R and B must be the same size." );
        }
        if (ndim == 3) {
            if (p != sr[2]) {
                mexErrMsgTxt( "The number of periods must be the same for arrays A and R." );
            }
        }
        Rpr = mxGetData(prhs[3]);        
        Rlpr = mxMalloc(m*m*p*sizeof(double));
        memcpy(Rlpr,Rpr,m*m*p*sizeof(double));
    }
    else {
        /* allocate matrix R=I if not provided */
        Rpr = mxCalloc(m*m*p,sizeof(double));
        for (k=0; k<p; k++) {
            for (i=0; i<m; i++) {
                Rpr[k*m*m+i*m+i] = 1;
            }
        }
        Rlpr = mxMalloc(m*m*p*sizeof(double));
        memcpy(Rlpr,Rpr,m*m*p*sizeof(double));
    }
    
    
    /* check fifth input argument */
    if (nrhs > 4) {
        if (!mxIsNumeric(prhs[4]) || mxIsSparse(prhs[4]) || !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])) {
            mexErrMsgTxt( "Input S must be a full double-precision real matrix." );
        }
        mwSize ds = mxGetNumberOfDimensions(prhs[4]);
        if (ds != ndim) {
            mexErrMsgTxt( "Input array A and S must have the same number of dimensions." );
        }
        const mwSize *ss = mxGetDimensions(prhs[4]);
        if (n != ss[0]) {
            mexErrMsgTxt( "The number of rows of arrays S and B must be the same size." );
        }
        if (m != ss[1]) {
            mexErrMsgTxt( "The number of columns of arrays S and B must be the same size." );
        }
        if (ndim == 3) {
            if (p != ss[2]) {
                mexErrMsgTxt( "The number of periods must be the same for arrays A and S." );
            }
        }
        Spr = mxGetData(prhs[4]);       
        jobs = 'N';
    }
    else {
        /* allocate matrix S=0 if not provided */
        Spr = mxCalloc(n*m*p,sizeof(double));
        jobs = 'Z';
    }
    
    
    /* check sixth input argument */
    if (nrhs > 5) {
        if (!mxIsNumeric(prhs[5]) || mxIsSparse(prhs[5]) || !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5])) {
            mexErrMsgTxt( "Input E must be a full double-precision real matrix." );
        }
        mwSize de = mxGetNumberOfDimensions(prhs[5]);
        if (de != ndim) {
            mexErrMsgTxt( "Input array A and E must have the same number of dimensions." );
        }
        const mwSize *se = mxGetDimensions(prhs[5]);
        if (n != se[0]) {
            mexErrMsgTxt( "The number of rows of arrays E and A must be the same size." );
        }
        if (n != se[1]) {
            mexErrMsgTxt( "The number of columns of arrays E and A must be the same size." );
        }
        if (ndim == 3) {
            if (p != se[2]) {
                mexErrMsgTxt( "The number of periods must be the same for arrays A and E." );
            }
        }
        Epr = mxGetData(prhs[3]);        
        jobe = 'N';
    }
    else {
        /* allocate matrix E=I if not provided */
        Epr = mxCalloc(m*m*p,sizeof(double));
        for (k=0; k<p; k++) {
            for (i=0; i<n; i++) {
                Epr[k*n*n+i*n+i] = 1;
            }
        }
        jobe = 'I';
    }
    
    
    /* allocate output arrays */
    dims[0] = n;
    dims[1] = n;
    dims[2] = p;
    plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    Xpr = mxGetData(plhs[0]);
    if (nlhs > 1) {
        dims[0] = m;
        dims[1] = n;
        dims[2] = p;
        plhs[1] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
        Kpr = mxGetData(plhs[1]);
    }
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(nn,1,mxCOMPLEX);
    }

    // return empty arrays if any size is zero
    if ((n == 0) || (m == 0)) {
        return;
    }
    
    
    // Perform SB02OY
    type = 'O';
    dico = 'D';
    jobb = 'B';
    fact = 'N';
    uplo = 'U';
    AFlpr = mxCalloc(nn*nn*p,sizeof(double));
    BFlpr = mxCalloc(nn*nn*p,sizeof(double));
    rnorm = mxCalloc(p,sizeof(double));
    iwork = mxCalloc(m,sizeof(mwSignedIndex));
    ldwork = max(3*m,2*n+m);
    dwork = mxMalloc(ldwork*sizeof(double));
    for (k=0; k<p; k++) {
        sb02oy( &type, &dico, &jobb, &fact, &uplo, &jobs, &jobe, &n, &m, &m,
                &Apr[k*n*n], &n, &Bpr[k*m*n], &n, &Qpr[k*n*n], &n, &Rpr[k*m*m], &m,
                &Spr[k*m*n], &n, &Epr[k*m*n], &n, &AFlpr[k*nn*nn], &nn, &BFlpr[k*nn*nn], &nn,
                &tol, iwork, dwork, &ldwork, &info );
        if (info != 0) {
            mxFree(Blpr);
            mxFree(Rlpr);
            mxFree(AFlpr);
            mxFree(BFlpr);
            mxFree(iwork);
            mxFree(dwork);
            mxFree(rnorm);
            if (info == 1) {
                mexErrMsgTxt("The computed extended matrix pencil is singular,\n possibly due to rounding errors.\n");
            }
            else {
                mexPrintf("SB02OY returned INFO=%d.\n",info);
                mexErrMsgTxt("SB02OY not successful.");
            }
        }
        rnorm[k] = dwork[1];
    }
    mxFree(dwork);
    mxFree(iwork);
    
    
    // Copy AF and BF results to complex H array
    Hlpr = mxMalloc(nn*nn*p*2*sizeof(double));
    Zlpr = mxMalloc(nn*nn*p*2*sizeof(double));
    Jlpr = mxMalloc(2*p*sizeof(mwSignedIndex));
    for (k=0; k<p; k++) {
        dlacpy( &full, &nn, &nn, &AFlpr[nn*nn*(p-k-1)], &nn, &Hlpr[nn*nn*2*k], &nn);
        Jlpr[2*k] = 1;
        dlacpy( &full, &nn, &nn, &BFlpr[nn*nn*((p-k) % p)], &nn, &Hlpr[nn*nn*(2*k+1)], &nn);
        Jlpr[2*k+1] = -1;
    }
    mxFree(AFlpr);
    mxFree(BFlpr);   
    
    
    // Perform ZPGHRD
    compq = 'I';
    ilo = 1;
    ihi = nn;
    ldwork = max(p,nn);
    dwork = mxMalloc(ldwork*sizeof(double));
    lzwork = max(p,2*nn);
    zwork = mxMalloc(2*ldwork*sizeof(double));  
    zpghrd( &compq, &p, &nn, &ilo, &ihi, Jlpr, (doublecomplex *)Hlpr, 
            &nn, &nn, (doublecomplex *)Zlpr, &nn, &nn,
            dwork, &ldwork, (doublecomplex *)zwork, &lzwork, &info );
    if (info != 0) {
        mxFree(Blpr);
        mxFree(Rlpr);
        mxFree(Jlpr);
        mxFree(Hlpr);
        mxFree(Zlpr);
        mxFree(dwork);
        mxFree(zwork);
        mxFree(rnorm);
        mexPrintf("ZPGHRD returned INFO=%d.\n",info);
        mexErrMsgTxt("ZPGHRD not successful.");
    }
    //mxFree(dwork);
    //mxFree(zwork);
    
    
    // Perform ZPGEQZ
    job = 'S';
    compq = 'V';
    alpha = mxCalloc(2*nn,sizeof(double));
    beta = mxCalloc(2*nn,sizeof(double));
    scal = mxCalloc(nn,sizeof(mwSignedIndex));
    zpgeqz( &job, &compq, &p, &nn, &ilo, &ihi, Jlpr, (doublecomplex *)Hlpr, 
            &nn, &nn, (doublecomplex *)alpha, (doublecomplex *)beta, scal, 
            (doublecomplex *)Zlpr, &nn, &nn,
            dwork, &ldwork, (doublecomplex *)zwork, &lzwork, &info );
    if (info != 0) {
        mxFree(Blpr);
        mxFree(Rlpr);
        mxFree(Jlpr);
        mxFree(Hlpr);
        mxFree(Zlpr);
        mxFree(alpha);
        mxFree(beta);
        mxFree(scal);
        mxFree(dwork);
        mxFree(zwork);
        mxFree(rnorm);
        mexPrintf("ZPGEQZ returned INFO=%d.\n",info);
        mexErrMsgTxt("ZPGEQZ not successful.");
    }
    mxFree(Jlpr);
    mxFree(Hlpr);
    mxFree(dwork);
    mxFree(zwork);
    
    if (nlhs > 2) {
        #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble *pc = mxGetComplexDoubles(plhs[2]);
        for (i=0; i<nn; i++) {
            pc[i].real = scalbn(alpha[2*i]/beta[2*i],scal[i]);
            pc[i].imag = scalbn(alpha[2*i+1]/beta[2*i+1],scal[i]);
        }
        #else
        double *pr = mxGetData(nlhs[2]);
        double *pi = mxGetImagData(nlhs[2]);
        for (i=0; i<nn; i++) {
            pr[i] = scalbn(alpha[2*i]/beta[2*i],scale[i]);
            pi[i] = scalbn(alpha[2*i+1]/beta[2*i+1],scale[i]);
        }
        #endif
    }
    mxFree(alpha);
    mxFree(beta); 
    mxFree(scal);
    

    // copy part Z to T and Z
    Tlpr = mxMalloc(p*n*n*sizeof(double));
    for (k=0; k<p; k++) {
        dlacpy( &full, &n, &n, &Zlpr[nn*nn*2*k], &nn, &Tlpr[n*n*(p-k-1)], &n );
        dlacpy( &full, &n, &n, &Zlpr[nn*nn*2*k+n], &nn, &Xpr[n*n*(p-k-1)], &n );
    }
    mxFree(Zlpr);

    
    // Perform MB02VD
    trans = 'N';
    ipiv = mxCalloc(n,sizeof(mwSignedIndex));
    for (k=0; k<p; k++) {
        mb02vd( &trans, &n, &n, &Tlpr[k*n*n], &n, ipiv, &Xpr[k*n*n], &n, &info );
        if (info != 0) {
            mxFree(Blpr);
            mxFree(Rlpr);
            mxFree(Tlpr);
            mxFree(ipiv);
            mxFree(rnorm);
            mexPrintf("MB02VD returned INFO=%d.\n",info);
            mexErrMsgTxt("MB02VD not successful.");
        }
    }
    mxFree(ipiv);
    mxFree(Tlpr);  
    
    // Make sure the solution matrix X is symmetric.
    for (k=0; k<p; k++) {
        for (i=0; i<n; i++) {
            nre = n-i;
            daxpy( &nre, &one,  &Xpr[k*n*n], &n, &Xpr[k*n*n], &ldone );
            dscal( &nre, &half, &Xpr[k*n*n], &ldone );
            dcopy( &nre, &Xpr[k*n*n], &ldone, &Xpr[k*n*n], &n );
        }
    }
    
    
    // perform SG02ND
    if (nlhs > 1) {
        job = 'K';
        jobx = 'N';
        Xlpr = mxCalloc(n*n,sizeof(double));
        ipiv = mxCalloc(n,sizeof(mwSignedIndex));
        oufact = mxCalloc(p*2,sizeof(mwSignedIndex));
        iwork = mxCalloc(m,sizeof(mwSignedIndex));
        ldwork = max(m*n,max(n+3*m+2,4*n+1));
        dwork = mxMalloc(ldwork*sizeof(double));
        for (k=0; k<p; k++) {
            if (k<(p-1)) {
                dlacpy( &full, &n, &n, &Xpr[(k+1)*n*n], &n, Xlpr, &n );
            }
            else {
                dlacpy( &full, &n, &n, &Xpr[0], &n, Xlpr, &n );
            }
            
            sg02nd( &dico, &jobe, &job, &jobx, &fact, &uplo, &jobs, &trans,
                    &n, &m, &m, &Apr[k*n*n], &n, &Epr[k*n*n], &n,
                    &Blpr[k*m*n], &n, &Rlpr[k*m*m], &m, ipiv, &Spr[k*m*n], &n,
                    Xlpr, &n, &rnorm[k], &Kpr[k*m*n], &m, NULL, &ldone, NULL, &ldone,
                    &oufact[2*k], iwork, twork, &ldmin, &info );
            
            if (info == 0) {
                ldwork = max(ldwork,max((mwSignedIndex)twork[0],max(n+3*m+2,4*n+1)));
                dwork = mxRealloc(dwork,ldwork*sizeof(double));
                sg02nd( &dico, &jobe, &job, &jobx, &fact, &uplo, &jobs, &trans,
                        &n, &m, &m, &Apr[k*n*n], &n, &Epr[k*n*n], &n, 
                        &Blpr[k*m*n], &n, &Rlpr[k*m*m], &m, ipiv, &Spr[k*m*n], &n,
                        Xlpr, &n, &rnorm[k], &Kpr[k*m*n], &m, NULL, &ldone, NULL, &ldone,
                        &oufact[2*k], iwork, dwork, &ldwork, &info );
            }
            
            if (info != 0) {
                mxFree(Blpr);
                mxFree(Rlpr);
                mxFree(Xlpr);
                mxFree(ipiv);
                mxFree(oufact);
                mxFree(iwork);
                mxFree(dwork);
                mxFree(rnorm);
                mexPrintf("SG02ND returned INFO=%d.\n",info);
                mexErrMsgTxt("SG02ND not successful.");
            }
        }
        mxFree(Xlpr);
        mxFree(ipiv);
        mxFree(oufact);
        mxFree(dwork);
        mxFree(iwork);
    }
    mxFree(Blpr);
    mxFree(Rlpr);
    mxFree(rnorm);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    if (nrhs < 3) {
        mexErrMsgTxt("DPREX requires at least three input arguments.");
    }
    if (nlhs < 1) {
        mexErrMsgTxt("DPREX requires at least one output arguments.");
    }
    if (nrhs > 7) {
        mexErrMsgTxt("Too many input arguments.");
    }
    if (nlhs > 3) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    /* select method */
    if ((nlhs == 7) && (strncmp(mxArrayToString(prhs[5]),"periodicqr",9) == 0)) {
        dprex_periodicqr(nlhs, plhs, nrhs, prhs);
    }
    else {
        dprex_complexqz(nlhs, plhs, nrhs, prhs);
    }

}
