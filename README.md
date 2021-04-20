These functions solve the Discrete-time Periodic Riccati Equation (DPRE) for periodic LQ state-feedback design. The functions compute the unique stabilizing solution X{k} of the discrete-time periodic Riccati equation and also returns the gain matrix K{k} in the state-feedback u{k} = -K{k}x{k}, where k = 1:P.

The m-file "dpre" solves the discrete-time periodic optimal control problem by a cyclic QZ or a Newton backward iteration method. These are not the fastest methods available, but work quite well. 

The mex-file "dprex" solves the discrete-time periodic optimal control problem by a periodic QR (using functions from matlab's internal slicot library) or a complex periodic QC method (using converted fortran to c code from the pqzschur library). The mex-file implementation is much faster, but requires compilation of the mex file which can be done by running make_dprex.m. 

[![View  Discrete-time Periodic Riccati Equation (DPRE) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/21379-discrete-time-periodic-riccati-equation-dpre)
