function [X,K] = dpre(A,B,Q,R,S,E,method)
%DPREX Discrete-time Periodic Riccati Equation (mex)
%  [X,K]=DPRE(A,B,Q,R,S,E) computes the unique stabilizing solution X{k},
%  k = 1:P, of the discrete-time periodic Riccati equation
%
%   E{k}'X{k}E{k} = A{k}'X{k+1}A{k} - (A{k}'X{k+1}B{k} + S{k})*...
%                 (B{k}'X{k+1}B{k} + R{k})\(A{k}'X{k+1}B{k} + S{k})' + Q{k}
%
%  When omitted, R, S and E are set to the default values R{k}=I, S{k}=0,
%  and E{k}=I. Beside the solution X{k}, DPRE also returns the gain matrix
%
%   K{k} = (B{k}'X{k+1}B{k} + R{k})\(B{k}'X{k+1}A{k} + S{k}'),
%
%  All input matrices have to be multidimensional arrays, like matrix 
%  A(N,N,P) and B(N,R,P). Output matrices are also multidimensional arrays
%  in the size of X(N,N,P) and K(R,N,P).
%
%  [X,K]=DPRE(A,B,Q,R,S,E,METHOD) specifies the method to be used.
%    'periodicqr' - Periodic QR algorithm method (mex using slicot)
%    'complexqz'  - Periodic complex QZ algorithm method (mex using pqzschur)
%
%  See also DPRE, DARE.

%  The mex-file code uses calls the periodic Hessenberg function MB03VD/MB03VY 
%  and the periodic Schur function MB03WD from Matlab's internal slicot library.
%  or the complex periodic Schur functions from PQZSCHUR.

 




