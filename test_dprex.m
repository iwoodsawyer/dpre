clear A B At Bt Q R L Lt

ns = 4;
nu = 3;
ny = 5;
p = 36;

sys = c2d(rss(ns,ny,nu,p),1/p);
for i=1:p
A(:,:,i) = sys.a(:,:,i);
B(:,:,i) = sys.b(:,:,i);
Q(:,:,i) = sys.c(:,:,i)'*sys.c(:,:,i);
end
tic
[X1,K1,E1] = dprex(A,B,Q,[],[],[],'periodicqr'); % periodic state feedback
toc
tic
[X2,K2,E2] = dprex(A,B,Q,[],[],[],'complexqz'); % periodic state feedback
toc

for i=1:p
At(:,:,p-i+1) = sys.a(:,:,i)';
Bt(:,:,p-i+1) = sys.c(:,:,i)';
R(:,:,p-i+1) = sys.b(:,:,i)*sys.b(:,:,i)';
end
tic
[X1t,L1t,Et1]  = dprex(At,Bt,R,[],[],[],'periodicqr'); % periodic state observer
toc
tic
[X2t,L2t,Et2] = dprex(At,Bt,R,[],[],[],'complexqz'); % periodic state observer
toc
for i=1:p
L(:,:,i) = L1t(:,:,p-i+1)';
end