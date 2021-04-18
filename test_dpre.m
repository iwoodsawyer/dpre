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
[X,K] = dpre(A,B,Q); % periodic state feedback
toc
 
for i=1:p
At(:,:,p-i+1) = sys.a(:,:,i)';
Bt(:,:,p-i+1) = sys.c(:,:,i)';
R(:,:,p-i+1) = sys.b(:,:,i)*sys.b(:,:,i)';
end
[X,Lt] = dpre(At,Bt,R); % periodic state observer
for i=1:p
L(:,:,i) = Lt(:,:,p-i+1)';
end