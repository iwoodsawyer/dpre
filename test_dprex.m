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
At(:,:,i) = sys.a(:,:,i)';
Bt(:,:,i) = sys.c(:,:,i)';
R(:,:,i) = sys.b(:,:,i)*sys.b(:,:,i)';
end
[X,Lt] = dprex(At,Bt,R); % periodic state observer
for i=1:p
L(:,:,i) = Lt(:,:,i)';
end