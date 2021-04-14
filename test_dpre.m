clear
sys = drss(50,11,3,36);
for i=1:36
A(:,:,i) = sys.a(:,:,i);
B(:,:,i) = sys.b(:,:,i);
Q(:,:,i) = sys.c(:,:,i)'*sys.c(:,:,i);
end
tic
[X,K] = dpre(A,B,Q); % periodic state feedback
toc
X(:,:,end)

tic
[X,K] = dpreb(A,B,Q); % periodic state feedback
toc
X(:,:,end)
 
for i=1:36
At(:,:,i) = sys.a(:,:,i)';
Bt(:,:,i) = sys.c(:,:,i)';
R(:,:,i) = sys.b(:,:,i)*sys.b(:,:,i)';
end
[X,Lt] = dpre(At,Bt,R); % periodic state observer
for i=1:36
L(:,:,i) = Lt(:,:,i)';
end