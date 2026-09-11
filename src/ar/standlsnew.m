function res=standlsnew(rniz,x)
% T=size(rniz,1);
p=size(x,1)-1;
M=tmatrix(rniz,p);
sa2=sigmaanew(rniz,x);
C=inv(M'*M)*sa2;
varijanse=diag(C);
res=sqrt(varijanse);