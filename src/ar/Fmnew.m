function res=Fmnew(x,r)
p=size(x,1)-1;
M=tmatrix(r,p);
res=M*x;