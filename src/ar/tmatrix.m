function res=tmatrix(r,p)
res=[];
T=size(r,1);
for i=1:p
    res=[res r(i:T-p+i-1)];
end
res=[res ones(T-p,1)];