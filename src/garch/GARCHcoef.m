function res=GARCHcoef(a,m,s)
%Feasible starting point
x0=[zeros(m,1);0.1;zeros(s,1)];
%Linear inequality constraints
eps=10^(-6);
A=[ones(1,m) 0 ones(1,s)];
b=1-eps;
%Box ogranicenja
lb=[zeros(m,1);eps;zeros(s,1)];

[x,f]=fmincon(@(x)fmlGARCH(a,x,m,s),x0,A,b,[],[],lb,[]);
res=x;