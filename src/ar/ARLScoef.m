function res=ARLScoef(r,p)
% x0=1*ones(p+1,1)/(p+1);
% % res=fminunc(@(x)flsARnew(x,r),x0);
% options = optimoptions('fminunc','OptimalityTolerance',10^(-6),'MaxIterations',1000,'StepTolerance',10^(-20),'MaxFunctionEvaluations',10000);
% [x,fval,exitflag,output] = fminunc(@(x)flsARnew(x,r),x0,options)

M=tmatrix(r,p);
y=r(p+1:end,1);
x=linsolve(M'*M,M'*y);



res=x;