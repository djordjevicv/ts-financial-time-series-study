function res=ARestimate(r,p)
x=ARLScoef(r,p);
se=standlsnew(r,x);
t=x./se;
% z=1.96;
% -norminv(0.05/2,0,1), written without the Statistics Toolbox.
z=sqrt(2)*erfcinv(0.05);
test=abs(t)>z;
tildex=x.*test;
res=[x';se';test';tildex'];