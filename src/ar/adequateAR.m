function [res at]=adequateAR(rniz,x)
T=size(rniz,1);
p=size(x,1)-1;
pred=Fmnew(x,rniz);
y=rniz(p+1:end,1);


% t=tmatrica(rniz,p);
% pred=Fm(x,t);
% y=rniz(p+1:T,1);
at=y-pred;
m=10;
% m=ceil(log(T));
g=sum(abs(sign(x)));
res=ljungbox3(at,m,g);


