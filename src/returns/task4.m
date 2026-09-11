function res=task4(P0,r)
T=size(r,1);
res=[P0];
% Pold=P0;
for i=1:T
%     Pt=Pold*exp(r(i));
    Pt=res(end)*exp(r(i));
    res=[res;Pt];
%     Pold=Pt;
end