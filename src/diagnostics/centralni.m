function res=centralni(niz,k)
T=size(niz,1);
m=mean(niz);
% res=sum((niz-m*ones(T,1)).^k)/(T-1);
res=sum((niz-m).^k)/(T-1);
