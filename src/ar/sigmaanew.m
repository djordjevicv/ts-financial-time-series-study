function res=sigmaanew(rniz,x)
T=size(rniz,1);
p=size(x,1)-1;
res=flsARnew(x,rniz)*2/(T-2*p-1);