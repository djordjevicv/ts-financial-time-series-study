function res=kpredictionsAR(rniz,x,k) 
p=size(x,1)-1;
T=size(rniz,1);
c=x(p+1);
if p>0
tt=rniz(T-p+1:T,1);
koef=x(1:p);
nizpred=[];
for i=1:k
    rtk=koef'*tt+c;
    nizpred=[nizpred;rtk];
    tt=[tt;rtk];
    tt(1)=[];
end

else
nizpred=c*ones(k,1);
end
res=nizpred;