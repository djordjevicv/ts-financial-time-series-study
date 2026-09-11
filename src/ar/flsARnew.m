function res=flsARnew(x,rniz)
% T=size(rniz,1);
p=size(x,1)-1;
pred=Fmnew(x,rniz);
y=rniz(p+1:end,1);
res=0.5*norm(pred-y)^2;
% plot(pred,'g')
% hold on
% plot(y,'r')
% plot(y-pred)