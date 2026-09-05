function res=acf(rniz,l)
r=mean(rniz);
% T=size(rniz,1);
d=sum((rniz-r).^2);
r1=rniz(l+1:end);
r2=rniz(1:end-l);
g=(r1-r)'*(r2-r);
res=g/d;