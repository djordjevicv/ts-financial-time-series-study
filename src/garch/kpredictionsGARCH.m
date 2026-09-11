function res=kpredictionsGARCH(a,x,m,s,k)
%x=(alfam,...,alfa1,alfa0,betas,...,beta1)

c=max(m,s);
alfa0=x(m+1);
T=size(a,1);
n=size(x,1);

s2=s2sequenceGARCH(a,x,m,s);
a2=a.^2;


if m>0
alfa=x(1:m);
poma=a2(T+1-m:T);
end
if s>0
beta=x(m+2:n);
poms=s2(T+1-s:T);
end

res=[];

if s==0
    if m==0
        %GARCH(0,0)
        res=[res;alfa0*ones(k,1)];
    else
        %GARCH(m,0)
        for i=1:k
            s2=alfa0+alfa'*poma;
            res=[res;s2];
            poma=[poma;s2];
            poma(1)=[];
        end
    end
else
    if m==0
        %GARC(0,s)
        for i=1:k
            s2=alfa0+beta'*poms;
            res=[res;s2];
            poms=[poms;s2];
            poms(1)=[];
        end
    else
        %GARCH(m,s)
        for i=1:k
            s2=alfa0+alfa'*poma+beta'*poms;
            res=[res;s2];
            poms=[poms;s2];
            poms(1)=[];
            poma=[poma;s2];
            poma(1)=[];            
        end
    end
end
        





