function res=s2sequenceGARCH(a,x,m,s)
%x=(alfam,...,alfa1,alfa0,betas,...,beta1)

c=max(m,s);
alfa0=x(m+1);
T=size(a,1);
n=size(x,1);

if m>0
alfa=x(1:m);
end
if s>0
beta=x(m+2:n);
end

pocetni=var(a);
poms=pocetni*ones(s,1);
res=[poms];

a2=a.^2;
if m<=s
    poma=a2(s+1-m:s);
else
    poma=[zeros(m-s,1);a2(1:s)];
end   

if s==0
    if m==0
        %GARCH(0,0)
        res=[res;alfa0*ones(T-s,1)];
    else
        %GARCH(m,0)
        for i=s+1:T
            s2=alfa0+alfa'*poma;
            res=[res;s2];
            poma=[poma;a2(i)];
            poma(1)=[];
        end
    end
else
    if m==0
        %GARCH(0,s)
        for i=s+1:T
            s2=alfa0+beta'*poms;
            res=[res;s2];
            poms=[poms;s2];
            poms(1)=[];
        end
    else
        %GARCH(m,s)
        for i=s+1:T
            s2=alfa0+alfa'*poma+beta'*poms;
            res=[res;s2];
            poms=[poms;s2];
            poms(1)=[];
            poma=[poma;a2(i)];
            poma(1)=[];            
        end
    end
end
        




