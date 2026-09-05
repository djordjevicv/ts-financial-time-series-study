function res=basicstat1(r)
T=size(r,1);
chis95=[3.84 5.99 7.81 9.48 11.07 12.59 14.07 15.51 16.92 18.31];
alfa=0.05;
% chis99=[6.63 9.21 11.34 13.28 15.09 16.81 18.48 20.09 21.66 23.21];
%2 degrees of freedom for chi^2 distribution
step=2;
%quantile chi^2_{0.95}
z=chis95(step);
zn=sqrt(2)*erfcinv(alfa);
    niz=r;
    m=momenat(niz,1);
    s2=centralni(niz,2);
    s=sqrt(s2);
    sem=s/sqrt(T);
    skew=centralni(niz,3)/s^3;
    seskew=sqrt(6/T);
    k=centralni(niz,4)/s^4;
    sek=sqrt(24/T);
    tskew=skew/seskew;
    tk=(k-3)/sek;
    tjb=tskew^2+tk^2;
    if abs(tskew)>abs(zn)
        test1=1;
        disp('test 1 - Reject H0');
    else
        test1=0;
        disp('test 1 - Do not reject H0');
    end
    if abs(tk)>abs(zn)
        test2=1;
        disp('test 2 - Reject H0');
    else
        test2=0;
        disp('test 2 - Do not reject H0');
    end
    if abs(tjb)>z
        test3=1;
        disp('test 3 - Reject H0');
    else
        test3=0;
        disp('test 3 - Do not reject H0');
    end
    
    res=[m sem s2 s skew seskew tskew zn k sek tk zn tjb z test1 test2 test3];
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
