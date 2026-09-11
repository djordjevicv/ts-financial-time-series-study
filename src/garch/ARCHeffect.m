function res=ARCHeffect(a)
m=ceil(log(size(a,1)));
alfa=0.05;
niz=a.^2;
res=ljungbox2(niz,m);
% res=[m;alfa;q;z];
q=res(3);
z=res(4);
if q>z
    disp('Reject H0')
%     disp('Moguce je da postoji linearna zavisnost sa nekom od prethodnih m vrednosti')
    disp('ARCH effect is statistically significant')
else
%     disp('Ne postoji linearna zavisnot')
    disp('ARCH is not statistically significant')
end
