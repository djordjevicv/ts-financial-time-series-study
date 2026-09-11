function res=ljungbox2(rniz,m)
disp('H0 (ro_1=...=ro_m=0)')
disp('We assume confidence level 0.95')
% disp('(Nezavisnost od prethodnih m vrednosti.)')
chis95=[3.84 5.99 7.81 9.48 11.07 12.59 14.07 15.51 16.92 18.31];
z=chis95(m);
% chis99=[6.63 9.21 11.34 13.28 15.09 16.81 18.48 20.09 21.66 23.21];
% if alfa==0.05
%     z=chis95(m-g);
% elseif alfa==0.01
%     z=chis99(m-g);
% else 
%     disp('Mozete izbrati nivo znacajnosti 0.05 ili 0.01.')
% end
T=size(rniz,1);
s=0;
for i=1:m
    s=s+((acf(rniz,i))^2)/(T-i);
end
q=T*(T+2)*s;
if q>z
    disp('Reject H0')
    disp('There may be some linear dependence with some of the previous values.')
else
    disp('Do not reject H0 - there is no significant linear dependence.')
end
alfa=0.05;
res=[m;alfa;q;z];

