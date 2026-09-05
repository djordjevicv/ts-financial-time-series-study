function res=acfgraf2(r,m)
res=[];
for l=1:m
    ri=acf(r,l);
    res=[res;ri];
end