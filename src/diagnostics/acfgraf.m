function res=acfgraf(rniz,m)
y=[];
for l=1:m
    rl=acf(rniz,l);
    y=[y;rl];
end
plot(y,'g*');
res=y;